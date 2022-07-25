#include <armadillo>
#include <gsl/gsl_sf_coupling.h>
#include "PhysicalConstants.hh"
#include "TwoBodyOperator.hh"
#include "Operator.hh"

Operator::~Operator()
{}

Operator::Operator()
{}

Operator::Operator(ModelSpace& ms, int rankJ, int rankP)
  : modelspace(&ms), rankJ(rankJ), rankP(rankP)
{
  ZeroBody = 0;
  int norbs = modelspace->GetNumberOrbits();
  OneBody = arma::mat(norbs, norbs, arma::fill::zeros);
  S = arma::mat(norbs, norbs, arma::fill::eye);
  TwoBody = TwoBodyOperator(*modelspace, rankJ, rankP);
  for (int i=0; i<norbs; i++){
    for (int j=i; j<norbs; j++){
      Orbit& oi = modelspace->GetOrbit(i);
      Orbit& oj = modelspace->GetOrbit(j);
      S(i,j) = MEOverlap(oi, oj, modelspace->GetZeta(), modelspace->GetProtonNumber());
      S(j,i) = S(i,j);
    }
  }
}

Operator::Operator(const Operator& op)
  : modelspace(op.modelspace), rankJ(op.rankJ), rankP(op.rankP), ZeroBody(op.ZeroBody), 
  S(op.S), OneBody(op.OneBody), TwoBody(op.TwoBody), orthonormalized(op.orthonormalized)
{}

Operator::Operator(Operator&& op)
  : modelspace(op.modelspace), rankJ(op.rankJ), rankP(op.rankP), ZeroBody(op.ZeroBody), 
  S(std::move(op.S)), OneBody(std::move(op.OneBody)), TwoBody(std::move(op.TwoBody)), orthonormalized(op.orthonormalized)
{}

Operator& Operator::operator=(const Operator& rhs) = default;
Operator& Operator::operator=(Operator&& rhs) = default;

Operator& Operator::operator*=(const double rhs)
{
  ZeroBody *= rhs;
  OneBody *= rhs;
  TwoBody *= rhs;
  return *this;
}

Operator Operator::operator*(const double rhs) const
{
  Operator opout = Operator(*this);
  opout *= rhs;
  return opout;
}

Operator operator*(const double lhs, const Operator& rhs)
{
   return rhs * lhs;
}
Operator operator*(const double lhs, const Operator&& rhs)
{
   return rhs * lhs;
}

Operator& Operator::operator/=(const double rhs)
{
  return *this *=(1.0/rhs);
}

Operator Operator::operator/(const double rhs) const
{
  Operator opout = Operator(*this);
  opout *= (1.0/rhs);
  return opout;
}

Operator& Operator::operator+=(const Operator& rhs)
{
  ZeroBody += rhs.ZeroBody;
  OneBody  += rhs.OneBody;
  TwoBody  += rhs.TwoBody;
  return *this;
}

Operator Operator::operator+(const Operator& rhs) const
{
  return ( Operator(*this) += rhs );
}

Operator& Operator::operator+=(const double& rhs)
{
  ZeroBody += rhs;
  return *this;
}

Operator Operator::operator+(const double& rhs) const
{
  return ( Operator(*this) += rhs );
}


Operator& Operator::operator-=(const Operator& rhs)
{
  ZeroBody -= rhs.ZeroBody;
  OneBody -= rhs.OneBody;
  TwoBody -= rhs.TwoBody;
  return *this;
}
  
Operator Operator::operator-(const Operator& rhs) const
{
  return ( Operator(*this) -= rhs );
}

Operator& Operator::operator-=(const double& rhs)
{
  ZeroBody -= rhs;
  return *this;
}

Operator Operator::operator-(const double& rhs) const
{
  return ( Operator(*this) -= rhs );
}

Operator Operator::operator-() const
{
  return (*this)*-1.0;
}


void Operator::SetDiracCoulombHamiltonian(bool OneBodyTerm, bool TwoBodyTerm)
{
  if(OneBodyTerm){ SetOneBodyDiracHamiltonian();}
  if(TwoBodyTerm){ SetTwoBodyCoulombInteraction();}
}

void Operator::SetCoulombHamiltonian(bool OneBodyTerm, bool TwoBodyTerm)
{
  if(OneBodyTerm){ SetOneBodyHamiltonian();}
  if(TwoBodyTerm){ SetTwoBodyCoulombInteraction();}
}

void Operator::SetOneBodyHamiltonian()
{
  Orbits orbits = modelspace->orbits;
  int norbs = orbits.GetNumberOrbits();
  double zeta = modelspace->GetZeta();
  double Z = modelspace->GetProtonNumber();
  for (int i1=0; i1<norbs; i1++){
    for (int i2=i1; i2<norbs; i2++){
      Orbit& o1 = orbits.GetOrbit(i1);
      Orbit& o2 = orbits.GetOrbit(i2);
      if(o1.kappa != o2.kappa) continue;
      OneBody(i1,i2) = MEKinetic(o1, o2, zeta, Z) - Z * MENuclPot(o1, o2, zeta, Z);
      OneBody(i2,i1) = OneBody(i1,i2);
    }
  }
}

void Operator::SetOneBodyDiracHamiltonian()
{
  Orbits orbits = modelspace->orbits;
  int norbs = orbits.GetNumberOrbits();
  double zeta = modelspace->GetZeta();
  double Z = modelspace->GetProtonNumber();
  for (int i1=0; i1<norbs; i1++){
    Orbit& o1 = orbits.GetOrbit(i1);
    for (int i2=i1; i2<norbs; i2++){
      Orbit& o2 = orbits.GetOrbit(i2);
      if(o1.kappa != o2.kappa) continue;
      double mass_term = 0;
      if(o1.ls ==-1) mass_term = -2*PhysConst::c * PhysConst::c*MEOverlap(o1, o2, zeta, Z);
      OneBody(i1,i2) = MEKinetic(o1, o2, zeta, Z) - Z * MENuclPot(o1, o2, zeta, Z) + mass_term;
      OneBody(i2,i1) = OneBody(i1,i2);
    }
  }
}

void Operator::SetTwoBodyCoulombInteraction()
{
  TwoBody.SetTwoBodyCoulombTerm();
}

/*
   (ab:J|U1 x U2|cd:J) = [ (a|U|c) (b|U|d)
   - (phase) (a|U|d) (b|U|c)
   - (phase) (a|U|d) (b|U|c)
   + (phase) (a|U|c) (b|U|d) ] / 2 sqrt((1+del_ab) (1+del_cd))
   */
arma::mat Operator::EmbedBasisTrans2(arma::mat T, TwoBodyChannel& tbc)
{
  Orbits& orbits = modelspace->GetOrbits();
  arma::mat U = arma::mat(tbc.GetNumberStates(), tbc.GetNumberStates(), arma::fill::zeros);
  for (int ileft=0; ileft < tbc.GetNumberStates(); ileft++){
    for (int iright=0; iright < tbc.GetNumberStates(); iright++){
      int i = tbc.GetOrbitIndex1(ileft);
      int j = tbc.GetOrbitIndex2(ileft);
      int k = tbc.GetOrbitIndex1(iright);
      int l = tbc.GetOrbitIndex2(iright);

      Orbit& oi = orbits.GetOrbit(i);
      Orbit& oj = orbits.GetOrbit(j);
      Orbit& ok = orbits.GetOrbit(k);
      Orbit& ol = orbits.GetOrbit(l);
      int phase_ij = pow( (-1), (oi.j2+oj.j2)/2-tbc.J );
      int phase_kl = pow( (-1), (ok.j2+ol.j2)/2-tbc.J );
      U(ileft, iright) = ((1+phase_ij*phase_kl) * T(i,k) * T(j,l) - (phase_ij+phase_kl) * T(i,l) * T(j,k)) * 0.5;
      if (i==j) {U(ileft, iright) /= sqrt(2);}
      if (k==l) {U(ileft, iright) /= sqrt(2);}
    }
  }
  return U;
}

void Operator::OrthoNormalize()
{
  if(orthonormalized){
    std::cout << "The operator is already orthonormalized; exiting..." << std::endl;
    return;
  }
  Orbits& orbits = modelspace->GetOrbits();
  TwoBodySpace& tbs = modelspace->GetTwoBodySpace();
  arma::mat L = arma::chol(S, "lower");
  arma::mat T = arma::inv(L);
  OneBody = T * OneBody * T.t();
  for (auto it : TwoBody.Channels){
    int ichbra = it.first[0];
    int ichket = it.first[1];
    TwoBodyOperatorChannel& opch = it.second;
    TwoBodyChannel * chbra = opch.chbra;
    TwoBodyChannel * chket = opch.chket;
    arma::mat Ubra = EmbedBasisTrans2(T, *chbra);
    arma::mat Uket;
    if(ichbra==ichket) {Uket = Ubra;}
    else {Uket = EmbedBasisTrans2(T, *chket);}
    auto& out = TwoBody.Channels[{ichbra,ichket}].MEs;
    out = Ubra * opch.MEs * Uket.t();
  }
  orthonormalized = true;
}

void Operator::Print()
{
  std::cout << "Zero-body part" << std::endl;
  std::cout << ZeroBody << std::endl;
  std::cout << "One-body part" << std::endl;
  std::cout << OneBody << std::endl;
  std::cout << "Two-body part" << std::endl;
  TwoBody.Print();
}

arma::vec Operator::DiagonalizeOneBody()
{
  int norbs = modelspace->GetNumberOrbits();
  arma::vec SPEs(norbs, arma::fill::zeros);
  OneBodySpace& obs = modelspace->GetOneBodySpace();
  for (int ich=0; ich< obs.GetNumberChannels(); ich++) {
    std::vector<unsigned long long> tmp(obs.channels[ich].begin(), obs.channels[ich].end());
    arma::uvec sub_idx(std::vector<unsigned long long>(tmp.begin(), tmp.end()));
    arma::mat Hch = OneBody(sub_idx, sub_idx);
    arma::vec eig;
    arma::mat vec;
    arma::eig_sym(eig, vec, Hch);
    SPEs(sub_idx) = eig;
  }
  return SPEs;
}

Operator Operator::DoNormalOrdering(int sign)
{
  std::cout << std::endl;
  std::cout << "Taking normal ordering..." << std::endl;
  modelspace->PrintHoleOrbits();
  Operator op(*this);
  if(op.rankJ==0 and op.rankP==1){
    for (auto hole_i : modelspace->hole_occ){
      int i = hole_i.first;
      Orbit& oi = modelspace->GetOrbit(i);
      double occ_i = hole_i.second;
      op.ZeroBody += (oi.j2+1) * sign * occ_i * OneBody(i,i);
      for (auto hole_j : modelspace->hole_occ){
        int j = hole_j.first;
        Orbit& oj = modelspace->GetOrbit(j);
        double occ_j = hole_j.second;
        if(i<j) continue;
        int Jmin = std::abs( oi.j2 - oj.j2)/2;
        int Jmax = (oi.j2 + oj.j2)/2;
        for (int J=Jmin; J<=Jmax; J++){
          if(i==j and J%2==1) continue;
          op.ZeroBody += (2*J+1) * occ_i * occ_j * TwoBody.Get2BME_J(J,J,i,j,i,j);
        }
      }
    }
  }

  arma::mat m = arma::mat(op.modelspace->GetNumberOrbits(), op.modelspace->GetNumberOrbits(), arma::fill::zeros);
  int norbs = modelspace->GetNumberOrbits();
  for (int a=0; a<norbs; a++){
    Orbit& oa = modelspace->GetOrbit(a);
    for (int b=0; b<norbs; b++){
      Orbit& ob = modelspace->GetOrbit(b);
      if(std::abs(oa.j2-ob.j2)/2 > op.rankJ or oa.j2+ob.j2 < op.rankJ) continue;
      if(pow(-1, oa.l+ob.l) != op.rankP) continue;

      for (auto hole_i : modelspace->hole_occ){
        int i = hole_i.first;
        Orbit& oi = modelspace->GetOrbit(i);
        double occ_i = hole_i.second;

        double norm = 1;
        if(a==i) norm *= sqrt(2);
        if(b==i) norm *= sqrt(2);
        int Jai_min = std::abs(oa.j2-oi.j2)/2;
        int Jai_max = (oa.j2+oi.j2)/2;
        int Jbi_min = std::abs(ob.j2-oi.j2)/2;
        int Jbi_max = (ob.j2+oi.j2)/2;
        for(int Jai=Jai_min; Jai<=Jai_max; Jai++){
          if(a==i and Jai%2==1) continue;
          for(int Jbi=Jbi_min; Jbi<=Jbi_max; Jbi++){
            if(b==i and Jbi%2==1) continue;
            if(std::abs(Jai-Jbi) > op.rankJ or Jai+Jbi < op.rankJ) continue;
            double Jfact = sqrt((2*Jai+1)*(2*Jbi+1));
            if(op.rankJ==0){
              op.OneBody(a,b) += Jfact * sign * occ_i * TwoBody.Get2BME_J(Jai,Jbi,a,i,b,i) * norm / (oa.j2+1);
              m(a,b) += Jfact * sign * occ_i * TwoBody.Get2BME_J(Jai,Jbi,a,i,b,i) * norm / (oa.j2+1);
            }
            else{
              double me = Jfact * sign * occ_i * pow(-1, (oa.j2+oi.j2)/2-Jbi+op.rankJ) *
                gsl_sf_coupling_6j(2*Jai, 2*Jbi, 2*op.rankJ, ob.j2, oa.j2, oi.j2) *
                TwoBody.Get2BME_J(Jai,Jbi,a,i,b,i) * norm;
              op.OneBody(a,b) += me;
            }
          }
        }
      }
    }
  }
  std::cout << std::endl;
  return op;
}
