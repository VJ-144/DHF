#include <iomanip>
#include <armadillo>
#include <gsl/gsl_sf_coupling.h>
#include "Orbits.hh"
#include "TwoBodySpace.hh"
#include "TwoBodyOperator.hh"

std::map<std::array<int,5>,double> CoulombIntegrals;

TwoBodyOperatorChannel::~TwoBodyOperatorChannel()
{}

TwoBodyOperatorChannel::TwoBodyOperatorChannel()
{}

TwoBodyOperatorChannel::TwoBodyOperatorChannel(TwoBodyChannel& chbra_in, TwoBodyChannel& chket_in)
  : chbra(&chbra_in), chket(&chket_in)
{
  MEs = arma::mat(chbra->GetNumberStates(), chket->GetNumberStates(), arma::fill::zeros);
}

TwoBodyOperator::~TwoBodyOperator()
{}

TwoBodyOperator::TwoBodyOperator()
{}

TwoBodyOperator::TwoBodyOperator(ModelSpace& ms, int rankJ, int rankP)
  : modelspace(&ms), rankJ(rankJ), rankP(rankP), NMesh(2000), Rmax(100), MeshType("Legendre")
{
  TwoBodySpace& tbs = modelspace->two;
  for (int ichbra=0; ichbra<tbs.GetNumberChannels(); ichbra++){
    for (int ichket=ichbra; ichket<tbs.GetNumberChannels(); ichket++){
      TwoBodyChannel& chbra = tbs.GetChannel(ichbra);
      TwoBodyChannel& chket = tbs.GetChannel(ichket);
      if(std::abs(chbra.J - chket.J) > rankJ or chbra.J+chket.J < rankJ) continue;
      if(chbra.Prty * chket.Prty * rankP == -1) continue;
      Channels[{ichbra,ichket}] = TwoBodyOperatorChannel(chbra, chket);
    }
  }
}

TwoBodyOperator& TwoBodyOperator::operator*=(const double rhs)
{
  for ( auto& itmat : Channels )
  {
    itmat.second.MEs *= rhs;
  }
  return *this;
}

TwoBodyOperator& TwoBodyOperator::operator+=(const TwoBodyOperator& rhs)
{
  for ( auto& itmat : rhs.Channels )
  {
    int ch_bra = itmat.first[0];
    int ch_ket = itmat.first[1];
    Channels[{ch_bra,ch_ket}].MEs += itmat.second.MEs;
  }
  return *this;
}

TwoBodyOperator& TwoBodyOperator::operator-=(const TwoBodyOperator& rhs)
{
  for ( auto& itmat : rhs.Channels )
  {
    auto ch_bra = itmat.first[0];
    auto ch_ket = itmat.first[1];
    Channels[{ch_bra,ch_ket}].MEs -= itmat.second.MEs;
  }
  return *this;
}



double TwoBodyOperator::Get2BME(int ichbra, int ichket, int i, int j, int k, int l)
{
  TwoBodySpace& tbs = modelspace->two;
  TwoBodyChannel& chbra = tbs.GetChannel(ichbra);
  TwoBodyChannel& chket = tbs.GetChannel(ichket);
  int idxbra = chbra.GetIndex(i,j);
  int idxket = chket.GetIndex(k,l);
  if(i==j and chbra.J%2==1) return 0;
  if(k==l and chket.J%2==1) return 0;
  if(chbra.Prty * chket.Prty * rankP == -1) return 0;
  if(std::abs(chbra.J - chket.J) > rankJ or chbra.J+chket.J < rankJ) return 0;
  double phase = 1;
  phase *= chbra.GetPhaseFactor(i,j);
  phase *= chket.GetPhaseFactor(k,l);
  if(ichbra>ichket){
    return Channels[{ichket,ichbra}].MEs(idxket,idxbra) * pow((-1), chbra.J-chket.J) * phase;
  }
  return Channels[{ichbra,ichket}].MEs(idxbra,idxket) * phase;
}

double TwoBodyOperator::Get2BME(int ichbra, int ichket, Orbit& oi, Orbit& oj, Orbit& ok, Orbit& ol)
{
  return Get2BME(ichbra, ichket, oi.idx, oj.idx, ok.idx, ol.idx);
}

double TwoBodyOperator::Get2BME_J(int Jbra, int Jket, Orbit& oi, Orbit& oj, Orbit& ok, Orbit& ol)
{
  TwoBodySpace & tbs = modelspace->two;
  int Prty_bra = pow((-1), oi.l+oj.l);
  int Prty_ket = pow((-1), ok.l+ol.l);
  int ichbra = tbs.GetChannelIndex(Jbra, Prty_bra);
  int ichket = tbs.GetChannelIndex(Jket, Prty_ket);
  return Get2BME(ichbra, ichket, oi, oj, ok, ol);
}

double TwoBodyOperator::Get2BME_J(int Jbra, int Jket, int i, int j, int k, int l)
{
  Orbits& orbits = modelspace->GetOrbits();
  Orbit& oi = orbits.GetOrbit(i);
  Orbit& oj = orbits.GetOrbit(j);
  Orbit& ok = orbits.GetOrbit(k);
  Orbit& ol = orbits.GetOrbit(l);
  return Get2BME_J(Jbra, Jket, oi, oj, ok, ol);
}

void TwoBodyOperator::Set2BME(int ichbra, int ichket, int i, int j, int k, int l, double me)
{
  TwoBodySpace& tbs = modelspace->GetTwoBodySpace();
  TwoBodyChannel& chbra = tbs.GetChannel(ichbra);
  TwoBodyChannel& chket = tbs.GetChannel(ichket);
  if(i==j and chbra.J%2==1) return;
  if(k==l and chket.J%2==1) return;
  if(chbra.Prty * chket.Prty * rankP == -1) return;
  if(std::abs(chbra.J - chket.J) > rankJ or chbra.J+chket.J < rankJ) return;
  int idxbra = chbra.GetIndex(i,j);
  int idxket = chket.GetIndex(k,l);
  double phase = 1;
  phase *= chbra.GetPhaseFactor(i,j);
  phase *= chket.GetPhaseFactor(k,l);
  if(ichbra>ichket){
    Channels[{ichket,ichbra}].MEs(idxket,idxbra)  = me * pow((-1), chbra.J-chket.J) * phase;
    return;
  }
  Channels[{ichbra,ichket}].MEs(idxbra,idxket) = me * phase;
  if(ichbra != ichket) return;
  Channels[{ichbra,ichket}].MEs(idxket,idxbra) = me * phase;
}

void TwoBodyOperator::Set2BME(int ichbra, int ichket, Orbit& oi, Orbit& oj, Orbit& ok, Orbit& ol, double me)
{
  Orbits& orbits = modelspace->orbits;
  int i = orbits.GetOrbitIndex(oi);
  int j = orbits.GetOrbitIndex(oj);
  int k = orbits.GetOrbitIndex(ok);
  int l = orbits.GetOrbitIndex(ol);
  Set2BME(ichbra, ichket, i, j, k, l, me);
}

void TwoBodyOperator::Set2BME_J(int Jbra, int Jket, Orbit& oi, Orbit& oj, Orbit& ok, Orbit& ol, double me)
{
  int Prty_bra = pow((-1), oi.l+oj.l);
  int Prty_ket = pow((-1), ok.l+ol.l);
  TwoBodySpace & tbs = modelspace->two;
  int ichbra = tbs.GetChannelIndex(Jbra, Prty_bra);
  int ichket = tbs.GetChannelIndex(Jket, Prty_ket);
  Set2BME(ichbra, ichket, oi, oj, ok, ol, me);
}

void TwoBodyOperator::Set2BME_J(int Jbra, int Jket, int i, int j, int k, int l, double me)
{
  Orbits& orbits = modelspace->orbits;
  Orbit& oi = orbits.GetOrbit(i);
  Orbit& oj = orbits.GetOrbit(j);
  Orbit& ok = orbits.GetOrbit(k);
  Orbit& ol = orbits.GetOrbit(l);
  Set2BME_J(Jbra, Jket, oi, oj, ok, ol, me);
}

void TwoBodyOperator::Print()
{
  for (auto it : Channels){
    int ichbra = it.first[0];
    int ichket = it.first[1];
    TwoBodyOperatorChannel& op_ch = it.second;
    TwoBodyChannel * chbra = op_ch.chbra;
    TwoBodyChannel * chket = op_ch.chket;
    std::cout
      << "Jbra=" << std::setw(4) << chbra->J
      << ", Prty bra=" << std::setw(4) << chbra->Prty
      << " | Jket=" << std::setw(4) << chket->J
      << ", Prty ket=" << std::setw(4) << chket->Prty << std::endl;
    arma::mat& m = op_ch.MEs;
    std::cout << m << std::endl;
  }
}

void TwoBodyOperator::SetTwoBodyCoulombTerm()
{
  if(std::ifstream(FileNameCoulomb).good()){
    ReadFileInteraction(FileNameCoulomb);
    return;
  }
  StoreCoulombIntegrals();
  TwoBodySpace& tbs = modelspace->GetTwoBodySpace();
  Orbits& orbits = modelspace->GetOrbits();
  for (int ich=0; ich<tbs.GetNumberChannels(); ich++){
    TwoBodyChannel& tbc = tbs.GetChannel(ich);
    #pragma omp parallel for
    for(int idxbra=0; idxbra<tbc.GetNumberStates(); idxbra++){
      for(int idxket=idxbra; idxket<tbc.GetNumberStates(); idxket++){
        int i = tbc.GetOrbitIndex1(idxbra);
        int j = tbc.GetOrbitIndex2(idxbra);
        int k = tbc.GetOrbitIndex1(idxket);
        int l = tbc.GetOrbitIndex2(idxket);
        Orbit& oi = orbits.GetOrbit(i);
        Orbit& oj = orbits.GetOrbit(j);
        Orbit& ok = orbits.GetOrbit(k);
        Orbit& ol = orbits.GetOrbit(l);
        double norm = 1.0;
        if (i==j) {norm /= sqrt(2);}
        if (k==l) {norm /= sqrt(2);}
        double v = 0;
        v = MECoulomb(oi, oj, ok, ol, tbc.J);
        v += MECoulomb(oi, oj, ol, ok, tbc.J) * pow((-1), ( (ok.j2+ol.j2)/2 - tbc.J + 1));
        //v = TestMECoulomb(oi, oj, ok, ol, tbc.J);
        //v += TestMECoulomb(oi, oj, ol, ok, tbc.J) * pow((-1), ( (ok.j2+ol.j2)/2 - tbc.J + 1));
        v *= norm;
        Set2BME(ich, ich, i, j, k, l, v);
      }
    }
  }
  if(not std::ifstream(FileNameCoulomb).good()){
    WriteFileInteraction(FileNameCoulomb);
  }
}

double TwoBodyOperator::MECoulomb(Orbit& o1, Orbit& o2, Orbit& o3, Orbit& o4, int J)
{
  if (o1.ls != o3.ls) return 0.0;
  if (o2.ls != o4.ls) return 0.0;
  double zeta = modelspace->GetZeta();
  double Z = modelspace->GetProtonNumber();
  int Lmin = std::max(std::abs(o1.j2-o3.j2), std::abs(o2.j2-o4.j2))/2;
  int Lmax = std::min(        (o1.j2+o3.j2),         (o2.j2+o4.j2))/2;
  double r = 0.0;
  for (int L=Lmin; L<=Lmax; L++) {
    if ( (o1.l+o3.l+L)%2 == 1 ) continue;
    if ( (o2.l+o4.l+L)%2 == 1 ) continue;
    if (abs(o1.l-o3.l) > L or o1.l+o3.l < L) continue;
    if (abs(o2.l-o4.l) > L or o2.l+o4.l < L) continue;
    double angular = gsl_sf_coupling_6j(o1.j2, o2.j2, 2*J, o4.j2, o3.j2, 2*L) *
      gsl_sf_coupling_3j(o1.j2, 2*L, o3.j2, -1, 0, 1) *
      gsl_sf_coupling_3j(o2.j2, 2*L, o4.j2, -1, 0, 1);
    if(std::abs(angular) < 1.e-8) continue;
    r += angular * GetCoulombIntegral(o1.idx,o2.idx,o3.idx,o4.idx,L);
  }
  r *= sqrt( (o1.j2+1) * (o2.j2+1) * (o3.j2+1) * (o4.j2+1)) * pow( (-1), (o1.j2+o3.j2)/2+J );
  //r *= sqrt( (o1.j2+1) * (o2.j2+1) * (o3.j2+1) * (o4.j2+1)) * pow( (-1), (o1.j2+o3.j2)/2+J ) * o1.ls * o2.ls;
  return r;
}

double TwoBodyOperator::TestMECoulomb(Orbit& o1, Orbit& o2, Orbit& o3, Orbit& o4, int J)
{
  if (o1.ls != o3.ls) return 0.0;
  if (o2.ls != o4.ls) return 0.0;
  double zeta = modelspace->GetZeta();
  double Z = modelspace->GetProtonNumber();
  int Lmin = std::max(std::abs(o1.j2-o3.j2), std::abs(o2.j2-o4.j2))/2;
  int Lmax = std::min(        (o1.j2+o3.j2),         (o2.j2+o4.j2))/2;
  gsl_integration_fixed_workspace *workspace;
  const gsl_integration_fixed_type *T = gsl_integration_fixed_legendre;
  bool norm_weight;
  if(MeshType=="Legendre"){
    norm_weight = false;
    T = gsl_integration_fixed_legendre;
    workspace = gsl_integration_fixed_alloc(T, NMesh, 0.0, Rmax, 0.0, 0.0);
  }
  else if(MeshType=="Laguerre"){
    T = gsl_integration_fixed_laguerre;
    bool norm_weight = true;
    workspace = gsl_integration_fixed_alloc(T, NMesh, 0.0, 2.0/zeta, 0.0, 0.0);
  }
  double r = 0.0;
  for (int L=Lmin; L<=Lmax; L++) {
    if ( (o1.l+o3.l+L)%2 == 1 ) continue;
    if ( (o2.l+o4.l+L)%2 == 1 ) continue;
    if (abs(o1.l-o3.l) > L or o1.l+o3.l < L) continue;
    if (abs(o2.l-o4.l) > L or o2.l+o4.l < L) continue;
    double angular = gsl_sf_coupling_6j(o1.j2, o2.j2, 2*J, o4.j2, o3.j2, 2*L) *
      gsl_sf_coupling_3j(o1.j2, 2*L, o3.j2, -1, 0, 1) *
      gsl_sf_coupling_3j(o2.j2, 2*L, o4.j2, -1, 0, 1);
    if(std::abs(angular) < 1.e-8) continue;
    double Integral = 0.0;
    for (int i=0; i<NMesh; i++){
      for (int j=0; j<NMesh; j++){
        double x = workspace->x[i];
        double wx = workspace->weights[i];
        double y = workspace->x[j];
        double wy = workspace->weights[j];
        Integral += wx * wy * 
          o1.RadialFunction(x, zeta, Z, norm_weight) *
          o2.RadialFunction(y, zeta, Z, norm_weight) *
          o3.RadialFunction(x, zeta, Z, norm_weight) *
          o4.RadialFunction(y, zeta, Z, norm_weight) *
          pow( std::min(x,y), L) / pow( std::max(x,y), (L+1) );
      }
    }
    r += angular * Integral;
  }
  r *= sqrt( (o1.j2+1) * (o2.j2+1) * (o3.j2+1) * (o4.j2+1)) * pow( (-1), (o1.j2+o3.j2)/2+J );
  return r;
}


void TwoBodyOperator::StoreCoulombIntegrals()
{
  Orbits& orbits = modelspace->GetOrbits();
  int nmax = orbits.GetNmax();
  int kappa_min = orbits.GetKappaMin();
  int kappa_max = orbits.GetKappaMax();
  double zeta = modelspace->GetZeta();
  double Z = modelspace->GetProtonNumber();


  gsl_integration_fixed_workspace *workspace;
  const gsl_integration_fixed_type *T_leg = gsl_integration_fixed_legendre;
  const gsl_integration_fixed_type *T_lag = gsl_integration_fixed_laguerre;
  bool norm_weight;
  if(MeshType=="Legendre"){
    int NMesh_lag = 100;
    if(NMesh < NMesh_lag) {
      std::cout << "Increase NMesh!" << std::endl;
      exit(0);
    }
    norm_weight = false;
    workspace = gsl_integration_fixed_alloc(T_leg, NMesh, 0.0, Rmax, 0.0, 0.0);
    gsl_integration_fixed_workspace *workspace_tmp1;
    gsl_integration_fixed_workspace *workspace_tmp2;
    workspace_tmp1 = gsl_integration_fixed_alloc(T_leg, NMesh-NMesh_lag, 0.0, Rmax, 0.0, 0.0);
    workspace_tmp2 = gsl_integration_fixed_alloc(T_lag, NMesh_lag, Rmax, 2.0/zeta, 0.0, 0.0);
    for (int i = 0; i<NMesh; i++){
      if(i < NMesh-NMesh_lag){
        workspace->x[i] = workspace_tmp1->x[i];
        workspace->weights[i] = workspace_tmp1->weights[i];
      }
      else{
        workspace->x[i] = workspace_tmp2->x[i-NMesh+NMesh_lag];
        workspace->weights[i] = workspace_tmp2->weights[i-NMesh+NMesh_lag] * exp(2.0*(workspace->x[i]-Rmax)/zeta);
      }
    }
  }
  else if(MeshType=="Laguerre"){
    norm_weight = true;
    workspace = gsl_integration_fixed_alloc(T_lag, NMesh, 0.0, 2.0/zeta, 0.0, 0.0);
  }

  //for (int i=0; i<NMesh; i++){
  //std::cout << std::setw(6) << i
  //  << std::setw(16) << std::setprecision(8) << workspace->x[i]
  //  << std::setw(16) << std::setprecision(8) << workspace->weights[i]
  //  << std::endl;
  //}

  std::vector<int> i1_list;
  std::vector<int> i3_list;
  std::map<std::array<int,2>,int> Index;
  int cnt = 0;
  for (int i1=0; i1<orbits.GetNumberOrbits(); i1++){
    for (int i3=0; i3<orbits.GetNumberOrbits(); i3++){
      Orbit & o1 = orbits.GetOrbit(i1);
      Orbit & o3 = orbits.GetOrbit(i3);
      if(o1.ls != o3.ls) continue;
      i1_list.push_back(i1);
      i3_list.push_back(i3);
      Index[{i1,i3}] = cnt;
      cnt += 1;
    }
  }

  std::map<std::array<int,3>,double> tmp;
  #pragma omp parallel for
  for (int i2=0; i2<NMesh; i2++){
    double r2 = workspace->x[i2];
    for (int iket=0; iket<Index.size(); iket++){
      int i1 = i1_list[iket];
      int i3 = i3_list[iket];
      Orbit & o1 = orbits.GetOrbit(i1);
      Orbit & o3 = orbits.GetOrbit(i3);
      std::vector<double> Rnl;
      for (int i1=0; i1<NMesh; i1++){
        double r = workspace->x[i1];
        double w = workspace->weights[i1];
        double rnl = w 
          * o1.RadialFunction(r, zeta, Z, norm_weight) 
          * o3.RadialFunction(r, zeta, Z, norm_weight);
        Rnl.push_back(rnl);
      }

      for (int L = std::max(0,std::abs(o1.l-o3.l)-1); L<=o1.l+o3.l+1; L++) {
        double Int = 0.0;
        for (int i1=0; i1<NMesh; i1++){
          double r1 = workspace->x[i1];
          Int += pow(std::min(r1,r2),L) / pow(std::max(r1,r2),L+1) * Rnl[i1];
        }
        #pragma omp critical
        tmp[{i2,L,iket}] = Int;
      }
    }
  }

  #pragma omp parallel for
  for (int ibra=0; ibra<Index.size(); ibra++){
    int i1 = i1_list[ibra];
    int i3 = i3_list[ibra];
    Orbit & o1 = orbits.GetOrbit(i1);
    Orbit & o3 = orbits.GetOrbit(i3);
    for (int iket=0; iket<Index.size(); iket++){
      int i2 = i1_list[iket];
      int i4 = i3_list[iket];
      Orbit & o2 = orbits.GetOrbit(i2);
      Orbit & o4 = orbits.GetOrbit(i4);

      int Lmin = std::max(0, std::max(std::abs(o1.l-o3.l)-1, std::abs(o2.l-o4.l)-1));
      int Lmax = std::min(o1.l+o3.l+1, o2.l+o4.l+1);
      for (int L=Lmin; L<=Lmax; L++){
        if((o1.l+o3.l+L)%2==1) continue;
        if((o2.l+o4.l+L)%2==1) continue;
        double Int = 0.0;
        for (int i=0; i<NMesh; i++){
          double r = workspace->x[i];
          double w = workspace->weights[i];
          Int += tmp[{i,L,iket}] * w * 
            o1.RadialFunction(r, zeta, Z, norm_weight) * 
            o3.RadialFunction(r, zeta, Z, norm_weight);
        }
        #pragma omp critical
        CoulombIntegrals[{i1,i2,i3,i4,L}] = Int;
      }
    }
  }
}

double TwoBodyOperator::GetCoulombIntegral(int i1, int i2, int i3, int i4, int L)
{
  auto it = CoulombIntegrals.find({i1,i2,i3,i4,L});
  if(it != CoulombIntegrals.end()){
    return it->second;
  }
  std::cout << "Something wrong..." << std::endl;
  return 0.0;
}

void TwoBodyOperator::ReadFileInteraction(std::string filename)
{
  if(filename=="") return;
  if(not std::ifstream(FileNameCoulomb).good()){
    std::cout << "File not found: " << FileNameCoulomb << std::endl;
  }
  long long int num_MEs = CountInteractionMEs();
  std::vector<double> MEs(num_MEs);
  std::ifstream input(filename, std::ios::in | std::ios::binary);
  input.read((char*) &MEs[0], num_MEs * sizeof(double));
  input.close();

  long long int cnt = 0;
  double zeta = modelspace->GetZeta();
  Orbits& orbits = modelspace->orbits;
  for (int i1=0; i1<orbits.GetNumberOrbits(); i1++){
    Orbit & o1 = orbits.GetOrbit(i1);
    for (int i2=i1; i2<orbits.GetNumberOrbits(); i2++){
      Orbit & o2 = orbits.GetOrbit(i2);

      for (int i3=0; i3<orbits.GetNumberOrbits(); i3++){
        Orbit & o3 = orbits.GetOrbit(i3);
        for (int i4=i3; i4<orbits.GetNumberOrbits(); i4++){
          Orbit & o4 = orbits.GetOrbit(i4);

          if(o1.ls + o2.ls != o3.ls + o4.ls) continue;
          if((o1.l + o2.l + o3.l + o4.l)%2 == 1) continue;
          int Jmin = std::max(std::abs(o1.j2-o2.j2), std::abs(o3.j2-o4.j2))/2;
          int Jmax = std::min(std::abs(o1.j2+o2.j2), std::abs(o3.j2+o4.j2))/2;
          for (int J=Jmin; J<=Jmax; J++){
            Set2BME_J(J,J,i1,i2,i3,i4,MEs[cnt]/zeta);
            cnt += 1;
          }
        }
      }
    }
  }
}

void TwoBodyOperator::WriteFileInteraction(std::string filename)
{
  if(filename=="") return;
  double zeta = modelspace->GetZeta();
  long long int num_MEs = CountInteractionMEs();
  std::vector<double> MEs(num_MEs,0);
  long long int cnt = 0;
  Orbits& orbits = modelspace->orbits;
  for (int i1=0; i1<orbits.GetNumberOrbits(); i1++){
    Orbit & o1 = orbits.GetOrbit(i1);
    for (int i2=i1; i2<orbits.GetNumberOrbits(); i2++){
      Orbit & o2 = orbits.GetOrbit(i2);

      for (int i3=0; i3<orbits.GetNumberOrbits(); i3++){
        Orbit & o3 = orbits.GetOrbit(i3);
        for (int i4=i3; i4<orbits.GetNumberOrbits(); i4++){
          Orbit & o4 = orbits.GetOrbit(i4);

          if(o1.ls + o2.ls != o3.ls + o4.ls) continue;
          if((o1.l + o2.l + o3.l + o4.l)%2 == 1) continue;
          int Jmin = std::max(std::abs(o1.j2-o2.j2), std::abs(o3.j2-o4.j2))/2;
          int Jmax = std::min(std::abs(o1.j2+o2.j2), std::abs(o3.j2+o4.j2))/2;
          for (int J=Jmin; J<=Jmax; J++){
            MEs[cnt] = Get2BME_J(J,J,i1,i2,i3,i4)*zeta;
            cnt += 1;
          }
        }
      }
    }
  }
  std::ofstream output(filename, std::ios::out | std::ios::binary);
  output.write((char*) &MEs[0], num_MEs * sizeof(double));
  output.close();
}

long long int TwoBodyOperator::CountInteractionMEs()
{
  long long int cnt=0;
  Orbits& orbits = modelspace->orbits;
  for (int i1=0; i1<orbits.GetNumberOrbits(); i1++){
    Orbit & o1 = orbits.GetOrbit(i1);
    for (int i2=i1; i2<orbits.GetNumberOrbits(); i2++){
      Orbit & o2 = orbits.GetOrbit(i2);

      for (int i3=0; i3<orbits.GetNumberOrbits(); i3++){
        Orbit & o3 = orbits.GetOrbit(i3);
        for (int i4=i3; i4<orbits.GetNumberOrbits(); i4++){
          Orbit & o4 = orbits.GetOrbit(i4);

          if(o1.ls + o2.ls != o3.ls + o4.ls) continue;
          if((o1.l + o2.l + o3.l + o4.l)%2 == 1) continue;
          int Jmin = std::max(std::abs(o1.j2-o2.j2), std::abs(o3.j2-o4.j2))/2;
          int Jmax = std::min(std::abs(o1.j2+o2.j2), std::abs(o3.j2+o4.j2))/2;
          for (int J=Jmin; J<=Jmax; J++){
            cnt += 1;
          }
        }
      }
    }
  }
  return cnt;
}
