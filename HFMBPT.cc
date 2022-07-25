#include "PhysicalConstants.hh"
#include "HFMBPT.hh"

HFMBPT::~HFMBPT()
{}

HFMBPT::HFMBPT(Operator & H, arma::mat C_CB2HF)
  : H(H), modelspace(H.GetModelSpace()), C_CB2HF(C_CB2HF)
{
  C_trans = C_CB2HF;
}

void HFMBPT::SetTransformationMatrix(std::string mode)
{
  if(mode=="CB2HF") C_trans = C_CB2HF;
  if(mode=="HF2NAT") C_trans = C_HF2NAT;
  if(mode=="CB2NAT") C_trans = C_CB2NAT;
}

Operator HFMBPT::TransformBasis(Operator& op)
{
  std::cout << "Operator basis transformation..." << std::endl;
  Operator opout = op;
  opout.OneBody = C_trans.t() * op.OneBody * C_trans;

  for(auto& it : op.TwoBody.Channels)
  {
    int ich_bra = it.first[0];
    int ich_ket = it.first[1];
    TwoBodyChannel& tbc_bra = op.modelspace->GetTwoBodyChannel(ich_bra);
    TwoBodyChannel& tbc_ket = op.modelspace->GetTwoBodyChannel(ich_ket);
    auto& IN = it.second.MEs;
    int nbras = IN.n_rows;
    int nkets = IN.n_cols;
    arma::mat Dbra(nbras,nbras);
    arma::mat Dket(nkets,nkets);

    for (int i_old = 0; i_old<nkets; ++i_old)
    {
      int i1 = tbc_ket.GetOrbitIndex1(i_old);
      int i2 = tbc_ket.GetOrbitIndex2(i_old);
      Orbit& o1 = op.modelspace->GetOrbit(i1);
      Orbit& o2 = op.modelspace->GetOrbit(i2);
      for (int i_new = 0; i_new<nkets; ++i_new)
      {
        int i3 = tbc_ket.GetOrbitIndex1(i_new);
        int i4 = tbc_ket.GetOrbitIndex2(i_new);
        Orbit& o3 = op.modelspace->GetOrbit(i3);
        Orbit& o4 = op.modelspace->GetOrbit(i4);
        Dket(i_old,i_new) = C_trans(i1,i3) * C_trans(i2,i4);
        if(i1 != i2)
        {
          Dket(i_old,i_new) += C_trans(i2,i3) * C_trans(i1,i4) 
            * tbc_ket.GetPhaseFactor(i2,i1); 
        }
        if (i1==i2)  Dket(i_old,i_new) *= PhysConst::SQRT2;
        if (i3==i4)  Dket(i_old,i_new) /= PhysConst::SQRT2;
      }
    }
    if (ich_bra == ich_ket) {
      Dbra = Dket.t();
    }
    else
    {
      for (int i_old=0; i_old<nbras; ++i_old)
      {
        int i1 = tbc_bra.GetOrbitIndex1(i_old);
        int i2 = tbc_bra.GetOrbitIndex2(i_old);
        Orbit& o1 = op.modelspace->GetOrbit(i1);
        Orbit& o2 = op.modelspace->GetOrbit(i2);
        for (int i_new=0; i_new<nbras; ++i_new)
        {
          int i3 = tbc_bra.GetOrbitIndex1(i_new);
          int i4 = tbc_bra.GetOrbitIndex2(i_new);
          Orbit& o3 = op.modelspace->GetOrbit(i3);
          Orbit& o4 = op.modelspace->GetOrbit(i4);
          Dbra(i_new,i_old) = C_trans(i3,i1) * C_trans(i4,i2);
          if (i3 != i4)
          {
            Dbra(i_new,i_old) += C_trans(i4,i1) * C_trans(i3,i2) 
              * tbc_bra.GetPhaseFactor(i4,i3);
          }
          if (i3==i4)  Dbra(i_new,i_old) *= PhysConst::SQRT2;
          if (i1==i2)  Dbra(i_new,i_old) /= PhysConst::SQRT2;
        }
      }
    }
    auto& OUT =  opout.TwoBody.Channels[{ich_bra,ich_ket}].MEs;
    OUT = Dbra * IN * Dket;
  }
  std::cout << "done" << std::endl;
  return opout;
}

double HFMBPT::GetMP2_Energy(Operator& op)
{
  double Emp2 = 0;
  for (int a : modelspace->particles){
    Orbit& oa = modelspace->GetOrbit(a);
    double ea = op.OneBody(a,a);
    for (int i : modelspace->holes){
      Orbit& oi = modelspace->GetOrbit(i);
      double ei = op.OneBody(i,i);
      if(std::abs(op.OneBody(i,a)) > 1.e-16){
        Emp2 += (oi.j2+1) * oi.occ * op.OneBody(i,a) * op.OneBody(a,i) / (ei-ea);
      }
      for (int b : modelspace->particles){
        if(a>b) continue;
        Orbit& ob = modelspace->GetOrbit(b);
        double eb = op.OneBody(b,b);
        for (int j : modelspace->holes){
          if(i>j) continue;
          Orbit& oj = modelspace->GetOrbit(j);
          double ej = op.OneBody(j,j);
          double denom = 1 / (ei + ej - ea - eb);
          
          int Jmin = std::max(std::abs(oi.j2-oj.j2), std::abs(oa.j2-ob.j2))/2;
          int Jmax = std::min(std::abs(oi.j2+oj.j2), std::abs(oa.j2+ob.j2))/2;
          for (int J=Jmin; J<=Jmax; J++){
            if(i==j and J%2) continue;
            if(a==b and J%2) continue;
            double me = op.TwoBody.Get2BME_J(J,J,i,j,a,b);
            Emp2 += (2*J+1) * oi.occ * oj.occ * me * me * denom;
          }
        }
      }
    }
  }
  return Emp2;
}
