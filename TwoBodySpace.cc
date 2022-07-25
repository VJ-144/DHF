#include <iomanip>
#include "Orbits.hh"
#include "TwoBodySpace.hh"

TwoBodyChannel::~TwoBodyChannel()
{}

TwoBodyChannel::TwoBodyChannel(int J, int Prty, Orbits& orbs)
  : J(J), Prty(Prty), orbits(&orbs)
{
  int num = 0;
  int norbs = orbits->GetNumberOrbits();
  for (int ip=0; ip<norbs; ip++){
    for (int iq=ip; iq<norbs; iq++){
      Orbit& op = orbits->GetOrbit(ip);
      Orbit& oq = orbits->GetOrbit(iq);
      if(pow((-1), (op.l+oq.l)) != Prty) continue;
      if(op.j2 + oq.j2 < 2*J) continue;
      if(std::abs(op.j2 - oq.j2) > 2*J) continue;
      if(ip==iq and J%2==1) continue;
      index1.push_back(ip);
      index2.push_back(iq);
      phase[{ip,iq}] = 1;
      index[{ip,iq}] = num;
      phase[{iq,ip}] = pow((-1), ((op.j2+oq.j2)/2-J+1));
      index[{iq,ip}] = num;
      num += 1;
    }
  }
  number_states = num;
}

void TwoBodyChannel::PrintChannel()
{
  std::cout << std::endl;
  int wint = 4;
  for (int idx=0; idx<GetNumberStates(); idx++)
  {
    int ip = GetOrbitIndex1(idx);
    int iq = GetOrbitIndex2(idx);
    Orbit& op = GetOrbit1(idx);
    Orbit& oq = GetOrbit2(idx);
    std::cout
      << "J =" << std::setw(wint) << J
      << ", Prty =" << std::setw(wint) << Prty
      << " | np =" << std::setw(wint) << op.n
      << ", lp =" << std::setw(wint) << op.l
      << ", j2p =" << std::setw(wint) << op.j2
      << ", lsp =" << std::setw(wint) << op.ls << " | "
      << "nq =" << std::setw(wint) << oq.n
      << ", lq =" << std::setw(wint) << oq.l
      << ", j2q =" << std::setw(wint) << oq.j2
      << ", lsq =" << std::setw(wint) << oq.ls << std::endl;
  }
}

TwoBodySpace::~TwoBodySpace()
{}

TwoBodySpace::TwoBodySpace()
{}

TwoBodySpace::TwoBodySpace(Orbits& orbs)
  : orbits(&orbs)
{
  int num = 0;
  for (int J=0; J<2*orbits->lmax+2; J++){
    for (int Prty : {1, -1}){
      TwoBodyChannel tbc = TwoBodyChannel(J, Prty, orbs);
      if(tbc.GetNumberStates() < 1) continue;
      index_from_JP[{J,Prty}] = num;
      channels.push_back(tbc);
      num += 1;
    }
  }
  number_channels = num;
}
