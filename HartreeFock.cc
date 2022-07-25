#include <iomanip>
#include <armadillo>

#include "TwoBodyOperator.hh"
#include "TwoBodySpace.hh"
#include "HartreeFock.hh"
#include "ModelSpace.hh"
#include "Operator.hh"
#include "Orbits.hh"

Monopole::~Monopole()
{}

Monopole::Monopole()
{}

Monopole::Monopole(Operator& H)
  : H(&H), modelspace(H.GetModelSpace())
{
  Orbits& orbits = modelspace->GetOrbits();
  int norbs = orbits.GetNumberOrbits();
  for (int i1=0; i1<norbs; i1++) {
    Orbit& o1 = orbits.GetOrbit(i1);
    for (int i3=0; i3<norbs; i3++) {
      Orbit& o3 = orbits.GetOrbit(i3);
      if (o1.kappa != o3.kappa) continue;
      for (int i2=0; i2<norbs; i2++) {
        for (int i4=0; i4<norbs; i4++) {
          Orbit& o2 = orbits.GetOrbit(i2);
          Orbit& o4 = orbits.GetOrbit(i4);
          if (o2.kappa != o4.kappa) continue;
          if (o1.ls + o2.ls != o3.ls + o4.ls) continue;
          idx_to_ijkl.push_back({i1,i2,i3,i4});
          double norm = 1.0;
          if (i1==i2) norm *= sqrt(2);
          if (i3==i4) norm *= sqrt(2);
          double v = 0.0;
          for (int J=std::abs(o1.j2-o2.j2)/2; J<=(o1.j2+o2.j2)/2; J++){
            if (i1==i2 and J%2==1) continue;
            if (i3==i4 and J%2==1) continue;
            v += (2*J+1) * H.Get2BME_J(J,J,o1,o2,o3,o4);
          }
          v *= norm / (o1.j2+1);
          vmon2.push_back(v);
        }
      }
    }
  }
}

void Monopole::Print() {
  for (int idx=0; idx<idx_to_ijkl.size(); idx++){
    std::cout << "i =" << std::setw(4) << idx_to_ijkl[idx][0]
      << ", j =" << std::setw(4) << idx_to_ijkl[idx][1]
      << ", k =" << std::setw(4) << idx_to_ijkl[idx][2]
      << ", l =" << std::setw(4) << idx_to_ijkl[idx][3]
      << ", Vmon =" << std::setw(12) << vmon2[idx] << std::endl;
  }
}

HartreeFock::~HartreeFock()
{}

HartreeFock::HartreeFock(Operator& H)
  : H(H), modelspace(H.GetModelSpace())
{
  if(not H.orthonormalized){
    std::cout << " OrthoNormalize the Hamiltonian first!" << std::endl;
    exit(0);
  }
  monopole = Monopole(H);
  Orbits orbits = modelspace->GetOrbits();
  int norbs = orbits.GetNumberOrbits();
  C = arma::mat(norbs, norbs, arma::fill::eye);
  rho = arma::mat(norbs, norbs, arma::fill::zeros);
  F = arma::mat(norbs, norbs, arma::fill::zeros);
  V = arma::mat(norbs, norbs, arma::fill::zeros);
  SPEs = arma::vec(norbs, arma::fill::zeros);
  S = H.S;

  r = UpdateFock();
  DiagonalizeFock();
  UpdateDensityMatrix();
  r = UpdateFock(0);
  CalcEnergy();
  PrintStatus(0);
}

void HartreeFock::Solve() {
  for (int n_iter=1; n_iter<1000; n_iter++) {
    DiagonalizeFock();
    UpdateDensityMatrix();
    r = UpdateFock(0);
    CalcEnergy();
    PrintStatus(n_iter);
    if (r < 1.e-8) break;
  }
}

double HartreeFock::UpdateFock(int n_itr)
{
  arma::mat Fock_old = F;
  Orbits& orbits = modelspace->GetOrbits();
  int norbs = orbits.GetNumberOrbits();
  V = arma::mat(norbs, norbs, arma::fill::zeros);
  if (n_itr != -1) {
    for (int idx=0; idx<monopole.idx_to_ijkl.size(); idx++) {
      int i = monopole.idx_to_ijkl[idx][0];
      int j = monopole.idx_to_ijkl[idx][1];
      int k = monopole.idx_to_ijkl[idx][2];
      int l = monopole.idx_to_ijkl[idx][3];
      V(i,k) += monopole.vmon2[idx] * rho(j,l);
    }
  }
  F = H.OneBody + V;
  double diff = 0.0;
  for (int i=0; i<norbs; i++) {
    for (int j=0; j<norbs; j++) {
      diff += sqrt( pow(( F(i,j) - Fock_old(i,j) ), 2) );
    }
  }
  return diff;
}

void HartreeFock::DiagonalizeFock() {
  OneBodySpace& obs = modelspace->GetOneBodySpace();
  for (int ich=0; ich< obs.GetNumberChannels(); ich++) {
    std::vector<unsigned long long> tmp(obs.channels[ich].begin(), obs.channels[ich].end());
    arma::uvec sub_idx(std::vector<unsigned long long>(tmp.begin(), tmp.end()));
    arma::mat Fch = F(sub_idx, sub_idx);
    arma::vec eig;
    arma::mat vec;
    arma::eig_sym(eig, vec, Fch);
    SPEs(sub_idx) = eig;
    C(sub_idx, sub_idx) = vec;
  }
}

void HartreeFock::UpdateDensityMatrix() {
  int norbs = modelspace->GetNumberOrbits();
  arma::mat tmp(norbs, norbs, arma::fill::zeros);
  modelspace->UpdateOrbitals(SPEs);
  for ( auto hole : modelspace->hole_occ) {
    tmp(hole.first,hole.first) = hole.second;
  }
  rho = C * tmp * C.t();
}

void HartreeFock::CalcEnergy() {
  Orbits& orbits = modelspace->GetOrbits();
  int norbs = orbits.GetNumberOrbits();
  E1 = 0.0;
  E2 = 0.0;
  for (int i=0; i<norbs; i++) {
    for (int j=0; j<norbs; j++) {
      Orbit& oi = orbits.GetOrbit(i);
      E1 += H.OneBody(i,j) * rho(i,j) * (oi.j2+1);
      E2 += V(i,j) * rho(i,j) * 0.5 * (oi.j2+1);
    }
  }
  EHF = E1+E2;
}

void HartreeFock::PrintStatus(int n_iter) {
  std::cout << "n iter:" << std::setw(4) << n_iter
    << ", OneBody: " << std::fixed << std::setw(12) << E1
    << ", TwoBody: " << std::fixed << std::setw(12) << E2
    << ", EHF: " << std::fixed << std::setw(12) << EHF
    << std::endl;
}

