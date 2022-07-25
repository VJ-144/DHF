#include <iomanip>
#include <set>
#include <unordered_set>
#include <iostream>
#include <armadillo>
#include "Orbits.hh"
#include "TwoBodySpace.hh"
#include "ModelSpace.hh"


OneBodySpace::~OneBodySpace()
{}

OneBodySpace::OneBodySpace()
{}

OneBodySpace::OneBodySpace(Orbits& orbs)
  : orbits(&orbs)
{
  std::set<int> tmp;
  for (Orbit& o : orbits->orbits) tmp.insert(o.kappa); 
  for (auto it : tmp) kappas.push_back(it);
  std::sort(kappas.begin(), kappas.end(), std::greater<int>());
  for (int channel_idx=0; channel_idx<kappas.size(); channel_idx++){
    std::vector <int> idxs;
    int n_large = 0;
    int n_small = 0;
    for (auto o : orbits->orbits){
      if (o.kappa != kappas[channel_idx]) continue;
      if (o.ls == 1) n_large += 1;
      if (o.ls ==-1) n_small += 1;
      idxs.push_back(orbits->GetOrbitIndex(o));
      orbit_index_to_channel_index[orbits->GetOrbitIndex(o)] = channel_idx;
    }
    channels.push_back(idxs);
    n_large_channel.push_back(n_large);
    n_small_channel.push_back(n_small);
  }
  number_channels = channels.size();
};

void OneBodySpace::PrintSpace()
{
  for (int idx=0; idx<GetNumberChannels(); idx++)
  {
    std::cout << std::endl;
    std::vector<int> Indices = channels[idx];
    for (int i : Indices){
      Orbit& o = orbits->GetOrbit(i);
      o.Print();
    }
  }
}

ModelSpace::~ModelSpace()
{}

ModelSpace::ModelSpace(int Ne, double Z, double zeta, Orbits orbs)
  : Ne(Ne), Z(Z), zeta(zeta), orbits(orbs)
{
  hole_occ = AssignHoles(Ne);
  one = OneBodySpace(orbits);
  two = TwoBodySpace(orbits);
  PrintHoleOrbits();
}

ModelSpace::ModelSpace(std::string atom, double zeta, Orbits orbs)
  : zeta(zeta), orbits(orbs)
{
  std::vector<std::string> periodic_table = {
    "NA",
    "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne",
    "Na", "Mg", "Al", "Si", "P",  "S",  "Cl", "Ar", "K",  "Ca",
    "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y",  "Zr",
    "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
    "Sb", "Te", "I",  "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
    "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
    "Lu", "Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
    "Pa", "U",  "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm",
    "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds",
    "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og" };
  
  int i = 0;
  while (not isdigit(atom[i]) and atom[i] != '+' and atom[i] != '-') i++;
  std::string element = atom.substr(0,i);
  auto it_elem = find(periodic_table.begin(),periodic_table.end(),element);
  Z = it_elem - periodic_table.begin();
  if(element==atom){
    Ne = Z;
  }
  else if (atom[i]=='+') {
    int n;
    std::stringstream( atom.substr(i+1,atom.size()-i-1)) >> n;
    Ne = Z - n;
  }
  else if (atom[i]=='-') {
    int n;
    std::stringstream( atom.substr(i+1,atom.size()-i-1)) >> n;
    Ne = Z + n;
  }
  else {
    std::cout << "Warning: Unknown format of atom" << std::endl;
  }
  hole_occ = AssignHoles(Ne);
  one = OneBodySpace(orbits);
  two = TwoBodySpace(orbits);
  PrintHoleOrbits();
}

std::map<int,double> ModelSpace::AssignHoles(int N_ele)
{
  //
  // K:  0s
  // L:  1s                                     1p 1p 1p
  // M:  2s                                     2p 2p 2p
  // N:  3s                      2d 2d 2d 2d 2d 3p 3p 3p
  // O:  4s                      3d 3d 3d 3d 3d 4p 4p 4p
  // P:  5s 3f 3f 3f 3f 3f 3f 3f 4d 4d 4d 4d 4d 5p 5p 5p
  // Q:  6s 4f 4f 4f 4f 4f 4f 4f 5d 5d 5d 5d 5d 6p 6p 6p
  //
  // e = n + l
  std::map<int,double> tmp_holes;
  int N = 0;
  for (int e=0; e<=orbits.GetNmax(); ++e) {
    int d = std::min(N_ele-N, 2); // n=e, l=0
    tmp_holes[orbits.GetOrbitIndex(e,0,1,1)] = d * 0.5;
    N += d;
    if(N==N_ele) return tmp_holes;
    for (int l=std::min(e,orbits.GetLmax()); l>=1; --l){
      int n = e - l - std::max(0,l-1);
      if(n > orbits.GetNmax() or n < 0) continue;
      for (int j2=std::abs(2*l-1); j2<=2*l+1; j2+=2)
      {
        d = std::min(N_ele-N, j2+1);
        tmp_holes[orbits.GetOrbitIndex(n,l,j2,1)] = d / (j2+1.0);
        N += d;
        if(N==N_ele) return tmp_holes;
      }
    }
  }
  std::cout << "Something seems wrong in AssignHoles, N=" << N << std::endl;
  std::cout << "Try again with a larger model space" << std::endl;
  exit(0);
  return tmp_holes;
}


void ModelSpace::PrintModelSpace(bool print_obs, bool print_tbs)
{
  std::cout << " Number of electrons: " << GetElectronNumber() << std::endl;
  std::cout << " Number of protons: " << GetProtonNumber() << std::endl;
  std::cout << " Basis parameter zeta: " << GetZeta() << std::endl;
  if(print_obs) one.PrintSpace();
  if(print_tbs) two.PrintSpace();
}

void ModelSpace::PrintHoleOrbits()
{
  std::cout << " Number of electrons: " << GetElectronNumber() << std::endl;
  std::cout << " List of hole orbits: " << std::endl;
 for (auto& it: hole_occ)
 {
   Orbit& o = orbits.GetOrbit(it.first);
   int width = 4;
   std::cout << " n =" << std::setw(width) << o.n << ", l =" << std::setw(width) << o.l << ", j2 =" << std::setw(width) << o.j2 << 
     ", prob = " << std::setw(width) << it.second << ", occ = " << std::setw(width) << it.second*(o.j2+1) << std::endl;
 }
}

//std::map<int,double> ModelSpace::GetElectronOccupation(arma::vec SPEs)
//{
//  std::map<int,double> tmp_holes;
//  OneBodySpace& obs = GetOneBodySpace();
//  for (int ich=0; ich< obs.GetNumberChannels(); ich++) {
//    std::vector<unsigned long long> tmp(obs.channels[ich].begin(), obs.channels[ich].end());
//    arma::uvec sub_idx(std::vector<unsigned long long>(tmp.begin(), tmp.end()));
//    arma::vec SPE_ch = SPEs(sub_idx);
//    arma::uvec sorted_idx = sub_idx(arma::sort_index(SPE_ch));
//    std::vector<double> tmp_hole;
//    for (auto & it : hole_occ){
//      Orbit& o_h = GetOrbit(it.first);
//      if (o_h.kappa != obs.kappas[ich]) continue;
//      tmp_hole.push_back(it.second);
//    }
//    int cnt = 0;
//    for (int i : sorted_idx){
//      cnt += 1;
//      if(cnt <= obs.n_small_channel[ich]) continue;
//      if(cnt-obs.n_small_channel[ich] > tmp_hole.size()) break;
//      tmp_holes[i] = tmp_hole[cnt-obs.n_small_channel[ich]-1];
//    }
//  }
//  return tmp_holes;
//}
//
void ModelSpace::UpdateOccupation(std::map<int,double> tmp_holes)
{
  hole_occ = tmp_holes;
  for (auto o : orbits.orbits) o.occ = 0;
  for (auto it : tmp_holes){
    Orbit & o = GetOrbit(it.first);
    o.occ = it.second;
  }
}

void ModelSpace::UpdateOrbitals(arma::vec SPEs)
{
  std::map<int,double> tmp_holes;
  OneBodySpace& obs = GetOneBodySpace();

  small_components.clear();
  large_components.clear();
  all_orbits.clear();
  holes.clear();
  particles.clear();
  for (int i=0; i<GetNumberOrbits(); i++) all_orbits.insert(i);
  for (int ich=0; ich< obs.GetNumberChannels(); ich++) {
    std::vector<double> tmp_holes_ch;
    for (auto & it : hole_occ){
      Orbit& o_h = GetOrbit(it.first);
      if (o_h.kappa != obs.kappas[ich]) continue;
      tmp_holes_ch.push_back(it.second);
    }

    std::vector<unsigned long long> tmp(obs.channels[ich].begin(), obs.channels[ich].end());
    arma::uvec sub_idx(std::vector<unsigned long long>(tmp.begin(), tmp.end()));
    arma::vec SPE_ch = SPEs(sub_idx);
    arma::uvec sorted_idx = sub_idx(arma::sort_index(SPE_ch));
    int cnt = 0;
    for (int i : sorted_idx){
      cnt += 1;
      if(cnt <= obs.n_small_channel[ich]){
        small_components.insert(i);
        all_orbits.insert(i);
      }
      else {
        large_components.insert(i);
        all_orbits.insert(i);
      }
      if(cnt <= obs.n_small_channel[ich]) continue;
      if(cnt-obs.n_small_channel[ich] > tmp_holes_ch.size()) break;
      tmp_holes[i] = tmp_holes_ch[cnt-obs.n_small_channel[ich]-1];
      holes.insert(i);
    }
  }
  hole_occ = tmp_holes;
  UpdateOccupation(tmp_holes);
  for (int i : all_orbits){
    auto it = holes.find(i);
    if(it != holes.end()) continue;
    particles.insert(i);
  }
  //std::cout << "all orbits: " ;
  //for (int i : all_orbits) std::cout << i << " ";
  //std::cout << std::endl;
  //std::cout << "holes: " ;
  //for (int i : holes) std::cout << i << " ";
  //std::cout << std::endl;
  //std::cout << "paraticles: " ;
  //for (int i : particles) std::cout << i << " ";
  //std::cout << std::endl;
  // TODO: core valence, qspace
}
