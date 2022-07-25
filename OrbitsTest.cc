#include <iostream>
#include <iomanip>
#include <armadillo>
#include "PhysicalConstants.hh"
#include "Parameters.hh"
#include "Orbits.hh"
#include "ModelSpace.hh"
#include "Operator.hh"

int main(int argc, char** argv)
{
  std::cout << " Orbit test " << std::endl;
  Parameters parameters(argc,argv);

  std::string orbitals = parameters.s("orbitals");
  std::string atom = parameters.s("atom");
  std::string radial_function_type = parameters.s("radial_function_type");
  double zeta_inv = parameters.d("zeta_inv");

  Orbits orbits = Orbits(orbitals, radial_function_type); 
  orbits.Print();
  ModelSpace ms = ModelSpace(atom, 1/zeta_inv, orbits);

  Operator H = Operator(ms);
  if(orbits.relativistic) H.SetDiracCoulombHamiltonian(true, false);
  else H.SetCoulombHamiltonian(true, false);
  H.OrthoNormalize();

  OneBodySpace& obs = ms.GetOneBodySpace();
  for (int ich=0; ich< obs.GetNumberChannels(); ich++) {
    std::vector<unsigned long long> tmp(obs.channels[ich].begin(), obs.channels[ich].end());
    arma::uvec sub_idx(std::vector<unsigned long long>(tmp.begin(), tmp.end()));
    arma::mat Hch = H.OneBody(sub_idx, sub_idx);
    arma::vec eig;
    arma::mat vec;
    arma::eig_sym(eig, vec, Hch);
    std::cout << " kappa: " << obs.GetKappa(ich) << std::endl;
    for (auto it : eig){
      std::cout << std::setw(16) << std::setprecision(8) << std::scientific << it << std::endl;
    }
  }
}

