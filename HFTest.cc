#include <armadillo>
#include "Parameters.hh"
#include "PhysicalConstants.hh"
#include "Orbits.hh"
#include "ModelSpace.hh"
#include "TwoBodyOperator.hh"
#include "Operator.hh"
#include "HartreeFock.hh"
#include "HFMBPT.hh"
int main(int argc, char** argv)
{
  std::cout << " HartreeFock test " << std::endl;
  std::cout << " Number of OpenMP threads: " <<  omp_get_max_threads() << std::endl;
  Parameters parameters(argc,argv);

  std::string orbitals = parameters.s("orbitals");
  std::string atom = parameters.s("atom");
  std::string radial_function_type = parameters.s("radial_function_type");
  std::string integral_mesh_type = parameters.s("eeintegral_mesh_type");
  std::string filename_coulomb = parameters.s("filename_coulomb");
  std::string filename_summary = parameters.s("filename_summary");
  int NMesh = parameters.i("NMesh");
  double integral_rmax = parameters.d("eeintegral_rmax");
  double zeta_inv = parameters.d("zeta_inv");
  Orbits orbits = Orbits(orbitals, radial_function_type); 
  orbits.Print();
  ModelSpace ms = ModelSpace(atom, 1/zeta_inv, orbits);

  Operator H = Operator(ms);
  H.TwoBody.SetNMesh(NMesh);
  H.TwoBody.SetRmax(integral_rmax);
  H.TwoBody.SetMeshType(integral_mesh_type);
  H.TwoBody.SetFileCoulomb(filename_coulomb);

  if(orbits.relativistic) H.SetDiracCoulombHamiltonian(true, true);
  else H.SetCoulombHamiltonian(true, true);
  H.OrthoNormalize();
  //H.Print();

  HartreeFock HF = HartreeFock(H);
  HF.Solve();
  std::cout << " HF energy: " << std::setw(16) << std::setprecision(8) << HF.EHF << std::endl;

  HFMBPT MBPT = HFMBPT(H, HF.C);
  Operator HNO = MBPT.TransformBasis(H);
  HNO.modelspace->UpdateOrbitals(HF.SPEs);
  HNO = HNO.DoNormalOrdering();
  std::cout << " EHF: " << std::setw(14) << std::setprecision(8) << HNO.ZeroBody << std::endl;
  std::cout << " MP2: " << std::setw(14) << std::setprecision(8) << HNO.ZeroBody + MBPT.GetMP2_Energy(HNO) << std::endl;

  if(filename_summary!=""){
    std::ofstream output;
    output.open(filename_summary,std::ofstream::out);
    output << " EHF: " << std::setw(14) << std::setprecision(8) << HNO.ZeroBody << std::endl;
    output << " MP2: " << std::setw(14) << std::setprecision(8) << HNO.ZeroBody + MBPT.GetMP2_Energy(HNO) << std::endl;
    output.close();
  }
}
