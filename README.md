
# Dirac-Hartree-Fock Implementation for C++

## Necessary Libraries:
GSL
Armadillo
Boost

## Instructions

Run makefile with "make" command to compile <br />
Specify parameters in run.py <br />
Execute run.py to run code


## parameters:

NMeshs - controls accuracy on 2-body integral/electron-electron interactions

orbitals_list - number and type of electron orbitals to be included in self-consistent field procedure of the DHF

atoms - Relevant atomic system

zetai_list - arbitary scaling parameter used to tweak convergence

radial_function_type - radial functions describing spatial dependence of Dirac Spinor wavefunction

eeintegral_mesh_type - Mesh type for 2-body integral

eeintegral_rmax - upper limit of 2-body integral 




