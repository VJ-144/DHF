
# Dirac-Hartree-Fock Implementation for C++
Includes a complete implementation and development of a working Dirac Hartree Fock framework for atomic systems, in which we hope we could apply our first-principles many-body method, the VS-IMSRG. The Dirac-Hartree-Fock is a gold standard computational approach for atomic systems and one of only a handful of such codes in the world, which opens an entirely new arena for ab initio calculations of potentially all atoms


While the quantity of results thus far obtained are limited to a few noble gases, the quality reflects the tremendous progress: we now are finding converged ground state energies for argon, in excellent agreement with experiment, which is the first time we have achieved this with the IMSRG. This in itself is a landmark accomplishment, and it only remains to test for open-shell atoms, higher electron systems, and implement isotope shift operators. If the results are of the expected quality, this will revolutionize the field of atomic theory as well nuclear physics where atomic input is needed.

## Necessary Libraries:
GSL
Armadillo
Boost

## Run Instructions
Compile program with the command:
"""
make clean; make
"""


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




