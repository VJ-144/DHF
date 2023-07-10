
# Dirac-Hartree-Fock Implementation for C++
Includes a complete implementation and development of a working Dirac Hartree Fock framework for atomic systems, in which we hope we could apply our first-principles many-body method, the VS-IMSRG. The Dirac-Hartree-Fock is a gold standard computational approach for atomic systems. The method utilises a self consistent field approach and the Dirac equation for the high accuracy results as intrinsic relativistic effects are intertwined into the calculation. This allows for post corrections of the calculation to be minimised (which is the current norm of these types of calculations) and is especially important for heavier atomic systems as they pertain stronger relativistic effects. Here we present one of only a handful of such codes in the world, which opens an entirely new arena for ab initio calculations of potentially all atoms


While the quantity of results thus far obtained are limited to a few noble gases, the quality reflects the tremendous progress: we now are finding converged ground state energies for argon, in excellent agreement with experiment, which is the first time we have achieved this with the IMSRG. This in itself is a landmark accomplishment, and it only remains to test for open-shell atoms, higher electron systems, and implement isotope shift operators. If the results are of the expected quality, this will revolutionize the field of atomic theory as well nuclear physics where atomic input is needed.

## Necessary Libraries:
GSL
Armadillo
Boost

## Run Instructions:
Compile program with the command:
```
make clean; make
```
A python wrapper is implemented to edit input parameters and adjust input parameters. After successful building only the execution of run.py is nessesary for calculations. 


## Parameters:

**NMeshs**                 - Controls the mesh grid for the integral evaluating the 2-body electron-electron (coulomb) interaction. This adjusts the                          fine width the integral is partitioned into. Reducing the NMesh to low will result in inaccurate results as the coulomb                          interaction ineffectivly calculated, this must be adusted relative to the atomic system of interest.

**orbitals_list**         - Number and type of electron orbitals to be included in self-consistent field calculation. This determines the model                              space in which the calculation are conducted. As such they must be considerable larger than the atomic configuration of                          the system of interest. A good general approach is to set the model space as x2+1 the atomic configuration. For example,                         Neon atomic config is 1s2 2s2 2p6 so a suitable model space would be 4s-4p or larger. This must be adjusted                                      approapriatly based on the atomic system of interest for an accurate calculation.

**atoms**                  - Atomic system of interest, a full list of atomic configurations can be found in the ModelSpace.cc

**zetai_list**             - Arbitary scaling parameter used to tweak convergence.

**radial_function_type**   - Radial functions describing spatial dependence of Dirac Spinor wavefunction. Refer to source code for alternate radial                           implementations.

**eeintegral_mesh_type**   - Mesh type for 2-body coulomb integral.

**eeintegral_rmax**        - upper limit of 2-body coulomb integral.




