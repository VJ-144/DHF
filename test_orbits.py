#!/usr/bin/env python3
import sys, os, itertools, subprocess
from sympy.physics.hydrogen import E_nl_dirac, E_nl

atoms = ["H"]
zetai_list = range(1,2)
zetai_list = [1]
#orbitals_list = [f"rel-s{n}-p{n}-d{n}" for n in range(4,5)]
orbitals_list = [f"rel-s{n}" for n in range(20,21)]
#orbitals_list = [f"rel-p{n}" for n in range(10,11)]
#orbitals_list = [f"rel-d{n}" for n in range(2,3)]
#orbitals_list = [f"nonrel-d{n}" for n in range(10,11)]
ARGS = {}

print("Non-rel hydrogen")
print(f"n=1, {E_nl(1)}")
print(f"n=2, {E_nl(2)}")
print(f"n=3, {E_nl(3)}")

print("Rel hydrogen")
print(f"n=0, l=0, j=1/2, {E_nl_dirac(1,0)}")
print(f"n=1, l=0, j=1/2, {E_nl_dirac(2,0)}")
print(f"n=0, l=1, j=1/2, {E_nl_dirac(2,1, spin_up=False)}")
print(f"n=0, l=1, j=3/2, {E_nl_dirac(2,1, spin_up=True)}")
print(f"n=0, l=2, j=3/2, {E_nl_dirac(3,2, spin_up=False)}")
print(f"n=0, l=2, j=5/2, {E_nl_dirac(3,2, spin_up=True)}")

for atom, orbitals, zeta_inv in itertools.product(atoms, orbitals_list, zetai_list):
    ARGS["atom"] = atom
    ARGS["orbitals"] = orbitals
    ARGS["zeta_inv"] = zeta_inv
    if(orbitals.find("nonrel") != -1): ARGS["radial_function_type"] = "NonRel_Laguerre"
    cmd = ' '.join(["time","./OrbitsTest.exe",] + [f"{x}={ARGS[x]}" for x in ARGS.keys()])
    subprocess.call(cmd, shell=True)    
