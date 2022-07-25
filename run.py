#!/usr/bin/env python3
import sys, os, itertools, subprocess

atoms = ["He"]
#zetai_list = range(2,10)
zetai_list = range(2,3)
#orbitals_list = [f"nonrel-s{n}-p{n}" for n in range(8,9)]
orbitals_list = [f"rel-s{n}-p{n}" for n in range(6,7)]
# orbitals_list = [f"nonrel-s{n}-p{n}-d{n}" for n in range(11,12)]
#orbitals_list = [f"rel-s{n}-p{n}" for n in range(2,3)]
#orbitals_list = [f"rel-s{n}" for n in range(8,9)]
ARGS = {}
# NMeshs = list(range(200,1000,100) ) + list(range(1000,11000,1000))
NMeshs = list(list(range(8000,9000,1000)))
# NMeshs = list(range(200, 400, 200))
for atom, orbitals, zeta_inv, NMesh in itertools.product(atoms, orbitals_list, zetai_list, NMeshs):
    ARGS["atom"] = atom
    ARGS["orbitals"] = orbitals
    ARGS["zeta_inv"] = zeta_inv
    ARGS["NMesh"] = NMesh
    ARGS["radial_function_type"] = "LSpinor"
    #ARGS["eeintegral_mesh_type"] = "Legendre"
    ARGS["eeintegral_mesh_type"] = "Laguerre"
    ARGS["eeintegral_rmax"] = 10
    if(orbitals.find("nonrel") != -1): ARGS["radial_function_type"] = "NonRel_Laguerre"
    #ARGS["filename_coulomb"] = f"eeCoulomb_{orbitals}_{ARGS['radial_function_type']}_MeshType{ARGS['eeintegral_mesh_type']}_NMesh{NMesh}.bin"
    ARGS["filename_summary"] = f"{atom}_{orbitals}_z{zeta_inv}_{ARGS['radial_function_type']}_NMesh{NMesh}.txt"

    cmd = ' '.join(["time","./HFTest.exe",] + [f"{x}={ARGS[x]}" for x in ARGS.keys()])
    subprocess.call(cmd, shell=True)    
