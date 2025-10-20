from ase.constraints import FixSymmetry
from ase.filters import FrechetCellFilter
#from ase.optimize import BFGS
from ase.io import read, write
import subprocess
import sys
import os

# ---------- input structure
# CrySPY outputs 'POSCAR' as an input file in work/xxxxxx directory
atoms = read('POSCAR', format='vasp')

# ---------- setting and run
#atoms.calc = XTB(method='GFN1-xTB') # GFN1-xTB, GFN2-xTB, GFN-FF
#atoms.set_constraint([FixSymmetry(atoms)])
#cell_filter = FrechetCellFilter(atoms, hydrostatic_strain=False)
#opt = BFGS(cell_filter)
if not os.path.exists("dftb_in.hsd"):
    subprocess.run(["cp", "./../../dftb_in.hsd", "./"], check=True)

# ---------- run
#converged = opt.run(fmax=0.01, steps=2000)
cpu_count = 1
subprocess.run(f"mpirun -np {cpu_count} dftb+ | tee dftb_out.log", shell=True)

# ---------- rule in ASE interface
# output file for energy: 'log.tote' in eV/cell
#                         CrySPY reads the last line of 'log.tote' file
# outimized structure: 'CONTCAR' file in vasp format
# check_opt: 'out_check_opt' file ('done' or 'not yet')
#                         CrySPY reads the last line of 'out_check_opt' file

# ------ energy
#e = cell_filter.atoms.get_total_energy()    # eV/cell
#with open('log.tote', mode='w') as f:
#    f.write(str(e))
energy = None
with open("detailed.out") as f:
    for line in f:
        if line.strip().startswith("Total energy:"):
            energy = float(line.split()[-2])
            break
with open('log.tote', 'w') as f:
    f.write(f"{energy:.6f}\n")

# ------ struc
#opt_atoms = cell_filter.atoms.copy()
#opt_atoms.set_constraint(None)    # remove constraint for pymatgen
#write('CONTCAR', opt_atoms, format='vasp', direct=True)
opt_atoms = read('geom.out.gen', format='gen')
write('CONTCAR', opt_atoms, format='vasp', direct=True)

# ------ check_opt
converged = os.path.exists('CONTCAR')
with open('out_check_opt', mode='w') as f:
    if converged:
        f.write('done\n')
    else:
        f.write('not yet\n')
