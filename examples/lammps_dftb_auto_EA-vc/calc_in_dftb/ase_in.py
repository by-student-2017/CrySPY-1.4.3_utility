from ase.constraints import FixSymmetry
from ase.filters import FrechetCellFilter
#from ase.optimize import BFGS
from ase.io import read, write
import subprocess
import sys
import os
import re
import subprocess

# ---------- input structure
# CrySPY outputs 'POSCAR' as an input file in work/xxxxxx directory
atoms = read('POSCAR', format='vasp')
cell = atoms.cell.lengths()  # [a, b, c]

# ---------- Extract kppvol from dftb_in.hsd (comment line #kppvol = value)
kppvol = None
if os.path.exists("dftb_in.hsd"):
    with open("dftb_in.hsd", "r") as f:
        for line in f:
            if line.strip().startswith("#kppvol"):
                match = re.search(r'#kppvol\\s*=\\s*(\\d+(\\.\\d+)?)', line)
                if match:
                    kppvol = float(match.group(1))
                break

# Default if not found
if kppvol is None:
    kppvol = 40.0

# ---------- calculate k-points based on cell lengths and kppvol
volume = cell[0] * cell[1] * cell[2]
total_kpoints = kppvol * volume
scale = total_kpoints ** (1/3)  # cubic root for 3D distribution
nkx = max(1, round(scale * cell[0] / sum(cell)))
nky = max(1, round(scale * cell[1] / sum(cell)))
nkz = max(1, round(scale * cell[2] / sum(cell)))

# ---------- setting and run
#atoms.calc = XTB(method='GFN1-xTB') # GFN1-xTB, GFN2-xTB, GFN-FF
#atoms.set_constraint([FixSymmetry(atoms)])
#cell_filter = FrechetCellFilter(atoms, hydrostatic_strain=False)
#opt = BFGS(cell_filter)
if not os.path.exists("dftb_in.hsd"):
    subprocess.run(["cp", "./../../dftb_in.hsd", "./"], check=True)

# Read dftb_in.hsd content
with open("dftb_in.hsd", "r") as f:
    content = f.read()

# Replace nkx, nky, nkz lines
content = re.sub(r'nkx', f'{nkx}', content)
content = re.sub(r'nky', f'{nky}', content)
content = re.sub(r'nkz', f'{nkz}', content)

# Write updated content back to dftb_in.hsd
with open("dftb_in.hsd", "w") as f:
    f.write(content)

# ---------- run
#converged = opt.run(fmax=0.01, steps=2000)
#cpu_count = os.cpu_count() or 1
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

# ------ check SCC
scc_failed = False
with open("dftb_out.log") as f:
    for line in f:
        if "SCC is NOT converged" in line:
            scc_failed = True
            break

energy = None
if not scc_failed:
    with open("detailed.out") as f:
        for line in f:
            if line.strip().startswith("Total energy:"):
                try:
                    energy = float(line.split()[-2])
                except ValueError:
                    energy = None
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
