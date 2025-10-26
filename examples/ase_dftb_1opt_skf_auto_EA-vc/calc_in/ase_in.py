from ase.constraints import FixSymmetry
from ase.filters import FrechetCellFilter
#from ase.optimize import BFGS
from ase.io import read, write
import subprocess
import sys
import os
import re
from collections import Counter
from ase.data import atomic_numbers

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

# --- MaxAngularMomentum block ---
def get_max_angular(element):
    Z = atomic_numbers[element]
    if Z <= 2:
        return "s"
    elif Z <= 18:
        return "p"
    elif Z <= 54:
        return "d"
    else:
        return "d"

elements = sorted(set(atoms.get_chemical_symbols()))
max_ang_block = "MaxAngularMomentum = {\n"
for el in elements:
    max_ang_block += f"    {el} = \"{get_max_angular(el)}\"\n"
max_ang_block += "  }"

# Add or replace in the content.
if re.search(r'MaxAngularMomentum\s*=\s*\{[^}]*\}', content, flags=re.DOTALL):
    content = re.sub(r'MaxAngularMomentum\s*=\s*\{[^}]*\}', max_ang_block, content, flags=re.DOTALL)
else:
    content += "\n" + max_ang_block

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

# ------ check ERROR!
error_detected = False
with open("dftb_out.log") as f:
    for line in f:
        if "ERROR!" in line:
            error_detected = True
            break

energy = None
if not error_detected:
    with open("detailed.out") as f:
        for line in f:
            if line.strip().startswith("Total energy:"):
                try:
                    energy = float(line.split()[-2])
                except ValueError:
                    energy = None
                break

# ---------- compute formation energy
symbols = atoms.get_chemical_symbols()
composition = Counter(symbols)
num_atoms = len(symbols)

# Get unique element order from POSCAR
unique_elements = []
for s in symbols:
    if s not in unique_elements:
        unique_elements.append(s)

# Read end_point from cryspy.in
end_point_values = []
if os.path.exists("./../../cryspy.in"):
    with open("./../../cryspy.in", "r") as f:
        for line in f:
            if line.strip().startswith("end_point"):
                nums = re.findall(r'[-+]?\d*\.\d+|\d+', line)
                if nums:
                    end_point_values = list(map(float, nums[1:]))  # skip 'end_point'
                break

# Compute reference energy
E_ref = 0.0
if end_point_values and len(end_point_values) == len(unique_elements):
    for i, elem in enumerate(unique_elements):
        E_ref += composition[elem] * end_point_values[i]

# Formation energy check
if energy is not None and E_ref != 0:
    E_form = (energy - E_ref) / num_atoms
    if E_form < -10.0*2 or E_form > 10.0*5:  # abnormal threshold
        energy = None

# ---------- write energy
with open('log.tote', 'w') as f:
    if energy is not None:
        f.write(f"{energy:.6f}\n")

# ------ struc
#opt_atoms = cell_filter.atoms.copy()
#opt_atoms.set_constraint(None)    # remove constraint for pymatgen
#write('CONTCAR', opt_atoms, format='vasp', direct=True)
opt_atoms = read('geom.out.gen', format='gen')
write('CONTCAR', opt_atoms, format='vasp', direct=True)

# ------ check_opt
#converged = os.path.exists('CONTCAR')
converged = os.path.exists('CONTCAR') and (energy is not None)
with open('out_check_opt', mode='w') as f:
    if converged:
        f.write('done\n')
    else:
        f.write('not yet\n')
