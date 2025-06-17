# ---------- import
from ase.constraints import FixSymmetry
from ase.filters import FrechetCellFilter
from ase.io import read, write
from ase.optimize import FIRE, BFGS, LBFGS
from chgnet.model import CHGNetCalculator

# ---------- input structure
# CrySPY outputs 'POSCAR' as an input file in work/xxxxxx directory
atoms = read('POSCAR')

# ---------- set up
atoms.calc = CHGNetCalculator()
atoms.set_constraint([FixSymmetry(atoms)])
cell_filter = FrechetCellFilter(atoms)
opt = BFGS(cell_filter, trajectory='opt.traj')

# ---------- run
opt.run(fmax=0.01, steps=2000)

# ---------- write structure
write('opt_struc.vasp', cell_filter.atoms, format='vasp', direct=True)

# ---------- opt. structure and energy
# [rule in ASE interface]
# output file for energy: 'log.tote' in eV/cell
#                         CrySPY reads the last line of 'log.tote'
# output file for structure: 'CONTCAR' in vasp format
# ------ energy
e = cell_filter.atoms.get_total_energy()    # eV/cell
with open('log.tote', mode='w') as f:
    f.write(str(e))

# -- for end_point
e_atom = e/len(atoms)
with open('end_point', mode='w') as f:
    f.write(str(e_atom))

# ------ struc
opt_atoms = cell_filter.atoms.copy()
opt_atoms.set_constraint(None)    # remove constraint for pymatgen
write('CONTCAR', opt_atoms, format='vasp', direct=True)


