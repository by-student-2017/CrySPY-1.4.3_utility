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
converged = opt.run(fmax=0.01, steps=2000)

# ---------- rule in ASE interface
# output file for energy: 'log.tote' in eV/cell
#                         CrySPY reads the last line of 'log.tote' file
# outimized structure: 'CONTCAR' file in vasp format
# check_opt: 'out_check_opt' file ('done' or 'not yet')
#                         CrySPY reads the last line of 'out_check_opt' file

# ------ energy
e = cell_filter.atoms.get_total_energy()    # eV/cell
with open('log.tote', mode='w') as f:
    f.write(str(e))

# ------ struc
opt_atoms = cell_filter.atoms.copy()
opt_atoms.set_constraint(None)    # remove constraint for pymatgen
write('CONTCAR', opt_atoms, format='vasp', direct=True)

# ------ check_opt
with open('out_check_opt', mode='w') as f:
    if converged:
        f.write('done\n')
    else:
        f.write('not yet\n')

