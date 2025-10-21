#!/usr/bin/env python
# -*- coding: utf-8 -*-

# ---------- import
from ase.constraints import FixSymmetry
from ase.filters import FrechetCellFilter
from ase.io import read, write
from ase.optimize import BFGS
from sevenn.calculator import SevenNetCalculator

# ---------- input structure
# CrySPY outputs 'POSCAR' as an input file in work/xxxxxx directory
atoms = read('POSCAR')

# ---------- set up SevenNet Calculator
# model_key は Matbench Discoveryで公開されている SevenNet モデル名
# 例: "7net-mf-ompa"
atoms.calc = SevenNetCalculator(model_key="7net-mf-ompa", device="cuda")  # or "cpu"

# ---------- apply symmetry constraint and cell filter
atoms.set_constraint([FixSymmetry(atoms)])
cell_filter = FrechetCellFilter(atoms)

# ---------- optimizer
opt = BFGS(cell_filter, trajectory='opt.traj')

# ---------- run optimization
#converged = opt.run(fmax=0.01, steps=2000)
converged = opt.run(fmax=0.05, steps=500)

# ---------- CrySPY interface rules
# output file for energy: 'log.tote' in eV/cell
# optimized structure: 'CONTCAR' in VASP format
# check_opt: 'out_check_opt' file ('done' or 'not yet')

# ------ energy
e = cell_filter.atoms.get_total_energy()  # eV/cell
with open('log.tote', mode='w') as f:
    f.write(str(e))

# ------ structure
opt_atoms = cell_filter.atoms.copy()
opt_atoms.set_constraint(None)  # remove constraint for pymatgen
write('CONTCAR', opt_atoms, format='vasp', direct=True)

# ------ check_opt
with open('out_check_opt', mode='w') as f:
    f.write('done\n' if converged else 'not yet\n')
