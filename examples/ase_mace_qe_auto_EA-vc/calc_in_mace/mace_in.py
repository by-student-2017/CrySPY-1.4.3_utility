#!/usr/bin/env python
# -*- coding: utf-8 -*-

# ---------- import
from ase.constraints import FixSymmetry
from ase.filters import FrechetCellFilter
from ase.io import read, write
from ase.optimize import BFGS
import mace_models

# ---------- input structure
atoms = read('POSCAR')

# ---------- set up MACE Calculator
# モデルはMatbench Discoveryで公開されているものを指定
# 例: "MACE-MP-0_small" または "MACE-MP-0"
model = mace_models.load("MACE-MP-0_small")
atoms.calc = model.get_calculator(dtype="float64")

# ---------- apply symmetry constraint and cell filter
atoms.set_constraint([FixSymmetry(atoms)])
cell_filter = FrechetCellFilter(atoms)

# ---------- optimizer
opt = BFGS(cell_filter, trajectory='opt.traj')

# ---------- run optimization
converged = opt.run(fmax=0.01, steps=2000)

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
opt_atoms.set_constraint(None)
write('CONTCAR', opt_atoms, format='vasp', direct=True)

# ------ check_opt
with open('out_check_opt', mode='w') as f:
    f.write('done\n' if converged else 'not yet\n')
