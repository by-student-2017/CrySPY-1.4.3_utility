#!/usr/bin/env python
# -*- coding: utf-8 -*-

# ---------- import
from ase.constraints import FixSymmetry
from ase.filters import FrechetCellFilter
from ase.io import read, write
from ase.optimize import BFGS
from orb_models.forcefield import pretrained
from orb_models.forcefield.calculator import ORBCalculator

# ---------- input structure
atoms = read('POSCAR')

# ---------- set up ORB v3 Calculator
# Conservativeモード推奨（力・応力がエネルギー勾配に基づく）
model = pretrained.orb_v3_conservative_inf_omat(device="cuda", precision="float32-high")
atoms.calc = ORBCalculator(model=model)

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

