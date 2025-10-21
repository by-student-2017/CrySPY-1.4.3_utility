#!/usr/bin/env python
# -*- coding: utf-8 -*-

# ---------- import
from ase.constraints import FixSymmetry
from ase.filters import FrechetCellFilter
from ase.io import read, write
from ase.optimize import BFGS
import mace_models
from ase.md.npt import NPT
from ase import units
import torch

# ---------- input structure
atoms = read('POSCAR')

# ---------- set up MACE Calculator
model = mace_models.load("MACE-MP-0_small")
atoms.calc = model.get_calculator(device="cuda", dtype="float32")  # GPU利用

# ---------- NPTダイナミクス（温度・圧力）
dyn = NPT(atoms,
          timestep=0.2 * units.fs,
          temperature_K=300,
          externalstress=0.624,  # eV/Å³ (100 GPa)
          ttime=100.0 * units.fs,
          pfactor=10.0)
dyn.run(steps=500)  # MDで条件反映

# ---------- apply symmetry constraint (必要ならここで適用)
# CrySPYが対称性保持を要求する場合のみ有効化
# atoms.set_constraint([FixSymmetry(atoms)])

# ---------- cell filter
cell_filter = FrechetCellFilter(atoms)

# ---------- optimizer
opt = BFGS(cell_filter, trajectory='opt.traj')

# ---------- run optimization
# 通常は十分なステップ数で実行
# converged = opt.run(fmax=0.01, steps=2000)
converged = opt.run(steps=1)  # CrySPY側で繰り返し制御する場合

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
