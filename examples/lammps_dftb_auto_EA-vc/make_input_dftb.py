#!/usr/bin/env python3
import sys
import os
import shutil
import math
from ase.build import bulk
from ase.io import write
from ase import Atoms

# ====== 1. 引数チェック ======
elements = sys.argv[1:]
if len(elements) < 1:
    print("Usage: python make_input_dftb.py Element1 [Element2 ...]")
    sys.exit(1)

# ====== 2. element_data読み込み ======
element_data = {}
if os.path.isfile('element_data_dftb.txt'):
    with open('element_data_dftb.txt') as f:
        for line in f:
            parts = line.split()
            if len(parts) >= 6 and not line.startswith('Element'):
                element_data[parts[0]] = {
                    'structure': parts[1],
                    'energy': float(parts[2]),
                    'a': float(parts[3]),
                    'b': float(parts[4]),
                    'c': float(parts[5])
                }

# ====== 3. cryspy.in生成 ======
max_len = max(len(el) for el in elements)
atype_line = "atype  = " + " ".join(el.ljust(max_len) for el in elements)
ll_nat_line = "ll_nat = " + " ".join("0".ljust(max_len) for _ in elements)
ul_nat_line = "ul_nat = " + " ".join("8".ljust(max_len) for _ in elements)
end_point_line = "end_point = " + " ".join("0.0".ljust(max_len) for _ in elements)

cryspy_lines = [
    "[basic]",
    "algo = EA-vc",
    "calc_code = ASE",
    "nstage = 1",
    "njob = 1",
    "jobcmd = bash",
    "jobfile = job_cryspy\n",
    "[structure]",
    atype_line,
    ll_nat_line,
    ul_nat_line + "\n",
    "[ASE]",
    "ase_python = ase_in.py\n",
    "[EA]",
    "n_pop = 20",
    "n_crsov = 5",
    "n_perm = 2",
    "n_strain = 2",
    "n_rand = 2",
    "n_add = 3",
    "n_elim = 3",
    "n_subs = 3",
    "target = random",
    "n_elite = 2",
    "n_fittest = 10",
    "slct_func = TNM",
    "t_size = 2",
    "maxgen_ea = 0",
    end_point_line + "\n",
    "[option]",
    "load_struc_flag = True",
    "symprec = 0.1",
    "out_cif = False\n"
]

with open('cryspy.in', 'w') as f:
    f.write("\n".join(cryspy_lines))
print("Created cryspy.in")

# ====== 4. calc_in置換 ======
src = 'calc_in_dftb'
dst = 'calc_in'
if os.path.exists(dst):
    shutil.rmtree(dst)
    print(f"Removed existing directory: {dst}")
if os.path.isdir(src):
    shutil.copytree(src, dst)
    print(f"Copied {src} to {dst}")
else:
    print(f"Error: Source directory {src} not found")

# ====== 5. 原子座標とセル生成関数 ======
def get_positions(struct, a, c):
    if struct == "bcc":
        return [(0,0,0), (0.5*a,0.5*a,0.5*a)]
    elif struct == "hcp":
        return [(0,0,0), (2/3*a,1/3*a,c/2)]
    elif struct == "dhcp":
        return [(0,0,0), (a/2,a*(3)**0.5/6,c/4), (0,0,c/2), (a/2,a*(3)**0.5/6,3*c/4)]
    elif struct == "dia":
        return [(0,0,0), (0.25*a,0.25*a,0.25*a)]
    else:
        return [(0,0,0)]

def get_cell_parameters(struct, a, b, c):
    if struct in ["fcc", "bcc", "sc", "dia"]:
        return [(a,0,0),(0,b,0),(0,0,c)]
    elif struct in ["hcp","dhcp"]:
        return [(a,0,0),(-a/2,a*math.sqrt(3)/2,0),(0,0,c)]
    elif struct == "ort":
        return [(a,0,0),(0,b,0),(0,0,c)]
    elif struct == "mon":
        beta = math.radians(100)
        return [(a,0,0),(0,b,0),(c*math.cos(beta),0,c*math.sin(beta))]
    else:
        return [(a,0,0),(0,b,0),(0,0,c)]

# ====== 6. 構造ディレクトリ生成 ======
xx_tmp_dir = 'Xx_tmp'
if not os.path.isdir(xx_tmp_dir):
    os.makedirs(xx_tmp_dir)

for el in elements:
    struct = element_data.get(el, {}).get('structure', 'sc')
    a_val = element_data.get(el, {}).get('a', 4.0)
    b_val = element_data.get(el, {}).get('b', a_val)
    c_val = element_data.get(el, {}).get('c', a_val)

    subdir_name = f"{el}_{struct}"
    os.makedirs(subdir_name, exist_ok=True)

    # ASEでgeometry.genとPOSCAR出力
    try:
        ase_struct_map = {'fcc':'fcc','bcc':'bcc','hcp':'hcp','dia':'diamond','sc':'sc'}
        if struct in ase_struct_map:
            atoms = bulk(el, ase_struct_map[struct], a=a_val)
        else:
            atoms = Atoms(el, positions=get_positions(struct,a_val,c_val),
                          cell=get_cell_parameters(struct,a_val,b_val,c_val), pbc=True)
    except Exception:
        atoms = bulk(el,'sc',a=4.0)

    write(os.path.join(subdir_name,'geometry.gen'), atoms, format='gen')
    write(os.path.join(subdir_name,'POSCAR'), atoms, format='vasp')

    print(f"Created {subdir_name} with geometry.gen and POSCAR.")

