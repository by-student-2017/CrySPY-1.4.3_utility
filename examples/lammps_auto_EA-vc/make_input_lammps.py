import sys
import os
import shutil
from ase.build import bulk
from ase.io import write
from ase import Atoms

elements = sys.argv[1:]
if len(elements) < 1:
    print("Usage: python make_input_lammps.py Element1 [Element2 ...]")
    sys.exit(1)

# ====== 1. element_data_lammps.txt読み込み ======
element_data = {}
with open('element_data_lammps.txt') as f:
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

# ====== 2. cryspy.in生成 ======
max_len = max(len(el) for el in elements)
atype_line = "atype  = " + " ".join(el.ljust(max_len) for el in elements)
ll_nat_line = "ll_nat = " + " ".join("0".ljust(max_len) for _ in elements)
ul_nat_line = "ul_nat = " + " ".join("12".ljust(max_len) for _ in elements)

end_points = [f"{element_data.get(el, {}).get('energy', -4.5000):.4f}" for el in elements]
end_point_line = "end_point  = " + " ".join(ep.ljust(max_len) for ep in end_points)

cryspy_lines = [
    "[basic]",
    "algo = EA-vc",
    "calc_code = LAMMPS",
    "nstage = 1",
    "njob = 5",
    "jobcmd = bash",
    "jobfile = job_cryspy\n",
    "[structure]",
    atype_line,
    ll_nat_line,
    ul_nat_line,
    "mindist_factor = 1.0\n",
    "[LAMMPS]",
    "lammps_infile = in.lmp",
    "lammps_outfile = out.lmp",
    "lammps_data = data.lmp",
    #"lammps_potential = CuAgAuNiPdPtAlPbFeMoTaWMgCoTiZr_Zhou04.eam.alloy\n",
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
    "[option]\n"
]

with open('cryspy.in', 'w') as f:
    f.write("\n".join(cryspy_lines))
print("Created cryspy.in")

# ====== 3. calc_inをコピー ======
calc_in_dir = 'calc_in'
calc_in_dst = 'calc_in'
if not os.path.isdir(calc_in_dir):
    print("Error: calc_in directory not found")
    sys.exit(1)
if os.path.exists(calc_in_dst):
    shutil.rmtree(calc_in_dst)
shutil.copytree(calc_in_dir, calc_in_dst)
print(f"Copied calc_in template to {calc_in_dst}")

# ====== 4. calc_in内のin_tmp.lmp編集 ======
in_tmp_file = os.path.join(calc_in_dst, 'in_tmp.lmp')
in_lmp_file = os.path.join(calc_in_dst, 'in.lmp')
if os.path.isfile(in_tmp_file):
    with open(in_tmp_file) as f:
        lines = f.readlines()
    elem_string = " ".join(elements)
    new_lines = [f'variable elem string "{elem_string}"\n' if line.strip().startswith("variable elem string") else line for line in lines]
    with open(in_lmp_file, 'w') as f:
        f.writelines(new_lines)
    os.remove(in_tmp_file)
    print(f"Created {in_lmp_file} with elements: {elem_string}")
else:
    print(f"Warning: {in_tmp_file} not found in calc_in")

# ====== 5. Xx_tmp確認 ======
xx_tmp_dir = 'Xx_tmp'
template_file = os.path.join(xx_tmp_dir, 'in_tmp.lmp')
if not os.path.isfile(template_file):
    print("Error: in_tmp.lmp not found in Xx_tmp")
    sys.exit(1)

# ====== 6. 各構造ディレクトリ処理 ======
for el in elements:
    struct = element_data.get(el, {}).get('structure', 'sc')
    a_val = element_data.get(el, {}).get('a', 4.0)
    b_val = element_data.get(el, {}).get('b', a_val)
    c_val = element_data.get(el, {}).get('c', a_val)

    subdir_name = f"{el}_{struct}"
    os.makedirs(subdir_name, exist_ok=True)

    # Copy Xx_tmp contents except in_tmp.lmp
    for item in os.listdir(xx_tmp_dir):
        if item == 'in_tmp.lmp':
            continue
        s = os.path.join(xx_tmp_dir, item)
        d = os.path.join(subdir_name, item)
        if os.path.isdir(s):
            shutil.copytree(s, d, dirs_exist_ok=True)
        else:
            shutil.copy2(s, d)

    # Edit in_tmp.lmp -> in.lmp
    with open(template_file) as f:
        lines = f.readlines()
    new_lines = [f'variable elem string "{el}"\n' if line.strip().startswith("variable elem string") else line for line in lines]
    with open(os.path.join(subdir_name, 'in.lmp'), 'w') as f:
        f.writelines(new_lines)

    # Generate data.lmp and POSCAR
    try:
        ase_struct_map = {
            'fcc': 'fcc', 'bcc': 'bcc', 'hcp': 'hcp',
            'dia': 'diamond', 'sc': 'sc',
            'rock': 'rocksalt', 'zb': 'zincblende', 'wz': 'wurtzite'
        }
        if struct in ase_struct_map:
            atoms = bulk(el, ase_struct_map[struct], a=a_val)
        elif struct == 'bct':
            cell = [[a_val, 0, 0], [0, a_val, 0], [0, 0, c_val]]
            atoms = Atoms(el, positions=[[0, 0, 0]], cell=cell, pbc=True)
        elif struct == 'dhcp':
            cell = [[a_val, 0, 0], [-a_val/2, (a_val*(3)**0.5)/2, 0], [0, 0, c_val]]
            positions = [
                [0, 0, 0], [a_val/2, (a_val*(3)**0.5)/6, c_val/4],
                [0, 0, c_val/2], [a_val/2, (a_val*(3)**0.5)/6, 3*c_val/4]
            ]
            atoms = Atoms(el*4, positions=positions, cell=cell, pbc=True)
        elif struct == 'ort':
            cell = [[a_val, 0, 0], [0, b_val, 0], [0, 0, c_val]]
            atoms = Atoms(el, positions=[[0, 0, 0]], cell=cell, pbc=True)
        elif struct == 'mon':
            cell = [[a_val, 0, 0], [0, b_val, 0], [0.1*a_val, 0, c_val]]
            atoms = Atoms(el, positions=[[0, 0, 0]], cell=cell, pbc=True)
        elif struct in ['αB12', 'βSn']:
            cell = [[a_val, 0, 0], [0, a_val, 0], [0, 0, c_val]]
            atoms = Atoms(el, positions=[[0, 0, 0]], cell=cell, pbc=True)
        else:
            atoms = bulk(el, 'sc', a=a_val)
    except Exception:
        atoms = bulk(el, 'sc', a=4.0)

    write(os.path.join(subdir_name, 'data.lmp'), atoms, format='lammps-data')
    write(os.path.join(subdir_name, 'POSCAR'), atoms, format='vasp')

    print(f"Created {subdir_name} with in.lmp, Xx_tmp contents, data.lmp, and POSCAR.")
