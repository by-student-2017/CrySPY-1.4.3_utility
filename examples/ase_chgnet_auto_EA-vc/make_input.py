import sys
import os
import shutil
from ase.build import bulk
from ase.io import write

# File containing element, structure, and energy data
energy_table_file = 'element_data_chgnet.txt'

# Read the energy table into a dictionary
element_data = {}
with open(energy_table_file, 'r') as f:
    for line in f:
        line = line.strip()
        # Skip header or empty lines
        if not line or line.startswith('Element') or line.startswith('='):
            continue
        parts = line.split()
        if len(parts) >= 6:
            element = parts[0]
            structure = parts[1]
            try:
                energy = float(parts[2])
            except ValueError:
                energy = -4.5000  # fallback energy
            element_data[element] = {
                'structure': structure,
                'energy': energy,
                'a': float(parts[3]),
                'b': float(parts[4]),
                'c': float(parts[5])
            }

# Get elements from command line arguments
elements = sys.argv[1:]
if len(elements) < 2:
    print("Usage: python make.py Element1 Element2 [Element3 ...]")
    sys.exit(1)

# Prepare cryspy.in content
cryspy_lines = []
cryspy_lines.append("[basic]")
cryspy_lines.append("algo = EA-vc")
cryspy_lines.append("calc_code = ASE")
cryspy_lines.append("nstage = 1")
cryspy_lines.append("njob = 5")
cryspy_lines.append("jobcmd = bash")
cryspy_lines.append("jobfile = job_cryspy\n")

# Determine column width
col_width = max(len(el) for el in elements) + 0  # 0 extra spaces for padding
cryspy_lines.append("[structure]")
cryspy_lines.append("atype  = " + " ".join(f"{el:<{col_width}}" for el in elements))
cryspy_lines.append("ll_nat = " + " ".join(f"{0:<{col_width}}" for _ in elements))
cryspy_lines.append("ul_nat = " + " ".join(f"{12:<{col_width}}" for _ in elements) + "\n")

cryspy_lines.append("[ASE]")
cryspy_lines.append("ase_python = chgnet_in.py\n")

cryspy_lines.append("[EA]")
cryspy_lines.append("n_pop = 20")
cryspy_lines.append("n_crsov = 5")
cryspy_lines.append("n_perm = 2")
cryspy_lines.append("n_strain = 2")
cryspy_lines.append("n_rand = 2")
cryspy_lines.append("n_add = 3")
cryspy_lines.append("n_elim = 3")
cryspy_lines.append("n_subs = 3")
cryspy_lines.append("target = random")
cryspy_lines.append("n_elite = 2")
cryspy_lines.append("n_fittest = 10")
cryspy_lines.append("slct_func = TNM")
cryspy_lines.append("t_size = 2")
cryspy_lines.append("maxgen_ea = 0")

# Generate end_point values from table or fallback
end_points = []
for el in elements:
    if el in element_data:
        end_points.append(f"{element_data[el]['energy']:.4f}")
    else:
        end_points.append("-4.5000")
cryspy_lines.append("end_point = " + " ".join(end_points) + "\n")

cryspy_lines.append("[option]\n")

# Write cryspy.in file
with open('cryspy.in', 'w') as f:
    f.write("\n".join(cryspy_lines))

print("cryspy.in generated with accurate end_point values.")

# Prepare POSCAR files
os.makedirs('poscars', exist_ok=True)

# Map abbreviations to ASE keywords
ase_struct_map = {
    'fcc': 'fcc',
    'bcc': 'bcc',
    'hcp': 'hcp',
    'dia': 'diamond',
    'bct': 'bcc',  # fallback to bcc for bct
    'dhcp': 'hcp', # fallback to hcp for dhcp
    'ort': 'sc',   # fallback
    'mon': 'sc',   # fallback
    'sc': 'sc'
}

for el in elements:
    struct = element_data.get(el, {}).get('structure', 'sc')
    ase_struct = ase_struct_map.get(struct, 'sc')
    a_val = element_data.get(el, {}).get('a', 4.0)
    try:
        atoms = bulk(el, ase_struct, a=a_val)
    except Exception:
        atoms = bulk(el, 'sc', a=4.0)
    poscar_path = os.path.join('poscars', f"POSCAR_{el}")
    write(poscar_path, atoms, format='vasp')
    print(f"Generated POSCAR for {el} ({struct} -> ASE: {ase_struct})")

print("POSCAR files generated in 'poscars' directory.")

# Additional step: create subdirectories in current working directory and copy files
xx_tmp_dir = 'Xx_tmp'
if not os.path.isdir(xx_tmp_dir):
    print(f"Warning: {xx_tmp_dir} directory does not exist. Skipping copy step.")
else:
    for el in elements:
        struct = element_data.get(el, {}).get('structure', 'sc')
        subdir_name = f"{el}_{struct}"
        os.makedirs(subdir_name, exist_ok=True)
        # Copy contents of Xx_tmp into subdir
        for item in os.listdir(xx_tmp_dir):
            s = os.path.join(xx_tmp_dir, item)
            d = os.path.join(subdir_name, item)
            if os.path.isdir(s):
                shutil.copytree(s, d, dirs_exist_ok=True)
            else:
                shutil.copy2(s, d)
        # Copy and rename POSCAR file
        poscar_src = os.path.join('poscars', f"POSCAR_{el}")
        poscar_dst = os.path.join(subdir_name, "POSCAR")
        if os.path.isfile(poscar_src):
            shutil.copy2(poscar_src, poscar_dst)
        print(f"Created directory {subdir_name} with POSCAR and copied Xx_tmp contents.")

# Delete the 'poscars' directory to clean up intermediate files
if os.path.isdir('poscars'):
    shutil.rmtree('poscars')
    print("Cleaned up: 'poscars' directory deleted.")
