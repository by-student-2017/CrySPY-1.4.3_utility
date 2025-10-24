import sys
import os
import shutil
import subprocess
import math

elements = sys.argv[1:]
if len(elements) < 2:
    print("Usage: python make_input_qe.py Element1 Element2 [Element3 ...]")
    sys.exit(1)

# ===== 擬ポテンシャルテーブル読み込み =====
pp_table = {}
base_url = "https://pseudopotentials.quantum-espresso.org/upf_files/"
with open('pp_table.txt') as f:
    for line in f:
        parts = line.split()
        if len(parts) >= 2 and not line.startswith('Element'):
            pp_table[parts[0]] = parts[1]

# ===== element_data_qe.txt読み込み =====
element_data = {}
with open('element_data_qe.txt') as f:
    for line in f:
        parts = line.split()
        if len(parts) >= 6 and not line.startswith('Element'):
            element_data[parts[0]] = {
                'structure': parts[1],
                'energy': float(parts[2]),
                'a': float(parts[3]),
                'b': float(parts[4]),
                'c': float(parts[5]),
                'maginit': float(parts[6])
            }

# ===== 擬ポテンシャルダウンロード =====
pp_dir = "pp"
os.makedirs(pp_dir, exist_ok=True)
for el in elements:
    pseudo = pp_table.get(el)
    if pseudo:
        file_path = os.path.join(pp_dir, pseudo)
        if not os.path.isfile(file_path):
            url = base_url + pseudo
            print(f"Downloading {pseudo} for {el}...")
            subprocess.run(["wget", "-q", "-O", file_path, url], check=True)
        else:
            print(f"{pseudo} already exists.")
    else:
        print(f"Warning: No pseudopotential entry for {el}")

# ===== CrySPY設定ファイル生成 =====
col_width = max(len(el) for el in elements)
cryspy_lines = [
    "[basic]",
    "algo = EA-vc",
    "calc_code = QE",
    "nstage = 2",
    "njob = 1",
    "jobcmd = bash",
    "jobfile = job_cryspy\n",
    "[structure]",
    "atype  = " + " ".join(f"{el:<{col_width}}" for el in elements),
    "ll_nat = " + " ".join(f"{0:<{col_width}}" for _ in elements),
    "ul_nat = " + " ".join(f"{12:<{col_width}}" for _ in elements) + "\n",
    "[QE]",
    "qe_infile = pwscf.in",
    "qe_outfile = pwscf.out",
    "kppvol = 25 40",
    "pv_term = True\n",
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
    "maxgen_ea = 5",
    "end_point = " + " ".join(f"{element_data.get(el, {}).get('energy', -4.5):.4f}" for el in elements) + "\n",
    "[option]\n"
]
with open('cryspy.in', 'w') as f:
    f.write("\n".join(cryspy_lines))
print("cryspy.in generated.")

# ===== calc_in生成 =====
src_dir = "calc_in_qe"
dst_dir = "calc_in"
if os.path.isdir(src_dir):
    shutil.copytree(src_dir, dst_dir, dirs_exist_ok=True)
else:
    os.makedirs(dst_dir, exist_ok=True)

template_files = [os.path.join(dst_dir, "1_pwscf_tmp.in"), os.path.join(dst_dir, "2_pwscf_tmp.in")]
target_files = [os.path.join(dst_dir, "1_pwscf.in"), os.path.join(dst_dir, "2_pwscf.in")]

for tmp_file, out_file in zip(template_files, target_files):
    if not os.path.isfile(tmp_file):
        continue
    with open(tmp_file) as f:
        lines = f.readlines()

    new_lines = []
    species_added = False
    skip_block = False
    nspin_value = None
    ntyp_value = len(elements)
    system_block = False

    for line in lines:
        if line.strip().startswith("&system"):
            system_block = True
        if system_block and "nspin" in line:
            nspin_value = int(line.split("=")[1].strip())
        if system_block and line.strip() == "/":
            if nspin_value == 2:
                for i in range(ntyp_value):
                    maginit = element_data.get(el, {}).get('maginit', 0.0)
                    new_lines.append(f"    starting_magnetization({i+1}) = {maginit:.2f}\n")
            system_block = False

        if "ntyp" in line:
            line = f"    ntyp = {ntyp_value}\n"

        if line.strip().startswith("ATOMIC_SPECIES") and not species_added:
            new_lines.append("ATOMIC_SPECIES\n")
            for el in elements:
                pseudo = pp_table.get(el, "UNKNOWN.UPF")
                new_lines.append(f"{el}  -1.0  {pseudo}\n")
            species_added = True
            skip_block = True
            continue

        if line.strip().startswith(("ATOMIC_POSITIONS", "CELL_PARAMETERS", "K_POINTS")):
            skip_block = True
            continue
        if skip_block and (line.strip() == "" or line.strip().startswith("/")):
            skip_block = False

        if not skip_block:
            new_lines.append(line)

    with open(out_file, 'w') as f:
        f.writelines(new_lines)

# 削除: *_tmp.in
for fname in os.listdir(dst_dir):
    if fname.endswith("_tmp.in"):
        os.remove(os.path.join(dst_dir, fname))

print("calc_in generated with 1_pwscf.in and 2_pwscf.in.")

# ===== 構造ディレクトリ生成 =====
xx_tmp_dir = 'Xx_tmp'
if not os.path.isdir(xx_tmp_dir):
    print("Warning: Xx_tmp not found.")
else:
    def get_positions(struct, a, c):
        if struct == "fcc":
            return [(0,0,0)]
        elif struct == "bcc":
            return [(0,0,0), (0.5*a,0.5*a,0.5*a)]
        elif struct == "hcp":
            return [(0,0,0), (2/3*a,1/3*a,c/2)]
        elif struct == "dhcp":
            return [
                (0,0,0),
                (a/2,a*(3)**0.5/6,c/4),
                (0,0,c/2),
                (a/2,a*(3)**0.5/6,3*c/4)
            ]
        elif struct == "dia":
            return [(0,0,0), (0.25*a,0.25*a,0.25*a)]
        else:
            return [(0,0,0)]

    def get_cell_parameters(struct, a, b, c):
        if struct in ["fcc", "bcc", "sc", "dia"]:
            return [(a,0,0),(0,b,0),(0,0,c)]
        elif struct in ["hcp","dhcp"]:
            return [(a,0,0),(-a/2,a*math.sqrt(3)/2,0),(0,0,c)]
        else:
            return [(a,0,0),(0,b,0),(0,0,c)]

    for el in elements:
        struct = element_data.get(el, {}).get('structure', 'fcc')
        a_val = element_data.get(el, {}).get('a', 4.0)
        b_val = element_data.get(el, {}).get('b', a_val)
        c_val = element_data.get(el, {}).get('c', a_val)

        subdir = f"{el}_{struct}"
        os.makedirs(subdir, exist_ok=True)

        # Xx_tmpコピー
        for item in os.listdir(xx_tmp_dir):
            s = os.path.join(xx_tmp_dir, item)
            d = os.path.join(subdir, item)
            if os.path.isdir(s):
                shutil.copytree(s, d, dirs_exist_ok=True)
            else:
                shutil.copy2(s, d)

        # pwscf_tmp.inを反映
        template_path = os.path.join(xx_tmp_dir, "pwscf_tmp.in")
        if os.path.isfile(template_path):
            with open(template_path) as f:
                lines = f.readlines()

            new_lines = []
            species_added = False
            skip_block = False
            positions = get_positions(struct, a_val, c_val)
            nat_value = len(positions)
            nspin_value = None
            system_block = False

            for line in lines:
                if line.strip().startswith("&system"):
                    system_block = True
                if system_block and "nspin" in line:
                    nspin_value = int(line.split("=")[1].strip())
                if system_block and line.strip() == "/":
                    if nspin_value == 2:
                        new_lines.append(f"    starting_magnetization(1) = {maginit:.2f}\n")
                    system_block = False

                if "nat" in line:
                    line = f"    nat = {nat_value}\n"
                if "ntyp" in line:
                    line = "    ntyp = 1\n"
                if line.strip().startswith("ATOMIC_SPECIES") and not species_added:
                    new_lines.append("ATOMIC_SPECIES\n")
                    pseudo = pp_table.get(el, "UNKNOWN.UPF")
                    new_lines.append(f"{el}  -1.0  {pseudo}\n")
                    species_added = True
                    skip_block = True
                    continue
                if line.strip().startswith(("ATOMIC_POSITIONS", "CELL_PARAMETERS", "K_POINTS")):
                    skip_block = True
                    continue
                if skip_block and (line.strip() == "" or line.strip().startswith("/")):
                    skip_block = False
                if not skip_block:
                    new_lines.append(line)

            # 追加: 座標・セル・k点
            new_lines.append("ATOMIC_POSITIONS (angstrom)\n")
            for pos in positions:
                new_lines.append(f"{el}  {pos[0]:.6f}  {pos[1]:.6f}  {pos[2]:.6f}\n")
            new_lines.append("CELL_PARAMETERS (angstrom)\n")
            for vec in get_cell_parameters(struct, a_val, b_val, c_val):
                new_lines.append(f"{vec[0]:.6f} {vec[1]:.6f} {vec[2]:.6f}\n")
            new_lines.append("K_POINTS automatic\n")
            nkx = max(1, int(round(40 / a_val)))
            nky = max(1, int(round(40 / b_val)))
            nkz = max(1, int(round(40 / c_val)))
            new_lines.append(f"{nkx} {nky} {nkz} 0 0 0\n")

            with open(os.path.join(subdir, "pwscf.in"), 'w') as f:
                f.writelines(new_lines)

        # 削除: pwscf_tmp.in
        tmp_in_path = os.path.join(subdir, "pwscf_tmp.in")
        if os.path.isfile(tmp_in_path):
            os.remove(tmp_in_path)
            print(f"Removed template file: {tmp_in_path}")

        print(f"Created {subdir} with pwscf.in and Xx_tmp contents.")
