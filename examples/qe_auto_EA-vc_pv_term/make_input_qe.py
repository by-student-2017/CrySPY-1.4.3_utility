import sys
import os
import shutil
import subprocess
import math

elements = sys.argv[1:]
if len(elements) < 1:
    print("Usage: python make_input_qe.py Element1 [Element2 ...]")
    sys.exit(1)

# 擬ポテンシャルテーブル読み込み
pp_table = {}
base_url = "https://pseudopotentials.quantum-espresso.org/upf_files/"
with open('pp_table.txt') as f:
    for line in f:
        parts = line.split()
        if len(parts) >= 2 and not line.startswith('Element'):
            pp_table[parts[0]] = parts[1]

# element_data_qe.txt読み込み（構造情報と格子定数）
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
                'c': float(parts[5])
            }

# ppディレクトリ作成
pp_dir = "pp"
os.makedirs(pp_dir, exist_ok=True)

# 必要な擬ポテンシャルのみダウンロード
for el in elements:
    pseudo = pp_table.get(el)
    if pseudo:
        file_path = os.path.join(pp_dir, pseudo)
        if not os.path.isfile(file_path):
            url = base_url + pseudo
            print(f"Downloading {pseudo} for {el} from {url}...")
            try:
                subprocess.run(["wget", "-q", "-O", file_path, url], check=True)
                print(f"Downloaded: {file_path}")
            except subprocess.CalledProcessError:
                print(f"Warning: Failed to download {pseudo} for {el}")
        else:
            print(f"{pseudo} already exists in {pp_dir}")
    else:
        print(f"Warning: No pseudopotential entry for {el} in pp_table.txt")


max_len = max(len(el) for el in elements)
atype_line = "atype  = " + " ".join(el.ljust(max_len) for el in elements)
ll_nat_line = "ll_nat = " + " ".join("0".ljust(max_len) for _ in elements)
ul_nat_line = "ul_nat = " + " ".join("12".ljust(max_len) for _ in elements)

end_points = []
for el in elements:
    if el in element_data and 'energy' in element_data[el]:
        end_points.append(f"{element_data[el]['energy']:.4f}")
    else:
        end_points.append("-4.5000")  # デフォルト値
end_point_line = "end_point  = " + " ".join(ep.ljust(max_len) for ep in end_points)

# CrySPY設定ファイル生成
cryspy_lines = [
    "[basic]",
    "algo = EA-vc",
    "calc_code = QE",
    "nstage = 2",
    "njob = 1",
    "jobcmd = bash",
    "jobfile = job_cryspy\n",
    "[structure]",
    atype_line,
    ll_nat_line,
    ul_nat_line + "\n",
    "[QE]",
    "qe_infile = pwscf.in",
    "qe_outfile = pwscf.out",
    "kppvol =  25  40",
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
    "maxgen_ea = 0",
    end_point_line + "\n",
    "[option]\n"
]

with open('cryspy.in', 'w') as f:
    f.write("\n".join(cryspy_lines))

# 原子座標生成
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

# 格子ベクトル生成
def get_cell_parameters(struct, a, b, c):
    if struct in ["fcc", "bcc", "sc", "dia"]:
        return [
            (a, 0.0, 0.0),
            (0.0, b, 0.0),
            (0.0, 0.0, c)
        ]
    elif struct in ["hcp", "dhcp"]:
        return [
            (a, 0.0, 0.0),
            (-a/2, a*math.sqrt(3)/2, 0.0),
            (0.0, 0.0, c)
        ]
    elif struct == "ort":  # orthorhombic
        return [
            (a, 0.0, 0.0),
            (0.0, b, 0.0),
            (0.0, 0.0, c)
        ]
    elif struct == "mon":  # monoclinic (β angle assumed ~90° for simplicity)
        beta = math.radians(100)  # example angle
        return [
            (a, 0.0, 0.0),
            (0.0, b, 0.0),
            (c*math.cos(beta), 0.0, c*math.sin(beta))
        ]
    else:
        return [
            (a, 0.0, 0.0),
            (0.0, b, 0.0),
            (0.0, 0.0, c)
        ]

# Xx_tmpからサブディレクトリ作成
xx_tmp_dir = 'Xx_tmp'
pwscf_template = os.path.join(xx_tmp_dir, 'pwscf_tmp.in')

if os.path.isfile(pwscf_template):
    for el in elements:
        struct = element_data.get(el, {}).get('structure', 'sc')
        a_val = element_data.get(el, {}).get('a', 4.0)
        b_val = element_data.get(el, {}).get('b', a_val)
        c_val = element_data.get(el, {}).get('c', a_val)

        positions = get_positions(struct, a_val, c_val)
        nat_value = len(positions)

        subdir_name = f"{el}_{struct}"
        os.makedirs(subdir_name, exist_ok=True)

        # Xx_tmpの内容コピー
        for item in os.listdir(xx_tmp_dir):
            if item == 'pwscf_tmp.in':
                continue
            s = os.path.join(xx_tmp_dir, item)
            d = os.path.join(subdir_name, item)
            if os.path.isdir(s):
                shutil.copytree(s, d, dirs_exist_ok=True)
            else:
                shutil.copy2(s, d)

        # pwscf.in再構築
        with open(pwscf_template) as f:
            lines = f.readlines()

        new_lines = []
        in_system_section = False
        for line in lines:
            if "starting_magnetization" in line or line.strip().startswith("ATOMIC_SPECIES") or line.strip().startswith("ATOMIC_POSITIONS"):
                continue

            if line.strip().lower().startswith("&system"):
                in_system_section = True

            if "nat" in line:
                line = f"    nat = {nat_value}\n"
            if "ntyp" in line:
                line = "    ntyp = 1\n"

            if in_system_section and line.strip() == "/":
                new_lines.append(f"    nspin = 2,\n")
                new_lines.append(f"    starting_magnetization(1) = 0.3\n")
                in_system_section = False

            new_lines.append(line)

        # ATOMIC_SPECIES
        pseudo = pp_table.get(el, "UNKNOWN.UPF")
        new_lines.append("ATOMIC_SPECIES\n")
        new_lines.append(f"{el.ljust(4)}  -1.0  {pseudo}\n")

        # ATOMIC_POSITIONS
        new_lines.append("ATOMIC_POSITIONS (angstrom)\n")
        for pos in positions:
            new_lines.append(f"{el}  {pos[0]:.6f}  {pos[1]:.6f}  {pos[2]:.6f}\n")

        # CELL_PARAMETERS
        new_lines.append("CELL_PARAMETERS (angstrom)\n")
        for vec in get_cell_parameters(struct, a_val, b_val, c_val):
            new_lines.append(f"{vec[0]:.6f} {vec[1]:.6f} {vec[2]:.6f}\n")

        # K_POINTS automatic
        k_density = 40
        nkx = max(1, int(round(k_density / a_val)))
        nky = max(1, int(round(k_density / b_val)))
        nkz = max(1, int(round(k_density / c_val)))
        new_lines.append("K_POINTS automatic\n")
        new_lines.append(f"{nkx} {nky} {nkz} 0 0 0\n")

        pwscf_dst = os.path.join(subdir_name, 'pwscf.in')
        with open(pwscf_dst, 'w') as f:
            f.writelines(new_lines)

        print(f"Created directory {subdir_name} with pwscf.in and Xx_tmp contents.")
else:
    print("Error: pwscf_tmp.in not found in Xx_tmp")
