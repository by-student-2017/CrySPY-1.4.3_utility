import sys
import os
import shutil
import subprocess

elements = sys.argv[1:]
if len(elements) < 2:
    print("Usage: python make_input_qe.py Element1 Element2 [Element3 ...]")
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

# 擬ポテンシャルダウンロード
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
    "kppvol =  25  40\n",
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
    "load_struc_flag = True\n"
]

with open('cryspy.in', 'w') as f:
    f.write("\n".join(cryspy_lines))

# calc_in_qe → calc_in にテンプレートコピー
os.makedirs('calc_in', exist_ok=True)

# カラム幅設定
col_width_el = max(len(el) for el in elements) + 2
col_width_mass = 6
col_width_pp = max(len(pp_table.get(el, "UNKNOWN.UPF")) for el in elements) + 2

ntyp = len(elements)

for i in [1, 2]:
    tmp_file = os.path.join('calc_in_qe', f"{i}_pwscf_tmp.in")
    out_file = os.path.join('calc_in', f"{i}_pwscf.in")

    if not os.path.isfile(tmp_file):
        print(f"Template {tmp_file} not found!")
        continue

    with open(tmp_file) as f:
        lines = f.readlines()

    new_lines = []
    in_system_section = False
    for line in lines:
        if line.strip().lower().startswith("&system"):
            in_system_section = True

        if "__NTYP__" in line:
            line = line.replace("__NTYP__", str(ntyp))

        if in_system_section and line.strip() == "/":
            new_lines.append(f"    nspin = 2,\n")
            for idx in range(ntyp):
                mag = 0.3 if idx == 0 else 0.0
                new_lines.append(f"    starting_magnetization({idx+1}) = {mag}\n")
            in_system_section = False

        new_lines.append(line)

        if line.strip().startswith("ATOMIC_SPECIES"):
            for el in elements:
                pseudo = pp_table.get(el, "UNKNOWN.UPF")
                formatted_line = (
                    f"{el.ljust(col_width_el)}"
                    f"{'-1.0'.ljust(col_width_mass)}"
                    f"{pseudo.ljust(col_width_pp)}\n"
                )
                new_lines.append(formatted_line)

    with open(out_file, 'w') as f:
        f.writelines(new_lines)

    print(f"Generated {out_file}")

# job_cryspyコピー
job_src = os.path.join('calc_in_qe', 'job_cryspy')
job_dst = os.path.join('calc_in', 'job_cryspy')
if os.path.isfile(job_src):
    shutil.copy2(job_src, job_dst)
    print(f"Copied job_cryspy to {job_dst}")
else:
    print("Warning: job_cryspy not found in calc_in_qe")
