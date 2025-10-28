import os
import re
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

cryspy_in_file = 'cryspy.in'
cryspy_rslt_file = 'data/cryspy_rslt'

# --- 元素と固定情報 ---
def extract_elements_and_fixed(file_path):
    elements, fixed_info = [], []
    ll_values, ul_values = [], []
    with open(file_path, 'r') as f:
        for line in f:
            if line.strip().startswith('atype'):
                elements = line.split('=')[1].strip().split()
            if line.strip().startswith('ll_nat'):
                ll_values = list(map(int, line.split('=')[1].strip().split()))
            if line.strip().startswith('ul_nat'):
                ul_values = list(map(int, line.split('=')[1].strip().split()))
    for el, ll, ul in zip(elements, ll_values, ul_values):
        if ll == ul:
            fixed_info.append(f"{el}={ll}")
    return elements, fixed_info

elements, fixed_info = extract_elements_and_fixed(cryspy_in_file)
fixed_elements = [fi.split('=')[0] for fi in fixed_info]

# --- 組成・エネルギー・追加情報 ---
pattern_comp = re.compile(r'\(([\d,\s]+)\)')

def read_rslt(file_path):
    ids, gens, compositions, energies = [], [], [], []
    spg_nums, spg_syms, spg_nums_opt, spg_syms_opt = [], [], [], []
    magmoms, opts = [], []

    header_map = {}
    with open(file_path, 'r') as f:
        lines = f.readlines()

    # ヘッダー行を解析
    for line in lines:
        if line.strip().startswith("Gen"):
            headers = line.strip().split()
            header_map = {name: idx for idx, name in enumerate(headers)}
            break

    # CrySPY形式の必須列
    required_cols = ["Gen", "Spg_num", "Spg_sym", "Spg_num_opt", "Spg_sym_opt", "Ef_eV_atom", "Magmom", "Opt"]
    for col in required_cols:
        if col not in header_map:
            raise ValueError(f"Missing column: {col}")

    # データ行を処理
    for line in lines:
        line = line.strip()
        if not line or line.startswith("Gen"):
            continue

        # 組成抽出
        match = pattern_comp.search(line)
        comp_nums = ()
        if match:
            comp_nums = tuple(map(int, match.group(1).split(',')))
            line_clean = pattern_comp.sub('', line).strip()
        else:
            line_clean = line

        parts = line_clean.split()
        if len(parts) < len(header_map):
            continue

        try:
            # IDはCrySPYの最初の列（0,1,2...）を使用
            id_val = parts[0]
            gen = parts[1]
            spg_num = parts[2]
            spg_sym = parts[3]
            spg_num_opt = parts[4]
            spg_sym_opt = parts[5]
            energy = float(parts[6])
            magmom = parts[-2]
            opt = parts[-1]
        except (ValueError, IndexError):
            continue

        ids.append(id_val)
        gens.append(gen)
        compositions.append(comp_nums)
        energies.append(energy)
        spg_nums.append(spg_num)
        spg_syms.append(spg_sym)
        spg_nums_opt.append(spg_num_opt)
        spg_syms_opt.append(spg_sym_opt)
        magmoms.append(magmom)
        opts.append(opt)

    return ids, gens, compositions, energies, spg_nums, spg_syms, spg_nums_opt, spg_syms_opt, magmoms, opts

(ids, gens, compositions, energies, spg_nums, spg_syms,
 spg_nums_opt, spg_syms_opt, magmoms, opts) = read_rslt(cryspy_rslt_file)

# --- 擬似三元組成 ---
pseudo_indices = [i for i, el in enumerate(elements) if el not in fixed_elements]
pseudo_compositions = []
for comp in compositions:
    selected = [comp[i] for i in pseudo_indices]
    total = sum(selected)
    pseudo_compositions.append(tuple(x/total if total > 0 else 0 for x in selected))

# --- 三角座標変換 ---
def barycentric_to_cartesian(a, b, c):
    x = b + c/2
    y = (np.sqrt(3)/2)*c
    return x, y

xy_points = [barycentric_to_cartesian(*pc) for pc in pseudo_compositions]
xs, ys = zip(*xy_points)

# --- 補間 ---
grid_x, grid_y = np.meshgrid(np.linspace(0, 1, 200), np.linspace(0, np.sqrt(3)/2, 200))
grid_z = griddata(xy_points, energies, (grid_x, grid_y), method='linear')

# --- プロット ---
fig, ax = plt.subplots(figsize=(8, 7))
contour = ax.contourf(grid_x, grid_y, grid_z, levels=20, cmap='coolwarm', alpha=0.8)
cbar = plt.colorbar(contour, ax=ax)
cbar.set_label("Ef_eV_atom")

scatter = ax.scatter(xs, ys, c=energies, cmap='coolwarm', s=60, edgecolors='white', linewidths=1,
                     vmin=min(energies), vmax=max(energies))

triangle = np.array([[0, 0], [1, 0], [0.5, np.sqrt(3)/2]])
ax.plot([triangle[0][0], triangle[1][0]], [triangle[0][1], triangle[1][1]], 'k-')
ax.plot([triangle[1][0], triangle[2][0]], [triangle[1][1], triangle[2][1]], 'k-')
ax.plot([triangle[2][0], triangle[0][0]], [triangle[2][1], triangle[0][1]], 'k-')

labels = [elements[i] for i in pseudo_indices]
ax.text(-0.05, -0.05, labels[0], ha='center', va='center', fontsize=12)
ax.text(1.05, -0.05, labels[1], ha='center', va='center', fontsize=12)
ax.text(0.5, np.sqrt(3)/2 + 0.05, labels[2], ha='center', va='center', fontsize=12)

fixed_text = ', '.join(fixed_info) if fixed_info else "None"
plt.title(f"Ternary Diagram with Energy Contour (Fixed: {fixed_text})", fontsize=14)

ax.set_xlim(-0.1, 1.1)
ax.set_ylim(-0.1, np.sqrt(3)/2 + 0.1)
ax.set_aspect('equal')
ax.axis('off')

plt.tight_layout()

# --- ホバー注釈 ---
annot = ax.text(0.02, 0.98, "", transform=ax.transAxes, fontsize=10, va='top')

def on_hover(event):
    if event.inaxes != ax:
        return
    closest_text = ""
    min_dist = float('inf')
    threshold = 15  # px以内ならヒット

    for i, (x, y) in enumerate(xy_points):
        xp, yp = ax.transData.transform((x, y))
        dist = np.hypot(event.x - xp, event.y - yp)
        if dist < min_dist and dist < threshold:
            min_dist = dist
            comp_text = ', '.join(f"{el}={count}" for el, count in zip(elements, compositions[i]))
            if fixed_info:
                comp_text += f" | Fixed: {', '.join(fixed_info)}"
            closest_text = (
                f"ID: {ids[i]}\nGen: {gens[i]}\nEf: {energies[i]:.3f} eV/atom\nComp: {comp_text}\n"
                f"Ini_Sym: {spg_nums[i]} ({spg_syms[i]})\n"
                f"Opt_sym: {spg_nums_opt[i]} ({spg_syms_opt[i]})\n"
                f"Magmom: {magmoms[i]}\nOpt: {opts[i]}"
            )

    annot.set_text(closest_text)
    fig.canvas.draw_idle()

fig.canvas.mpl_connect('motion_notify_event', on_hover)

plt.show()
