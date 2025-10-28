import re
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from matplotlib.path import Path

cryspy_in_file = 'cryspy.in'
cryspy_rslt_file = 'data/cryspy_rslt'

# --- 元素と固定情報抽出 ---
def extract_elements_and_fixed(file_path):
    elements, ll_values, ul_values = [], [], []
    with open(file_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith('atype'):
                elements = line.split('=')[1].split()
            elif line.startswith('ll_nat'):
                ll_values = list(map(int, line.split('=')[1].split()))
            elif line.startswith('ul_nat'):
                ul_values = list(map(int, line.split('=')[1].split()))
    fixed_info = [f"{el}={ll}" for el, ll, ul in zip(elements, ll_values, ul_values) if ll == ul]
    return elements, fixed_info

elements, fixed_info = extract_elements_and_fixed(cryspy_in_file)
fixed_elements = [fi.split('=')[0] for fi in fixed_info]

# --- 結果ファイル読み込み ---
pattern_comp = re.compile(r'\(([\d,\s]+)\)')

def read_rslt(file_path):
    ids, gens, compositions, energies = [], [], [], []
    spg_nums, spg_syms, spg_nums_opt, spg_syms_opt, magmoms, opts = [], [], [], [], [], []

    with open(file_path) as f:
        lines = f.readlines()

    try:
        headers = next(line.strip().split() for line in lines if line.strip().startswith("Gen"))
    except StopIteration:
        raise ValueError("ヘッダー行 (Gen ...) が見つかりません。ファイル内容を確認してください。")

    header_map = {name: idx for idx, name in enumerate(headers)}

    for line in lines:
        line = line.strip()
        if not line or line.startswith("Gen"):
            continue

        match = pattern_comp.search(line)
        comp_nums = tuple(map(int, match.group(1).split(','))) if match else ()
        line_clean = pattern_comp.sub('', line).strip()
        parts = line_clean.split()
        if len(parts) < len(header_map):
            continue

        try:
            ids.append(parts[0])
            gens.append(parts[1])
            spg_nums.append(parts[2])
            spg_syms.append(parts[3])
            spg_nums_opt.append(parts[4])
            spg_syms_opt.append(parts[5])
            energies.append(float(parts[6]) if parts[6] != 'NaN' else np.nan)
            magmoms.append(parts[-2])
            opts.append(parts[-1])
            compositions.append(comp_nums)
        except ValueError:
            continue

    return ids, gens, compositions, energies, spg_nums, spg_syms, spg_nums_opt, spg_syms_opt, magmoms, opts

(ids, gens, compositions, energies, spg_nums, spg_syms,
 spg_nums_opt, spg_syms_opt, magmoms, opts) = read_rslt(cryspy_rslt_file)

# --- 擬似三元組成 ---
pseudo_indices = [i for i, el in enumerate(elements) if el not in fixed_elements]
pseudo_compositions = [
    tuple(x / sum(comp[i] for i in pseudo_indices) if sum(comp[i] for i in pseudo_indices) > 0 else 0
          for x in [comp[i] for i in pseudo_indices])
    for comp in compositions
]

# --- 三角座標変換 ---
def barycentric_to_cartesian(a, b, c):
    return b + c / 2, (np.sqrt(3) / 2) * c

xy_points = [barycentric_to_cartesian(*pc) for pc in pseudo_compositions]

# --- NaN 除外 ---
energies = np.array(energies)
xy_points = np.array(xy_points)
valid_mask = ~np.isnan(energies)
xy_points_valid = xy_points[valid_mask]
energies_valid = energies[valid_mask]

# --- 補間 ---
grid_x, grid_y = np.meshgrid(np.linspace(0, 1, 200), np.linspace(0, np.sqrt(3)/2, 200))
grid_z = griddata(xy_points_valid, energies_valid, (grid_x, grid_y), method='linear')

# --- 三角形外をマスク ---
triangle = np.array([[0, 0], [1, 0], [0.5, np.sqrt(3)/2]])
path = Path(triangle)
mask = ~path.contains_points(np.c_[grid_x.ravel(), grid_y.ravel()])
grid_z = grid_z.copy()
grid_z.ravel()[mask] = np.nan

# --- プロット ---
fig, ax = plt.subplots(figsize=(8, 7))
contour = ax.contourf(grid_x, grid_y, grid_z, levels=20, cmap='coolwarm', alpha=0.8)
cbar = plt.colorbar(contour, ax=ax)
cbar.set_label("Formation energy, Ef (eV/atom)")

ax.scatter(xy_points_valid[:, 0], xy_points_valid[:, 1], c=energies_valid, cmap='coolwarm', s=60,
           edgecolors='white', linewidths=1, vmin=min(energies_valid), vmax=max(energies_valid))

# 三角形枠
ax.plot(*triangle[[0, 1]].T, 'k-')
ax.plot(*triangle[[1, 2]].T, 'k-')
ax.plot(*triangle[[2, 0]].T, 'k-')

labels = [elements[i] for i in pseudo_indices]
ax.text(-0.05, -0.05, labels[0], ha='center', va='center', fontsize=12)
ax.text(1.05, -0.05, labels[1], ha='center', va='center', fontsize=12)
ax.text(0.5, np.sqrt(3)/2 + 0.05, labels[2], ha='center', va='center', fontsize=12)

plt.title(f"Ternary Diagram with Energy Contour (Fixed: {', '.join(fixed_info) or 'None'})", fontsize=14)
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
    closest_text, min_dist = "", float('inf')
    for i, (x, y) in enumerate(xy_points_valid):
        xp, yp = ax.transData.transform((x, y))
        dist = np.hypot(event.x - xp, event.y - yp)
        if dist < min_dist and dist < 15:
            min_dist = dist
            comp_text = ', '.join(f"{el}={count}" for el, count in zip(elements, compositions[i]))
            if fixed_info:
                comp_text += f" | Fixed: {', '.join(fixed_info)}"
            closest_text = (
                f"ID: {ids[i]}\nGen: {gens[i]}\nEf: {energies_valid[i]:.3f} eV/atom\nComp: {comp_text}\n"
                f"Ini_Sym: {spg_nums[i]} ({spg_syms[i]})\n"
                f"Opt_sym: {spg_nums_opt[i]} ({spg_syms_opt[i]})\n"
                f"Magmom: {magmoms[i]}\nOpt: {opts[i]}"
            )
    annot.set_text(closest_text)
    fig.canvas.draw_idle()

fig.canvas.mpl_connect('motion_notify_event', on_hover)
plt.show()
