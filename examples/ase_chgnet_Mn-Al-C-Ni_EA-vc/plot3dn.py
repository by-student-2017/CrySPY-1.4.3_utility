import os
import re
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from mpl_toolkits.mplot3d import proj3d

# ディレクトリとファイル選択
hull_dir = 'data/convex_hull'
cryspy_file = 'cryspy.in'

# 最新世代のファイルを選択
pattern = re.compile(r'hull_dist_all_gen_(\d+)')
max_gen = -1
hull_file = None
for fname in os.listdir(hull_dir):
    match = pattern.match(fname)
    if match:
        gen_num = int(match.group(1))
        if gen_num > max_gen:
            max_gen = gen_num
            hull_file = os.path.join(hull_dir, fname)

if hull_file is None:
    raise FileNotFoundError("No hull_dist_all_gen_* file found in data/convex_hull")

# 元素名をcryspy.inから取得
elements = []
if os.path.exists(cryspy_file):
    with open(cryspy_file, 'r') as f:
        for line in f:
            if line.strip().startswith('atype'):
                parts = line.split('=')
                if len(parts) > 1:
                    elements = parts[1].strip().split()
                break
if not elements:
    elements = ['Elem1', 'Elem2', 'Elem3', 'Elem4']

# データ読み込み
data = []
with open(hull_file, 'r') as f:
    for line in f:
        if line.strip() == '' or line.strip().startswith('ID'):
            continue
        line_clean = line.replace('(', '').replace(')', '').replace(',', ' ')
        parts = line_clean.split()
        if len(parts) >= 6:
            struct_id = int(parts[0])
            hull_dist = float(parts[1])
            composition = tuple(map(int, parts[2:6]))
            total_atoms = sum(composition)
            if total_atoms == 0:
                continue
            fractions = [c / total_atoms for c in composition]
            data.append({
                'ID': struct_id,
                'HullDistance': hull_dist,
                'Fractions': fractions,
                'Composition': composition
            })

if not data:
    raise ValueError("No valid data found in the hull file.")

# 四面体の頂点
v0 = np.array([0, 0, 0])
v1 = np.array([1, 0, 0])
v2 = np.array([0.5, np.sqrt(3)/2, 0])
v3 = np.array([0.5, np.sqrt(3)/6, np.sqrt(6)/3])
vertices = [v0, v1, v2, v3]

# バリセントリック→デカルト変換
def barycentric_to_cartesian(frac):
    return frac[0]*v0 + frac[1]*v1 + frac[2]*v2 + frac[3]*v3

coords = np.array([barycentric_to_cartesian(d['Fractions']) for d in data])
hull_distances = np.array([d['HullDistance'] for d in data])
ids = np.array([d['ID'] for d in data])
compositions = [d['Composition'] for d in data]

# Hull上とその他を分離
hull_mask = hull_distances == 0
hull_points = coords[hull_mask]
hull_ids = ids[hull_mask]
hull_comps = [compositions[i] for i in np.where(hull_mask)[0]]
other_points = coords[~hull_mask]
other_distances = hull_distances[~hull_mask]
other_ids = ids[~hull_mask]
other_comps = [compositions[i] for i in np.where(~hull_mask)[0]]

# Hull点をID順に並べ替え
sorted_indices = np.argsort(hull_ids)
sorted_hull_points = hull_points[sorted_indices]
sorted_hull_comps = [hull_comps[i] for i in sorted_indices]

# 組成ラベルを生成する関数
def format_composition(comp):
    return "-".join(f"{elements[i]}{comp[i]}" for i in range(len(elements)))

# プロット
fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111, projection='3d')

# 四面体の面
faces = [[v0, v1, v2], [v0, v1, v3], [v0, v2, v3], [v1, v2, v3]]
ax.add_collection3d(Poly3DCollection(faces, facecolors='lightgray', alpha=0.1, edgecolors='k'))

# 軸範囲を調整
all_coords = np.vstack(vertices)
x_min, x_max = all_coords[:, 0].min(), all_coords[:, 0].max()
y_min, y_max = all_coords[:, 1].min(), all_coords[:, 1].max()
z_min, z_max = all_coords[:, 2].min(), all_coords[:, 2].max()
margin = -0.1
ax.set_xlim(x_min - margin, x_max + margin)
ax.set_ylim(y_min - margin, y_max + margin)
ax.set_zlim(z_min - margin, z_max + margin)

# その他の点
sc = ax.scatter(other_points[:,0], other_points[:,1], other_points[:,2],
                c=other_distances, cmap='viridis', s=60, depthshade=True, label='Off Hull')

# Hull点と線
if hull_points.size > 0:
    hull_sc = ax.scatter(sorted_hull_points[:,0], sorted_hull_points[:,1], sorted_hull_points[:,2],
                          color='red', s=80, marker='^', label='On Hull')
    for i in range(len(sorted_hull_points)-1):
        p1 = sorted_hull_points[i]
        p2 = sorted_hull_points[i+1]
        ax.plot([p1[0], p2[0]], [p1[1], p2[1]], [p1[2], p2[2]], color='red', linewidth=2)

# カラーバー
cbar = plt.colorbar(sc, ax=ax, shrink=0.6)
cbar.set_label('Hull Distance (eV/atom)', fontsize=12)

# 頂点ラベル
for i, v in enumerate(vertices):
    ax.text(v[0], v[1], v[2], elements[i], fontsize=14, fontweight='bold')

ax.set_title(f'Quaternary Tetrahedral Plot (Gen {max_gen})', fontsize=16)
ax.set_axis_off()
ax.legend()

# 注釈
annot = ax.text2D(0.05, 0.95, "", transform=ax.transAxes)

# ホバーイベント
def on_hover(event):
    if event.inaxes != ax:
        return
    closest_text = ""
    min_dist = float('inf')
    threshold = 15  # px以内ならヒット
    # 3D座標をスクリーン座標に変換
    def project_to_2d(point):
        x2, y2, _ = proj3d.proj_transform(point[0], point[1], point[2], ax.get_proj())
        return ax.transData.transform((x2, y2))

    # Off Hull
    for i, point in enumerate(other_points):
        xp, yp = project_to_2d(point)
        dist = np.hypot(event.x - xp, event.y - yp)
        if dist < min_dist and dist < threshold:
            min_dist = dist
            closest_text = (f"Off Hull -> ID: {other_ids[i]}, "
                            f"HullDist: {other_distances[i]:.3f} eV/atom\n"
                            f"Comp: {format_composition(other_comps[i])}")

    # On Hull
    for i, point in enumerate(sorted_hull_points):
        xp, yp = project_to_2d(point)
        dist = np.hypot(event.x - xp, event.y - yp)
        if dist < min_dist and dist < threshold:
            min_dist = dist
            closest_text = (f"On Hull -> ID: {hull_ids[sorted_indices][i]}, "
                            f"HullDist: 0.000 eV/atom\n"
                            f"Comp: {format_composition(sorted_hull_comps[i])}")

    annot.set_text(closest_text)
    fig.canvas.draw_idle()

fig.canvas.mpl_connect('motion_notify_event', on_hover)

plt.show()
