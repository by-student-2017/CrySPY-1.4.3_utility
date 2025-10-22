import pandas as pd
import matplotlib.pyplot as plt

file_path = 'data/cryspy_rslt'

# ヘッダーを自動認識して読み込み
df = pd.read_csv(file_path, sep=r'\s+')

# Optが'done'の行のみ抽出
df_done = df[df['Opt'] == 'done']

# 外れ値を除外: E_eV_atom が -15 より大きいものだけ残す
df_filtered = df_done[df_done['E_eV_atom'] > -15]

# 平均と標準偏差を計算
mean_val = df_filtered['E_eV_atom'].mean()
std_val = df_filtered['E_eV_atom'].std()

ymin = mean_val - 1.0 * std_val
ymax = mean_val + 1.0 * std_val
df_final = df_filtered[(df_filtered['E_eV_atom'] >= ymin) & (df_filtered['E_eV_atom'] <= ymax)]

# インデックスをIDとして使用
ids = df_filtered.index
energy_values = df_filtered['E_eV_atom']

# 最小エネルギーの情報（外れ値除外後）
min_idx = energy_values.idxmin()
min_id = min_idx
min_energy = energy_values.min()

# プロット設定
#plt.rcParams.update({'font.family': 'Times New Roman', 'font.size': 14})
#plt.rcParams.update({'font.family': 'Arial', 'font.size': 14})
#plt.rcParams.update({'font.size': 14})
plt.rcParams.update({'font.family': 'Liberation Serif', 'font.size': 14}) # Times New Roman互換
#plt.rcParams.update({'font.family': 'Liberation Sans', 'font.size': 14}) # Arial互換
fig, ax = plt.subplots(figsize=(10, 6))

# 散布図（CrySPYの指定を再現）
ax.scatter(
    ids,
    energy_values,
    #color='blue',
    color='#ADD8E6',
    s=12**2,                # markersize=12 → sは面積なので12^2
    edgecolors='black',     # marker_edge_color
    linewidths=1.0,         # marker_edge_width
    alpha=1.0,
    label='Energy',
    picker=True             # クリックイベント: On
)

# 最小値の点
ax.scatter(
    min_id,
    min_energy,
    color='red',
    s=12**2,
    edgecolors='black',
    linewidths=1.0,
    alpha=1.0,
    zorder=5,
    label=f'Min Energy (ID={min_id})'
)

# 破線を追加（最小エネルギーの位置）
ax.axhline(y=min_energy, color='black', linestyle='--', linewidth=1.5, alpha=0.7)
ax.axvline(x=min_id, color='black', linestyle='--', linewidth=1.5, alpha=0.7)

# 最小エネルギーの行データを取得
min_row = df_filtered.loc[min_id]
spg_num = min_row['Spg_num']
spg_sym = min_row['Spg_sym']
spg_num_opt = min_row['Spg_num_opt']
spg_sym_opt = min_row['Spg_sym_opt']
#magmom = min_row['Magmom']

# タイトル用の文字列を整形
title_text = (
    f"Lowest Energy: ID={min_id}, Energy={min_energy:.3f} eV/atom, Magmom={min_row['Magmom']}, Opt={min_row['Opt']}\n"
    f"Init-Spg: {spg_num} ({spg_sym}), OPT-Spg: {spg_num_opt} ({spg_sym_opt})"
)

ymin = min_energy - 0.2 * std_val
ymax = mean_val + 1.0 * std_val

# 軸・タイトル設定
ax.set_xlabel('Structure ID')
ax.set_ylabel('Energy (eV/atom)')
ax.set_ylim(ymin, ymax)      # ymin, ymax
#ax.set_title(None)          # タイトル非表示
ax.set_title(title_text, fontsize=12, loc='center')
ax.tick_params(direction='in', length=6, width=1)
ax.grid(True)
ax.legend()
fig.tight_layout()

# アノテーション（最初は非表示）
annot = ax.annotate("", xy=(0,0), xytext=(15,15), textcoords="offset points",
                    bbox=dict(boxstyle="round", fc="w"),
                    arrowprops=dict(arrowstyle="->"))
annot.set_visible(False)

# クリックイベント
def on_pick(event):
    ind = event.ind[0]
    row = df_filtered.iloc[ind]
    x = ids[ind]
    y = energy_values.iloc[ind]
    text = (
        f"ID={x}\n"
        f"E={y:.3f} eV/atom\n"
        f"Init-Spg: {row['Spg_num']} ({row['Spg_sym']})\n"
        f"OPT-Spg: {row['Spg_num_opt']} ({row['Spg_sym_opt']})\n"
        f"Magmom: {row['Magmom']}, Opt: {row['Opt']}"
    )
    annot.xy = (x, y)
    annot.set_text(text)
    annot.set_visible(True)
    fig.canvas.draw_idle()

fig.canvas.mpl_connect('pick_event', on_pick)

# SVG形式で保存
fig.savefig('energy_scatter_plot.svg', format='svg')
plt.show()
