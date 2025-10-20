import matplotlib.pyplot as plt
import numpy as np
import re

# File paths
cryspy_file = 'cryspy.in'
hull_file = 'data/convex_hull/hull_dist_all_gen_1'

# Step 1: Extract elements and fixed info from cryspy.in
def extract_elements_and_fixed(file_path):
    elements = []
    fixed_info = []
    ll_values = []
    ul_values = []
    with open(file_path, 'r') as f:
        for line in f:
            if line.strip().startswith('atype'):
                parts = line.split('=')
                if len(parts) > 1:
                    elements = parts[1].strip().split()
            if line.strip().startswith('ll_nat'):
                ll_parts = line.split('=')
                if len(ll_parts) > 1:
                    ll_values = list(map(int, ll_parts[1].strip().split()))
            if line.strip().startswith('ul_nat'):
                ul_parts = line.split('=')
                if len(ul_parts) > 1:
                    ul_values = list(map(int, ul_parts[1].strip().split()))
    # Determine fixed elements
    for el, ll, ul in zip(elements, ll_values, ul_values):
        if ll == ul:
            fixed_info.append(f"{el}={ll}")
    return elements, fixed_info

elements, fixed_info = extract_elements_and_fixed(cryspy_file)

# Step 2: Read compositions and hull distances
def read_compositions_with_hull(file_path):
    data = []
    pattern_comp = re.compile(r'\(([\d,\s]+)\)')
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('ID'):
                continue
            parts = line.split()
            try:
                hull_dist = float(parts[1])
            except ValueError:
                continue
            match = pattern_comp.search(line)
            if match:
                nums = tuple(map(int, match.group(1).split(',')))
                data.append((nums, hull_dist))
    return data

comp_hull_data = read_compositions_with_hull(hull_file)

# Separate compositions and hull distances
compositions = [item[0] for item in comp_hull_data]
hull_dists = [item[1] for item in comp_hull_data]

# Step 3: Normalize compositions
normalized_compositions = []
for comp in compositions:
    total = sum(comp)
    normalized_compositions.append(tuple(x/total for x in comp))

# Step 4: Convert to pseudo-ternary ignoring fixed elements
# Assuming elements order: Al, C, Mn, Ni (example)
# We'll exclude any fixed element from pseudo-ternary
# For simplicity, exclude elements that appear in fixed_info
fixed_elements = [fi.split('=')[0] for fi in fixed_info]
pseudo_indices = [i for i, el in enumerate(elements) if el not in fixed_elements]

pseudo_compositions = []
for comp in normalized_compositions:
    selected = [comp[i] for i in pseudo_indices]
    total = sum(selected)
    pseudo_compositions.append(tuple(x/total for x in selected))

# Separate on-hull and off-hull based on hull distance
on_hull_data = [(pseudo, original_comp) for pseudo, original_comp, dist in zip(pseudo_compositions, compositions, hull_dists) if abs(dist) < 1e-8]
off_hull_data = [pseudo for pseudo, dist in zip(pseudo_compositions, hull_dists) if abs(dist) >= 1e-8]

# Step 5: Plot ternary diagram with annotations and fixed info in title
def plot_ternary(on_data, off_data, labels, fixed_info):
    def convert_to_xy(data):
        xs, ys = [], []
        for a, b, c in data:
            x = b + c/2
            y = (np.sqrt(3)/2)*c
            xs.append(x)
            ys.append(y)
        return xs, ys

    xs_on, ys_on = convert_to_xy([item[0] for item in on_data])
    xs_off, ys_off = convert_to_xy(off_data)

    fig, ax = plt.subplots(figsize=(8, 7))
    ax.scatter(xs_on, ys_on, color='red', s=50, alpha=0.7, label='On Hull')
    ax.scatter(xs_off, ys_off, color='blue', s=50, alpha=0.7, label='Off Hull')

    # Annotate on-hull points with full composition
    for (x, y), (_, comp) in zip(zip(xs_on, ys_on), on_data):
        annotation = ''.join(f"{el}{count}" for el, count in zip(elements, comp))
        ax.text(x + 0.02, y, annotation, fontsize=8, color='black')

    # Draw triangle
    triangle = np.array([[0, 0], [1, 0], [0.5, np.sqrt(3)/2]])
    ax.plot([triangle[0][0], triangle[1][0]], [triangle[0][1], triangle[1][1]], 'k-')
    ax.plot([triangle[1][0], triangle[2][0]], [triangle[1][1], triangle[2][1]], 'k-')
    ax.plot([triangle[2][0], triangle[0][0]], [triangle[2][1], triangle[0][1]], 'k-')

    # Axis labels from extracted elements (pseudo_indices)
    ax.text(-0.05, -0.05, labels[0], ha='center', va='center', fontsize=12)
    ax.text(1.05, -0.05, labels[1], ha='center', va='center', fontsize=12)
    ax.text(0.5, np.sqrt(3)/2 + 0.05, labels[2], ha='center', va='center', fontsize=12)

    # Add fixed info to title
    fixed_text = ', '.join(fixed_info)
    plt.title(f"CrySPY-style Pseudo-Ternary Diagram (Fixed: {fixed_text})", fontsize=14)

    ax.set_xlim(-0.1, 1.1)
    ax.set_ylim(-0.1, np.sqrt(3)/2 + 0.1)
    ax.set_aspect('equal')
    ax.axis('off')
    plt.legend(loc='upper right')
    plt.savefig('pseudo_ternary_filtered.svg', format='svg')
    plt.show()

# Automatically derive labels from elements excluding fixed ones
labels = [elements[i] for i in pseudo_indices]
plot_ternary(on_hull_data, off_hull_data, labels, fixed_info)

print("SVG file pseudo_ternary_filtered.svg generated with fixed element info in title and labels auto-derived.")
