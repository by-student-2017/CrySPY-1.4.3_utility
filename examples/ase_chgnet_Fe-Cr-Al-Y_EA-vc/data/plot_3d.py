import re
import pandas as pd
import plotly.express as px

# File path
file_path = 'data/convex_hull/hull_dist_all_gen_1'

# Read and parse the file
data = []
with open(file_path, 'r') as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith('ID'):
            continue
        # Match pattern: ID hull_distance (composition)
        match = re.match(r'(\d+)\s+([\d\.]+)\s+\((.*?)\)', line)
        if match:
            struct_id = int(match.group(1))
            hull_dist = float(match.group(2))
            composition = tuple(map(int, match.group(3).split(',')))
            data.append({'ID': struct_id, 'HullDistance': hull_dist, 'Composition': composition})

# Determine number of elements
if not data:
    raise ValueError("No data found in the file.")
num_elements = len(data[0]['Composition'])

# Try to extract element names from header line if available
element_names = []
with open(file_path, 'r') as f:
    header_line = f.readline()
    # If header contains element names, parse them (not guaranteed)
    # We'll fallback to generic names if not found
    # For now, assume generic names
element_names = [f'Elem{i+1}' for i in range(num_elements)]

# Compute atomic fractions and prepare dataframe
rows = []
for entry in data:
    comp = entry['Composition']
    total_atoms = sum(comp)
    fractions = [c / total_atoms for c in comp]
    row = {'ID': entry['ID'], 'HullDistance': entry['HullDistance']}
    for i, frac in enumerate(fractions):
        row[element_names[i]] = frac
    rows.append(row)

df = pd.DataFrame(rows)

# Create 3D scatter plot
x_elem, y_elem, z_elem = element_names[0], element_names[1], element_names[2]
color_elem = element_names[3] if num_elements > 3 else element_names[2]

fig = px.scatter_3d(
    df, x=x_elem, y=y_elem, z=z_elem,
    color=color_elem, size='HullDistance',
    hover_data=['ID', 'HullDistance'],
    title='CrySPY Convex Hull 3D Plot'
)

# Highlight hull points (HullDistance == 0)
hull_points = df[df['HullDistance'] == 0]
if not hull_points.empty:
    fig.add_scatter3d(
        x=hull_points[x_elem],
        y=hull_points[y_elem],
        z=hull_points[z_elem],
        mode='markers',
        marker=dict(size=6, color='red', symbol='circle'),
        name='On Hull'
    )

# Save the plot as HTML and PNG
fig.write_html('convex_hull_3d.html')
fig.write_image('convex_hull_3d.png')

print("3D convex hull plot saved as convex_hull_3d.html and convex_hull_3d.png")
