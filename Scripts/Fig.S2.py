import numpy as np
from ase.io import read
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import LinearSegmentedColormap

def get_bond_pairs(positions, max_bond_length=1.6):
    bonds = []
    n_atoms = len(positions)
    for i in range(n_atoms):
        for j in range(i + 1, n_atoms):
            dist = np.linalg.norm(positions[i] - positions[j])
            if dist < max_bond_length:
                bonds.append((i, j))
    return bonds

def plot_structure_comparison_3d(poscar1_path, poscar2_path, labels=('DFT', 'BCM')):
    structure1 = read(poscar1_path, format='vasp')
    structure2 = read(poscar2_path, format='vasp')
    
    positions1 = structure1.get_positions()
    positions2 = structure2.get_positions()
    
    fig = plt.figure(figsize=(18, 10))
    ax = fig.add_subplot(111, projection='3d')
    
    bonds1 = get_bond_pairs(positions1)
    bonds2 = get_bond_pairs(positions2)
    
    # Custom colormaps
    cmap1 = LinearSegmentedColormap.from_list("cmap1", ["#f28e2b", "#f28e2b"])
    cmap2 = LinearSegmentedColormap.from_list("cmap2", ["#76b7b2", "#76b7b2"])
    
    # Draw bonds
    for bonds, positions, cmap in [(bonds1, positions1, cmap1), (bonds2, positions2, cmap2)]:
        for bond in bonds:
            i, j = bond
            z = (positions[i, 2] + positions[j, 2]) / 2
            color = cmap((z - positions[:, 2].min()) / (positions[:, 2].max() - positions[:, 2].min()))
            ax.plot([positions[i, 0], positions[j, 0]], 
                    [positions[i, 1], positions[j, 1]], 
                    [positions[i, 2], positions[j, 2]], 
                    color=color, linewidth=2, alpha=0.8)
    
    # Plot atoms
    for positions, cmap, label in [(positions1, cmap1, labels[0]), (positions2, cmap2, labels[1])]:
        colors = cmap((positions[:, 2] - positions[:, 2].min()) / (positions[:, 2].max() - positions[:, 2].min()))
        scatter = ax.scatter(positions[:, 0], positions[:, 1], positions[:, 2], 
                             c=colors, s=400, edgecolors='k', linewidth=1.5, alpha=0.9, label=label)

    # Remove background and grid lines
    ax.set_facecolor('white')  # Set background to white
    ax.grid(False)  # Disable grid lines

    # Set axis pane fill colors to none and edge colors to black
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False
    
    # Change the edge color of the panes
    #ax.xaxis.pane.set_edgecolor('black')
    #ax.yaxis.pane.set_edgecolor('black')
    #ax.zaxis.pane.set_edgecolor('black')

    # Set axis labels
    ax.set_xlabel('X-axis (Å)', fontsize=24, labelpad=15)
    ax.set_ylabel('Y-axis (Å)', fontsize=24, labelpad=15)
    ax.set_zlabel('Z-axis (Å)', fontsize=24, labelpad=15)
    ax.legend(fontsize=20, loc='upper right')
    # Increase the font size of the tick labels for x, y, and z axes
    ax.tick_params(axis='x', labelsize=24)  # Adjust the size as needed
    ax.tick_params(axis='y', labelsize=24)  # Adjust the size as needed
    ax.tick_params(axis='z', labelsize=24)  # Adjust the size as needed

    # Set equal aspect ratio for the 3D plot
    max_range = np.array([positions1[:, 0].ptp(), positions1[:, 1].ptp(), positions1[:, 2].ptp(),
                          positions2[:, 0].ptp(), positions2[:, 1].ptp(), positions2[:, 2].ptp()]).max() / 2.0
    mid_x = (positions1[:, 0].mean() + positions2[:, 0].mean()) / 2
    mid_y = (positions1[:, 1].mean() + positions2[:, 1].mean()) / 2
    mid_z = (positions1[:, 2].mean() + positions2[:, 2].mean()) / 2
    ax.set_xlim(mid_x - max_range, mid_x + max_range)
    ax.set_ylim(mid_y - max_range, mid_y + max_range)
    ax.set_zlim(mid_z - max_range, mid_z + max_range)

    # Add gray border around the plot
    bbox_props = dict(boxstyle="round,pad=0.3", edgecolor='gray', facecolor='none', lw=2)
    ax.text2D(0.5, 0.5, '', ha='center', va='center', fontsize=28, bbox=bbox_props)

    plt.tight_layout()
    
    # Save the 3D structure comparison figure
    plt.savefig('3d_structure_comparison.png', dpi=900, bbox_inches='tight')
    plt.close()

# Example usage
plot_structure_comparison_3d('/home/hgpark/Desktop/test/1_4-cyanocyclohexene/dft/CONTCAR', 
                             '/home/hgpark/Desktop/test/1_4-cyanocyclohexene/bcm/POSCAR')

