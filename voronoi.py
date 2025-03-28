import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat
from scipy.spatial import Voronoi, ConvexHull, voronoi_plot_2d, distance_matrix

# Load particle data from .mat file
data = loadmat("./particle.mat")
particle_info = data['particle']

area_threshold = 1

# Extract x, y coordinates and g^2 values
x_coords = np.array([entry['x'][0][0] for entry in particle_info])
y_coords = np.array([entry['y'][0][0] for entry in particle_info])
g2_values = np.array([entry['g2'][0][0] for entry in particle_info])
points = np.column_stack((x_coords, y_coords))

# Plotting function
def plot_voronoi(vor, finite_indices, g2_values, area_threshold):
    fig, ax = plt.subplots(figsize=(12, 12))
    colormap = plt.cm.get_cmap('viridis')
    norm = plt.Normalize(g2_values.min(), g2_values.max())
    for idx, index in enumerate(finite_indices):
        region = vor.regions[vor.point_region[index]]
        if ConvexHull(vor.vertices[region]).volume <= average_area * area_threshold:
            polygon = [vor.vertices[i] for i in region]
            color = colormap(norm(g2_values[index]))
            plt.fill(*zip(*polygon), color=color, edgecolor='k', alpha=0.6)
    plt.plot(vor.points[:, 0], vor.points[:, 1], 'ro', markersize=2)
    num_particles = len(finite_indices)
    volumefraction = round(volume_fraction_uniform * 100, 3)
    plt.annotate(f"Including {num_particles} particles", (0.05, 0.95), xycoords='axes fraction', fontsize=12, color='black')
    plt.annotate(f"Volumefraction: {volumefraction}%", (0.05, 0.90), xycoords='axes fraction', fontsize=12, color='black')
    # plt.annotate(f"Cell Area Threshold: {cell_area_threshold}", (0.05, 0.85), xycoords='axes fraction', fontsize=12, color='black')

    cbar = plt.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=colormap))
    cbar.set_label('$g^2$ value', rotation=270, labelpad=15)
    ax.set_aspect('equal')
    ax.set_title(f"Voronoi Cells with Area <= {area_threshold} * Average Area")
    plt.show()



# Remove acute angles
def acute_angle(polygon):
    """
    Checks if a polygon contains an acute angle
    """
    for i in range(len(polygon)):
        p1 = np.array(polygon[i - 1])
        p2 = np.array(polygon[i])
        p3 = np.array(polygon[(i + 1) % len(polygon)])
        
        v1 = p1 - p2
        v2 = p3 - p2
        
        cosine_angle = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
        angle = np.arccos(cosine_angle)
        
        if angle < np.pi / 3:  # angle < x degrees
            return True
    return False


# Calculate Voronoi diagram
vor_original = Voronoi(points)

# Identify finite regions
finite_regions = [region for region in vor_original.regions if -1 not in region and len(region) > 0]
finite_particle_indices = [i for i in range(len(vor_original.points)) if vor_original.regions[vor_original.point_region[i]] in finite_regions]

# Calculate areas of finite Voronoi cells
areas = [ConvexHull(vor_original.vertices[region]).volume for region in finite_regions]
average_area = np.mean(areas)

# Identify the Voronoi cells that are less than or equal to cell_area_threshold times the average area
valid_indices = [index for index in finite_particle_indices if ConvexHull(vor_original.vertices[vor_original.regions[vor_original.point_region[index]]]).volume <= area_threshold * average_area]

# Using the minimum distance between two particles (which is equivalent to the diameter) to get a uniform radius
dist_matrix = distance_matrix(points[valid_indices], points[valid_indices])
np.fill_diagonal(dist_matrix, np.inf)
min_distance = np.min(dist_matrix)
radius_uniform = min_distance / 2

# Calculate the total area occupied by all valid particles based on this uniform radius
particles_area = len(valid_indices) * np.pi * radius_uniform**2

# Calculate the area of the valid Voronoi cells
valid_area = [ConvexHull(vor_original.vertices[vor_original.regions[vor_original.point_region[index]]]).volume for index in valid_indices]
total_area = sum(valid_area)

# Calculate the ratio of the total valid particle area to the total valid Voronoi area using the uniform radius
volume_fraction_uniform = particles_area / total_area


# Plotting the Voronoi cells with colors based on g^2 values and annotation
plot_voronoi(vor_original, valid_indices, g2_values, area_threshold)


# Filter out Voronoi cells that contain an acute angle
no_acute_cell = [index for index in valid_indices if not acute_angle(vor_original.vertices[vor_original.regions[vor_original.point_region[index]]])]

# 重新计算有效面积
# Calculate the area of the valid Voronoi cells
valid_area = [ConvexHull(vor_original.vertices[vor_original.regions[vor_original.point_region[index]]]).volume for index in no_acute_cell]
total_area = sum(valid_area)

# Calculate the total area of valid particles based on this uniform radius
particles_area = len(no_acute_cell) * np.pi * radius_uniform**2

# Calculate the ratio of the total valid particle area to the total valid Voronoi area using the uniform radius
volume_fraction_uniform = particles_area / total_area

# Plot the Voronoi diagram without cells that have an acute angle using the correct function

plot_voronoi(vor_original, no_acute_cell, g2_values, area_threshold)

