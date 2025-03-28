

% Clear workspace and close figures
clear all;
close all;
clc;

% Load particle data from .mat file
load('particle.mat');

% Define area threshold for filtering voronoi cells
area_threshold = 1;

% Extract x, y coordinates and g^2 values
x_coords = zeros(length(particle), 1);
y_coords = zeros(length(particle), 1);
g2_values = zeros(length(particle), 1);

for i = 1:length(particle)
    x_coords(i) = particle(i).x;
    y_coords(i) = particle(i).y;
    g2_values(i) = particle(i).g2;
end

points = [x_coords, y_coords];

% Calculate Voronoi diagram
[vor_vertices, vor_cells] = voronoin(points);

% Identify finite regions
finite_regions = [];
finite_particle_indices = [];

for i = 1:length(vor_cells)
    if ~any(ismember(vor_cells{i}, 1)) && ~isempty(vor_cells{i})
        finite_regions = [finite_regions; i];
        finite_particle_indices = [finite_particle_indices; i];
    end
end

% Calculate areas of finite Voronoi cells
areas = zeros(length(finite_regions), 1);
for i = 1:length(finite_regions)
    region_idx = finite_regions(i);
    vertices = vor_vertices(vor_cells{region_idx}, :);
    areas(i) = polyarea(vertices(:,1), vertices(:,2));
end

average_area = mean(areas);

% Identify the Voronoi cells that are less than or equal to area_threshold times the average area
valid_indices = [];
for i = 1:length(finite_particle_indices)
    index = finite_particle_indices(i);
    vertices = vor_vertices(vor_cells{index}, :);
    area = polyarea(vertices(:,1), vertices(:,2));
    if area <= area_threshold * average_area
        valid_indices = [valid_indices; index];
    end
end

% Using the minimum distance between two particles to get a uniform radius
dist_matrix = pdist2(points(valid_indices,:), points(valid_indices,:));
dist_matrix(eye(size(dist_matrix))==1) = Inf;
min_distance = min(dist_matrix(:));
radius_uniform = min_distance / 2;

% Calculate the total area occupied by all valid particles based on this uniform radius
particles_area = length(valid_indices) * pi * radius_uniform^2;

% Calculate the area of the valid Voronoi cells
valid_area = zeros(length(valid_indices), 1);
for i = 1:length(valid_indices)
    index = valid_indices(i);
    vertices = vor_vertices(vor_cells{index}, :);
    valid_area(i) = polyarea(vertices(:,1), vertices(:,2));
end
total_area = sum(valid_area);

% Calculate the volume fraction
volume_fraction_uniform = particles_area / total_area;
volumefraction_percent = round(volume_fraction_uniform * 100, 3);

% Plot the Voronoi diagram with the valid cells
plot_voronoi_diagram(vor_vertices, vor_cells, valid_indices, g2_values, average_area, area_threshold, volume_fraction_uniform, x_coords, y_coords);

% Filter out Voronoi cells that contain an acute angle
no_acute_cell = [];
for i = 1:length(valid_indices)
    index = valid_indices(i);
    vertices = vor_vertices(vor_cells{index}, :);
    if ~has_acute_angle(vertices)
        no_acute_cell = [no_acute_cell; index];
    end
end

% Recalculate valid area without acute angles
valid_area_no_acute = zeros(length(no_acute_cell), 1);
for i = 1:length(no_acute_cell)
    index = no_acute_cell(i);
    vertices = vor_vertices(vor_cells{index}, :);
    valid_area_no_acute(i) = polyarea(vertices(:,1), vertices(:,2));
end
total_area_no_acute = sum(valid_area_no_acute);

% Recalculate particle area and volume fraction
particles_area_no_acute = length(no_acute_cell) * pi * radius_uniform^2;
volume_fraction_uniform_no_acute = particles_area_no_acute / total_area_no_acute;

% Plot the Voronoi diagram without cells that have an acute angle
figure;
plot_voronoi_diagram(vor_vertices, vor_cells, no_acute_cell, g2_values, average_area, area_threshold, volume_fraction_uniform_no_acute, x_coords, y_coords);
title('Voronoi Cells with No Acute Angles');

% Save results to a file
result_data = struct('volume_fraction', volume_fraction_uniform, ...
                     'volume_fraction_no_acute', volume_fraction_uniform_no_acute, ...
                     'num_particles', length(valid_indices), ...
                     'num_particles_no_acute', length(no_acute_cell), ...
                     'radius_uniform', radius_uniform);
save('voronoi_results.mat', 'result_data');

% Display results
fprintf('Analysis completed:\n');
fprintf('  Number of particles: %d\n', length(valid_indices));
fprintf('  Number of particles without acute angles: %d\n', length(no_acute_cell));
fprintf('  Volume fraction: %.2f%%\n', volume_fraction_uniform * 100);
fprintf('  Volume fraction (no acute): %.2f%%\n', volume_fraction_uniform_no_acute * 100);
fprintf('  Results saved to voronoi_results.mat\n');

% Helper functions

function result = has_acute_angle(polygon)
    % Checks if a polygon contains an acute angle
    n = size(polygon, 1);
    result = false;
    
    for i = 1:n
        p1 = polygon(mod(i-2, n)+1, :);
        p2 = polygon(i, :);
        p3 = polygon(mod(i, n)+1, :);
        
        v1 = p1 - p2;
        v2 = p3 - p2;
        
        cosine_angle = dot(v1, v2) / (norm(v1) * norm(v2));
        angle = acos(cosine_angle);
        
        if angle < pi / 3  % angle < 60 degrees
            result = true;
            return;
        end
    end
end

function plot_voronoi_diagram(vor_vertices, vor_cells, indices, g2_values, avg_area, area_threshold, volume_fraction, x_coords, y_coords)
    figure('Position', [100, 100, 800, 800]);
    
    % Create colormap
    colormap(viridis(256));
    cmap = colormap;
    
    % Normalize g2 values
    min_g2 = min(g2_values);
    max_g2 = max(g2_values);
    norm_factor = size(cmap, 1);
    
    hold on;
    for i = 1:length(indices)
        index = indices(i);
        vertices = vor_vertices(vor_cells{index}, :);
        
        % Skip if vertices include point at infinity
        if any(vertices(:,1) > 1e10) || any(vertices(:,2) > 1e10)
            continue;
        end
        
        % Calculate area of this cell to check threshold
        cell_area = polyarea(vertices(:,1), vertices(:,2));
        if cell_area <= avg_area * area_threshold
            % Get color based on g2 value
            g2_val = g2_values(index);
            color_idx = min(max(1, round((g2_val - min_g2) / (max_g2 - min_g2) * norm_factor)), norm_factor);
            color = cmap(color_idx, :);
            
            % Fill polygon with this color
            fill(vertices(:,1), vertices(:,2), color, 'EdgeColor', 'k', 'FaceAlpha', 0.6);
        end
    end
    
    % Plot all points
    plot(x_coords, y_coords, 'ro', 'MarkerSize', 2);
    
    % Add annotations
    text(min(x_coords) + 50, max(y_coords) - 50, ['Including ' num2str(length(indices)) ' particles'], 'FontSize', 12);
    text(min(x_coords) + 50, max(y_coords) - 100, ['Volumefraction: ' num2str(volume_fraction * 100) '%'], 'FontSize', 12);
    
    % Add colorbar
    c = colorbar;
    ylabel(c, 'g^2 value', 'Rotation', 270, 'FontSize', 12);
    c.Label.Position(1) = 3.5;
    
    axis equal;
    title(['Voronoi Cells with Area <= ' num2str(area_threshold) ' * Average Area']);
    hold off;
end

function cmap = viridis(m)
    % Implementation of viridis colormap (similar to Python's viridis)
    if nargin < 1
       m = 256;
    end
    
    % Base colors from Python's viridis
    base_colors = [
        0.267004, 0.004874, 0.329415;
        0.282623, 0.140926, 0.457517;
        0.253935, 0.265254, 0.529983;
        0.202430, 0.367780, 0.553248;
        0.146972, 0.439703, 0.537059;
        0.101579, 0.494877, 0.505496;
        0.116257, 0.546610, 0.456395;
        0.199742, 0.589765, 0.392134;
        0.320620, 0.622924, 0.301775;
        0.454454, 0.643837, 0.212687;
        0.581269, 0.653565, 0.134670;
        0.697123, 0.653075, 0.074708;
        0.800692, 0.640797, 0.028838;
        0.887749, 0.618842, 0.017182;
        0.962581, 0.585847, 0.045005;
        0.998271, 0.546236, 0.120419;
        0.994141, 0.476813, 0.176758;
        0.974139, 0.393609, 0.207318;
        0.939215, 0.299271, 0.204227;
        0.880267, 0.196756, 0.168791;
    ];
    
    % Interpolate to desired number of colors
    n_colors = size(base_colors, 1);
    interp_vals = linspace(1, n_colors, m);
    
    % Initialize output colormap
    cmap = zeros(m, 3);
    
    % Interpolate each RGB channel
    for i = 1:3
        cmap(:, i) = interp1(1:n_colors, base_colors(:, i), interp_vals, 'pchip');
    end
end 