%Updated to first release version of PeGS2 for 05/2024
%Adapted from voronoi.m script to fit the PeGS2 modular framework
%This module performs Voronoi analysis on particle data to calculate volume fractions and analyze local structure
% 
% **I**The input of this module is the standardised MATLAB **particle** structure array. Each particle in the image has the following fields:
% 
% - `id` : assigned ID number of the particle
% - `x`: x coordinate of the particle centre, in *pixels*
% - `y` : y coordinate of the particle centre, in *pixels*
% - `r` : radius of the particle, in *pixels*
% - other fields from the particle structure (can be used if available)
% 
% In addition to this, the following user inputs are required:
% 
% - `area_threshold` : the maximum area multiplier relative to average area for valid Voronoi cells (default: 1)
% - `acute_angle_filter` : whether to filter out Voronoi cells with acute angles (default: true)
% - `acute_angle_threshold` : angle in degrees below which a cell corner is considered acute (default: 60)
% - `save_figures` : whether to save the Voronoi diagram figures (default: true)
% 
% **O**
% The output of this module is voronoi analysis data including:
% - Voronoi diagrams colored by g2 values (saved as images)
% - Volume fraction calculations
% - Voronoi cell area statistics
% - Results saved in voronoi_results.mat
%
% ---
% 
% ## Usage
% This module follows the PeGS2 format of the following main blocks
% 
% - `File management and structure initialization` : house keeping
% - `Loading in data and calculation of Voronoi tessellation`
% - `Area analysis and volume fraction calculation`
% - `Displaying data and saving results` : if verbose
% - `saving parameters`
% 
% Below the main function is a subfunction used to set defaults

function out = voronoi(fileParams, vaParams, verbose)

%% FILE MANAGEMENT
if not(isfolder(fullfile(fileParams.topDir, fileParams.voronoiDir))) %make a new folder for voronoi analysis
    mkdir(fullfile(fileParams.topDir, fileParams.voronoiDir));
end

% Find the solved files to use for Voronoi analysis
files = dir(fullfile(fileParams.topDir, fileParams.solvedDir, '*solved.mat')); %which files are we processing?
nFrames = length(files); %how many files are we processing?
if nFrames == 0
    % Try using particle data if solved data doesn't exist
    files = dir(fullfile(fileParams.topDir, fileParams.particleDir, '*_centers.txt')); 
    nFrames = length(files);
    if nFrames == 0
        error(['No data found in: ', fullfile(fileParams.topDir, fileParams.solvedDir, '*solved.mat'), ' or ', ...
               fullfile(fileParams.topDir, fileParams.particleDir, '*_centers.txt'), ' --check path']);
    else
        data_source = 'centers';
        disp(['Using particle center data from ', num2str(nFrames), ' files for Voronoi analysis']);
    end
else
    data_source = 'solved';
    disp(['Using solved data from ', num2str(nFrames), ' files for Voronoi analysis']);
end

%% DEFAULT PARAMETERS NEEDED TO RUN THIS SCRIPT ARE SET BELOW
vaParams = setupParams(vaParams);

%% PROCESS EACH FRAME
for frame = 1:nFrames
    % Clear variables for this frame
    clearvars particle points vor_vertices vor_cells;
    
    % Load particle data
    if strcmp(data_source, 'solved')
        % Load from solved data
        pres = load(fullfile(files(frame).folder, files(frame).name)); 
        particle = pres.particle;
    else
        % Load from centers file
        centers_data = readmatrix(fullfile(files(frame).folder, files(frame).name));
        % Create basic particle structure from centers data
        for i = 1:size(centers_data, 1)
            particle(i).id = i;
            particle(i).x = centers_data(i, 1);
            particle(i).y = centers_data(i, 2);
            particle(i).r = centers_data(i, 3);
            particle(i).g2 = 0; % Default g2 value if not available
        end
    end
    
    % Extract coordinates for Voronoi analysis
    NN = length(particle);
    x_coords = zeros(NN, 1);
    y_coords = zeros(NN, 1);
    g2_values = zeros(NN, 1);
    radii = zeros(NN, 1);
    
    for i = 1:NN
        x_coords(i) = particle(i).x;
        y_coords(i) = particle(i).y;
        radii(i) = particle(i).r;
        if isfield(particle, 'g2')
            g2_values(i) = particle(i).g2;
        end
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
        if area <= vaParams.area_threshold * average_area
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
    
    % Filter out Voronoi cells that contain an acute angle if requested
    if vaParams.acute_angle_filter
        no_acute_cell = [];
        for i = 1:length(valid_indices)
            index = valid_indices(i);
            vertices = vor_vertices(vor_cells{index}, :);
            if ~has_acute_angle(vertices, vaParams.acute_angle_threshold)
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
    else
        no_acute_cell = valid_indices;
        volume_fraction_uniform_no_acute = volume_fraction_uniform;
    end
    
    % Generate filename for this frame
    if isfield(fileParams, 'frameIdInd')
        current_file = files(frame).name;
        frameid = str2double(current_file(fileParams.frameIdInd:fileParams.frameIdInd+3));
    else
        frameid = frame;
    end
    
    % Create base name for output files
    if strcmp(data_source, 'solved')
        base_name = strrep(files(frame).name, '_solved.mat', '');
    else
        base_name = strrep(files(frame).name, '_centers.txt', '');
    end
    
    % Save results for this frame
    result_data = struct('volume_fraction', volume_fraction_uniform, ...
                         'volume_fraction_no_acute', volume_fraction_uniform_no_acute, ...
                         'num_particles', length(valid_indices), ...
                         'num_particles_no_acute', length(no_acute_cell), ...
                         'radius_uniform', radius_uniform, ...
                         'average_area', average_area);
     
    % Save a single result file with all necessary information
    save(fullfile(fileParams.topDir, fileParams.voronoiDir, [base_name, '_voronoi.mat']), 'result_data', 'valid_indices', 'no_acute_cell', 'vor_vertices', 'vor_cells');
    
    % Plot the Voronoi diagrams if verbose or save_figures is true
    if verbose || vaParams.save_figures
        % Create only one final version of the Voronoi diagram (with acute angles removed)
        h = figure('Visible', verbose);
        
        % Directly draw the Voronoi diagram with acute angles removed
        plot_voronoi_diagram(vor_vertices, vor_cells, no_acute_cell, g2_values, average_area, vaParams.area_threshold, volume_fraction_uniform_no_acute, x_coords, y_coords);
        title(['Voronoi Analysis - Volume Fraction: ' num2str(volume_fraction_uniform_no_acute * 100, '%.2f') '%']);
        
        % Save the image
        if vaParams.save_figures
            saveas(h, fullfile(fileParams.topDir, fileParams.voronoiDir, [base_name, '_voronoi.png']));
        end
    end
    
    % Display results for this frame
    if verbose
        fprintf('Frame %d analysis:\n', frameid);
        fprintf('  Number of particles: %d\n', length(valid_indices));
        if vaParams.acute_angle_filter
            fprintf('  Number of particles without acute angles: %d\n', length(no_acute_cell));
        end
        fprintf('  Volume fraction: %.2f%%\n', volume_fraction_uniform * 100);
        if vaParams.acute_angle_filter
            fprintf('  Volume fraction (no acute): %.2f%%\n', volume_fraction_uniform_no_acute * 100);
        end
        fprintf('  Average Voronoi cell area: %.2f pixels²\n', average_area);
        fprintf('  Uniform particle radius: %.2f pixels\n', radius_uniform);
    end
end

% Create a summary file with results from all frames
summarize_voronoi_results(fileParams, nFrames, verbose);

% Save parameters used in this run
fields = fieldnames(vaParams);
for i = 1:length(fields)
    fileParams.(fields{i}) = vaParams.(fields{i});
end

fileParams.time = datetime("now");
fields = fieldnames(fileParams);
C = struct2cell(fileParams);
vaParamsCell = [fields C];

writecell(vaParamsCell, fullfile(fileParams.topDir, fileParams.voronoiDir, 'voronoi_params.txt'), 'Delimiter', 'tab');

if verbose
    disp('Completed voronoi()');
end

out = true;
end

function result = has_acute_angle(polygon, threshold_degrees)
    % Checks if a polygon contains an angle below the threshold (in degrees)
    threshold_radians = threshold_degrees * (pi/180);
    n = size(polygon, 1);
    result = false;
    
    for i = 1:n
        p1 = polygon(mod(i-2, n)+1, :);
        p2 = polygon(i, :);
        p3 = polygon(mod(i, n)+1, :);
        
        v1 = p1 - p2;
        v2 = p3 - p2;
        
        cosine_angle = dot(v1, v2) / (norm(v1) * norm(v2));
        angle = acos(min(max(cosine_angle, -1), 1)); % Ensure value is within valid range
        
        if angle < threshold_radians
            result = true;
            return;
        end
    end
end

function plot_voronoi_diagram(vor_vertices, vor_cells, indices, g2_values, avg_area, area_threshold, volume_fraction, x_coords, y_coords)
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
    text(min(x_coords) + 50, max(y_coords) - 100, ['Volume fraction: ' num2str(volume_fraction * 100, '%.2f') '%'], 'FontSize', 12);
    
    % Add colorbar
    c = colorbar;
    ylabel(c, 'g^2 value', 'Rotation', 270, 'FontSize', 12);
    c.Label.Position(1) = 3.5;
    
    axis equal;
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

function summarize_voronoi_results(fileParams, nFrames, verbose)
     % Create a summary of all voronoi results
     voronoi_files = dir(fullfile(fileParams.topDir, fileParams.voronoiDir, '*_voronoi.mat'));
     
     if length(voronoi_files) ~= nFrames
         warning('Number of voronoi results files does not match number of processed frames');
     end
     
     % Initialize arrays to store frame-by-frame results
     frame_ids = zeros(length(voronoi_files), 1);
     volume_fractions = zeros(length(voronoi_files), 1);
     volume_fractions_no_acute = zeros(length(voronoi_files), 1);
     num_particles = zeros(length(voronoi_files), 1);
     num_particles_no_acute = zeros(length(voronoi_files), 1);
     
     % Process each voronoi result file
     for i = 1:length(voronoi_files)
         % Load the file
         data = load(fullfile(voronoi_files(i).folder, voronoi_files(i).name));
         
         % Extract frame ID from filename if possible
         try
             if isfield(fileParams, 'frameIdInd')
                 frame_ids(i) = str2double(voronoi_files(i).name(fileParams.frameIdInd:fileParams.frameIdInd+3));
             else
                 frame_ids(i) = i;
             end
         catch
             frame_ids(i) = i;
         end
         
         % Store results
         volume_fractions(i) = data.result_data.volume_fraction;
         volume_fractions_no_acute(i) = data.result_data.volume_fraction_no_acute;
         num_particles(i) = data.result_data.num_particles;
         num_particles_no_acute(i) = data.result_data.num_particles_no_acute;
     end
     
     % Create summary structure
     summary = struct();
     summary.frame_ids = frame_ids;
     summary.volume_fractions = volume_fractions;
     summary.volume_fractions_no_acute = volume_fractions_no_acute;
     summary.num_particles = num_particles;
     summary.num_particles_no_acute = num_particles_no_acute;
     summary.mean_volume_fraction = mean(volume_fractions);
     summary.mean_volume_fraction_no_acute = mean(volume_fractions_no_acute);
     
     % Save only one summary data file
     save(fullfile(fileParams.topDir, fileParams.voronoiDir, 'voronoi_summary.mat'), 'summary');
     
     % Display summary if verbose
     if verbose
         fprintf('\nSummary of voronoi analysis across %d frames:\n', length(voronoi_files));
         fprintf('  Average volume fraction: %.2f%%\n', summary.mean_volume_fraction_no_acute * 100);
         fprintf('  Average number of particles: %.1f\n', mean(num_particles_no_acute));
         fprintf('  Results saved to voronoi_summary.mat\n');
     end
end

function params = setupParams(params)
    % Set default parameters if not provided by the user
    
    if isfield(params, 'area_threshold') == 0
        params.area_threshold = 1; % Default area threshold multiplier
    end
    
    if isfield(params, 'acute_angle_filter') == 0
        params.acute_angle_filter = true; % Filter out cells with acute angles by default
    end
    
    if isfield(params, 'acute_angle_threshold') == 0
        params.acute_angle_threshold = 60; % Default acute angle threshold in degrees
    end
    
    if isfield(params, 'save_figures') == 0
        params.save_figures = true; % Save figures by default
    end
end 