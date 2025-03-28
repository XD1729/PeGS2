% Updated by Xie Dong 27/03/2025
% Adapted from Katarzyna Maksymiuk's adaptation of Jonathan Kollmer's PeGS 1.0
% Extracted from PeGSDiskPrep_static.m
%
% # imagePreprocess
% 
% This module performs lighting correction and contrast enhancement on photoelastic images.
% It evens out the lighting in the image to provide better quality images for analysis.
%
% **i/o**
% The input of this module are experimental images. It processes RGB images to prepare them
% for better particle detection and force estimation in subsequent modules.
% 
% **i** 
% N input images
% 
% **o** 
% The output are N processed images saved as:
% * Enhanced images: with even lighting and improved contrast
% * Parameters file: settings used in preprocessing
% 
% **parameters**
% The user set parameters are held inside of the structure <ipParams> in the main function (PeGSModular)
% and as <p> in the imagePreprocess function. Parameters include:
% - `light_correction_coefficients`: coefficients for correcting lighting based on distance from center [a, b, c]
% - `I_low`, `I_high`: intensity adjustment limits for contrast enhancement
% - `auto_detect`: automatically detect image type and apply appropriate parameters
% - `obstacle_mask`: apply obstacle masking to blank out specific areas
% - `region_specific`: apply region-specific contrast enhancement
%
% **notes**
% This module is designed to be inserted before particleDetect in the PeGS workflow.
% The lighting correction is particularly useful for uneven illumination.


function out = imagePreprocess(fileParams, ipParams, verbose)
%fileParams = directories and image name pattern
%ipParams = parameters for image preprocessing

%% FILE MANAGEMENT
% Set up output directory for preprocessed images
if isfield(fileParams,'processedImgDir') == 0
    fileParams.processedImgDir = 'processed_images';
end

% Create output directory if it doesn't exist
processedImgDirPath = fullfile(fileParams.topDir, fileParams.processedImgDir);
if ~exist(processedImgDirPath, 'dir')
    mkdir(processedImgDirPath);
end

% Get image files
if isfield(fileParams,'imgReg') == 0
    fileParams.imgReg = '*.jpg'; %image format and regex
end

imgPath = fullfile(fileParams.topDir, fileParams.imgDir);
imgFiles = dir(fullfile(imgPath, fileParams.imgReg));
nFrames = length(imgFiles);

if nFrames == 0
    error(['No image files found in ', imgPath, ' with pattern ', fileParams.imgReg]);
end

% Get parameters if not set
p = paramsSetUp(ipParams);

% Save the parameters used
paramFile = fullfile(processedImgDirPath, 'imagePreprocess_params.txt');
paramFileID = fopen(paramFile, 'w');
fprintf(paramFileID, 'light_correction_coefficients: [%g, %g, %g]\n', p.light_correction_coefficients);
fprintf(paramFileID, 'I_low: %g\n', p.I_low);
fprintf(paramFileID, 'I_high: %g\n', p.I_high);
fprintf(paramFileID, 'auto_detect: %d\n', p.auto_detect);
fprintf(paramFileID, 'obstacle_mask: %d\n', p.obstacle_mask);
fprintf(paramFileID, 'region_specific: %d\n', p.region_specific);
fclose(paramFileID);

if verbose
    disp('Starting image preprocessing...');
    disp(['Processing ', num2str(nFrames), ' images']);
end

%% MAIN PROCESSING LOOP
for frame = 1:nFrames
    % Read the image
    imageFile = fullfile(imgPath, imgFiles(frame).name);
    img = imread(imageFile);
    info = imfinfo(imageFile);
    Width = double(info.Width);
    Height = double(info.Height);
    
    % Auto-detect image type and set appropriate parameters if enabled
    if p.auto_detect
        [img_type, img_params] = detectImageType(img, imgFiles(frame).name);
        
        % Update parameters based on detected image type
        if strcmp(img_type, 'test')
            p.light_correction_coefficients = img_params.light_correction_coefficients;
            p.I_low = img_params.I_low;
            p.I_high = img_params.I_high;
            p.obstacle_mask = img_params.obstacle_mask;
            p.region_specific = img_params.region_specific;
        elseif strcmp(img_type, 'step')
            p.light_correction_coefficients = img_params.light_correction_coefficients;
            p.I_low = img_params.I_low;
            p.I_high = img_params.I_high;
            p.obstacle_mask = img_params.obstacle_mask;
            p.region_specific = img_params.region_specific;
        end
        
        if verbose
            disp(['Detected image type: ', img_type]);
        end
    end
    
    if verbose
        disp(['Processing image ', num2str(frame), ' of ', num2str(nFrames), ': ', imgFiles(frame).name]);
    end
    
    %% STEP 1: Convert to grayscale
    gray_img = rgb2gray(img);
    gray_img = im2double(gray_img);
    
    %% STEP 2: Evening out the lighting
    % Calculate distance from center for each pixel
    coefficients = p.light_correction_coefficients;
    
    % Skip lighting correction if coefficients are all zero
    if sum(abs(coefficients)) > 0
        coefficients_corr = [-coefficients(1) coefficients(2) 0];
        centre = [0.5*Height 0.5*Width];
        pixel_y = linspace(1, Height, Height);
        pixel_x = linspace(1, Width, Width);
        
        dist_centre = zeros(Height, Width);
        for i = 1:Width
            for n = 1:Height
                dist_centre(n, i) = sqrt((pixel_x(i) - centre(2))^2 + (pixel_y(n) - centre(1))^2);
            end
        end
        
        % Correct intensities in grayscale image
        I_corr_photo = polyval(coefficients_corr, dist_centre);
        corrected_img = I_corr_photo + gray_img;
        
        % Ensure values are in valid range [0, 1]
        corrected_img = max(0, min(1, corrected_img));
    else
        % Skip lighting correction - use grayscale directly
        corrected_img = gray_img;
    end
    
    %% STEP 3: Apply obstacle mask if needed
    if p.obstacle_mask
        % Define obstacle region
        obstacle_xs = [Width/2-550 Width/2+100 Width/2+750];
        obstacle_ys = [Height Height-500 Height];
        obstacle_blank = roipoly(corrected_img, obstacle_xs, obstacle_ys);
        
        % Apply mask
        corrected_img = double(~obstacle_blank).*corrected_img;
        
        % Blank out bottom region (similar to original PeGSDiskPrep)
        if strcmp(img_type, 'test')
            corrected_img(Height-250:Height, :) = 0;
        end
    end
    
    %% STEP 4: Contrast enhancement (global and region-specific)
    if p.region_specific && strcmp(img_type, 'test')
        % Global adjustment first
        photo_corr_3 = imadjust(corrected_img, [p.I_low p.I_high]);
        
        % Define parameters for region-specific adjustments
        I_low_tl = 0.09; I_high_tl = 1;    % Top-left corner
        I_low_tr = 0.09; I_high_tr = 1;    % Top-right corner
        I_low_el = 0;    I_high_el = 1;    % Left edge
        I_low_er = 0;    I_high_er = 1;    % Right edge
        
        % Region 1: Top-left corner adjustment
        xs_tl = [1 650 1300 650];
        ys_tl = [Height/5 Height/5-600 Height/5 Height/5+600];
        tl = roipoly(photo_corr_3, xs_tl, ys_tl);
        photo_corr_4 = double(tl).*corrected_img;
        photo_corr_4 = imadjust(photo_corr_4, [I_low_tl I_high_tl]);
        photo_corr_5 = double(~tl).*photo_corr_3;
        photo_corr_5 = imadd(photo_corr_4, photo_corr_5);
        
        % Region 2: Top-right corner adjustment
        xs_tr = [Width-1500 Width-800 Width Width-650];
        ys_tr = [Height/5-200 Height/5-1000 Height/5-170 Height/5+900];
        tr = roipoly(photo_corr_5, xs_tr, ys_tr);
        photo_corr_6 = double(tr).*corrected_img;
        photo_corr_6 = imadjust(photo_corr_6, [I_low_tr I_high_tr]);
        photo_corr_7 = double(~tr).*photo_corr_5;
        photo_corr_7 = imadd(photo_corr_7, photo_corr_6);
        
        % Region 3: Left edge adjustment
        xs_l_edge = [600 1000 Width/3+600 Width/3];
        ys_l_edge = [Height/3+600 Height/3+220 Height-200 Height];
        l_edge = roipoly(photo_corr_7, xs_l_edge, ys_l_edge);
        photo_corr_8 = double(l_edge).*photo_corr_7;
        photo_corr_8 = imadjust(photo_corr_8, [I_low_el I_high_el]);
        photo_corr_9 = double(~l_edge).*photo_corr_7;
        photo_corr_9 = imadd(photo_corr_8, photo_corr_9);
        
        % Region 4: Right edge adjustment
        xs_r_edge = [Width*2/3-1000 Width-1000 Width Width*2/3];
        ys_r_edge = [Height Height/3 Height/3 Height];
        r_edge = roipoly(photo_corr_9, xs_r_edge, ys_r_edge);
        photo_corr_10 = double(r_edge).*photo_corr_9;
        photo_corr_10 = imadjust(photo_corr_10, [I_low_er I_high_er]);
        photo_corr_11 = double(~r_edge).*photo_corr_9;
        enhanced_img = imadd(photo_corr_10, photo_corr_11);
        
        % Final brightness adjustment
        enhanced_img = imadjust(enhanced_img);
    else
        % Standard global contrast enhancement
        enhanced_img = imadjust(corrected_img, [p.I_low p.I_high]);
    end
    
    %% Save processed images
    [~, baseName, ext] = fileparts(imgFiles(frame).name);
    
    % Save the enhanced image
    enhancedImgFile = fullfile(processedImgDirPath, [baseName, '_enhanced', ext]);
    imwrite(enhanced_img, enhancedImgFile);
    
    % Display results if verbose
    if verbose
        figure(1);
        imshow(gray_img);
        title('Original Grayscale Image');
        
        figure(2);
        imshow(corrected_img);
        title('Light-Corrected Image');
        
        figure(3);
        imshow(enhanced_img);
        title('Enhanced Image');
        
        drawnow;
    end
end

if verbose
    disp('Image preprocessing complete.');
end

out = true;
end

%% Image type detection function
function [img_type, img_params] = detectImageType(img, filename)
    % Detect image type based on filename and image characteristics
    % Returns the image type and appropriate parameters
    
    img_type = 'unknown';
    
    % Check filename for classification hints
    if contains(lower(filename), 'test')
        img_type = 'test';
    elseif contains(lower(filename), 'step')
        img_type = 'step';
    end
    
    % If filename doesn't provide a clear indication, analyze image properties
    if strcmp(img_type, 'unknown')
        % Convert to grayscale for analysis
        gray_img = rgb2gray(img);
        
        % Calculate average brightness and standard deviation
        avg_brightness = mean(gray_img(:));
        std_brightness = std(double(gray_img(:)));
        
        % Check RGB channel distribution for typical photoelastic patterns
        red_channel = img(:,:,1);
        green_channel = img(:,:,2);
        
        red_mean = mean(red_channel(:));
        green_mean = mean(green_channel(:));
        
        % Feature-based classification
        if red_mean > 1.5 * green_mean
            img_type = 'test'; % Typical for test.JPG where red channel is much stronger
        else
            img_type = 'step'; % Default to step type
        end
    end
    
    % Set appropriate parameters based on image type
    if strcmp(img_type, 'test')
        % Parameters optimized for test.JPG type
        img_params.light_correction_coefficients = [-9.268015495383530e-09, 1.192483586444351e-06, 0.076128557100614];
        img_params.I_low = 0.09;
        img_params.I_high = 1;
        img_params.obstacle_mask = true;
        img_params.region_specific = true;
    elseif strcmp(img_type, 'step')
        % Parameters optimized for Step09 series
        img_params.light_correction_coefficients = [0, 0, 0]; % Minimal or no lighting correction
        img_params.I_low = 0.05;
        img_params.I_high = 0.95;
        img_params.obstacle_mask = false;
        img_params.region_specific = false;
    else
        % Default parameters for unknown image types
        img_params.light_correction_coefficients = [0, 0, 0];
        img_params.I_low = 0.05;
        img_params.I_high = 0.95;
        img_params.obstacle_mask = false;
        img_params.region_specific = false;
    end
end

%% Parameter default function
function p = paramsSetUp(ipParams)
    % Set default parameters if not provided by user
    p = ipParams;
    
    % Enable auto-detect by default
    if ~isfield(p, 'auto_detect')
        p.auto_detect = true;
    end
    
    % Lighting correction coefficients
    if ~isfield(p, 'light_correction_coefficients')
        p.light_correction_coefficients = [-9.268015495383530e-09, 1.192483586444351e-06, 0.076128557100614];
    end
    
    % Global intensity adjustment
    if ~isfield(p, 'I_low')
        p.I_low = 0.09;
    end
    if ~isfield(p, 'I_high')
        p.I_high = 1;
    end
    
    % Obstacle masking (for blocking out areas near edges/obstacle)
    if ~isfield(p, 'obstacle_mask')
        p.obstacle_mask = true; % Default to true
    end
    
    % Region-specific processing (for test.JPG)
    if ~isfield(p, 'region_specific')
        p.region_specific = true; % Default to true
    end
end 
