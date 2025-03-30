%Updated to first release version of PeGS2 by Carmen Lee 29/9/24
% Written for PeGS 2.0 by Kerstin Nordstrom 8/24 
%Adapted from Carmen Lee's adaptation of Jonathan Kollmer's PeGS 1.0
%

% # particleDetect
% 
% 
% **i/o**
% The input of this module are experimental images. It is assumed that the image is RGB, the particles are detected using the red channel and the green channel contains the force information. 
% 
% **i** 
% N input images
% 
% 
% **o** 
% The output are N txt files in the format:
% *Centers text files: X Y radius edge
% *Params text file: parameters used in circle finding, et al
% 
% **parameters**
% The user set parameters are held inside of the structure <pdParams> in the main function (PeGSModular) and as <p> in the particleDetect function. The fields self populate if not set by the user with default values (done in the function paramsSetUp). The parameters are
% - `boundaryType` : shape of the boundary holding the particles, this is for identifying edge particles
% - `radiusRange`: an array consisting of the lower and upper radius bounds for the particles (pixels)
% - `dtol` : the distance tolerance (pixels) between the edge of the particle and the wall to be assigned as as edge
% - `sensitivity` : the sensitivity to for circle finding. Higher sensitivity == find more circles
% - `edgeThresh` : to assign what is the edge of an object in the image (note: not the same edge as the edge flag, and this is set inside the main function after the image has been loaded)
% - `imageType` : automatically detected image type (e.g., "test" or "step09") to adjust parameters accordingly
% 
% **notes**
% If the polariscope is set up for being in transmission rather than reflection, imfindcircles might work better if you set 'ObjectPolarity' to 'dark' rather than 'bright'
% 


function out = particleDetect(fileParams, pdParams, verbose)
%fileParams = directories and image name pattern
%pdParams = parameters for particle detect

%% FILE MANAGEMENT


    if ~exist(fullfile(fileParams.topDir, fileParams.particleDir) , 'dir')
        mkdir(fullfile(fileParams.topDir, fileParams.particleDir))
    end

    
    if verbose
        disp('starting particleDetect() to find all particle centroids and save results in particleDir')
    end

%% Set up for parameters is found in the function below
    pdParams = paramsSetUp(pdParams); %see below for defaults for params




%% IMPORT DATA & BEGIN ANALYSIS

    images = dir(fullfile(fileParams.topDir, fileParams.imgDir,fileParams.imgReg));
    nFrames = length(images);
    

    for frame = 1:nFrames

        im = imread(fullfile(images(frame).folder,images(frame).name)); %read image
        
        %% Check image type for parameter adjustment
        currentImage = images(frame).name;
        currentParams = pdParams;
        
        % Check if current image is test.jpg or similar
        if contains(lower(currentImage), 'test')
            % Parameters for test.jpg type images (larger particles)
            currentParams.imageType = "test";
            currentParams.radiusRange = [110 160]; % Larger particles in test.jpg
            currentParams.sensitivity = 0.95; % May need different sensitivity
            
            if verbose
                disp(['Detected test.jpg image type - using larger radius range: [', ...
                     num2str(currentParams.radiusRange(1)), ' ', ...
                     num2str(currentParams.radiusRange(2)), ']']);
            end
        else
            % Default parameters for Step09 series
            currentParams.imageType = "step09";
            % Keep existing parameters from pdParams
        end
        
        if verbose
            disp(['Processing image: ', currentImage, ' with image type: ', currentParams.imageType]);
        end
        
        %get red channel for particle detection
        red = im(:,:,1); %particles
        green = im(:,:,2); %photoelastic signal
        red = imsubtract(red, green*0.05); %some green bleeds through, makes sharper red

        %circle detection
        [centers, radii, ~] = imfindcircles(red,currentParams.radiusRange,'ObjectPolarity','bright','Method','TwoStage','Sensitivity',currentParams.sensitivity);%, 'EdgeThreshold',p.edgeThresh);

        %classify edge particles
        if currentParams.boundaryType == "rectangle"
            
            lpos = min(centers(:,1)-radii);
            rpos = max(centers(:,1)+radii);
            upos = max(centers(:,2)+radii);
            bpos = min(centers(:,2)-radii);
            lwi = centers(:,1)-radii <= lpos+currentParams.dtol;
            rwi = centers(:,1)+radii >= rpos-currentParams.dtol;
            uwi = centers(:,2)+radii >= upos-currentParams.dtol;
            bwi = centers(:,2)-radii <= bpos+currentParams.dtol; %need to add edge case of corner particle

            edges = zeros(length(radii), 1);
            edges(rwi) = 1; %right
            edges(lwi) = -1; %left
            edges(uwi) = 2;  %"upper" - as vertical pixels are backwards from cartesian, actually bottom of image
            edges(bwi) = -2; %"bottom" - see above comment
            %interior particles are 0
        end %for edge detection
        
        %display circles
        if verbose == true
            imshow(red);
            viscircles(centers, radii);
            hold on
            str=string(edges);
            text(centers(:,1), centers(:,2),str,'Color','red','FontSize',14)
            drawnow;
            hold off
        end


        %make matrix of positions x ,y, radii, edge classification
        particle = [centers(:,1), centers(:,2), radii, edges];

        %save to text file

        txtfilename=[images(frame).name(1:end-4), '_centers.txt'];
        writematrix(particle, fullfile(fileParams.topDir, fileParams.particleDir,txtfilename), 'delimiter',',')

    end %end for loop over images




%save sample image if verbose
    if verbose == true
        imfilename=['Centers_',images(frame).name];
        saveas(gcf,fullfile(fileParams.topDir, fileParams.particleDir,imfilename))
    end

%% saving parameters in and finishing module

fields = fieldnames(pdParams);
for i = 1:length(fields)
    fileParams.(fields{i}) = pdParams.(fields{i});
end

fileParams.lastimagename=images(frame).name;
fileParams.time = datetime("now");
fields = fieldnames(fileParams);
C=struct2cell(fileParams);
pdParams = [fields C];
writecell(pdParams,fullfile(fileParams.topDir, fileParams.particleDir,'particleDetect_params.txt'),'Delimiter','tab')


if verbose 
        disp('done with particleDetect()');
end

out = true;

end

%% defaults for particleDetectModule
function p = paramsSetUp(p)
%classify boundary type
if isfield(p,'boundaryType') == 0
    p.boundaryType = "rectangle";
end

%set radius range
if isfield(p,'radiusRange') == 0
    p.radiusRange = [45 80];  % Default radius range for Step09 series
end

%classify edge particles with tolerance
if isfield(p,'dtol') == 0
    p.dtol = 10;
end


%sensitivity for Hough
if isfield(p,'sensitivity') == 0
    p.sensitivity = 0.945;
end


if isfield(p,'edgeThresh') == 0
    p.edgeThresh = 0.02;
end

% Set default image type
if isfield(p,'imageType') == 0
    p.imageType = "step09";  % Default image type
end
end
