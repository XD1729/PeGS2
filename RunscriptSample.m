  %% a sample runscript showing all of the module parameters

%% file parameters
fileParams.topDir = './testdata/'; % where the top directory of the data will be stored
fileParams.imgDir = 'images'; % where the images are saved
fileParams.particleDir = 'particles'; %output directory for particle information
fileParams.contactDir = 'contacts'; % output directory for contact information
fileParams.solvedDir = 'solved'; % output directory for solved force information
fileParams.adjacencyDir = 'adjacency'; % adjacency list directory
fileParams.imgReg = '*.jpg'; %image format and regex
fileParams.frameIdInd = 8; %the index in the file names where sequential numbering starts, optional, remove if unwanted or if naming conventions are irregular

%% verbose
verbose = true; % do you want plots?
 
%% particleDetect parameters

pdParams.boundaryType = "rectangle";
pdParams.radiusRange = [45 80]; %range in pixels to look for particles
pdParams.dtol = 10; %classify edge particles with tolerance
pdParams.sensitivity = 0.945; %sensitivity for Hough
pdParams.edgeThresh = 0.02; %senstitivity for Hough

%% particleTrack parameters
ptParams.skipValue = 20; % padding for data array, I recommend ~10% of packing size, but it will give warning if it needs to be larger

%% contactDetect parameters

cdParams.metersperpixel = .007/160; %meters/pixel
cdParams.fsigma = 140;  %photoelastic stress coefficient
cdParams.g2cal = 100; %calibration Value for the g^2 method, can be computed by joG2cal.m (PEGS 1.0 version)
cdParams.dtol = 10; %how far away can the outlines of 2 particles be to still be considered Neighbors
cdParams.contactG2Threshold = 0.5; %sum of g2 in a contact area larger than this determines a valid contact
cdParams.CR = 10; %contact radius over which contact gradient is calculated
cdParams.imadjust_limits = [0,0.65]; %adjust contrast in green channel
cdParams.rednormal = 2; %fractional amount to subtract the red channel from the green channel (Rimg/rednormal)
cdParams.figverbose = true; %show figures

%% diskSolve parameters

dsParams.algorithm = 'levenberg-marquardt'; %Algorithm to use. Other options: 'trust-region-reflective', 'interior-point'
dsParams.maxIterations = 200; % maximum iterations allowed for each fit
dsParams.maxFunctionEvaluations = 400; %maximum number of function evaluations for each fit
dsParams.functionTolerance = 0.01; % if subsequent fits change by less than this amount, stop
dsParams.scaling = 0.5; %scale image by this much
dsParams.maskradius = 0.96; %percent of disk to use for fit (cuts outer edge)
dsParams.original = 1; % Set 1Run original unvectorised version of disk solver
dsParams.vectorise = 0; %Run vectorised version of disk solver (coming soon)

%% adjacencyMatrix parameters


amParams.go = true; %compile adjaceceny lists for each image, set to false if you have already done this and only want to have the big list
amParams.fmin = 0.000001;%minimum force (in Newtons) to consider a contact a valid contact
amParams.fmax = 1000; %maximum force (in Newtons) to consider a contact a valid contact
amParams.emax = 2800; %maximum fit error/residual to consider a contact a valid contact
amParams.skipvalue = 20; %I chose this as a result of my system size, could and should be altered based on your specific system and variability in finding particles


%% compile module parameters into their own structure to pass into PeGSModular

moduleParams.pdParams = pdParams;
moduleParams.ptParams = ptParams;
moduleParams.cdParams = cdParams;
moduleParams.dsParams = dsParams;
moduleParams.amParams = amParams;

%% run
PeGSModular(fileParams, moduleParams, verbose)