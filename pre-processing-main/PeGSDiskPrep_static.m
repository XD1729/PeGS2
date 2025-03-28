% Particle Detector, Neighbour finder and Contact validator that creates input files for J. E. Kollmer's Photoelastic Disk Solver
% Particle Detection and Neigbour Finding Adapted from J. E. Kollmer's Earlier Script (joCentersNonMonodisperse.m as of 2016/05/03)
% Photoelastic Disk Solver inspired from peDiskSolve by James Puckett (Phd-Thesis 2012) http://nile.physics.ncsu.edu

% If you use this please cite the follwoing paper
% K.E. Daniels, J. E. Kollmer & J. G. Puckett, "Photoelastic force measurements in granular materials", Rev. Sci. Inst. (2017)
% DOI: 10.1063/1.4983049

% edited on 2018/08/09 by Joshua Miller (jsmille9@ncsu.edu)
%           2022/12/25 by Katarzyna Maksymiuk
%           2023 by Xie Dong

close all
clear all
clc


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           User defined values                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Image information
verbose = true; %Generates lots of plots showing results
directory = './';
% files = dir([directory, '2.JPG']); %Which files are we processing?

files = dir([directory, 'IMG_0078.JPG']);
%files = dir('Centers*0001.txt'); %Alternatively, centers files can be loaded. This requires that both particle detections be flagged false however.
nFrames = length(files); %How many files are we processing ?
light_correction_coefficients = [-9.268015495383530e-09,1.192483586444351e-06,0.076128557100614];  %(file PP_light_correction.mlx can be used to obtain these coefficients)

%% Uncomment the code below if the diameter of the circle is unknown
% imageFile = [directory,files(1).name];
% img = imread(imageFile);
% figure(101), imshow(img(2000:3000,2500:3500))  %Show excerpt of the image
% d = drawline;  %Draw a line across a particle's diameter
% pos = d.Position;
% diffPos = diff(pos);
% radius_est = hypot(diffPos(1), diffPos(2))/2  %Estimate the radius [px]
radius_desired = 122;  %[px]
pxPerMeter = (14.8e-3)/radius_desired;  %[m/px]

%% Hough Transform Values
doParticleDetectionH = true; %Detect particles using Hough Transform?

RlargeH = [110 160]; %115 150What radius (in pixels) range do we expect for the large discs?
SL = 0.97; %Sensitivity of the Hough Transform large disc detetcor, exact value is Voodo magic...
RsmallH = [45 55]; %What radius (in pixels) range do we expect for the small discs?可以有大小两种圆盘
SS = 0.0; %Sensitivity of the Hough Transform small disc detetcor, set to zero in monodisperse system

HoughDebug = false; %Debugs Hough Sensitivities so particles are found "better" -- does not really work well here, finds overlapping circles
DS = 0.0025; % How much should we adjust sensitivity if wrong number of particles are found
NsmallH = 0; %Number of small discs. Only used in Hough Debug.
NlargeH = 190; %Number of large discs. Only used in Hough Debug.

%% Convolution Method Values -- unused in this project
doParticleDetectionC = false; %Detect particles using convolution method?
ConvDebug = false;

RlargeC = 155; %What radius (in pixels) do we expect for the large discs?
RsmallC = 47; %What radius (in pixels) do we expect for the small discs? 
%Note: The above can be input in a range ONLY if the ConvDebug is set to
%Note: true. Otherwise, the program needs the radius that works.
NsmallC = 0; %Number of small discs. Needed for Convolution.
NlargeC = 33; %Number of large discs. Needed for Convolution.

%% Neighbour Finding Values
findNeighbours = true; 

fsigma = 2416; %Photoelastic stress coefficient [Br] (can be derived in CV_stress_optic_coefficient.mlx)
g2cal = 0.4890; %Calibration Value for the g^2 method [1/N] (can be derived in CV_Force_G2_calibration.mlx)
dtol = 25; %How far away can the outlines of 2 particles be to still be considered Neighbours [px]

% contactG2Threshold = 0.0036; %Sum of g2 in a contact area larger than this determines a valid contact  <-- choice to customise this in lines 455-545
contactG2Threshold = 0.0015;
CR = 30; %Radius around a contactact point that is checked for contact validation [px]

%% Values for intensity adjustment
% Set values for the whole image (best not go below the maximum intensity in the corrected photo for the upper bound)
I_low_all = 0.09; 
I_high_all = 1;

% Set values for adjusting the top left and right corners of the image, if not needed set to the same values as the whole image
I_low_tl = 0.09; 
I_high_tl = 1;

I_low_tr = 0.09; 
I_high_tr = 1;

% Set values for adjusting the left and right edges of the image (near the walls), if not needed set to 0 and 1
I_low_el = 0; 
I_high_el = 1;


I_low_er = 0; 

I_high_er = 1;

%% The rheomether dimensions
% (1) The right wall
% %Uncomment to establish points of reference the position of system boundary on the right (need to comment lines 101-103 then)
%     imageFile = [directory,files(1).name];
%     img = imread(imageFile);
%     figure(102), imshow(img)
%     line_right = drawline;  %Draw along the right boundary, towards the top of image
%     pos_right = line_right.Position
%     slope_right = (pos_right(2,2) - pos_right(1,2))/(pos_right(2,1) - pos_right(1,1))
%     b_right = pos_right(1,2) - slope_right*pos_right(1,1)
%     x_right = pos_right(1,1):1:pos_right(2,1);
    slope_right = -1.062027231467473;
    b_right = 7.223605453733785e+03;
    x_right = 3276:1:5616.5;
    for i = 1:length(x_right)
        y_right(i) = slope_right.*x_right(i) + b_right;
    end
    r_wall_lim = 1.6*radius_desired;   %Distance which limits particle recognition at the wall
% (2) The left wall
% % Uncomment to establish points of reference for the position of system boundary on the left (need to comment lines 118-120 then)
%     imageFile = [directory,files(1).name];
%     img = imread(imageFile);
%     figure(102), imshow(img)
%     line_left = drawline;  %Draw along the left boundary, towards the top of image
%     pos_left = line_left.Position
%     slope_left = (pos_left(2,2) - pos_left(1,2))/(pos_left(2,1) - pos_left(1,1))
%     b_left = pos_left(1,2) - slope_left*pos_left(1,1)
%     x_left = pos_left(2,1):1:pos_left(1,1);
    slope_left = 0.907699665231946;
    b_left = 1.503783855085417e+03;
    x_left = 1:1:2469;
    for i = 1:length(x_left)
        y_left(i) = slope_left.*x_left(i) + b_left;
    end
    l_wall_lim = 1.6*radius_desired;   %Distance which limits particle recognition at the wall
 % (3) The upper, curved wall
    r_rheo = (0.5-0.04)/pxPerMeter;  
    x_cent_rheo = 2808 + 120;
    y_cent_rheo = 3744 + 240;
    th = 0:pi/100:pi;
    x_wall_up = zeros(1,100);
    y_wall_up = zeros(1,100);
    for i=1:100
        x_wall_up(i) = -r_rheo*cos(th(i))+(x_cent_rheo);
        y_wall_up(i) = -r_rheo*sin(th(i))+(y_cent_rheo);
    end
    up_wall_lim = 1.6*radius_desired;   %Distance which limits particle recognition at the wall
 % (4) The area where the force could have been applied
    x_area = [2100 3700];
    y_area = [2850 3744];
    
% Uncomment below to see the positions of rheometer dimensions derived above
%    imageFile = [directory,files(1).name];
%    img = imread(imageFile);
%    figure(103), imshow(img)
%    hold on 
%  %(1) The right wall
%    plot(x_right, y_right, 'b', 'LineWidth', 4)
%    xs_lim_r = [(min(x_right)-r_wall_lim); max(x_right); max(x_right)];
%    ys_lim_r = [max(y_right); min(y_right)-r_wall_lim; max(y_right)];
%    plot(xs_lim_r, ys_lim_r, 'c:', 'LineWidth', 1)
%  %(2) The left wall
%    plot(x_left, y_left, 'b', 'LineWidth', 4)
%    xs_lim_l = [(max(x_left)+l_wall_lim); min(x_left); min(x_left)];
%    ys_lim_l = [max(y_left); min(y_left)-l_wall_lim; max(y_left)];
%    plot(xs_lim_l, ys_lim_l, 'c:', 'LineWidth', 1)
%  %(3) The upper wall
%    plot(x_wall_up, y_wall_up, 'b', 'LineWidth', 4)
%    plot(x_wall_up, y_wall_up+up_wall_lim, 'c:', 'LineWidth', 1)
 %(4) The area where the force could have been applied
%    plot([x_area(1) x_area(1) x_area(2) x_area(2)], [y_area(2) y_area(1) y_area(1) y_area(2)], 'y', 'LineWidth', 4) 
%      hold off   
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 User Input Not Required Below This Line                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if HoughDebug
    imageFile = [directory,files(1).name]; %input filename
    img = imread(imageFile);
    Rimg = img(:,:,1);
    [SL, SS] = PeGSHoughDebug(Rimg, RlargeH, SL, RsmallH, SS, DS, NlargeH, NsmallH);    
elseif ConvDebug
    imageFile = [directory,files(1).name]; %input filename
    img = imread(imageFile);
    [RsmallC, RlargeC, NlargeC, NsmallC] = PeGSConvDebug(img, RsmallC, NsmallC, RlargeC, NlargeC);    
end


for frame = 1:nFrames %Loops for total number of images
    
    if doParticleDetectionH || doParticleDetectionC
        imageFile = [directory,files(frame).name];  %Input filename
        img = imread(imageFile);  %Read an experimental color image
        info = imfinfo(imageFile);  %Save information about the image
        Width = double(info.Width);
        Height = double(info.Height); 
%         Rimg = img(:,:,1);  %Particle image -- different definition in this project
%         Gimg = img(:,:,2);  %Force image -- different definition in this project
    else
        centersfile = [directory, files(frame).name]; %input filename
        gImgFile = [directory, 'Img',centersfile(8:end-3),'jpg'];  %adjusted force image filename
        rImgFile = [directory, 'Img',centersfile(8:end-3),'jpg'];  %adjusted force image filename
        Rimg = imread(rImgFile); %particle image
        Gimg = imread(gImgFile); %force image
    end
    


    %% Photo pre-processing -- Step 1: Channel distinction

    % 'Red channel' image is used for the detection of particles, so visible edges are beneficial
%     photo = img;    
%     photo_red_ch = photo(:,:,1);  %Save the red channel of the image
%     photo_BW = im2double(photo_red_ch); 
%     Rimg = imadjust(photo_BW);  % Adjust the top and bottom 1% intensities, for better circle finding
% 
%     % Black out bright areas near the obstacle, to avoid particle detection there
%     obstacle_xs = [Width/2-550 Width/2+100 Width/2+750];
%     obstacle_ys = [Height Height-500 Height];
%     obstacle_blank = roipoly(Rimg, obstacle_xs, obstacle_ys);
% 
%     Rimg = double(~obstacle_blank).*Rimg;
%     Rimg(Height-250:Height, :) = 0;


% % Error in PeGSDiskPrep_static (line 588)
% %     particle(p_no(i)).z = particle(p_no(i)).z - 1;  %Reduce the coordination number


    photo = img;  
    photo_gray = rgb2gray(photo);  
    photo_BWall = im2double(photo_gray);  
    New_img = imadjust(photo_BWall);  % Adjust the top and bottom 1% intensities, for better circle finding

        % Black out bright areas near the obstacle, to avoid particle detection there
    obstacle_xs = [Width/2-550 Width/2+100 Width/2+750];
    obstacle_ys = [Height Height-500 Height];
    obstacle_blank = roipoly(New_img, obstacle_xs, obstacle_ys);
    New_img = double(~obstacle_blank).*New_img;
    New_img(Height-250:Height, :) = 0;


    %% Photo pre-processing - Step 2: Evening out the lighting in the image
    % Adjusting the light intensity relative to the distance from centre of the image for even distribution
    coefficients = light_correction_coefficients;
    coefficients_corr = [-coefficients(1) coefficients(2) 0];
    centre = [0.5*Height 0.5*Width];
    pixel_y = linspace(1, Height, Height);
    pixel_x = linspace(1, Width, Width);
    % Find distance from image centre of each pixel
    for i = 1:Width
        for n = 1:Height
        dist_centre(n, i) = sqrt((abs(pixel_x(i)) - centre(2)).^2 + (abs(pixel_y(n)) - centre(1)).^2);
        end
    end
    % Correct the intensities in grey-scale image according to the intensity correction function
    photo = rgb2gray(img);
    photo = im2double(photo);
    I_corr_photo = polyval(coefficients_corr, dist_centre);
    photo_corr = I_corr_photo + photo;
   
    %% Photo pre-processing -- Step 3: Adjustment of the grey scale image
    % Adjusted grey-scale image can be used to estimate the forces on the particle
    % Elimate bright areas near right and left walls      
    for i = 1:length(x_right)
        x_right(i) = round(x_right(i));
        y_right(i) = round(y_right(i));
        A = x_right(i):1:Width;
        B = zeros(1,length(A));
        B(1,:) = y_right(i);
        photo_corr(B, A, :) = 0;
    end
    for i = 1:length(x_left)
        x_left(i) = round(x_left(i));
        y_left(i) = round(y_left(i));
        A = 1:1:x_left(i);
        B = zeros(1,length(A));
        B(1,:) = y_left(i);
        photo_corr(B, A, :) = 0;
    end  
    % Eliminate bright area near the obstacle
    photo_corr_2 = photo_corr(1:Height, 1:Width);
    obstacle_xs = [Width/2-550 Width/2+100 Width/2+750];
    obstacle_ys = [Height Height-500 Height];
    obstacle_blank = roipoly(photo_corr_2, obstacle_xs, obstacle_ys);
    photo_corr_2 = double(~obstacle_blank).*photo_corr_2;
    photo_corr_2(Height-250:Height, :) = 0;
    
    % Adjusts the image, saturating all intensities below a set amount to zero
    photo_corr_3 = imadjust(photo_corr_2, [I_low_all I_high_all]);
    %figure(104), imshow(photo_corr_3)  %Uncomment to see the adjusted image

    % Code below can further adjust areas near walls and in upper corners 
    % (This can eliminate the bright particle edges, but requires additional force calibration for accurate force estimation)
    % (1) The correction at the upper left corner
    xs_tl = [1 650 1300 650];
    ys_tl = [min(y_left) min(y_left)-600 min(y_left) min(y_left)+600];
    tl = roipoly(photo_corr_3,xs_tl,ys_tl);
    photo_corr_4 = double(tl).*photo_corr_2;
    photo_corr_4 = imadjust(photo_corr_4, [I_low_tl I_high_tl]);
    photo_corr_5 = double(~tl).*photo_corr_3;
    photo_corr_5 = imadd(photo_corr_4, photo_corr_5);
    % (2) The correction at the upper right corner
    xs_tr = [Width-1500 Width-800 Width Width-650];
    ys_tr = [min(y_left)-200 min(y_left)-1000 min(y_left)-170 min(y_left)+900];
    tr = roipoly(photo_corr_5,xs_tr,ys_tr);
    photo_corr_6 = double(tr).*photo_corr_2;
    photo_corr_6 = imadjust(photo_corr_6, [I_low_tr I_high_tr]);
    photo_corr_7 = double(~tr).*photo_corr_5;
    photo_corr_7 = imadd(photo_corr_7, photo_corr_6);

    % (3) The correction near the left wall
    xs_l_edge = [600 1000 max(x_left)+600 max(x_left)];
    ys_l_edge = [min(y_left)+600 min(y_left)+220 max(y_left)-200 max(y_left)];
    l_edge = roipoly(photo_corr_7,xs_l_edge,ys_l_edge);
    photo_corr_8 = double(l_edge).*photo_corr_7;
    photo_corr_8 = imadjust(photo_corr_8, [I_low_el I_high_el]);
    photo_corr_9 = double(~l_edge).*photo_corr_7;
    photo_corr_9 = imadd(photo_corr_8, photo_corr_9);
    % (4) The correction near the right wall
    xs_r_edge = [(min(x_right)-1000) Width-1000 Width min(x_right)];
    ys_r_edge = [max(y_right) min(y_right) min(y_right) max(y_right)];
    r_edge = roipoly(photo_corr_9,xs_r_edge,ys_r_edge);
    photo_corr_10 = double(r_edge).*photo_corr_9;
    photo_corr_10 = imadjust(photo_corr_10, [I_low_er I_high_er]);
    photo_corr_11 = double(~r_edge).*photo_corr_9;
    photo_corr_11 = imadd(photo_corr_10, photo_corr_11);   
    
    % Image used for force estimation, adjusted same as the calibration ones
    photo_for_forces = photo_corr_11;
    Gimg = photo_corr_11;  %Equivalent of Gimg used in original PeGS algorithm
    % Brightened image used for pseudo image creation
    photo_corr_12 = imadjust(photo_corr_11);  %Top and botttom 1% intensities saturated
    photo_for_pseudo = photo_corr_12;
%     figure(105), imshow(photo_corr_12)  %Uncomment to see the adjusted image
        
 %%  Results section
    if (verbose)
        figure(1); %Draw the particle Image
%         imshow(Rimg);
          imshow(New_img);
        figure(2); %Draw the Force Image
        imshow(photo_for_pseudo); 
    end
    
    if doParticleDetectionH 
%         particle = PeGSDiskFindH(Rimg, RlargeH, SL, RsmallH, SS, pxPerMeter, fsigma);
        particle = PeGSDiskFindH(New_img, RlargeH, SL, RsmallH, SS, pxPerMeter, fsigma);  

    elseif doParticleDetectionC
        particle = PeGSDiskFindC(img, RsmallC, NsmallC, RlargeC, NlargeC);
    else
        pData = dlmread(centersfile); %Read Position data from centers file
        N = size(pData,1);
        particle(1:N) = struct('id',0,'x',0,'y',0,'r',0,'rm',0,'color','','fsigma',0,'z',0,'f',0,'g2',0,'forces',[],'betas',[],'alphas',[],'neighbours',[],'contactG2s',[],'forceImage',[]);
        for n=1:N %Bookkeeping
            particle(n).id= n;
            particle(n).x = pData(n,2); %-xoffset;
            particle(n).y = pData(n,3); %-yoffset;
            particle(n).r = pData(n,4);
        end        
    end
    
    N0 = length(particle);  %Number of all found circles
    % Set the recognised particle radius to the desired value uniform for all circles
    for n = 1:N0
        particle(n).r = radius_desired;
        particle(n).rm = particle(n).r * pxPerMeter;
    end
    % Ensure that there is no overlapping circles recognised
    for a = 1:N0
        for n = 1:N0
            dist_parts_cent = sqrt((particle(n).x-particle(a).x)^2 + (particle(n).y-particle(a).y)^2);
            if particle(n).id > particle(a).id && dist_parts_cent < 1.8*particle(n).r
                particle(n).id = 1000;
            end
        end
    end
    particle([particle.id]==1000) = [];  %Remove all particles recognised as overlapping
   
    % Sort the particle IDs with increasing y cooridnate of the centre
    T_particle = struct2table(particle);
    sorted_T = sortrows(T_particle, 'y');
    particle = table2struct(sorted_T);
    particles = 1:1:length(particle);
    for p = 1:length(particle)
        particle(p).id = p;
    end
    
    N = length(particle);  %Actual number of particles on the image
    if(verbose)
        %add some information about the particles to the plots
        figure(1)
        for n=1:N
            viscircles([particle(n).x; particle(n).y]', particle(n).r,'EdgeColor',particle(n).color, 'LineStyle', 'none');  %Draw particle outline
            hold on
            plot(particle(n).x,particle(n).y,'rx');  %Mark particle centers
            text(particle(n).x,particle(n).y,num2str(particle(n).id),'Color','g');  %Write particle IDs
        end
        figure(2)
        for n=1:N
            viscircles([particle(n).x; particle(n).y]', particle(n).r,'EdgeColor',particle(n).color);  %Draw particle outline
            hold on
            plot(particle(n).x,particle(n).y,'rx');  %Mark particle centers
            text(particle(n).x,particle(n).y,num2str(particle(n).id),'Color','b');  %Write particle IDs
        end
        drawnow;
    end
    
%% Force estimation for each particle
for n=1:N
        %Create a circular mask
        r = particle(n).r;
        mask = abs(-r:r);
        mask = mask.^2 + mask.^2';
        mask1 = double(sqrt(mask) <= r);

        %This crops out a particle
        cropXstart = round(particle(n).x-r);
        cropXstop = round(particle(n).x-r)+ size(mask1,1)-1;
        cropYstart = round(particle(n).y-r);
        cropYstop = round(particle(n).y-r)+ size(mask1,2)-1;
        
        % Code below allows for cropping a cut-off particle that exceeds the
        % image bounds, however ammending still required for PeGSDiskSolve
        if cropYstop > info.Height
            yoverhang = abs(cropYstop - info.Height);
            yoverhang1 = length(mask1)-yoverhang;
            mask1n = mask1(1:yoverhang1, :);
            cropYstop = info.Height;
        elseif cropYstart <= 0
            yoverhang = abs(cropYstart);
            mask1n = mask1(yoverhang+2:length(mask1), :);   
            cropYstart = 1;
        else mask1n = mask1;
        end
        if cropXstop > info.Width
            xoverhang = cropXstop - info.Width;
            xoverhang1 = abs(length(mask1)-xoverhang);
            if cropYstop == info.Height
                mask1n = mask1(1:yoverhang1, 1:xoverhang1);
            elseif cropYstart == 1
                mask1n = mask1(yoverhang+2:length(mask1), 1:xoverhang1);
            else  mask1n = mask1(:, 1:xoverhang1);
            end
            cropXstop = info.Width; 
        elseif cropXstart <= 0
            xoverhang = abs(cropXstart);
            if cropYstop == info.Height
                mask1n = mask1(1:yoverhang1, xoverhang+2:length(mask1));
            elseif cropYstart == 1
                mask1n = mask1(yoverhang+2:length(mask1), xoverhang+2:length(mask1));
            else  mask1n = mask1(:, xoverhang+2:length(mask1));
            end
            cropXstart = 1;
        end
        
        % Use the brighter image to create the pseudo image
        cimg = photo_for_pseudo(cropYstart:cropYstop, cropXstart:cropXstop);
        particleImg = cimg.*mask1n;
        particle(n).forceImage=particleImg;
        
        % Use the darker image to estimate the forces on the particle
        fimg = photo_for_forces(cropYstart:cropYstop, cropXstart:cropXstop);
        particleFImg = fimg.*mask1n;
        
        %Create a circular mask with a radius that is one pixel smaller
        %for cropping out the relevant gradient
        mask2 = double(sqrt(mask) <= r-1);

        if cropXstart == 1
           mask2n = mask2(1:size(mask1n, 1), (size(mask2, 2)-size(mask1n, 2)+1):size(mask2, 2));
        elseif cropXstop == info.Width
           mask2n = mask2(1:size(mask1n, 1), 1:size(mask1n, 2));
        elseif cropYstart == 1
           mask2n = mask2((size(mask2, 1)-size(mask1n, 1)+1):size(mask2, 1), 1:size(mask1n, 2));
        elseif cropYstop == info.Height
           mask2n = mask2(1:size(mask1n, 1), 1:size(mask1n, 2));
        else mask2n = mask2;
        end
        
        %Compute G^2 for each particle
        [gx,gy] = gradient(particleFImg);
        g2 = (gx.^2 + gy.^2).*mask2n;
        particle(n).g2 = sum(sum(g2));  %G2 value for the particle
        particle(n).f = particle(n).g2/g2cal;  %Force estimate for the particle
        particle(n).CR = CR;  %CR value for the particle
end
    

%% Flexible selection for the contactG2Threshold value 
% Code below allows for flexible adjustment of the contactG2Threshold value
% (This can be used to increase the threshold near the walls rather than 
%  increasing the cut-off at these sections and customising calibration 
%  procedure, or to consider other g2 qualities in the threshold selection)
for p = 1:length(particle)
    % (1) If the contactG2Threshold value should be uniform for all particles:
         particle(p).contactG2Threshold = contactG2Threshold;
    % (2) If the contactG2Threshold value should depend on particle's position on the image
    % (a) Set rectangular bounds:
    %     if particle(p).y <= 1520 && particle(p).x <= 1260
    %        particle(p).contactG2Threshold = 8;
    %     elseif particle(p).y <= 1520 && particle(p).x >= 4025
    %        particle(p).contactG2Threshold = 8;
    %     elseif particle(p).y <= 600
    %        particle(p).contactG2Threshold = 8;
    %     else 
    %        particle(p).contactG2Threshold = 8;
    %     end
    % (b) A particluar distance from a particluar particle:
    %     d = sqrt((abs(particle(p).x) - abs(particle(139).x))^2 + (abs(particle(p).y) - abs(particle(139).y))^2);
    %     if d >= sqrt((abs(particle(128).x) - abs(particle(139).x))^2 + (abs(particle(128).y) - abs(particle(139).y))^2)
    %        particle(p).contactG2Threshold = 2;
    %     else particle(p).contactG2Threshold = 9;
    %     end
    % (c) Being at the right or left wall:
        % First row at the wall:
        xs_tri_l1 = [(max(x_left)+2*radius_desired); min(x_left); min(x_left)];
        ys_tri_l1 = [max(y_left); min(y_left)-2*radius_desired; max(y_left)];
        in_l1 = inpolygon(particle(p).x, particle(p).y, xs_tri_l1, ys_tri_l1);
        xs_tri_r1 = [(min(x_right)-2*radius_desired); max(x_right); max(x_right)];
        ys_tri_r1 = [max(y_right); min(y_right)-2*radius_desired; max(y_right)];
        in_r1 = inpolygon(particle(p).x, particle(p).y, xs_tri_r1, ys_tri_r1);
        % Second row at the wall:
        xs_tri_l2 = [(max(x_left)+2*radius_desired); (max(x_left)+4*radius_desired); (min(x_left)+3*radius_desired); (min(x_left)+1.5*radius_desired)];
        ys_tri_l2 = [max(y_left); max(y_left); min(y_left)-4*radius_desired; min(y_left)-2*radius_desired];
        in_l2 = inpolygon(particle(p).x, particle(p).y, xs_tri_l2, ys_tri_l2);
        xs_tri_r2 = [(min(x_right)-2*radius_desired); (min(x_right)-4*radius_desired); (max(x_right)-2.5*radius_desired); (max(x_right)-1.5*radius_desired)];
        ys_tri_r2 = [max(y_right); max(y_right); min(y_right)-4*radius_desired; min(y_right)-2*radius_desired];
        in_r2 = inpolygon(particle(p).x, particle(p).y, xs_tri_r2, ys_tri_r2);
        % Third/fourth row at the wall:
        xs_tri_l3 = [(max(x_left)+4*radius_desired); (max(x_left)+6*radius_desired); (min(x_left)+3.5*radius_desired); (min(x_left)+2*radius_desired)];
        ys_tri_l3 = [max(y_left); max(y_left); min(y_left)-6*radius_desired; min(y_left)-4*radius_desired];
        in_l3 = inpolygon(particle(p).x, particle(p).y, xs_tri_l3, ys_tri_l3);
        xs_tri_r3 = [(max(x_right)-12*radius_desired);(max(x_right)-9*radius_desired); (max(x_right)-3*radius_desired); (max(x_right)-5*radius_desired)];
        ys_tri_r3 = [0.5*max(y_right)-3*radius_desired; 0.5*max(y_right)+0*radius_desired; min(y_right)-3*radius_desired; min(y_right)-5*radius_desired];
        in_r3 = inpolygon(particle(p).x, particle(p).y, xs_tri_r3, ys_tri_r3);
        % Force application area:
        in_force_app = inpolygon(particle(p).x, particle(p).y, x_area, y_area);
        % First row at the upper wall:
        in_u1 = inpolygon(particle(p).x, particle(p).y-1.6*particle(1).r, x_wall_up, y_wall_up);
        if in_force_app == 1  %Limit within the force application area
            particle(p).contactG2Threshold = 0.0500;
%         elseif in_u1 == 0  %Limit within the first particle row at upper wall
% %             particle(p).CR = 20;
%             particle(p).contactG2Threshold = 0.0130;
        elseif in_r1 == 1  %Limit within the first particle row at right wall
            particle(p).contactG2Threshold = 0.0280;   
        elseif in_l1 == 1  %Limit within the first particle row at left wall
            particle(p).contactG2Threshold = 0.0080;
        elseif in_r2 == 1  %Limit within the second particle row at right wall
            particle(p).contactG2Threshold = 0.0290;  
        elseif in_l2 == 1  %Limit within the second particle row at right wall
            particle(p).contactG2Threshold = 0.0110;
%         elseif in_l3 == 1  %Limit within the third particle row at right wall
%             particle(p).contactG2Threshold = 0.003;
        elseif in_r3 == 1  %Limit within the third particle row at right wall
            particle(p).contactG2Threshold = 0.0088;
        end
    
    % (3) If the contactG2Threshold should be unique to each particle, for example based on the mean G2 around its perimeter
    % Find information about all possible contacts around the particle's circumference
%         mask = abs(-CR:CR);
%         mask = mask.^2 + mask.^2';
%         mask = double(sqrt(mask) <= CR-1);         
%         x_cent = particle(p).x;
%         y_cent = particle(p).y;
%         rad = particle(p).r;
%         theta_bl = linspace(0,2*pi,100); 
%     
%         for n = 1:length(theta_bl)
%             y = y_cent-(rad-CR)*sin(theta_bl(n));
%             x = x_cent-(rad-CR)*cos(theta_bl(n));
%             contactImg = (imcrop(photo_for_forces,[x-CR y-CR CR*2 CR*2]));
%             contactImg = contactImg.*mask; 
%             [gx,gy] = gradient(contactImg);
%             g2 = (gx.^2 + gy.^2).*mask;
%             contactG2(n) = sum(sum(g2));
%         end
%         particle(p).contactG2Threshold = (mean(contactG2));
end

%% Neighbour finding and plotting the results
    
    if findNeighbours    
        particle = PeGSNeighbourFind(Gimg, contactG2Threshold, dtol, CR, verbose, particle, slope_right, x_right, y_right, r_wall_lim, slope_left, x_left, y_left, l_wall_lim, r_rheo, x_cent_rheo, y_cent_rheo, th, x_wall_up, y_wall_up, up_wall_lim, x_area, y_area);    
    end
     
    % Uncomment below to remove contacts manually 
    %%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    %%

%     p_no = [1,2，98]; %IDs of particles with false contacts
    p_no = [];  

    c_no = [1,1,4];  
    for i = 1:length(p_no)
    particle(p_no(i)).z = particle(p_no(i)).z - 1;  %Reduce the coordination number
    particle(p_no(i)).betas(c_no(i)) = [];  %Remove the contact angle info
    particle(p_no(i)).neighbours(c_no(i)) = [];  %Remove the neighbour info
    particle(p_no(i)).contactG2s(c_no(i)) = [];  %Remove the contact G2 info
    end    
       
    figure(3);
    %%%%%%%%%%%%%%%%%%% force chain %%%%%%%%%%%%%%%%%%%%%%%%
    imshow(img); hold on
    for n = 1:N
        z = particle(n).z; %get particle coordination number
        if (z>0) %if the particle does have contacts
            for m = 1:z %for each contact
                %draw contact lines
                lineX(1)=particle(n).x;
                lineY(1)=particle(n).y;
                lineX(2) = lineX(1) + particle(n).r * cos(particle(n).betas(m));
                lineY(2) = lineY(1) + particle(n).r * sin(particle(n).betas(m));
                cX = lineX(1) + (particle(n).r-CR) * cos(particle(n).betas(m));
                cY = lineY(1) + (particle(n).r-CR) * sin(particle(n).betas(m));
                hold on; % Don't blow away the image.
                plot(lineX, lineY,':y','LineWidth',2);hold on;
                text(particle(n).x,particle(n).y,num2str(particle(n).id),'Color','g');  %Write particle IDs
            end
        else
            text(particle(n).x,particle(n).y,num2str(particle(n).id),'Color','g');  %Write particle IDs
        end
    end


    figure(4);
    imshow(Gimg); 
    hold on;

    % Determine the maximum force magnitude for normalization
    max_force = max([particle.f]);

    % Set the range for line widths
    min_line_width = 1;
    max_line_width = 8;

    % Define color extremes for minimum and maximum forces
    color_min = [0, 0, 1];  % Blue for minimum force
    color_max = [1, 0, 0];  % Red for maximum force

    for n = 1:N
        z = particle(n).z; %get particle coordination number
        if (z>0) %if the particle does have contacts
            for m = 1:z %for each contact

                % Draw contact lines
                lineX(1) = particle(n).x;
                lineY(1) = particle(n).y;
                lineX(2) = lineX(1) + particle(n).r * cos(particle(n).betas(m));
                lineY(2) = lineY(1) + particle(n).r * sin(particle(n).betas(m));

                normalized_force = particle(n).f / max_force;
                line_width = min_line_width + normalized_force * (max_line_width - min_line_width); % 线宽
                
                % Calculate color based on force magnitude using linear interpolation
                line_color = color_min + normalized_force * (color_max - color_min);
                
                plot(lineX, lineY, '-', 'Color', line_color, 'LineWidth', line_width); hold on;
            end
        end
    end

    hold off;
    
    for n = 1:N
        % Crop out a pixel wide cross-section area of each particle's centre
        cropXstart = round(particle(n).x-r);
        cropXstop = round(particle(n).x+r);
        cropYstart = round(particle(n).y-1);
        cropYstop = round(particle(n).y+1);
        fimg = photo_for_pseudo(cropYstart:cropYstop, cropXstart:cropXstop);

        % Find the number of bright fringes from pixel intensities
        n_ones = find(fimg(2,:) > 0.50);
        if n_ones > 0
            n_p_peaks = zeros(1,length(n_ones));
            for i=1:length(n_ones)-1
                if n_ones(i+1)-n_ones(i) > 20  %Find number of peaks by considering the number of bright pixel clusters
                    n_p_peaks(i) = i;
                end
            end
            particle(n).peaks = nnz(n_p_peaks)+1;
        else particle(n).peaks = 0;
        end
        
        if particle(n).peaks == 2 && particle(n).f < 2.0  %Restrict the force on particle for formation of the first fringe
            particle(n).peaks = 1;
        end
    end
    
    %Save what we got so far
    save([directory, files(frame).name(1:end-4),'_contact.mat'],'particle');
    
end

