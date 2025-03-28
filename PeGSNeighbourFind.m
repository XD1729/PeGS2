%{
概括：
1.函数定义与初始化:
函数的输出是更新后的 particle。
函数内部定义了数组 xmat、ymat 和 rmat 来存储 particle 结构中的 x、y 和 r 字段。
2.提取 particle 数据:
使用循环遍历 particle 结构，提取每个粒子的 x、y 和 r 值，并将它们存储在先前初始化的 xmat、ymat 和 rmat 数组中。
计算距离矩阵:
3.使用 pdist2 函数计算粒子中心位置的距离矩阵 dmat。
通过将 rmat 与其转置相加，为每个粒子创建一个组合半径矩阵。
============================
粒子处理：
确定接触的粒子:
使用 friendmat 变量来确定哪些粒子是相互接触的。
使用 find 函数从 friendmat 中提取相互接触的粒子索引 f1 和 f2。

对每一对接触的粒子进行处理:
对于每一对接触的粒子，函数提取了它们的 x、y 和 r 值。
函数计算了接触区域的位置，并使用 imcrop 函数从原始图像 Gimg 中裁剪出接触区域的图像。
函数计算了接触区域的梯度平方 g2。
如果接触区域的梯度平方超过了设定的阈值 contactG2Threshold，则将这两个粒子视为接触的，并更新了 particle 结构中的相关字段，例如 z、contactG2s、contactIs、neighbours 和 betas。

考虑左墙的影响:
函数定义了一个与左墙接触的三角形区域。
对于每一个粒子，函数检查了它是否与左墙接触，并计算了接触区域的梯度平方 g2。
如果粒子与左墙接触并且接触区域的梯度平方超过了设定的阈值 contactG2Threshold，则更新了 particle 结构中的相关字段。
总体来说，函数 PeGSNeighbourFind.m 的主要目的是检测粒子之间以及粒子与墙之间的接触。这是通过计算接触区域的梯度平方来完成的，如果梯度平方超过了设定的阈值，则两个粒子或粒子与墙之间被视为接触的。函数还更新了 particle 结构中的相关字段，以记录这些接触信息。
%}


function [particle] = PeGSNeighbourFind(Gimg, contactG2Threshold, dtol, CR, verbose, particle, slope_right, x_right, y_right, r_wall_lim, slope_left, x_left, y_left, l_wall_lim, r_rheo, x_cent_rheo, y_cent_rheo, th, x_wall_up, y_wall_up, up_wall_lim, x_area, y_area)
%% Original code for contact detection based on physical particle touch
N = length(particle);
xmat = zeros([N,1]);
ymat = zeros([N,1]); 
rmat = zeros([N,1]);

for l = 1:N
    xmat(l) = particle(l).x;
    ymat(l) = particle(l).y; %Pulls data from particle structure
    rmat(l) = particle(l).r;
end
rmats = rmat; %Saves our radius matrix for later

dmat = pdist2([xmat,ymat],[xmat,ymat]); %Creates a distance matrix for particle center locations
rmat = rmat + rmat'; %Makes a combination of radii for each particle

friendmat = dmat < (rmat + dtol) & dmat~=0; %Logical "friend" matrix

friendmat = triu(friendmat); %Only examine the upper triangle portion (no repeats)
[f1, f2] = find(friendmat == 1); %Creates an index of particles that are considered touching

for l = 1:length(f1)
    x1 = particle(f1(l)).x;
    y1 = particle(f1(l)).y;
    r1 = particle(f1(l)).r;
    x2 = particle(f2(l)).x;
    y2 = particle(f2(l)).y;
    r2 = particle(f2(l)).r;
    
    mask_f1 = abs(-particle(f1(l)).CR:particle(f1(l)).CR);
    mask_f1 = mask_f1.^2 + mask_f1.^2';
    mask_f1 = double(sqrt(mask_f1) <= particle(f1(l)).CR-1);
    
    contactXp1 = x1 + (r1 - particle(f1(l)).CR) * cos(atan2(y2-y1,x2-x1));
    contactYp1 = y1 + (r1 - particle(f1(l)).CR) * sin(atan2(y2-y1,x2-x1));
    
    contactXp2 = x1 + (r1 + particle(f2(l)).CR + dmat(f1(l),f2(l)) - rmat(f1(l),f2(l))) * cos(atan2(y2-y1,x2-x1));
    contactYp2 = y1 + (r1 + particle(f2(l)).CR + dmat(f1(l),f2(l)) - rmat(f1(l),f2(l))) * sin(atan2(y2-y1,x2-x1));
    
    contactImg = (imcrop(Gimg,[contactXp1-particle(f1(l)).CR contactYp1-particle(f1(l)).CR particle(f1(l)).CR*2 particle(f1(l)).CR*2]));
    contactImg = contactImg.*mask_f1;
    
    [gx,gy] = gradient(contactImg);
    g2 = (gx.^2 + gy.^2).*mask_f1;
    contactG2p1 = sum(sum(g2));
    contactIp1 = sum(sum(contactImg));
    
    mask_f2 = abs(-particle(f2(l)).CR:particle(f2(l)).CR);
    mask_f2 = mask_f2.^2 + mask_f2.^2';
    mask_f2 = double(sqrt(mask_f2) <= particle(f2(l)).CR-1);
    
    contactImg = (imcrop(Gimg,[contactXp2-particle(f2(l)).CR contactYp2-particle(f2(l)).CR particle(f2(l)).CR*2 particle(f2(l)).CR*2]));
    contactImg = contactImg.*mask_f2;
    
    [gx,gy] = gradient(contactImg);
    g2 = (gx.^2 + gy.^2).*mask_f2;
    contactG2p2 = sum(sum(g2));
    contactIp2 = sum(sum(contactImg));
    
    %if we declare our contact valid
    if(contactG2p1 > particle(f1(l)).contactG2Threshold && contactG2p2 > particle(f2(l)).contactG2Threshold)
        %cI = sum(sum(contactImg)); %Use integrated intensity instead of g2
        %if(cI > contactIThreshold) %Use integrated intensity instead of g2
        %Plot contact area
        if (verbose)
            display(['contact found between particles ',num2str(f1(l)),' and ',num2str(f2(l))]);
            viscircles([contactXp1; contactYp1]', particle(f1(l)).CR,'EdgeColor','g');
            viscircles([contactXp2; contactYp2]', particle(f2(l)).CR,'EdgeColor','g');
            %text(contactXp1, contactYp1,num2str(contactG2p1),'Color','w');
            %drawnow;
        end
        %this is a valid contact, remember it
        particle(f1(l)).z= particle(f1(l)).z +1; %increase coordination number
        particle(f1(l)).contactG2s(particle(f1(l)).z)=contactG2p1; %remember the g2 value of the current contact area
        particle(f1(l)).contactIs(particle(f1(l)).z)=contactIp1;
        particle(f1(l)).neighbours(particle(f1(l)).z) = f2(l); %particle m is now noted as a neigbour in the particle l datastructure
        particle(f1(l)).betas(particle(f1(l)).z) = atan2(y2-y1,x2-x1); %the contact angle to particle m is now noted in the particle l datastructure
        particle(f2(l)).z= particle(f2(l)).z +1; %increase coordination number
        particle(f2(l)).contactG2s(particle(f2(l)).z)=contactG2p2; %remember the g2 value of the current contact area
        particle(f2(l)).contactIs(particle(f2(l)).z)=contactIp2;
        particle(f2(l)).neighbours(particle(f2(l)).z) = f1(l); %particle m is now noted as a neigbour in the particle l datastructure
        particle(f2(l)).betas(particle(f2(l)).z) = atan2(y1-y2,x1-x2);
    end
end

%% New code for custom rheometer boundaries:

% Accounting for the left wall:
xs_tri_l = [(max(x_left)+l_wall_lim); min(x_left); min(x_left)];
ys_tri_l = [max(y_left); min(y_left)-l_wall_lim; max(y_left)];
for p = 1:length(particle)
    % Get the necessary particle info
    y_cent = particle(p).y;
    x_cent = particle(p).x;
    rad = particle(p).r;
    CR = particle(p).CR;
    mask = abs(-CR:CR);
    mask = mask.^2 + mask.^2';
    mask = double(sqrt(mask) <= CR-1);
    % Consider contact perpendicular to wall
    y = y_cent+(rad-CR)*sin(atan(slope_left));
    x = x_cent-(rad-CR)*cos(atan(slope_left));
    contactImg = im2double(imcrop(Gimg,[x-CR y-CR CR*2 CR*2]));
    contactImg = contactImg.*mask; 
    [gx,gy] = gradient(contactImg);
    g2 = (gx.^2 + gy.^2).*mask;
    contactG2 = sum(sum(g2));
    in_wall_l = inpolygon(particle(p).x, particle(p).y, xs_tri_l, ys_tri_l);
    % For particles located within a set distance from the wall, check the G2 threshold criterium
    if in_wall_l == 1 && contactG2 >= particle(p).contactG2Threshold
        if(verbose)
            %text(x,y,num2str(contactG2),'Color','w');
            viscircles([x; y]', CR,'EdgeColor','y');
        end
        % Remember the parameters of a valid contact
        particle(p).z = particle(p).z +1;
        particle(p).betas(particle(p).z) = atan(slope_left)+pi/2;
        particle(p).contactG2s = [particle(p).contactG2s contactG2];
        particle(p).neighbours(particle(p).z) = -1;  %Neighbour ID set to an arbitratry negative number
    end
end

% Accounting for the right wall:
xs_tri = [(min(x_right)-r_wall_lim); max(x_right); max(x_right)];
ys_tri = [max(y_right); min(y_right)-r_wall_lim; max(y_right)];
for p = 1:length(particle)
    % Get the necessary particle info
    y_cent = particle(p).y;
    x_cent = particle(p).x;
    rad = particle(p).r;
    y = (y_cent-(rad-CR)*sin(atan(slope_right)));
    x = x_cent+(rad-CR)*cos(atan(slope_right));
    CR = particle(p).CR;
    mask = abs(-CR:CR);
    mask = mask.^2 + mask.^2';
    mask = double(sqrt(mask) <= CR-1);
    % Consider contact perpendicular to wall
    contactImg = im2double(imcrop(Gimg,[x-CR y-CR CR*2 CR*2]));
    contactImg = contactImg.*mask; 
    [gx,gy] = gradient(contactImg);
    g2 = (gx.^2 + gy.^2).*mask;
    contactG2 = sum(sum(g2));
    in_wall_r = inpolygon(particle(p).x, particle(p).y, xs_tri, ys_tri);
    % For particles located within a set distance from the wall, check the G2 threshold criterium
    if in_wall_r == 1 && contactG2 >= particle(p).contactG2Threshold
        if(verbose)
            %text(x,y,num2str(contactG2),'Color','w');
            viscircles([x; y]', CR,'EdgeColor','b');
        end
        % Remember parameters of a valid contact
        particle(p).z = particle(p).z +1;
        particle(p).betas(particle(p).z) = atan(slope_right)+pi/2;
        particle(p).contactG2s = [particle(p).contactG2s contactG2];
        particle(p).neighbours(particle(p).z) = -1;  %Neighbour ID set to an arbitratry negative number
    end
end

% Accounting for the upper wall
for p = 1:length(particle)
    % Consider particles located within a set distance from the wall
    in_wall_up = inpolygon(particle(p).x, particle(p).y-up_wall_lim, x_wall_up, y_wall_up);
    if in_wall_up == 0
        for n = 1:length(th)
            % Get the necessary particle parameters
            y_cent = particle(p).y;
            x_cent = particle(p).x;
            rad = particle(p).r;
            y(n) = y_cent-(rad-CR)*sin(th(n));
            x(n) = x_cent-(rad-CR)*cos(th(n));
            CR = particle(p).CR;
            mask = abs(-CR:CR);
            mask = mask.^2 + mask.^2';
            mask = double(sqrt(mask) <= CR-1);
            % Find distance of each potential contact from the centre of rheometer
            dist_rheo_cent(n) = sqrt((x_cent_rheo - abs(x(n)))^2+(y_cent_rheo - abs(y(n)))^2);
        end
        % Use the contact with maximum distance to check against the G2 criterium
        [max_dist, id_dist] = max(dist_rheo_cent);
        if max_dist >= 0.90*r_rheo  %Check that the particle is indeed near the upper wall
            contactImg = im2double(imcrop(Gimg,[x(id_dist)-CR y(id_dist)-CR CR*2 CR*2]));
            contactImg = contactImg.*mask;
            [gx,gy] = gradient(contactImg);
            g2 = (gx.^2 + gy.^2).*mask;
            contactG2 = sum(sum(g2));
            if contactG2 >= particle(p).contactG2Threshold
                if(verbose)
%                     text(x(id_dist),y(id_dist),num2str(contactG2),'Color','w');
                    viscircles([x(id_dist); y(id_dist)]', CR,'EdgeColor','g');
                end
                % Remember parameters of a valid contact
                particle(p).z = particle(p).z +1;
                particle(p).betas(particle(p).z) = th(id_dist)-pi;
                particle(p).contactG2s = [particle(p).contactG2s contactG2];
                particle(p).neighbours(particle(p).z) = -2;  %Neighbour ID set to an arbitratry negative number
            end
        end
    end
end

% Acount for applied force
for p = 1:length(particle)
    % Gather the necessary particle parameters
    CR = particle(p).CR;
    mask = abs(-CR:CR);
    mask = mask.^2 + mask.^2';
    mask = double(sqrt(mask) <= CR-1);
    x_cent = particle(p).x;
    y_cent = particle(p).y;
    rad = particle(p).r;
    % Consider particles found within the specified area of possible contact application
    if (particle(p).x > x_area(1) && particle(p).x < x_area(2)) && (particle(p).y > y_area(1) && particle(p).x < y_area(2))
        theta_bl = linspace(0,2*pi,100); %% create arbitrary angles along the circle
        % For all angles find potential contact locations
        for n = 1:length(theta_bl)
            y = y_cent-(rad-CR)*sin(theta_bl(n));
            x = x_cent-(rad-CR)*cos(theta_bl(n));
            contactImg = im2double(imcrop(Gimg,[x-CR y-CR CR*2 CR*2]));
            contactImg = contactImg.*mask;
            %viscircles([x; y]', CR,'EdgeColor','w');  %Uncomment to see the little circular masks
            [gx,gy] = gradient(contactImg);
            g2 = (gx.^2 + gy.^2).*mask;
            contactG2 = sum(sum(g2));
            
            % Coordinates of the contact points found originally
            y_old_ct = y_cent+(rad-CR).*sin(particle(p).betas(:));
            x_old_ct = x_cent+(rad-CR).*cos(particle(p).betas(:));
            
            angle_th = 50*(pi/180);  %Threshold angle limiting the distance between the contact points,
            % (From theory 60 degrees but 50 degrees used here to account for any error during image catpure)
            
            for i = 1:length(particle(p).betas)
                %           % The minimum and maximum angle around the contact point found originally
                beta_min = particle(p).betas(i)-angle_th;
                beta_max = particle(p).betas(i)+angle_th;
                
                % Code below specifies when a potential contact should be excluded based on distance from the already found contact
                % (Conditions below are required to account for the angle sign change at the mid-point of left hemisphere.)
                if particle(p).betas(i) > 0 || (particle(p).betas(i) == pi || particle(p).betas(i) == -pi)
                    if beta_max > pi
                        beta_max = beta_max - 2*pi;
                        if theta_bl(n)-pi < 0
                            if theta_bl(n)-pi <= beta_max && theta_bl(n)-pi >= -pi
                                y = NaN;
                                x = NaN;
                            end
                        elseif theta_bl(n)-pi >= 0
                            if theta_bl(n)-pi >= beta_min && theta_bl(n)-pi <= pi
                                y = NaN;
                                x = NaN;
                            end
                        end
                    else
                        if theta_bl(n)-pi <= beta_max && theta_bl(n)-pi >= beta_min
                            y = NaN;
                            x = NaN;
                        end
                    end
                elseif particle(p).betas(i) <= 0
                    if beta_min < -pi
                        beta_min = beta_min + 2*pi;
                        if theta_bl(n)-pi < 0
                            if theta_bl(n)-pi <= beta_max && theta_bl(n)-pi >= -pi
                                y = NaN;
                                x = NaN;
                            end
                        elseif theta_bl(n)-pi >= 0
                            if theta_bl(n)-pi >= beta_min && theta_bl(n)-pi <= pi
                                y = NaN;
                                x = NaN;
                            end
                        end
                    else
                        if theta_bl(n)-pi <= beta_max && theta_bl(n)-pi >= beta_min
                            y = NaN;
                            x = NaN;
                        end
                    end
                end
            end
            
            % Find potential contacts not excluded by their proximity to the originally found contacts
            if(contactG2 > particle(p).contactG2Threshold) && isnan(x) == 0
                %Plot contact area
                particle(p).contactG2(n) = contactG2;
                particle(p).theta = theta_bl;
                particle(p).x_new(n) = x;
                particle(p).y_new(n) = y;
                if (verbose)
                    viscircles([x; y]', CR,'EdgeColor','b');  %New potential contacts in blue
                    for n = 1:length(x_old_ct)
                        viscircles([x_old_ct(n); y_old_ct(n)]', CR,'EdgeColor','r');  %Original contacts in red
                    end
                    %%text(x, y,num2str(contactG2),'Color','w');
                end
            else particle(p).x_new(n) = 0;
                particle(p).y_new(n) = 0;
            end
        end
        
        % Now choose the point with the greatest intensity gradient and exclude any others in close proximity to it
        % (not the optimal approach at greater applied forces, due to multiple fringes)
        if (particle(p).z) > 0  %If a particle already has a connection, it may be part of a contact network
            if sum((particle(p).x_new)) > 0  %Change 'if' to 'while' to find multiple connections
                [max_cont, idx] = max(particle(p).contactG2);
                y_diff = particle(p).y_new(idx)-y_cent;
                x_diff = particle(p).x_new(idx)-x_cent;
                if y_diff > 0 && x_diff > 0
                    atan_max = atan(y_diff/x_diff);
                elseif y_diff < 0 && x_diff < 0
                    atan_max = atan(y_diff/x_diff) - pi;
                elseif y_diff > 0 && x_diff < 0
                    atan_max = atan(y_diff/x_diff) + pi;
                elseif y_diff < 0 && x_diff > 0
                    atan_max = atan(y_diff/x_diff);
                elseif y_diff == 0 && x_diff > 0
                    atan_max = 0;
                elseif y_diff == 0 && x_diff < 0
                    atan_max = -pi;
                elseif y_diff > 0 && x_diff == 0
                    atan_max = pi/2;
                elseif y_diff < 0 && x_diff == 0
                    atan_max = -pi/2;
                end
                theta_max = atan_max+angle_th;
                theta_min = atan_max-angle_th;
                for a = 1:length(particle(p).x_new)
                    if a == idx
                        % Remember parameters of the contact point with greatest G2 value
                        particle(p).z = particle(p).z +1;
                        particle(p).betas(particle(p).z) = atan_max;
                        particle(p).y_new(a) = 0;
                        particle(p).x_new(a) = 0;
                        particle(p).contactG2s = [particle(p).contactG2s particle(p).contactG2(a)];
                        particle(p).contactG2(a) = 0;
                        particle(p).neighbours(particle(p).z) = -1; %Neighbour ID set to an arbitratry negative number
                        % Conditions below used when detection of multiple contacts is required
                    else
                        if atan_max > 0
                            theta_max = atan_max + angle_th;
                            if theta_max > pi
                                theta_max = theta_max - 2*pi;
                                if theta_bl(a)-pi < 0
                                    if theta_bl(a)-pi >= -pi && theta_bl(a)-pi <= theta_max
                                        particle(p).y_new(a) = 0;
                                        particle(p).x_new(a) = 0;
                                        particle(p).contactG2(a) = 0;
                                    end
                                elseif theta_bl(a)-pi >= 0
                                    if theta_bl(a)-pi >= theta_min && theta_bl(a)-pi <= pi
                                        particle(p).y_new(a) = 0;
                                        particle(p).x_new(a) = 0;
                                        particle(p).contactG2(a) = 0;
                                    end
                                end
                            else
                                if theta_bl(a)-pi <= theta_max && theta_bl(a)-pi >= theta_min
                                    particle(p).y_new(a) = 0;
                                    particle(p).x_new(a) = 0;
                                    particle(p).contactG2(a) = 0;
                                end
                            end
                        elseif atan_max <= 0
                            theta_min = atan_max - angle_th;
                            if theta_min < -pi
                                theta_min = theta_min + 2*pi;
                                if theta_bl(a)-pi < 0
                                    if theta_bl(a)-pi >= -pi && theta_bl(a)-pi <= theta_max
                                        particle(p).y_new(a) = 0;
                                        particle(p).x_new(a) = 0;
                                        particle(p).contactG2(a) = 0;
                                    end
                                elseif theta_bl(a)-pi >= 0
                                    if theta_bl(a)-pi >= theta_min && theta_bl(a)-pi <= pi
                                        particle(p).y_new(a) = 0;
                                        particle(p).x_new(a) = 0;
                                        particle(p).contactG2(a) = 0;
                                    end
                                end
                            else
                                if theta_bl(a)-pi <= theta_max && theta_bl(a)-pi >= theta_min
                                    particle(p).y_new(a) = 0;
                                    particle(p).x_new(a) = 0;
                                    particle(p).contactG2(a) = 0;
                                end
                            end
                        end
                    end
                end
            end
        elseif (particle(p).z) == 0
            display(['Particle ', num2str(particle(p).id), ' has no connections.'])
        end
    end
end
        
%% Original code below - can find contacts at rectangular walls

% %Check if any of the walls is a neighbour as well
% %TODO: Do this in a less verbose manner, I remember there was a function
% %somewhere that tells if somethig is inside a polygon or not
% 
% circs = [ymat, xmat, rmats]; %Makes a circs matrix from old matrices
% 
% rightwall = max(circs(:,2) + circs(:,3));
% leftwall = min(circs(:,2) - circs(:,3)); %Finds our theorhetical wall locations
% topwall = min(circs(:,1) - circs(:,3));
% bottomwall = max(circs(:,1) + circs(:,3));
% 
% rwi = find(circs(:,2) + circs(:,3) + dtol*1.5 >= rightwall);
% lwi = find(circs(:,2) - circs(:,3) - dtol*1.5 <= leftwall) %Indexes based on particles that would be considered to be touching the wall
% bwi = find(circs(:,1) + circs(:,3) + dtol*25 >= bottomwall)
% twi = find(circs(:,1) - circs(:,3) - dtol*1.5 <= topwall);
% 
% for l = 1:length(lwi) %Runs through each index to check for contacts via gradients
%     x = circs(lwi(l),2);
%     y = circs(lwi(l),1);
%     r = round(circs(lwi(l),3));
%     
%     contactX = x-(r-CR)
%     contactY = y
%     
%     contactImg = (imcrop(Gimg,[contactX-CR contactY-CR CR*2 CR*2]));
%     contactImg = contactImg.*mask;
%     imshow(contactImg)
%     
%     [gx,gy] = gradient(contactImg);
%     g2 = (gx.^2 + gy.^2);
%     contactG2 = sum(sum(g2))
%     
%     if(contactG2 > contactG2Threshold)
%         cI = sum(sum(contactImg));
%         %if(cI > contactIThreshold)
%         %this is a valid contact, remember it
%         if(verbose)
%             text(contactX,contactY,num2str(contactG2),'Color','w');
%             viscircles([contactX; contactY]', CR,'EdgeColor','w');
%         end
%         particle(lwi(l)).z= particle(lwi(l)).z +1; %increase coordination number
%         particle(lwi(l)).contactG2s(particle(lwi(l)).z)=contactG2;
%         particle(lwi(l)).contactIs(particle(lwi(l)).z)=cI;
%         particle(lwi(l)).neighbours(particle(lwi(l)).z) = -1; %the wall is now noted as a neigbour in the particle l datastructure
%         particle(lwi(l)).betas(particle(lwi(l)).z) = pi; %the contact angle to the wall is now noted in the particle l datastructure
%     end
% end
% 
% for l = 1:length(rwi)
%     x = circs(rwi(l),2);
%     y = circs(rwi(l),1);
%     r = round(circs(rwi(l),3));
%     
%     contactX = x+(r-CR);
%     contactY = y;
%     
%     contactImg = im2double(imcrop(Gimg,[contactX-CR contactY-CR CR*2 CR*2]));
%     contactImg = contactImg.*mask;
%     
%     [gx,gy] = gradient(contactImg);
%     g2 = (gx.^2 + gy.^2);
%     contactG2 = sum(sum(g2));
%     
%     if(contactG2 > contactG2Threshold)
%         cI = sum(sum(contactImg));
%         %if(cI > contactIThreshold)
%         %this is a valid contact, remember it
%         if(verbose)
%             text(contactX,contactY,num2str(contactG2),'Color','w');
%             viscircles([contactX; contactY]', CR,'EdgeColor','w');
%         end
%         particle(rwi(l)).z= particle(rwi(l)).z +1; %increase coordination number
%         particle(rwi(l)).contactG2s(particle(rwi(l)).z)=contactG2;
%         particle(rwi(l)).contactIs(particle(rwi(l)).z)=cI;
%         particle(rwi(l)).neighbours(particle(rwi(l)).z) = -2; %the wall is now noted as a neigbour in the particle l datastructure
%         particle(rwi(l)).betas(particle(rwi(l)).z) = 0; %the contact angle to the wall is now noted in the particle l datastructure
%     end
% end
% 
% for l = 1:length(twi)
%     x = circs(twi(l),2);
%     y = circs(twi(l),1);
%     r = round(circs(twi(l),3));
%     
%     contactX = x;
%     contactY = y-(r-CR);
%     
%     contactImg = im2double(imcrop(Gimg,[contactX-CR contactY-CR CR*2 CR*2]));
%     contactImg = contactImg.*mask;
%     
%     [gx,gy] = gradient(contactImg);
%     g2 = (gx.^2 + gy.^2);
%     contactG2 = sum(sum(g2));
%     
%     if(contactG2 > contactG2Threshold)
%         cI = sum(sum(contactImg));
%         %if(cI > contactIThreshold)
%         %this is a valid contact, remember it
%         if(verbose)
%             text(contactX,contactY,num2str(contactG2),'Color','w');
%             viscircles([contactX; contactY]', CR,'EdgeColor','w');
%         end
%         particle(twi(l)).z= particle(twi(l)).z +1; %increase coordination number
%         particle(twi(l)).contactG2s(particle(twi(l)).z)=contactG2;
%         particle(twi(l)).contactIs(particle(twi(l)).z)=cI;
%         particle(twi(l)).neighbours(particle(twi(l)).z) = -3; %the wall is now noted as a neigbour in the particle l datastructure
%         particle(twi(l)).betas(particle(twi(l)).z) = -pi/2; %the contact angle to the wall is now noted in the particle l datastructure
%     end
% end
% 
% for l = 1:length(bwi)
%     x = circs(bwi(l),2);
%     y = circs(bwi(l),1);
%     r = round(circs(bwi(l),3));
%     
%     contactX = x
%     contactY = y+(r-CR)
%     
%     contactImg = (imcrop(Gimg,[contactX-CR contactY-CR CR*2 CR*2]));
%     contactImg = contactImg.*mask;
%     %imshow(contactImg)
%     
%     [gx,gy] = gradient(contactImg);
%     g2 = (gx.^2 + gy.^2);
%     contactG2 = sum(sum(g2))
%     
%     if(contactG2 > contactG2Threshold)
%         cI = sum(sum(contactImg));
%         %if(cI > contactIThreshold)
%         %this is a valid contact, remember it
%         if(verbose)
%             text(contactX,contactY,num2str(contactG2),'Color','w');
%             viscircles([contactX; contactY]', CR,'EdgeColor','w');
%         end
%         particle(bwi(l)).z= particle(bwi(l)).z +1; %increase coordination number
%         particle(bwi(l)).contactG2s(particle(bwi(l)).z)=contactG2;
%         particle(bwi(l)).contactIs(particle(bwi(l)).z)=cI;
%         particle(bwi(l)).neighbours(particle(bwi(l)).z) = -4; %the wall is now noted as a neigbour in the particle l datastructure
%         particle(bwi(l)).betas(particle(bwi(l)).z) = pi/2; %the contact angle to the wall is now noted in the particle l datastructure
%     end
% end

end