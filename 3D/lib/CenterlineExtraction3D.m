function [Q,A] = CenterlineExtraction3D(I,t1,t2,epsilon,Threshold)
%% Extraction of the centerline in 3D images
%
%
% Details   This function contains the two stages to extract the centerline 
%           Stage1:enhancement of the centerline.
%           Stage2: Linked the point and extract the centerline
%
%           as implemented for:
%
%          'Sequential Graph-based Extraction of Curvilinear Structures'
%           by Shuaa S. Alharbi, Chris G. Willcocks, Philip T. G. Jackson, Haifa F. Alhasson
%           and Boguslaw Obara in Signal, Image and Video Processing Journal (2019).
% Usage     [Q,A] = CenterlineExtraction2D(I,t1,t2,epsilon,0)
%           [Q,A] = CenterlineExtraction2D(I,t1,t2,epsilon,1)
% Inputs    I - orginal 3D image
%           t1 - optimisation parameter t1 for step1
%           t2 - optimisation parameter t2 for step2
%           epsilon - tunable parameter
%           Threshold - 0 (without threshold 'slower result', need to increase t2 to very lage value i.e (50000)) 
%                     - 1 (with threshold to reduce the set of nodes 'faster results')
%                      
% Outputs   Q - a 3D matrix that define the detected curve
%           A - a 3D accumulative map calculate by ccumulating every prominent peak
%               position and increment it location each time detected.
%
% Copyright 2019 Shuaa S. Alharbi, Durham University, UK
% 
% License   See included <a href="./LICENSE/">file</a> or visit
%           <a href="">The GitHub
%              Repository</a>
%
% See also main.m
%
%
%
%
%********************************************** Stage 1 **********************************************%
%
%
%
%
% Display parameters
kd = 200; k = kd; 
count = 0; % Counter for display
% Input parameters
[d1,d2,d3] = size(I);
A = zeros([d1,d2,d3]); % Output image - Accumulative map
temp = max(d1,d2); % Maximum length
rmax = max(temp,d3);
% First display
figure; subplot(1,2,1); imagesc(max(I,[],3)); colormap gray; axis equal; axis tight; axis off;
title('Orginal Image');
ha = subplot(1,2,2); imagesc(max(A,[],3)); colormap gray; axis equal; axis tight; axis off; 

%% Loop to sample image
for f = 1 : t1  
    count = count + 1; 
    %% Random Point on the image 
    p = [1 + (d1-1).*rand(1), 1 + (d2-1).*rand(1),1 + (d3-1).*rand(1)];

    %% Correct normal distribution because of randn
    d = randn(1,3); % directional vector (plane direction)
    d = bsxfun(@rdivide,d,sqrt(sum(d.^2,2)));
    len = round(rand(1)*rmax); % random length

%% Sample point
i = [1 0 0]; % Base vector
a = cross(i,d);
b = cross(a,d);
hiProm = 0; % Highest prominence
hiLoc = 0; % Location of the highest prominence
hiTheta = 1; % Theta value that used to find the maximum prominent peak 
    for theta = 1 : 360
         disp(theta)
         dir = len*(sin(theta)*a+cos(theta)*b);
         dir = bsxfun(@rdivide,dir,sqrt(sum(dir.^2,2))); %Normlize vector 
  % Samples points
         irow3 = interp3(I, linspace(p(:,2), p(:,2)+dir(:,2)*len, len), ...
                            linspace(p(:,1), p(:,1)+dir(:,1)*len, len), ...
                            linspace(p(:,3), p(:,3)+dir(:,3)*len, len));
        irow3(isnan(irow3)) = [];
        if isempty(irow3) || length(irow3)<3 ; continue; end
        
    %% Peak finding based on the prominence
     [pks, locs, ~, proms] = findpeaks(irow3);
      if ~isempty(pks)
        [prom, ind] = max(proms); 
        locs = locs(ind); % Location of the peak on irow3 (sample points)
        if prom > hiProm
            hiProm = prom; % Update the prom with the maximum one
            hiLoc = locs;
            hiTheta = theta;
        end
      end
    end % End loop of (theta)
    %% Calculate the accumulative map A
    dir = (a*sin(hiTheta) + b*cos(hiTheta));
    dir = bsxfun(@rdivide,dir,sqrt(sum(dir.^2,2)));
    pos3d = p + hiLoc*dir; % Position of the maximum prominent peak
    
     if  pos3d(1)>size(I,1) || pos3d(2)>size(I,2) || pos3d(3)>size(I,3) || ...
         pos3d(1)<1 || pos3d(2)<1 || pos3d(3)<1
         continue
     end
    pos3d = round(pos3d);
    idx2 = sub2ind(size(A),pos3d(1),pos3d(2),pos3d(3));
    A(idx2) = A(idx2) + 1;
    
    %% Plot the accumulative map
    if (k-count)<0 % displays every k-iter
        k = count+kd;
        imagesc(max(A,[],3),'Parent',ha); colormap gray; axis equal; axis tight; axis off;
        drawnow; title('Probability Map');fprintf('number of iteration\n');
        num2str(count)
        pause(0.2)
    end % End(display)
end % End loop
A = (A - min(A(:)))/(max(A(:))-min(A(:))); % Normalise accumulative map
%
%
%
%
%********************************************** Stage 2 **********************************************%
%
%
%
%
% Input parameters
Q = zeros(size(I)); % matrix to store the best path each time extract
im = max(I,[],3);
% First display
figure,imshow(im);colormap gray; axis equal; axis tight; axis off; hold on;
%% Convert the orginal image to grayscale image
I = mat2gray(I);

%% Convert 3D accumulation map A to a graph (Create sparse matrix G - graph representation) 
fprintf('**************** Convert 3D accumulation map to graph ****************\n');
C = 1./log(A+1); % Cost function
G = im2graph3D(C); % Convet to graph

%% Caculate a set of the shortest path and choose the best one
M = graythresh(A); % Threshold by Otsu's method
ind = find(A>M);
for j = 1 : epsilon
    bestPath = [];
    bestPathSum = 0;
    for i = 1 : t2
        disp(i)
       % Choose points randomly with range not outside(ind)
        if Threshold == 0 % without threshold
            endpoint1 = randi(numel(A));
            endpoint2 = randi(numel(A));
        elseif Threshold == 1 % with threshold
            endpoint1 = ind(randi(numel(ind))); 
            endpoint2 = ind(randi(numel(ind)));
        else
            fprintf('**** Choose (1) to get result faster OR choose (0) to get result without threshold ****\n')
        end
        % Find shortest path
        [~, path, ~] = graphshortestpath(G, endpoint1, endpoint2); % Default uses Dijkstra
        
        %% Chosen of the best path
        path(Q(path) == 1) = inf; 
        endp = find(path == inf, 1);
        path(endp:length(path)) = [];
       % Linear search about the largest intensity path(sum along its profile)
        pathSum = sum(A(path));
        if (pathSum > bestPathSum)
            bestPathSum = pathSum;
            bestPath = path;
        end
    end
    
   %% Display the sequential centerline extraction
    Q(bestPath) = 1;
    [x,y,z] = ind2sub(size(I), bestPath);
    plot3(y(:), x(:), z(:),'LineWidth', 1);
    drawnow;
end % End loop
end %End function