function [Q,A] = CenterlineExtraction2D(I,t1,t2,epsilon,Threshold)
%% Extraction of the centerline based on an iterative graph-based optimization
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
% Usage     [Q,A] = CenterlineExtraction2D(I,0)
%           [Q,A] = CenterlineExtraction2D(I,1)
% Inputs    I - orginal 2D image
%           Threshold - 0 (without threshold 'slower result', need to increase t2 to very lage value i.e (50000)) 
%                     - 1 (with threshold to reduce the set of nodes 'faster results')
%                      
% Outputs   Q - a 2D matrix that define the detected curve
%           A - a 2D accumulative map calculate by ccumulating every prominent peak
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
count = 0; 
% Input parameters
s = 1.0; % pixel sampling (step size)
[m,n] = size(I);
tmax = max(m,n); % Maximum line length
A = zeros([m,n]);% Output image - Accumulative map
% First display
figure; subplot(1,2,1); imagesc(I); colormap gray; axis equal; axis tight; axis off;
title('Orginal Image');
ha = subplot(1,2,2); imagesc(A); colormap gray; axis equal; axis tight; axis off;

%% Loop to sample image using a line
for i = 1 : t1 %loop
    count = count + 1;  
    % integrate over all possible line segment "positions"
    p = [1 + (m-1).*rand(1),1 + (n-1).*rand(1)];
    
    %% Correct normal distribution because of randn
    d = randn(1,2); % directional vector (line direction) 
    d = bsxfun(@rdivide,d,sqrt(sum(d.^2,2))); 
    len = 0:s:rand(1)*tmax; % line length ( l in the paper)
    
    %% Samples points define by line equation r = p + dt
    r = [];
    r(1,:) = repmat(p(1),[length(len),1]) + repmat(d(1),[length(len),1]).*len';
    r(2,:) = repmat(p(2),[length(len),1]) + repmat(d(2),[length(len),1]).*len';
    
    %% Remove sample points that go outside the image
    keep = (r(1,:) >= 1) & (r(1,:) <= m);
    keep = keep & (r(2,:) >= 1) & (r(2,:) <= n);
    r = r(:, keep);

    %% Plot sample points
%     figure; imagesc(I); colormap gray; axis equal; axis tight; axis off;
%     hold on;
%     plot(r(2,:),r(1,:),'g-','LineWidth',2.7);
%     hold on;
%     plot(p(2),p(1),'r*','LineWidth',2.0);
    %% Interpolation for 2-D gridded data in meshgrid format
    if length(r)<3; continue; end
    irow = interp2(I,r(2,:),r(1,:));

    %% Peak finding in the line 1D
    [pks, locs, ~, proms] = findpeaks(irow);
    
    %% Calculate the accumulative map A
    if ~isempty(pks)
        [~, idx] = max(proms); % find peak based on prominent 
        locs = locs(idx);
        x = round(r(1,locs)); % double pixel position to integer position
        y = round(r(2,locs));
        idx = sub2ind([m,n],x,y);
        A(idx) = A(idx) + 1;
    else
        continue
    end
    %% Plot the accumulative map A
    if (k-count)<0 % displays every k-iterations
        k = count + kd;
        imagesc(A,'Parent',ha); colormap gray; axis equal; axis tight; axis off;
        drawnow; title('Accumulative Map'); fprintf('number of iteration\n');
        num2str(count)
        pause(0.2)
    end % End of peak detecting
end % End of loop
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
figure, imshow(I); hold on; % First display

%% Convert 2D accumulation map A to a graph (Create sparse matrix G - graph representation) 
fprintf('**************** Convert 2D accumulation map to graph ****************\n');
C = 1./A; % Cost function
G = im2graph2D(C); % Convert to graph

%% Caculate a set of the shortest path and choose the best one
M = graythresh(A); % Threshold by Otsu's method
ind = find(A>M);
for j = 1 : epsilon % Loop for sequential extraction of the centerline 
    bestPath = [];
    bestPathSum = 0; % Sum of the intensity (Wight to choose the best path)
    for i = 1 : t2 % Loop to produce a set of the shortest paths
        disp(i)
        % Choose the endpoints randomly with range not outside(ind)
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
        
%% plot the shortest path and endpoints
%         [x1,y1]=ind2sub(size(A),endpoint1);
%         [x2,y2]=ind2sub(size(A),endpoint2);
%         [xp,yp] = ind2sub(size(A),path);
%         figure,imshow(I)
%         hold on
%         plot(yp,xp,'-', 'Color', [228,26,28]/255, 'LineWidth',2.0)
%         plot(y1,x1,'x', 'Color', [55,126,184]/255, 'LineWidth',2.0)
%         plot(y2,x2,'x', 'Color', [55,126,184]/255, 'LineWidth',2.0)

        %% Chosen of the best path
        path(Q(path) == 1) = inf; % To avoide choose the same path againe
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
      [y,x] = ind2sub(size(I), bestPath);
      plot(x(:), y(:),'LineWidth', 1);
      drawnow;                                                                           
end % End loop
end % End function 