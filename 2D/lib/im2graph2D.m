function G = im2graph2D(C)
% PURPOSE: 
% To represent a 2D image as undirect connected weighted graph G(N,E).
% Each pixel corresponds to a node N.
% Graph edges E related to the node N are defined by 8-connected
% neighbourhood of pixels.
% The weight of the edges are estemated using cost function C.
%
% Usage:
% G = im2graph(A)
%
% INPUT: 
% A - Accumulative matrix of dimensions MxN
%       
% OUTPUT: 
% C -  Cost matrix.
% Nonzero entries in matrix G represent the weights of the edges.
% 
% AUTHOR:
% Shuaa S. Alharbi
%
% VERSION:
% 0.1 - 14/08/2016 First implementation

%% Image indexes
idx = find(ones(size(C)));
[m,n] = size(C);

%% Index offsets
% East         M
% Southeast    M + 1
% South        1
% Southwest   -M + 1
% West        -M
% Northwest   -M - 1
% North       -1
% Northeast    M - 1
E   =  m;
SE  =  m + 1;
S   =  1;
SW  = -m + 1;
W   = -m;
NW  = -m - 1;
N   = -1;
NE  =  m - 1;
%% Indexes
idxE  = idx + E;
idxSE = idx + SE;
idxS  = idx + S;
idxSW = idx + SW;
idxW  = idx + W;
idxNW = idx + NW;
idxN  = idx + N;
idxNE = idx + NE;
%% Compute the values of all the north, south, ..., neighbors
pE  = [idx, idxE ]; % pairs
pSE = [idx, idxSE];
pS  = [idx, idxS ];
pSW = [idx, idxSW];
pW  = [idx, idxW ];
pNW = [idx, idxNW];
pN  = [idx, idxN ];
pNE = [idx, idxNE];

%% Edges
%East neighbors of the current pixels
pE  = pE( ismember(idxE, idx),:);
%SouthEast neighbors of the current pixels
dSE = ismember(idxSE,idx);
%disconnect the last row(No southeast)
for i = 1 : n
    didx = sub2ind([m,n],m,i);
    dSE(didx) = 0;
end
pSE = pSE(dSE,:);
%South neighbors of the current pixels
dS = ismember(idxS, idx);
%disconnect the last row(No south)
for i = 1 : n
    didx = sub2ind([m,n],m,i);
    dS(didx) = 0;
end
pS  = pS(dS,:);
%SouthWest neighbors of the current pixels
dSW = ismember(idxSW,idx);
%disconnect the last row(No southwest)
for i = 1 : n
    didx = sub2ind([m,n],m,i);
    dSW(didx) = 0;
end
pSW = pSW(dSW,:);
%West neighbors of the current pixels
pW  = pW( ismember(idxW, idx),:);
%NorthWest neighbors of the current pixels
dNW = ismember(idxNW,idx);
%disconnect the first row(No Northwest)
for i = 3 : n
    didx = sub2ind([m,n],1,i);
    dNW(didx) = 0;
end
pNW = pNW(dNW,:);
%North neighbors of the current pixels
dN = ismember(idxN, idx);
%disconnect the first row(No North)
for i = 2 : n
    didx = sub2ind([m,n],1,i);
    dN(didx) = 0;
end
pN  = pN(dN,:);
%NorthEast neighbors of the current pixels
dNE = ismember(idxNE,idx);
%disconnect the first row(No Northeast)
for i = 1 : n
    didx = sub2ind([m,n],1,i);
    dNE(didx) = 0;
end
pNE = pNE(dNE,:);

%% Weights
nE  = C(pE( :,2));
nSE = C(pSE(:,2));
nS  = C(pS( :,2));
nSW = C(pSW(:,2));
nW  = C(pW( :,2));
nNW = C(pNW(:,2));
nN  = C(pN( :,2)); 
nNE = C(pNE(:,2));
%-----------------Graph representation-------------------%
Weights = [nE; nSE; nS; nSW; nW; nNW; nN; nNE]';
Edges = [pE; pSE; pS; pSW; pW; pNW; pN; pNE]';
G = sparse(Edges(1,:),Edges(2,:),double(Weights));%sparse matrix 
end % End function