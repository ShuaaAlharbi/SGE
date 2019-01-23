function G = im2graph3D(C)
% PURPOSE: 
% To represent a 3D image as  undirect connected weighted graph G(N,E).
% Each Voxel corresponds to a node N. 
% Graph edges E related to the node N are defined by 26 neighbors(6 faces, 12 edges and 8 courner).
% The weight of the edges are estemated using cost function C.
%
% Usage:
% G = im2graph3D(A)
%
%INPUT: 
% C - Cost matrix MxNxL       
%
%OUTPUT: 
% G - A sparse matrix that represents a graph.
% Nonzero entries in matrix G represent the weights of the edges.
% 
% AUTHOR:
% Shuaa S. Alharbi
%
% VERSION:
% 0.1 - 14/08/2016 First implementation

%% Image indexes
idx = find(ones(size(C)));

%% Index offsets (6 faces & 12 edges & 8 corners)
[m,n,L] = size(C);
%---6 faces--- %
E   =  m; % East  
B   =  1; % Bottom(South)
W   = -m; % West
T   = -1; % Top(North)
Page = m * n; % third dimention
%---12 edges--- %
SE  =  m + 1; % Southeast of the current idx
SW  = -m + 1; % Southwest of the current idx
NW  = -m - 1; % Northwest of the current idx
NE  =  m - 1; % Northeast of the current idx
EF  = -Page + m; % East(east of the front)
BF  =  Page - 1; % Bottom(south of the front)
WF  = -Page - m; % West(of the front)
TF  = -Page - 1;  % Top(North of the front)
EBa = Page + m; % East(east of the back)
BBa = Page + 1; % Bottom(south of the back)
WBa = Page - m; % West(of the back)
TBa = Page -1;  % Top(North of the back)
%---8 corners--- %
SEF = -Page + SE; % Southeast(of the front)
SWF = -Page + SW; % Southwest(of the front)
NWF = -Page + NW; % Northwest(of the front)
NEF = -Page + NE; % Notheast(of the front)
SEBa = Page + SE; % Southeast(of the back)
SWBa = Page + SW; % Southwest(of the back)
NWBa = Page + NW; % Northwest(of the back)
NEBa = Page + NE; % Notheast(of the back)

%% Indexes
% 6 faces
idxE  = idx + E;
idxB  = idx + B;
idxW  = idx + W;
idxT  = idx + T;
idxF  = idx - Page;
idxBa = idx + Page;
% 12 edges
idxSE  = idx + SE;
idxSW  = idx + SW;
idxNW  = idx + NW;
idxNE  = idx + NE;
idxEF = idx + EF;
idxBF = idx - BF;
idxWF  = idx + WF;
idxTF  = idx + TF;
idxEBa  = idx + EBa;
idxBBa = idx + BBa;
idxWBa  = idx + WBa;
idxTBa = idx + TBa;
% 8 courner
idxSEF = idx + SEF;
idxSWF = idx + SWF;
idxNWF = idx + NWF;
idxNEF = idx + NEF;
idxSEBa = idx + SEBa;
idxSWBa = idx + SWBa;
idxNWBa = idx + NWBa;
idxNEBa = idx + NEBa;

%% Compute the values of all the north, south, ..., neighbors
pE  = [idx, idxE]; % pairs 6 faces
pB = [idx, idxB];
pW  = [idx, idxW];
pT = [idx, idxT];
pF  = [idx, idxF];
pBa = [idx, idxBa];
%
pSE  = [idx, idxSE]; % pairs 12 edges
pSW = [idx, idxSW];
pNW  = [idx, idxNW];
pNE = [idx, idxNE];
pEF  = [idx, idxEF];
pBF = [idx, idxBF];
pWF = [idx, idxWF];
pTF = [idx, idxTF];
pEBa  = [idx, idxEBa];
pBBa = [idx, idxBBa];
pWBa = [idx, idxWBa];
pTBa = [idx, idxTBa];
%
pSEF = [idx, idxSEF]; % pairs 8 corners
pSWF = [idx, idxSWF];
pNWF = [idx, idxNWF];
pNEF = [idx, idxNEF];
pSEBa = [idx, idxSEBa];
pSWBa = [idx, idxSWBa];
pNWBa = [idx, idxNWBa];
pNEBa = [idx, idxNEBa];

%% Edages
%-------6 faces-------%
%East neighbors of the current pixels
dE = ismember(idxE, idx);
%disconnect the last cols(No east neighbors)
for j = 1 : L 
    for i = 1 : m
        didx = sub2ind([m,n,L],i,n,j);
        dE(didx) = 0; 
    end
end
pE  = pE(dE,:);
%Bottom neighbors of the current pixels
dB = ismember(idxB, idx);
%disconnect the last rows(No bottom(south)neighbors)
for j = 1 : L
    for i = 1 : n
        didx = sub2ind([m,n,L],m,i,j);
        dB(didx) = 0;
    end
end
pB  = pB(dB,:);
%West neighbors of the current pixels
dW = ismember(idxW, idx);
%disconnect the first cols(No west neighbors)
for j = 1 : L
    for i = 1 : m
        didx = sub2ind([m,n,L],i,1,j);
        dW(didx) = 0;
    end
end
pW  = pW(dW,:);
%Top neighbors of the current pixels
dT = ismember(idxT, idx);
%disconnect the first rows(No top neighbors)
for j = 1 : L
    for i = 1 : n
        didx = sub2ind([m,n,L],1,i,j);
        dT(didx) = 0;
    end
end
pT  = pT(dT,:);
%Front neighbors of the current pixels
dF = ismember(idxF, idx);
%disconnect the first page(No front neighbors)
for j = 1 : m
    for i = 1 : n
        didx = sub2ind([m,n,L],j,i,1);
        dF(didx) = 0;
    end
end
pF  = pF(dF,:);
%Back neighbors of the current pixels
dBa = ismember(idxBa, idx);
%disconnect the last page(No back neighbors)
for j = 1 : m
    for i = 1 : n
        didx = sub2ind([m,n,L],j,i,L);
        dBa(didx) = 0;
    end
end
pBa  = pBa(dBa,:);

%%
%-------12 edges-------%
%SouthEast neighbors of the current pixels
dSE = ismember(idxSE,idx);
%disconnect the last rows & cols(No southeast)
for j = 1 : L 
    for i = 1 : n
        didx = sub2ind([m,n,L],m,i,j);
        dSE(didx) = 0;
    end
end
for j = 1 : L 
    for i = 1 : m
    didx = sub2ind([m,n,L],i,n,j);
    dSE(didx) = 0;
    end
end
pSE = pSE(dSE,:);
%SouthWest neighbors of the current pixels
dSW = ismember(idxSW,idx);
%disconnect the last rows & frist cols(No southwest)
for j = 1 : L 
    for i = 1 : n
        didx = sub2ind([m,n,L],m,i,j);
        dSW(didx) = 0;
    end
end
for j = 1 : L 
    for i = 1 : m
        didx = sub2ind([m,n,L],i,1,j);
        dSW(didx) = 0;
    end
end
pSW = pSW(dSW,:);
%NorthWest neighbors of the current pixels
dNW = ismember(idxNW,idx);
%disconnect the first rows & first cols(No northwest)
for j = 1 : L 
    for i = 1 : n
        didx = sub2ind([m,n,L],1,i,j);
        dNW(didx) = 0;
    end
end
for j = 1 : L 
    for i = 1 : m
        didx = sub2ind([m,n,L],i,1,j);
        dNW(didx) = 0;
    end
end
pNW = pNW(dNW,:);
%NorthEast neighbors of the current pixels
dNE = ismember(idxNE,idx);
%disconnect the first rows & last cols(No northeast)
for j = 1 : L 
    for i = 1 : n
    didx = sub2ind([m,n,L],1,i,j);
    dNE(didx) = 0;
    end
end
for j = 1 : L
    for i = 1 : m
        didx = sub2ind([m,n,L],i,n,j);
        dNE(didx) = 0;
    end
end
pNE = pNE(dNE,:);
%East neighbors of the front pixels
dEF = ismember(idxEF, idx); 
for j = 1 : L 
    for i = 1 : m
        didx = sub2ind([m,n,L],i,n,j);
        dEF(didx) = 0;
    end
end
pEF = pEF(dEF,:);
%Bottom neighbors of the front pixels
dBF = ismember(idxBF, idx); 
for j = 1 : L
    for i = 1 : n
        didx = sub2ind([m,n,L],m,i,j);
        dBF(didx) = 0;
    end
end
pBF = pBF(dBF,:);
%West neighbors of the front pixels
dWF = ismember(idxWF, idx); 
for j = 1 : L
    for i = 1 : m
        didx = sub2ind([m,n,L],i,1,j);
        dWF(didx) = 0;
    end
end
pWF = pWF(dWF,:);
%Top neighbors of the front pixels
dTF = ismember(idxTF, idx); 
for j = 1 : L 
    for i = 1 : n
        didx = sub2ind([m,n,L],1,i,j);
        dTF(didx) = 0;
    end
end
pTF = pTF(dTF,:);
%East neighbors of the back pixels
dEBa = ismember(idxEBa, idx); 
for j = 1 : L 
    for i = 1 : m
        didx = sub2ind([m,n,L],i,n,j);
        dEBa(didx) = 0;
    end
end
pEBa = pEBa(dEBa,:);
%Bottom neighbors of the back pixels
dBBa = ismember(idxBBa, idx); 
for j = 1 : L
    for i = 1 : n
        didx = sub2ind([m,n,L],m,i,j);
        dBBa(didx) = 0;
    end
end
pBBa = pBBa(dBBa,:);
%West neighbors of the back pixels
dWBa = ismember(idxWBa, idx); 
for j = 1 : L 
    for i = 1 : m
        didx = sub2ind([m,n,L],i,1,j);
        dWBa(didx) = 0;
    end
end
pWBa = pWBa(dWBa,:);
%Top neighbors of the back pixels
dTBa = ismember(idxTBa, idx); 
for j = 1 : L
    for i = 1 : n
        didx = sub2ind([m,n,L],1,i,j);
        dTBa(didx) = 0;
    end
end
pTBa = pTBa(dTBa,:);
%%
%-------8 corners-------%
%SouthEast neighbors of the front pixels
dSEF = ismember(idxSEF, idx); 
for j = 1 : L 
    didx = sub2ind([m,n,L],m,n,j);
    dSEF(didx) = 0;
end
pSEF = pSEF(dSEF,:);
%SouthWest neighbors of the front pixels
dSWF = ismember(idxSWF, idx); 
for j = 1 : L 
    didx = sub2ind([m,n,L],m,1,j);
    dSWF(didx) = 0;
end
pSWF = pSWF(dSWF,:);
%NorthWest neighbors of the front pixels
dNWF = ismember(idxNWF, idx); 
for j = 1 : L
    didx = sub2ind([m,n,L],1,1,j);
    dNWF(didx) = 0;
end
pNWF = pNWF(dNWF,:);
%NorthEast neighbors of the front pixels
dNEF = ismember(idxNEF, idx); 
for j = 1 : L 
    didx = sub2ind([m,n,L],1,n,j);
    dNEF(didx) = 0;
end
pNEF = pNEF(dNEF,:);
%SouthEast neighbors of the back pixels
dSEBa = ismember(idxSEBa, idx); 
for j = 1 : L 
    didx = sub2ind([m,n,L],m,n,j);
    dSEBa(didx) = 0;
end
pSEBa = pSEBa(dSEBa,:);
%SouthWest neighbors of the back pixels
dSWBa = ismember(idxSWBa, idx); 
for j = 1 : L 
    didx = sub2ind([m,n,L],m,1,j);
    dSWBa(didx) = 0;
end
pSWBa = pSWBa(dSWBa,:);
%NorthWest neighbors of the back pixels
dNWBa = ismember(idxNWBa, idx); 
for j = 1 : L
    didx = sub2ind([m,n,L],1,1,j);
    dNWBa(didx) = 0;
end
pNWBa = pNWBa(dNWBa,:);
%%NorthEast neighbors of the back pixels
dNEBa = ismember(idxNEBa, idx); 
for j = 1 : L
    didx = sub2ind([m,n,L],1,n,j);
    dNEBa(didx) = 0;
end
pNEBa = pNEBa(dNEBa,:);
%% Weights
%-------6 faces-------%
nE = C(pE( :,2));
nB = C(pB(:,2));
nW = C(pW( :,2));
nT = C(pT(:,2));
nF = C(pF( :,2));
nBa = C(pBa(:,2));
%-------12 edges-------%
nSE = C(pSE( :,2));
nSW = C(pSW(:,2));
nNW = C(pNW( :,2));
nNE = C(pNE(:,2));
nEE1 = C(pEF( :,2));
nBF = C(pBF(:,2));
nWF = C(pWF(:,2));
nTF = C(pTF(:,2));
nEBa = C(pEBa( :,2));
nBBa = C(pBBa(:,2));
nWBa = C(pWBa(:,2));
nTBa = C(pTBa(:,2));
%-------8 corners-------%
nSEF = C(pSEF(:,2));
nSWF = C(pSWF(:,2));
nNWF = C(pNWF( :,2));
nNEF = C(pNEF(:,2));
nSEBa = C(pSEBa( :,2));
nSWBa = C(pSWBa(:,2));
nNWBa = C(pNWBa(:,2));
nNEBa = C(pNEBa(:,2));
%-----------------Graph representation-------------------%
Weights = [nE; nB; nW; nT; nF; nBa; nSE; nSW; nNW; nNE; nEE1; nBF; nWF; nTF; nEBa; nBBa; nWBa; nTBa; nSEF; nSWF; nNWF; nNEF; nSEBa; nSWBa; nNWBa; nNEBa]';
Edges = [pE; pB; pW; pT; pF; pBa; pSE; pSW; pNW; pNE; pEF; pBF; pWF; pTF; pEBa; pBBa; pWBa; pTBa; pSEF; pSWF; pNWF; pNEF; pSEBa; pSWBa; pNWBa; pNEBa]';
G = sparse(Edges(1,:),Edges(2,:),double(Weights));%sparse matrix
%% End
end