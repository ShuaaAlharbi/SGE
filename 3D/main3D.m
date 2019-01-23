%% 3D image %%
%% Clear ALL
clc; clear all; close all;
%% Add path
addpath('../3D/exampleImages/');
addpath('../3D/lib/');
addpath('../3D/plots/');

%% Load images
foldername = './exampleImages/';
foldername_plots1 = 'plots/accumulativeMap/mat/';
foldername_plots2 = 'plots/accumulativeMap/png/';
foldername_plots3 = 'plots/extractCenterline/mat/';
foldername_plots4 = 'plots/extractCenterline/png/';
foldername_plots5 = 'plots/extractCenterline/fig/';
filename = 'neuralNetwork';
namefig = 'neuralNetworkFig';
ext = '.tif'; 
ext2 = '.mat';
extR = '.png';
extfig = '.png';
I = double(load3D(fullfile(foldername,[filename ext])));

%% Normalize & Complement
I = (I - min(I(:)))/(max(I(:))-min(I(:))); % Normalise orginal image
%I = imcomplement(I); % for some images, when they have white background and black object

%% Parameters setting
t1 = 20000; % optimization parameter
t2 = 3000; % optimization parameter
epsilon = 6; % tunable parameter

%% Call function
[Q,A] = CenterlineExtraction3D(I,t1,t2,epsilon,1); 

%% Save and Write
%Save as 'png' image (Top-down View)
a = max(A,[],3); 
b = max(Q,[],3); 
imwrite(a,(fullfile(foldername_plots2,[filename extR])),'XResolution',300);
imwrite(b,(fullfile(foldername_plots4,[filename extR])),'XResolution',300);
%Save as 'mat' file
save((fullfile(foldername_plots1,[filename ext2])),'A');
save((fullfile(foldername_plots3,[filename ext2])),'Q');
%Print figure (3D view)
figure(2),imshow(I),hold on,
[xi,yi,zi] = ind2sub(size(Q),find(Q>0));
plot(yi,xi,zi,'.r','color',[77 175 74]/255,'MarkerSize',3);
saveas(gcf,fullfile(foldername_plots5,[namefig extfig]));
print(fullfile(foldername_plots5,namefig),'-r300','-dpng');