%% 2D image %%
%% Clear ALL
clc; clear all; close all;
%% Add path
addpath('../2D/exampleImages/'); 
addpath('../2D/lib/');  
addpath('../2D/plots/');

%% Load images
foldername = './exampleImages/';
foldername_plots1 = 'plots/accumulativeMap/mat/';
foldername_plots2 = 'plots/accumulativeMap/png/';
foldername_plots3 = 'plots/extractCenterline/mat/';
foldername_plots4 = 'plots/extractCenterline/png/';
foldername_plots5 = 'plots/extractCenterline/fig/';
filename = 'retina';
namefig = 'retinaFig';
ext = '.png';
ext2 = '.mat';
extfig = '.fig';
I = double(imread(fullfile(foldername,[filename ext])));

%% Normalize & Complement
I = (I - min(I(:)))/(max(I(:))-min(I(:))); % Normalise orginal image
%I = imcomplement(I); % for some images, when they have white background and black object

%% Parameters setting
t1 = 15000; % optimization parameter
t2 = 2000; % optimization parameter
epsilon = 3; % tunable parameter

%% Call function
[Q,A] = CenterlineExtraction2D(I,t1,t2,epsilon,1); 

%% Save and Write
%Save as 'png' image
imwrite(A,(fullfile(foldername_plots2,[filename ext])),'XResolution',300);
imwrite(Q,(fullfile(foldername_plots4,[filename ext])),'XResolution',300);
%Save as 'mat' file
save((fullfile(foldername_plots1,[filename ext2])),'A');
save((fullfile(foldername_plots3,[filename ext2])),'Q');
%Print figure
saveas(gcf,fullfile(foldername_plots5,[namefig extfig]));
print(fullfile(foldername_plots5,namefig),'-r300','-dpng');