% makefile for the complete GSH circle for a particular model
clear;
close all;
clc;

HOME = pwd;
addpath([HOME '/Data']);
addpath([HOME '/Tools']);

% Model
% Load previous saved model

%model_name = 'Crust01_crust';
%load(model_name);

% Construct new model
inputModel    

%%%%%%%%%%%%%%%%%%% Computation area %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Part that can be modified %%%%%%%%%%%%%%%%%%%%%%%

latLim =    [-89.5 89.5 0.1];  % [deg] min latitude, max latitude, resolution latitude (preferable similar to latitude)
lonLim =    [0 360 0.1];%[-180 180 1] [deg] min longitude, max longitude, resolution longitude (preferable similar to latitude)
height =    225000.0; % height of computation above spheroid 225000
SHbounds =  [3 110]; % Truncation settings: lower limit, upper limit SH-coefficients used

%%%%%%%%%%%%%% Part that can be modified %%%%%%%%%%%%%%%%%%%%%%%

%% Global Spherical Harmonic Analysis 

% % tic;
% % [V] = model_SH_analysis(Model);
% % toc
% % 
% save(['Results/' Model.name '.mat'],'V')

%% Global Spherical Harmonic Synthesis

V1=readtable("V1.xlsx");
V1=table2array(V1);
V1=V1(:,1:4);
V1 = sortrows(V1,2);

V1=[0 0 1 0; V1];
%save("V1.mat", "V1");
%V1=load("V1.mat");

tic;
[data] = model_SH_synthesis(lonLim,latLim,height,SHbounds,V1,Model);
toc

%% Save data  

DATE = datestr(now);
save(['Results/data_' Model.name '_' num2str(SHbounds(1)) '_' num2str(SHbounds(2)) '_1_highresolution'  '.mat'],'data')


