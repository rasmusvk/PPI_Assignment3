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
%Q3_inputModel    

%%%%%%%%%%%%%%%%%%% Computation area %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Part that can be modified %%%%%%%%%%%%%%%%%%%%%%%

latLim =    [-89.5 89.5 1];  % [deg] min latitude, max latitude, resolution latitude (preferable similar to latitude)
lonLim =    [0.5 359.5 1];% [deg] min longitude, max longitude, resolution longitude (preferable similar to latitude)
height =    0.0; % height of computation above spheroid
SHbounds =  [3 110]; % Truncation settings: lower limit, upper limit SH-coefficients used

%%%%%%%%%%%%%% Part that can be modified %%%%%%%%%%%%%%%%%%%%%%%
%% read the gravity observation
load('Results/data_Crust10_crust_3_110_rightsize_0km.mat')
g_obs=data.vec.R; %180x360 data of observed 
% ft="C:\Skolan\AATM\Planetary sciences\Assignment 3\GSH-main\Data\Q3.bd2.gmt";%chaning the second boundary
% A=dlmread(ft);
% dlmwrite(ft, [A(:,1), A(:,2), -50.*ones(180.*360, 1)], 'delimiter', ' ');
epsilon=[];
alpha=500;
e=1;
while abs(e)>1e-6 && length(epsilon)<50
    tic;
    SHbounds=[3 110]; %low-pass filter
    
    %Q3_inputModel
    %
    [V] = model_SH_analysis(Model); %initialized
    %V(1000:end, 3:4)=V(1000:end, 3:4).*0.01;
    %[n, dV]=degreeVariance(V);
    %figure;
    %loglog(n(3:end), dV(3:end), "*"); xlabel("Spherical harmonic degree"); ylabel("Variance"); hold on;
    
    
    [data] = model_SH_synthesis(lonLim,latLim,height,SHbounds,V,Model);
    g_model=data.vec.R;
    residual=g_obs-g_model;
    figure;
    imagesc(g_model); c=colorbar;
    residualvector=trans(residual);
%     if abs(mean(residualvector))>abs(mean(e))
%        break
%     end
    A=dlmread(ft);
    bound2=A(:,3); %read the boundary
    bound2_new=bound2+residualvector.*alpha;
    B=[A(:,1), A(:,2), bound2_new];
    %Model.l2.bound=B;
    dlmwrite(ft, B, 'delimiter', ' ');
    max_elevationchange=max(residualvector*alpha)
    e=mean(residualvector)
    epsilon=[epsilon, e];
    toc
end



%% functions
function vector=trans(A)
    v=zeros(180.*360,1);
    i=1;
    for lats=1:1:180
        for lons=1:1:360
            v(i)=A(lats, lons);%in km
            i=i+1;
        end
    end
   vector=v;
end

