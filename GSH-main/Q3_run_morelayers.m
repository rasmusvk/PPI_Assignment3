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
Q3_inputModel_morelayers 

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
alpha=100e3;
e=1;
%init
[V] = model_SH_analysis(Model); %initialized
[data] = model_SH_synthesis(lonLim,latLim,height,SHbounds,V,Model);
g_model=data.vec.R;
residual=g_obs-g_model; 

figure;
lat=-89.5:1:89.5;
lon=0.5:1:359.5;
figure;
imagesc(lon, lat, residual.*1e5); c=colorbar;
set(gca,'YDir','normal','Fontsize',11)
xlabel('Longitude (\circ)','Fontsize',22)
ylabel('Latitude (\circ)','Fontsize',22)
ylabel(c,'Residual (mGal)','Fontsize',22)

%%
%matrix for if we add or subtract mass. If 1, we have too much mass in our model whch means that we 
%need to lower the crustal depth. This means that bd5 needs to be lowered,
%allowing more crustal density. if 0, bd3 needs to be changed.
A=residual<0; 

while abs(e)>1e-6 && length(epsilon)<20
    tic;
    %SHbounds=[3 110]; %low-pass filter
    
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
    %residualvector=trans(residual);
%     if abs(mean(residualvector))>abs(mean(e))
%        break
%     end
    %A=dlmread(ft);
    
    
    for lats=1:1:180
        for lons=1:1:360
            if A(lats, lons)%change bd5
               Model.l5.bound(lats, lons)=Model.l5.bound(lats, lons)+residual(lats, lons).*alpha;
            else %change bd3
               Model.l3.bound(lats, lons)=Model.l3.bound(lats, lons)+residual(lats, lons).*alpha;
            end
        end
    end   
    
    
%     bound_u=Model.l2.bound; %read the upper boundary, in matrix form
%     og_u=Model.l2.og;
%     bound_l=Model.l3.bound; %read the lower boundary, in matrix form
%     og_l=Model.l3.og;
%     newbound=bound_u+residual.*alpha;
%     for lats=1:1:180
%         for lons=1:1:360
%             i=0;
%             while i<11
%                if newbound(lats, lons)<=bound_l
%                     
%                i=i+1;
%                end
%             end
%         end
%     end
    
    %B=[A(:,1), A(:,2), bound2_new];
    %Model.l2.bound=B;
    %dlmwrite(ft, B, 'delimiter', ' ');
    max_elevationchange=max(residual(:).*alpha)
    e=mean(residual(:))
    epsilon=[epsilon, e];
    toc
end

%% plot

depth=ones(180,360);
for lats=1:1:180
    for lons=1:1:360
        if A(lats, lons)% bd5
           depth(lats, lons)=Model.l5.bound(lats, lons);
        else % bd3
           depth(lats, lons)=Model.l3.bound(lats, lons);
        end
    end
end   
%%
lat=-89.5:1:89.5;
lon=0.5:1:359.5;
figure;
imagesc(lon, lat, depth); c=colorbar;
set(gca,'YDir','normal','Fontsize',11)
xlabel('Longitude (\circ)','Fontsize',22)
ylabel('Latitude (\circ)','Fontsize',22)
ylabel(c,'Depth (km)','Fontsize',22)

%%
figure;
imagesc(lon, lat, residual.*1e5); c=colorbar;
set(gca,'YDir','normal','Fontsize',11)
xlabel('Longitude (\circ)','Fontsize',22)
ylabel('Latitude (\circ)','Fontsize',22)
ylabel(c,'Residual (mGal)','Fontsize',22)

%% save for later use

A=dlmread("C:\Skolan\AATM\Planetary sciences\Assignment 3\GSH-main\Data\crust1.bd1.gmt");
lon=A(:,1);
lat=A(:,2);

crustalthickness=zeros(180.*360,1);
i=1;
for lats=1:1:180
    for lons=1:1:360
        crustalthickness(i)=depth(lats, lons)./1e3;%in km
        i=i+1;
    end
end

ft="C:\Skolan\AATM\Planetary sciences\Assignment 3\GSH-main\Data\Q6_3.bd2.gmt";

At=[lon, lat, crustalthickness];

dlmwrite(ft, At, 'delimiter', ' ');





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



