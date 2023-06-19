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

filename = "C:\Skolan\AATM\Planetary sciences\Assignment 3\shape+topography\topography.img";
resolution = 4;
% Read in the file.
f = fopen(filename,'r','ieee-be');
topography = fread(f,[360*resolution Inf],'int16')';
%topography=flip(topography);
fclose(f);

topography = imresize(topography, [180, 360], 'bilinear');

%%%%%%%%%%%%%%%%%%% Computation area %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Part that can be modified %%%%%%%%%%%%%%%%%%%%%%%

latLim =    [-89.5 89.5 1];  % [deg] min latitude, max latitude, resolution latitude (preferable similar to latitude)
lonLim =    [0.5 359.5 1];% [deg] min longitude, max longitude, resolution longitude (preferable similar to latitude)
height =    0.0; % height of computation above spheroid
SHbounds =  [3 179]; % Truncation settings: lower limit, upper limit SH-coefficients used

%%%%%%%%%%%%%% Part that can be modified %%%%%%%%%%%%%%%%%%%%%%%

%% Global Spherical Harmonic Analysis 
Q6_input_q3
tic;
%[V3] = model_SH_analysis(Model);
%fil="C:\Skolan\AATM\Planetary sciences\Assignment 3\GSH-main\Data\V_q3_forQ6.txt";
%V3=dlmread(fil, ' ');
%V3=model_SH_analysis(Model);
moho3=dlmread(Model.l2.bound, ' ');
moho3=gmt2matrix(moho3).*1e3;
topo3=dlmread(Model.l1.bound, ' ');
topo3=gmt2matrix(topo3).*1e3;
[V3] = calcV(Model,moho3, topography);
[data3] = model_SH_synthesis(lonLim,latLim,height,SHbounds,V3,Model);
toc
%%
Q6_input_q4
tic;
moho4=dlmread(Model.l2.bound, ' ');
moho4=gmt2matrix(moho4).*1e3;
topo4=dlmread(Model.l1.bound, ' ');
topo4=gmt2matrix(topo4).*1e3;
[V4] = calcV(Model,moho4, topography);%model_SH_analysis(Model);
[data4] = model_SH_synthesis(lonLim,latLim,height,SHbounds,V4,Model);
toc
%%
Q6_input_q5
tic;
moho5=dlmread(Model.l2.bound, ' ');
moho5=gmt2matrix(moho5).*1e3;
topo5=dlmread(Model.l1.bound, ' ');
topo5=gmt2matrix(topo5).*1e3;
[V5] = calcV(Model,moho5, topography);%model_SH_analysis(Model);
[data5] = model_SH_synthesis(lonLim,latLim,height,SHbounds,V5,Model);
toc
%save(['Results/' Model.name '.mat'],'V')


V1=readtable("Q2_V1.xlsx");
V1=table2array(V1);
V1=V1(:,1:4);
V1 = sortrows(V1,2);

V1=[0 0 1 0; V1];


%% plot
lat=-89.5:1:89.5;
lon=0.5:1:359.5;
load('Results/data_Crust10_crust_3_110_rightsize_0km.mat')
g_obs1=data.vec.R; %180x360 data of observed 
g_3=data3.vec.R;
g_4=data4.vec.R;
g_5=data5.vec.R;
%% plot gravity
figure;
imagesc(lon, lat, g_3.*1e5); c=colorbar;
set(gca,'YDir','normal','Fontsize',11)
xlabel('Longitude (\circ)','Fontsize',22)
ylabel('Latitude (\circ)','Fontsize',22)
ylabel(c,'mGal','Fontsize',22)
title("Q3")


figure;
imagesc(lon, lat, g_4.*1e5); c=colorbar;
set(gca,'YDir','normal','Fontsize',11)
xlabel('Longitude (\circ)','Fontsize',22)
ylabel('Latitude (\circ)','Fontsize',22)
ylabel(c,'mGal','Fontsize',22)
title("Q4")


figure;
imagesc(lon, lat, g_5.*1e5); c=colorbar;
set(gca,'YDir','normal','Fontsize',11)
xlabel('Longitude (\circ)','Fontsize',22)
ylabel('Latitude (\circ)','Fontsize',22)
ylabel(c,'mGal','Fontsize',22)
title("Q5")

figure;
imagesc(lon, lat, g_obs1.*1e5); c=colorbar;
set(gca,'YDir','normal','Fontsize',11)
xlabel('Longitude (\circ)','Fontsize',22)
ylabel('Latitude (\circ)','Fontsize',22)
ylabel(c,'mGal','Fontsize',22)
title("Observed")


%% plot residual


figure;
imagesc(lon, lat, (g_obs1-g_3).*1e5); c=colorbar;
set(gca,'YDir','normal','Fontsize',11)
xlabel('Longitude (\circ)','Fontsize',22)
ylabel('Latitude (\circ)','Fontsize',22)
ylabel(c,'mGal','Fontsize',22)
title("Q3")
error3=sqrt(sum(sum((g_obs1-g_3).^2)))

figure;
imagesc(lon, lat, (g_obs1-g_4).*1e5); c=colorbar;
set(gca,'YDir','normal','Fontsize',11)
xlabel('Longitude (\circ)','Fontsize',22)
ylabel('Latitude (\circ)','Fontsize',22)
ylabel(c,'mGal','Fontsize',22)
title("Q4")
error4=sqrt(sum(sum((g_obs1-g_4).^2)))

figure;
imagesc(lon, lat, (g_obs1-g_5.*1e5)); c=colorbar;
set(gca,'YDir','normal','Fontsize',11)
xlabel('Longitude (\circ)','Fontsize',22)
ylabel('Latitude (\circ)','Fontsize',22)
ylabel(c,'mGal','Fontsize',22)
title("Q5")
error5=sqrt(sum(sum((g_obs1-g_5).^2)))


%% plot degree variance
figure
[nV, dV]=degreeVariance(V1);
loglog(nV(1:110), dV(1:110).*1e10, "*"); hold on;
[nV3, dV3]=degreeVariance(V3);
loglog(nV3(1:110), dV3(1:110).*1e10, "*"); hold on;
[nV4, dV4]=degreeVariance(V4);
loglog(nV4(1:110), dV4(1:110).*1e10, "*"); hold on;
[nV5, dV5]=degreeVariance(V5);
loglog(nV5(1:110), dV5(1:110).*1e10, "*"); hold on; legend("Observed","Q3-inverse", "Q4-Airy", "Q5-flexure");
xlabel("Spherical harmonic degree"); ylabel("Variance (mGal^2)");

%% plot crustal thickness

figure;
imagesc(lon, lat, flip(moho3./1e3)); c=colorbar;
set(gca,'YDir','normal','Fontsize',11)
xlabel('Longitude (\circ)','Fontsize',22)
ylabel('Latitude (\circ)','Fontsize',22)
ylabel(c,'km','Fontsize',22)
title("Q3")


figure;
imagesc(lon, lat, flip(moho4)./1e3); c=colorbar;
set(gca,'YDir','normal','Fontsize',11)
xlabel('Longitude (\circ)','Fontsize',22)
ylabel('Latitude (\circ)','Fontsize',22)
ylabel(c,'km','Fontsize',22)
title("Q4")


figure;
imagesc(lon, lat, flip(moho5)./1e3); c=colorbar;
set(gca,'YDir','normal','Fontsize',11)
xlabel('Longitude (\circ)','Fontsize',22)
ylabel('Latitude (\circ)','Fontsize',22)
ylabel(c,'km','Fontsize',22)
title("Q5")


%%
function V =calcV(Model, moho, topo)
top = 25000;

bot = -300000; min(min(moho));

thick_lay = 25000;

top_layer = top:-thick_lay:bot;

bot_layer = [top_layer(2:end) bot];

for numl = 1:length(top_layer)

    disp(['Constructing layer number ' num2str(numl) '...'])

    ubound = top_layer(numl);

    lbound = bot_layer(numl);

    % LAB is shallower than 100 km

    upper_LAY = topo; 

    upper_LAY(topo>ubound) = ubound;                                                                                                                                       

    upper_LAY(topo<lbound) = lbound;

    low_LAY = moho; 

    low_LAY(moho>ubound) = ubound;                                                                                                                                       

    low_LAY(moho<lbound) = lbound; 

    % crust

    Model.l1.bound = upper_LAY;        

    Model.l1.dens  = Model.rho_c;

    % ice

    % crust

    Model.l2.bound = low_LAY; 

    Model.l2.dens = Model.rho_m;

    % mantle

    Model.l3.bound = lbound; 

    % perform spherical harmonic analyses and synthesis       

    [Vlay] = model_SH_analysis(Model);

    % add to previous coefficients           

    if numl == 1

        V_Model = Vlay;

    else

        V_Model(:,3) = V_Model(:,3) + Vlay(:,3);

        V_Model(:,4) = V_Model(:,4) + Vlay(:,4);   

    end

end 
V=V_Model;
end


