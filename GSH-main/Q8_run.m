% makefile for the complete GSH circle for a particular model
% clear;
% close all;
% clc;

HOME = pwd;
addpath([HOME '/Data']);
addpath([HOME '/Tools']);

% Model
% Load previous saved model

%model_name = 'Crust01_crust';
%load(model_name);

% Construct new model
Q6_input_q5

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
SHbounds =  [3 110]; % Truncation settings: lower limit, upper limit SH-coefficients used

%%%%%%%%%%%%%% Part that can be modified %%%%%%%%%%%%%%%%%%%%%%%
%% read the gravity observation
V1=readtable("Q2_V1.xlsx");
V1=table2array(V1);
V1=V1(:,1:4);
V1 = sortrows(V1,2);

V1=[0 0 1 0; V1];
[data] = model_SH_synthesis(lonLim,latLim,height,SHbounds,V1,Model);
g_obs=data.vec.R;
% ft="C:\Skolan\AATM\Planetary sciences\Assignment 3\GSH-main\Data\Q3.bd2.gmt";%chaning the second boundary
% A=dlmread(ft);
% dlmwrite(ft, [A(:,1), A(:,2), -50.*ones(180.*360, 1)], 'delimiter', ' ');
epsilon=[];
alpha=1000;
e=1;
%init
moho=dlmread(Model.l2.bound, ' ');
moho=gmt2matrix(moho).*1e3;
[V] = calcV(Model, moho, topography);%model_SH_analysis(Model); %initialized
[data] = model_SH_synthesis(lonLim,latLim,height,SHbounds,V,Model);
g_model=data.vec.R;
residual=g_obs-g_model; 

%figure;
lat=-89.5:1:89.5;
lon=0.5:1:359.5;
% figure;
% imagesc(lon, lat, residual.*1e5); c=colorbar;
% set(gca,'YDir','normal','Fontsize',11)
% xlabel('Longitude (\circ)','Fontsize',22)
% ylabel('Latitude (\circ)','Fontsize',22)
% ylabel(c,'Residual (mGal)','Fontsize',22)

%%

%moho=Model.l2.bound; %initialize the boundary between crust and mantle
dens=Model.l1.dens;
while abs(e)>5e-7 && length(epsilon)<20
    tic;
    %SHbounds=[3 110]; %low-pass filter
    
    %Q3_inputModel
    %
    
    
    %V(1000:end, 3:4)=V(1000:end, 3:4).*0.01;
    %[n, dV]=degreeVariance(V);
    %figure;
    %loglog(n(3:end), dV(3:end), "*"); xlabel("Spherical harmonic degree"); ylabel("Variance"); hold on;
    
    
    [data] = model_SH_synthesis(lonLim,latLim,height,SHbounds,V,Model);
    g_model=data.vec.R;
    residual=g_obs-g_model;
%     figure;
%     imagesc(g_model); c=colorbar;

    Model.l1.dens=Model.l1.dens+alpha.*residual;
    V=calcV(Model, moho, topography);
    
    
    
    max_densitychange=max(abs(residual(:).*alpha))
    e=mean(residual(:))
    epsilon=[epsilon, e];
    toc
end

%% plot

% depth=ones(180,360);
% for lats=1:1:180
%     for lons=1:1:360
%         if A(lats, lons)% bd5
%            depth(lats, lons)=Model.l5.bound(lats, lons);
%         else % bd3
%            depth(lats, lons)=Model.l3.bound(lats, lons);
%         end
%     end
% end   


%%
lat=-89.5:1:89.5;
lon=0.5:1:359.5;
figure;
imagesc(lon, lat, Model.l1.dens); c=colorbar;
set(gca,'YDir','normal','Fontsize',11)
xlabel('Longitude (\circ)','Fontsize',22)
ylabel('Latitude (\circ)','Fontsize',22)
ylabel(c,'Density (kg/m^3)','Fontsize',22)

%%
figure;
imagesc(lon, lat, residual.*1e5); c=colorbar;
set(gca,'YDir','normal','Fontsize',11)
xlabel('Longitude (\circ)','Fontsize',22)
ylabel('Latitude (\circ)','Fontsize',22)
ylabel(c,'Residual (mGal)','Fontsize',22)

%% save for later use
% 
% A=dlmread("C:\Skolan\AATM\Planetary sciences\Assignment 3\GSH-main\Data\crust1.bd1.gmt");
% lon=A(:,1);
% lat=A(:,2);
% 
% crustalthickness=zeros(180.*360,1);
% i=1;
% for lats=1:1:180
%     for lons=1:1:360
%         crustalthickness(i)=moho(lats, lons)./1e3;%in km
%         i=i+1;
%     end
% end
% 
% ft="C:\Skolan\AATM\Planetary sciences\Assignment 3\GSH-main\Data\Q8.bd2.gmt";
% 
% At=[lon, lat, crustalthickness];
% 
% dlmwrite(ft, At, 'delimiter', ' ');





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

    %Model.l1.dens  = Model.rho_c;

    % ice

    % crust

    Model.l2.bound = low_LAY; 

    %Model.l2.dens = Model.rho_m;

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