% clear all;
% close all;
% clc;

load('Results/data_Crust10_crust_3_110_rightsize_0km.mat')
lon_o = data.grd.lon(1,:);
lats_o = data.grd.lat(:,1);


filename = "C:\Skolan\AATM\Planetary sciences\Assignment 3\shape+topography\topography.img";
resolution = 4;
% Read in the file.
f = fopen(filename,'r','ieee-be');
topography = fread(f,[360*resolution Inf],'int16')';
topography=flip(topography);
fclose(f);

lon_t=0:0.25:359.75;
lat_t=-90:0.25:89.75;

%% bouguercorrection
g_observed=data.vec.R;

G=6.6743e-11;
rho_topo=2960;
h=topography;
deltag_b=2*pi*G*rho_topo.*h;

resized_deltag_b = imresize(deltag_b, [length(lats_o), length(lon_o)], 'bilinear');
resized_deltag_b=flip(resized_deltag_b);
%imagesc(lon_o,lats_o,(resized_deltag_b)); c=colorbar;



%%
A_b=flip(g_observed)-resized_deltag_b; %calculating the anomaly

figure;
imagesc(lon_o,lats_o,A_b.*1e5); c=colorbar;
hold on
%plot(long,lat,'k','LineWidth',1.5);%coastline
xlim([min(lon_o) max(lon_o)])
ylim([min(lats_o) max(lats_o)])
hold off
xlabel('Longitude [^o]')
ylabel('Latitude [^o]')
%title(['Bouguer anomaly'])
ylabel(c,'mGal') 
set(gca,'YDir','normal')









