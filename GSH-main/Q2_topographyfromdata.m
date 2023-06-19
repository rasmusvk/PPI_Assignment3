%image=imread("C:\Skolan\AATM\Planetary sciences\Assignment 3\shape+topography\areoid.img")
filename = "C:\Skolan\AATM\Planetary sciences\Assignment 3\shape+topography\topography.img";
resolution = 4;
% Read in the file.
f = fopen(filename,'r','ieee-be');
topography = fread(f,[360*resolution Inf],'int16')';
topography=flip(topography);
fclose(f);

%% Topography is the planetary radius minus the areoid radius at a given longitude and latitude."
lon=0:0.25:359.75;
lat=-90:0.25:89.75; %data is flipped

imagesc(lon,lat,topography);c=colorbar; 
hold on
%plot(long,lat,'k','LineWidth',1.5);%coastline
xlim([min(lon) max(lon)])
ylim([min(lat) max(lat)])
hold off
xlabel('Longitude [^o]')
ylabel('Latitude [^o]')
title(['Topography'])
ylabel(c,'m') 
set(gca,'YDir','normal')


%% bouguer
G=6.6743e-11;
rho_topo=2960;
h=topography;
deltag_b=2*pi*G*rho_topo.*h;

imagesc(lon,lat,deltag_b.*1e5);c=colorbar; 
hold on
%plot(long,lat,'k','LineWidth',1.5);%coastline
xlim([min(lon) max(lon)])
ylim([min(lat) max(lat)])
hold off
xlabel('Longitude [^o]')
ylabel('Latitude [^o]')
title(['Bouguer correction'])
ylabel(c,'mGal') 
set(gca,'YDir','normal')



%% 







