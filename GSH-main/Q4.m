% close all;
% clear all;
% clc;

filename = "C:\Skolan\AATM\Planetary sciences\Assignment 3\shape+topography\topography.img";
resolution = 4;
% Read in the file.
f = fopen(filename,'r','ieee-be');
topography = fread(f,[360*resolution Inf],'int16')';
%topography=flip(topography);
fclose(f);

topography = imresize(topography, [180, 360], 'bilinear').*1e-3;

D=50;
rho_crust=2960;
rho_mantle=3440;

r=topography.*(rho_crust/(rho_mantle-rho_crust));

depth=D+r;

lat=-89.5:1:89.5;
lon=0.5:1:359.5;

% Plot the surface
plot(1:1:length(depth(100,:)),depth(100, :)); xlim([0 360]); %latitude 100
xlabel('Longitude (\circ)','Fontsize',22); ylabel("Depth (km)", "Fontsize", 22);

figure;
imagesc(lon, lat, flip(depth)); c=colorbar;
set(gca,'YDir','normal','Fontsize',11)
xlabel('Longitude (\circ)','Fontsize',22)
ylabel('Latitude (\circ)','Fontsize',22)
ylabel(c,'Depth (km)','Fontsize',22)
%s=surf(X, Y, r); s.EdgeColor="none";
%% save for later use


A=dlmread("C:\Skolan\AATM\Planetary sciences\Assignment 3\GSH-main\Data\crust1.bd1.gmt");
lon=A(:,1);
lat=A(:,2);
depth=depth;
imagesc(depth)
crustalthickness=zeros(180.*360,1);
i=1;
for lats=1:1:180
    for lons=1:1:360
        crustalthickness(i)=-depth(lats, lons);%already in km
        i=i+1;
    end
end
%imagesc(crustalthickness)


%%
ft="C:\Skolan\AATM\Planetary sciences\Assignment 3\GSH-main\Data\Q6_4.bd2.gmt";

At=[lon, lat, crustalthickness];

dlmwrite(ft, At, 'delimiter', ' ');






