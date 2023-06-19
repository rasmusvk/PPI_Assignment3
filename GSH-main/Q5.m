clear;
close all;
clc;

addpath("C:\Skolan\AATM\Planetary sciences\Assignment 3\GSH-main\Tools")

%%
filename = "C:\Skolan\AATM\Planetary sciences\Assignment 3\shape+topography\topography.img";
resolution = 4;

% Read in the file.
f = fopen(filename,'r','ieee-be');
topography = fread(f,[360*resolution Inf],'int16')';
fclose(f); 

%%

latLimT = [-90+(1/resolution/2) 90-(1/resolution/2) 1/resolution]; 
lonLimT = [1/resolution/2  360-(1/resolution/2) 1/resolution]; 

lonT = lonLimT(1):lonLimT(3):lonLimT(2);
latT = fliplr(latLimT(1):latLimT(3):latLimT(2));
LonT = repmat(lonT,length(latT),1);
LatT = repmat(latT',1,length(lonT));

aa = 18;
figure
imagesc(lonT,latT,topography./1e3);cc=colorbar;
xlabel('Longitude (\circ)','Fontsize',aa)
ylabel('Latitude (\circ)','Fontsize',aa)
ylabel(cc,'Topography (km)','Fontsize',aa)
set(gca,'YDir','normal','Fontsize',aa)

%% airy model, here everything is in km
topography_s = imresize(topography, [180, 360], 'bilinear').*1e-3;

D=50;
rho_crust=2960;
rho_mantle=3440;

r=topography_s.*(rho_crust/(rho_mantle-rho_crust));
imagesc(r)
depth=D+r;



%% GSHA analysis

cs = GSHA(r.*1e3,179);
sc = cs2sc(cs);
% [Clm,Slm,llvec,mmvec]=sc2vecml(sc, 110);
% V=[llvec', mmvec', Clm, Slm];
% [nV, dV]=degreeVariance(V);
% loglog(nV, dV);


n = 1:size(sc,1);

T_e=90e3;
E=120e9;
sigma=0.25; %poisson
D = E*T_e^3/(12*(1-sigma^2));

g=3.76;
rho_mantle=3440;
rho_crust=2960;
R=3396000;

PHI = (1 + (D)/((rho_mantle-rho_crust)*g).*(2.*(n+1)./(2*R)).^4).^(-1);


%multiplying with flecure response
sc_flex = zeros(size(sc));

for m = 1:size(sc,2)
    sc_flex(:,m) = sc(:,m).*PHI';
end

% [Clm,Slm,llvec,mmvec]=sc2vecml(sc_flex, 110);
% V=[llvec', mmvec', Clm, Slm];
% [nV, dV]=degreeVariance(V);
% loglog(nV(3:end), dV(3:end), "*");


%% GSHS

mapf = GSHS(sc_flex,lonT,90-latT,110);
depth=50e3+mapf; %adding back the same 
    
%%
figure
imagesc(lonT,latT,depth./1e3);cc=colorbar;
xlabel('Longitude (\circ)','Fontsize',aa)
ylabel('Latitude (\circ)','Fontsize',aa)
ylabel(cc,'Crustal thickness (km)','Fontsize',aa)
set(gca,'YDir','normal','Fontsize',11)

% figure
% imagesc(lonT,latT,(topography-mapf)./1e3);cc=colorbar;
% xlabel('Longitude (\circ)','Fontsize',aa)
% ylabel('Latitude (\circ)','Fontsize',aa)
% ylabel(cc,'Difference (km)','Fontsize',aa)
% set(gca,'YDir','normal','Fontsize',aa)


%% save for later use
A=dlmread("C:\Skolan\AATM\Planetary sciences\Assignment 3\GSH-main\Data\crust1.bd1.gmt");
lon=A(:,1);
lat=A(:,2);

crustalthickness=zeros(180.*360,1);
i=1;
depth=imresize(depth, [180, 360], 'bilinear');
for lats=1:1:180
    for lons=1:1:360
        crustalthickness(i)=depth(lats, lons)./1e3;%in km
        i=i+1;
    end
end

ft="C:\Skolan\AATM\Planetary sciences\Assignment 3\GSH-main\Data\Q6_5.bd2.gmt";

At=[lon, lat, -1.*crustalthickness];

dlmwrite(ft, At, 'delimiter', ' ');


