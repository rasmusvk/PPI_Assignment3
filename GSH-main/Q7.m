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
latLim =    [-89.5 89.5 1];  % [deg] min latitude, max latitude, resolution latitude (preferable similar to latitude)
lonLim =    [0.5 359.5 1];% [deg] min longitude, max longitude, resolution longitude (preferable similar to latitude)
height =    0.0; % height of computation above spheroid
SHbounds =  [3 179]; % Truncation settings: lower limit, upper limit SH-coefficients used

latLimT = [-90+(1/resolution/2) 90-(1/resolution/2) 1/resolution]; 
lonLimT = [1/resolution/2  360-(1/resolution/2) 1/resolution]; 

lonT = lonLimT(1):lonLimT(3):lonLimT(2);
latT = fliplr(latLimT(1):latLimT(3):latLimT(2));
LonT = repmat(lonT,length(latT),1);
LatT = repmat(latT',1,length(lonT));
% 
aa = 18;
% figure
% imagesc(lonT,latT,topography./1e3);cc=colorbar;
% xlabel('Longitude (\circ)','Fontsize',aa)
% ylabel('Latitude (\circ)','Fontsize',aa)
% ylabel(cc,'Topography (km)','Fontsize',aa)
% set(gca,'YDir','normal','Fontsize',aa)

%% airy model, here everything is in km
topography_s = imresize(topography, [180, 360], 'bilinear').*1e-3;

D=50;
rho_crust=2960;
rho_mantle=3440;

r=topography_s.*(rho_crust/(rho_mantle-rho_crust));

depth=D+r;



%% GSHA analysis

cs = GSHA(r.*1e3,179);
sc = cs2sc(cs);
% [Clm,Slm,llvec,mmvec]=sc2vecml(sc, 110);
% V=[llvec', mmvec', Clm, Slm];
% [nV, dV]=degreeVariance(V);
% loglog(nV, dV);


n = 1:size(sc,1);

T_e=[80, 85, 90, 95, 100, 105, 110].*1e3;
E=120e9;
sigma=0.25; %poisson


g=3.76;
rho_mantle=3440;
rho_crust=2960;
R=3396000;
var=zeros(180,length(T_e));
% A=cell(1, length(T_e));
ii=1;
figure
for i=T_e
    tic;
    D = E*i^3/(12*(1-sigma^2));
    PHI = (1 + (D)/((rho_mantle-rho_crust)*g).*(2.*(n+1)./(2*R)).^4).^(-1);
    semilogx(n, PHI); hold on;

    %multiplying with flecure response
    sc_flex = zeros(size(sc));

    for m = 1:size(sc,2)
        sc_flex(:,m) = sc(:,m).*PHI';
    end
    mapf = 50e3+GSHS(sc_flex,lonT,90-latT,179);
%     figure;
    %imagesc(mapf)
    Q7_input
    
    Model.l2.bound=-1.*imresize(mapf, [180, 360], 'bilinear');
    
    [V] = model_SH_analysis(Model);
    %[data] = model_SH_synthesis(lonLim,latLim,height,SHbounds,V,Model);
    [nV, dV]=degreeVariance(V);
    var(:,ii)=dV;
    %loglog(nV(1:110), dV(1:110).*1e10, "*"); hold on;
    toc
    ii=ii+1;
    %A{1,i}=50e3+mapf;
end



%% GSHS
V=readtable("Q2_V1.xlsx");
V=table2array(V);
V=V(:,1:4);
V = sortrows(V,2);

V=[0 0 1 0; V];
[nV, dV]=degreeVariance(V);
lat=-89.5:1:89.5;
lon=0.5:1:359.5;

%[data] = model_SH_synthesis(lonLim,latLim,height,[3 110],V,Model);
%%
figure;
for x =1:length(T_e)
    loglog(nV(5:30), var(5:30, x).*1e10, "*"); hold on;
end

loglog(nV(5:30), dV(5:30).*1e10, "*r"); hold on;

legend("T_e=80", "T_e=85", "T_e=90", "T_e=95", "T_e=100", "T_e=105", "T_e=110", "Observed");
xlabel("Spherical harmonic degree (n)"); ylabel("Power spectrum (mGal)"); xlim([5 30]);
%% 
diff=(var(5:30,:)-dV(5:30)).*1e5;
val=zeros(1,length(T_e));
for i =1:length(T_e)
    val(i)=sqrt(sum(diff(:,i).^2));
end
%%
% figure
% imagesc(lonT,latT,A{1}./1e3);cc=colorbar;
% xlabel('Longitude (\circ)','Fontsize',aa)
% ylabel('Latitude (\circ)','Fontsize',aa)
% ylabel(cc,'Crustal thickness (km)','Fontsize',aa)
% set(gca,'YDir','normal','Fontsize',11)




