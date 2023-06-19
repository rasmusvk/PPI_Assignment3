% This is an input file for the GSHA procedure
%
% It should contain a list of names and location of geometry boundaries followed by a
% list of names for density values

filename = "C:\Skolan\AATM\Planetary sciences\Assignment 3\shape+topography\topography.img";
resolution = 4;
% Read in the file.
f = fopen(filename,'r','ieee-be');
topography = fread(f,[360*resolution Inf],'int16')';
%topography=flip(topography);
fclose(f);

topography = imresize(topography, [180, 360], 'bilinear');


% Single values for density files are allowed.
HOME = pwd;
Model = struct();

Model.number_of_layers = 2;
Model.name = 'Q3_morelayers';

%in meters! https://pds-geosciences.wustl.edu/mro/mro-m-rss-5-sdp-v1/mrors_1xxx/data/shadr/jgmro_110c_sha.tab
% Additional variables, mars in km: 0.3396000000000000E+04, 0.4282837564100000E+05
Model.GM = 0.42828375641E14; %earth=3.9860004415E14
Model.Re_analyse = 3396000;
Model.Re = 3396000; %earth=6378136.30
Model.geoid = 'none';
Model.nmax = 179; %180/resolution-1     
%Model.correct_depth = 0;

% % Topography layer
Model.l1.bound = topography;%if character-->km, 
Model.l1.dens  = cb(2960);%if char in g/cm3, if numeric-->kg/m^3
% %Model.l1.alpha = 

% 0 km, this is what we are changing at first. 
Model.l2.bound = cb(-50e3);%in m
Model.l2.dens  = cb(2960); %in kg/m^3
Model.l2.og=-25e3;
% %Model.l2.alpha = 
% 

% -25 km
Model.l3.bound = cb(-300e3);
% Save model in .mat file for use of the new software

Model.rho_m=3440;
Model.rho_c=2960; %2960



save([HOME '/Data/' Model.name '.mat'],'Model')


function vector=cb(depth)
   v=depth.*ones(180,360);
   vector=v;
end






