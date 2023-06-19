% This is an input file for the GSHA procedure
%
% It should contain a list of names and location of geometry boundaries followed by a
% list of names for density values
% 
% Single values for density files are allowed.
HOME = pwd;
Model = struct();

Model.number_of_layers = 2;
Model.name = 'Q6_4crust';

%in meters! https://pds-geosciences.wustl.edu/mro/mro-m-rss-5-sdp-v1/mrors_1xxx/data/shadr/jgmro_110c_sha.tab
% Additional variables, mars in km: 0.3396000000000000E+04, 0.4282837564100000E+05
Model.GM = 0.42828375641E14; %earth=3.9860004415E14
Model.Re_analyse = 3396000;
Model.Re = 3396000; %earth=6378136.30
Model.geoid = 'none';
Model.nmax = 179; %180/resolution-1     
%Model.correct_depth = 0;

% % Topography layer
Model.l1.bound = [HOME '/Data/Q6_4.bd1.gmt'];%if character-->km, 
Model.l1.dens  = 2960;%in g/cm3 if numeric-->kg/m^3
% %Model.l1.alpha = 

% Crust layer
Model.l2.bound = [HOME '/Data/Q6_4.bd2.gmt'];%Q3
Model.l2.dens  = 3440; %in kg/m^3
% %Model.l2.alpha = 
% 
% % Mantle layer
Model.l3.bound = [HOME '/Data/Q6_4.bd3.gmt'];
%Model.l3.dens  = 3.450;
Model.rho_m=3440;
Model.rho_c=2960;
% Save model in .mat file for use of the new software

save([HOME '/Data/' Model.name '.mat'],'Model')