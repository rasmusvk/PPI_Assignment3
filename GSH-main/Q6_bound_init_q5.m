filename = "C:\Skolan\AATM\Planetary sciences\Assignment 3\shape+topography\topography.img";
resolution = 4;
% Read in the file.
f = fopen(filename,'r','ieee-be');
topography = fread(f,[360*resolution Inf],'int16')';
%topography=flip(topography);
fclose(f);

topography = imresize(topography, [180, 360], 'bilinear');
%x=0:1:179;
%y=0:1:359;
%imagesc(x,y,(topography)); c=colorbar;
%A=dlmread("C:\Skolan\AATM\Planetary sciences\Assignment 3\GSH-main\Data\crust1.bd1.gmt");
A=dlmread("C:\Skolan\AATM\Planetary sciences\Assignment 3\GSH-main\Data\crust1.bd1.gmt");
lon=A(:,1);
lat=A(:,2);

topographyvector=zeros(180.*360,1);
i=1;
for lats=1:1:180
    for lons=1:1:360
        topographyvector(i)=topography(lats, lons)./1e3;%in km
        i=i+1;
    end
end


crustvector=-50.*ones(180.*360, 1); 
mantlevector=-500.*ones(180.*360, 1);


ft="C:\Skolan\AATM\Planetary sciences\Assignment 3\GSH-main\Data\Q3.bd1.gmt";
fc="C:\Skolan\AATM\Planetary sciences\Assignment 3\GSH-main\Data\Q3.bd2.gmt";
fm="C:\Skolan\AATM\Planetary sciences\Assignment 3\GSH-main\Data\Q3.bd3.gmt";

At=[lon, lat, topographyvector];
Ac=[lon, lat, crustvector];
Am=[lon, lat, mantlevector];

dlmwrite(ft, At, 'delimiter', ' ');
dlmwrite(fc, Ac, 'delimiter', ' ');
dlmwrite(fm, Am, 'delimiter', ' ');




%%
function vector=trans(A)
    v=zeros(180.*360,1);
    i=1;
    for lons=1:1:360
        for lats=1:1:180
            v(i)=A(lats, lons);%in km
            i=i+1;
        end
    end
   vector=v;
end





