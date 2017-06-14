function [ h, zice, mask ] = ROMSgrdthin_sm( h, zice, mask_rho, hmin, zmin, land )
%ROMSGRDTHIN_SM Checks for bathymetry points that are thinner than a
%minimum allowed thickness and fixes them
%   h - ROMS bathymetry
%   zice - ROMS ice shelf
%   mask_rho - ROMS land mask for rho points
%   hmin - Minimum water column thickness
%   zmin - Minimum ice shelf thickness
%   land - Binary option to set too thin points to land (for h only) - 1
%   sets it to land, 0 keeps it water

% Set zice at mask = 0 to bathy

zice(mask_rho==0) = -h(mask_rho==0);

%Check for nans before continuing!  Must have nans in zice to not mess with
%non-ice ocean points!!!
zice(zice==0)=NaN;

%h(h<hmin) = hmin;
zice(zice>-zmin) = -zmin;

% Adjust thin water layers under ice shelf - initial pass

wc_thick = h+zice;
wc_thick(mask_rho==0)=hmin; %Arbitrarily set water column thickness to not modify land points
thin = wc_thick<hmin;  %Points that are too thin (459 total, 252 @ 1/2)
pts = length(find(thin==1));

if(land)
    land = wc_thick<hmin/2; %Set points that are less than hmin to land
    mask_rho(land)=0; %Readjust mask_rho

end

%% Deepen thin points still left:
counter = 0;

while(pts>0)
    h(thin) = h(thin) + .5; %Deepen bathy by .5m at thin points
    zice(thin) = zice(thin) + .5; %Shallow ice by .5m at thin points
    
    %Reset h & zice mins:
    h(h<hmin) = hmin;
    zice(zice>-zmin) = -zmin;
    
    %Re-calculate wc thickness & thin points
    
    wc_thick = h+zice;
    wc_thick(mask_rho==0)=hmin;
    thin = wc_thick<hmin;
    pts = length(find(thin==1));
    
    counter = counter + 1;  
    
end

%Reset water zice points to 0 ice thickness
zice(isnan(zice))=0;
mask = mask_rho;

end

