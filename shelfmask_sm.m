function [ mask ] = shelfmask_sm( run, onshelf, iceshelf, shelfbreak )
%SHELFMASK_SM - Creates a mask for the Ross Sea grid for on/off shelf
%measurements based on a certain depth set as the shelf break.
%   INPUT:  run - name of folder for a given ROMS run, directory is already
%                 specified in function; e.g. '016'
%           onshelf - binary option for specifying only volumes on the
%                 shelf (true - 1) or off the shelf (false - 0)
%           iceshelf - binary option for including (true-1) or not including
%                 (false-0) volumes under the ice shelf
%           shelfbreak - positive integer value of depth of shelf break;
%                 e.g. 1500
%           To include entire model domain, set both onshelf and iceshelf
%           to true, and set shelfbreak to a number larger than the deepest
%           depth in the domain, say 10000 for example.
%           To include only under iceshelf values, set onshelf to false and
%           iceshelf to true.
%           To include only offshelf values, sent both to false.
%   OUTPUT: mask - Grid the size of the domain where all areas but the
%   specified have been masked out.  Takes into account deep areas on the
%   shelf for relatively shallow values of shelfbreak

filelist = romsinitialize_sm('his',run);
mask_rho = nc_varget(filelist(1,:),'mask_rho');
zice = nc_varget(filelist(1,:),'zice');
h = nc_varget(filelist(1,:), 'h'); %meters

%mask_deep covers deep on shelf points
mask_deep = zeros(size(h));
if(size(h,1)<800)  %If 5km grid
    mask_deep(1:180,:)=1;
    mask_deep(1:250,1:100)=1;
else %If 1.5km grid
    mask_deep(1:600,:)=1;
    mask_deep(1:800,1:350)=1;
end

mask_shelf = zeros(size(h));
mask_ice = ones(size(h));

if(onshelf) %set on or off shelf to 1
    shelf = h>40 & h<shelfbreak;
elseif(iceshelf) %If not 'onshelf' but 'iceshelf' - do iceshelf only
    shelf = zice<=-10;
else
    shelf = h>=shelfbreak;
end
mask_shelf(shelf)=1;

if(onshelf) %If onshelf, include deep values
    mask_shelf = ceil((mask_shelf+mask_deep)./2);
end

if(~onshelf && iceshelf) %If iceshelf only
    if(size(h,1)<800) %If 5km grid
        mask_shelf(165:370,:)=0; %Glaciers north of IS
        mask_shelf(140:370,225:279)=0; %Glaciers to the east
    else
        mask_shelf(:,775:930)=0; %Glaciers east
        mask_shelf(560:1230,:)=0; %Glaciers north
    end
end

ice = zice<=-10;
if(~iceshelf)
    mask_ice(ice)=0;
end
mask = mask_rho.*mask_ice.*mask_shelf;

end

