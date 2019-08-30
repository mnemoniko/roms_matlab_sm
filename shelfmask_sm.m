function [ mask ] = shelfmask_sm( in_dir, onshelf, iceshelf, shelfbreak )
%SHELFMASK_SM - Creates a mask for the Ross Sea grid for on/off shelf
%measurements based on a certain depth set as the shelf break.
%   INPUT:  in_dir - directory for ROMS output
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

filelist = romsinitialize_sm(in_dir,'his');
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

%Set on shelf mask:
mask_shelf = zeros(size(h));
shelf = h>40 & h<shelfbreak;
mask_shelf(shelf)=1;
mask_shelf = ceil((mask_shelf+mask_deep)./2);

if(~onshelf)
    mask_shelf = abs(mask_shelf-1);
end

%Set ice mask:
mask_ice = zeros(size(h));
ice = zice<=-10;
mask_ice(ice)=1;


%Iceshelf only special option:
if(~onshelf && iceshelf) %If iceshelf only
    if(size(h,1)<800) %If 5km grid
        mask_ice(165:370,:)=0; %Glaciers north of IS
        mask_ice(140:370,225:279)=0; %Glaciers to the east
    else
        mask_ice(:,775:930)=0; %Glaciers east
        mask_ice(560:1230,:)=0; %Glaciers north
    end
    mask = mask_ice.*mask_rho;
elseif(onshelf && iceshelf) %ice shelf + on shelf
    mask = ceil((mask_ice+mask_shelf)./2).*mask_rho;
else %No ice shelf cases
    mask = mask_shelf.*mask_rho;
end


end

