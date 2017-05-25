function [ volume ] = volume_sm(run)
%VOLUME_SM Calculates cell by cell volume for the ROSS grid, given specific
%parameters as input. Uses romsinitialize_sm & shelfmask_sm.
%   INPUT:  run - name of folder for a given ROMS run, directory is already
%                 specified in function; e.g. '016'
%   OUTPUT: volume - volume of every cell in model domain. Given in km^3
%           Also saves a .mat file sm_tools/{run}volume.mat.


filelist = romsinitialize_sm('avg',run);  %Whatever is in the output folder
nfiles = size(filelist,1);
 
Vtransform = nc_varget(filelist(1,:),'Vtransform');
Vstretching = nc_varget(filelist(1,:),'Vstretching');
theta_s = nc_varget(filelist(1,:),'theta_s');
theta_b = nc_varget(filelist(1,:),'theta_b');
h = nc_varget(filelist(1,:), 'h'); %meters
zice = nc_varget(filelist(1,:),'zice'); %Ice shelf thickness
hc = nc_varget(filelist(1,:),'hc'); %critical depth
s_rho = nc_varget(filelist(1,:),'s_rho');
N = length(s_rho);
igrid = 5;  %w velocity points, top & bottom of each layer

zeta = [];
for i =1:nfiles
    temp = nc_varget(filelist(i,:),'zeta');
    if(ndims(temp)<3)
        temp = shiftdim(temp,-1);
    end
    zeta=cat(1,zeta,temp);
end
clear temp

zeta(find(isnan(zeta)))=0; %#ok<*FNDSB>
steps = size(zeta,1);
X = size(zeta,2);  Y = size(zeta,3);
z = zeros(steps,N+1,X,Y);
for i=1:steps
    temp = set_depth(Vtransform, Vstretching, theta_s, theta_b, hc, N, igrid, h, squeeze(zeta(i,:,:))+zice);
    z(i,:,:,:)=shiftdim(temp,2);
    display(['Finished depth i = ' num2str(i)]);
end
clear filelist nfiles Vtransform Vstretching theta_s theta_b h zice hc s_rho
clear igrid zeta temp

if(X>500)
    grdsize = 1.5;
else
    grdsize = 5.0;
end
volume = zeros(steps,N,X,Y);
volume(1:steps,1:N,1:X,1:Y) = abs(z(1:steps,1:N,1:X,1:Y)-z(1:steps,2:N+1,1:X,1:Y));
volume(find(isnan(volume)))=0;
volume = grdsize.*grdsize.*volume./1e3; %m to km3

save(['sm_tools/' run 'Volume.mat'],'volume','-v7.3');

end

