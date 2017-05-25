function [ mld ] = mld_sm( run )
%MLD_SM Takes model output from Ross Sea ROMS simulations and calculates
%mixed layer depth using de Boyer Montegut et al 2004
%   INPUT:  run - run number for ROMS simulation
%   OUTPUT: mld - overall mixed layer depth
%   This function also saves a .mat file [run 'MLD.mat'] format with mixed
%   layer depth for density and temperature, and time.

files = romsinitialize_sm('avg',run);
nfiles = size(files,1);

time = [];
for i=1:nfiles
    time = cat(1,time,nc_varget(files(i,:),'ocean_time'));
end
time = time./3600; time = time./24; %days
time = time -2190; %set 0 to sept 15 2010
time = datenum(2010,9,15+time,0,0,0);

temp = [];
for i=1:nfiles
    temp2 = nc_varget(files(i,:),'temp');
    if (ndims(temp2)<4)
        temp2 = shiftdim(temp2,-1);
    end
    temp = cat(1,temp,temp2);
end

salt = [];
for i=1:nfiles
    temp2 = nc_varget(files(i,:),'salt');
    if (ndims(temp2)<4)
        temp2 = shiftdim(temp2,-1);
    end
   salt = cat(1,salt,temp2);
end

zeta = []; %Needed for depth calculation
for i=1:nfiles
    temp2 = nc_varget(files(i,:),'zeta');
    if (ndims(temp2)<3)
        temp2 = shiftdim(temp2,-1);
    end
    zeta = cat(1,zeta,temp2);
end

%Set depths of layers
Vtransform = nc_varget(files(1,:),'Vtransform');
Vstretching = nc_varget(files(1,:),'Vstretching');
theta_s = nc_varget(files(1,:),'theta_s');
theta_b = nc_varget(files(1,:),'theta_b');
h = nc_varget(files(1,:), 'h'); %meters
hc = nc_varget(files(1,:),'hc'); %critical depth
s_rho = nc_varget(files(1,:),'s_rho');
zice = nc_varget(files(1,:),'zice');
lon = nc_varget(files(1,:),'lon_rho');
lat = nc_varget(files(1,:),'lat_rho');
mask_rho=nc_varget(files(1,:),'mask_rho');
N = length(s_rho);
igrid = 1;

zeta(find(isnan(zeta)))=0; %#ok<*FNDSB>
steps = size(zeta,1);
X = size(zeta,2);  Y = size(zeta,3);
z = zeros(steps,N,X,Y);
for i=1:steps
    temp2 = set_depth(Vtransform, Vstretching, theta_s, theta_b, hc, N, igrid, h, squeeze(zeta(i,:,:))+zice);
    z(i,:,:,:)=shiftdim(temp2,2);
    display(['Finished depth i = ' num2str(i)]);
end
clear Vtransform Vstretching theta_s theta_b hc N igrid s_rho temp2;

% Convert rho to potential rho

rho = zeros(size(z));
for n=1:steps
    for i=1:24
        pres = gsw_p_from_z(squeeze(z(n,i,:,:)),lat);
        SA = gsw_SA_from_SP(squeeze(salt(n,i,:,:)),pres,lon,lat);
        CT = gsw_CT_from_pt(SA,squeeze(temp(n,i,:,:)));
        rho(n,i,:,:) = gsw_rho_CT_exact(SA,CT,10.1325);
    end
end

display('Finished calculating potential density');

%Set surface layer as reference layer - NO, change to 10m layer - check
%this, need more changes to make it work
% tempRef = zeros(steps,X,Y); rhoRef = tempRef;
% for n=1:steps
%     for i=1:X
%         for j=1:Y
%             temp2 = abs(squeeze(z(n,:,i,j)) + 10); %Offset depths by 10m
%             ref = find(min(temp2)); %Find value closest to 0 for ref depth
%             tempRef(n,i,j)=squeeze(temp(n,ref,i,j));
%             rhoRef(n,i,j)=squeeze(rho(n,ref,i,j));
%         end
%     end
% end

%Set reference layer as 1st layer down
tempRef = temp(:,24,:,:);
rhoRef = rho(:,24,:,:);

%Calculate |deltaT| from surface layer and find first point greater than
%0.2 deg C
%Calculate deltaRho from surface and find first point greater than 0.03
%kg/m^3
deltaT = zeros(size(temp)); deltaRho = deltaT;
for i=1:24 %Don't calculate surface - set to 0
    deltaT(:,i,:,:) = abs(temp(:,i,:,:)-tempRef);
    deltaRho(:,i,:,:) = rho(:,i,:,:)-rhoRef;
end
%deltaT(:,24,:,:)=0; deltaRho(:,24,:,:)=0; %Ignore surface layer when calculating MLD

mldR = zeros(steps,X,Y); mldT = mldR;
for n = 1:steps
    for i=1:X
        for j=1:Y
            if(mask_rho(i,j)==0) %If land, skip
                continue;
            end
            dT = squeeze(deltaT(n,:,i,j));
            layer = find(dT>.2,1,'last');
            if(isempty(layer))
                mldT(n,i,j)=squeeze(z(n,1,i,j)); %Set to bottom depth
            else
                mldT(n,i,j)=squeeze(z(n,layer,i,j));
            end
            
            dR = squeeze(deltaRho(n,:,i,j));
            layer = find(dR>0.03,1,'last');
            if(isempty(layer))
                mldR(n,i,j)=squeeze(z(n,1,i,j));
            else
                mldR(n,i,j)=squeeze(z(n,layer,i,j));
            end
        end
    end
end

display('Finished calculating mld. Saving variables...');

%Compare two depths, and choose shallower of two for fully mixed layer
mld = max(mldT,mldR); %max, because z is negative

save([run 'MLD.mat'],'mld', 'mldR','mldT','run','time','-v7.3');

figure;
pcolor(squeeze(mld(1,:,:)));
shading flat
colorbar
title('MLD at time 1');

end

