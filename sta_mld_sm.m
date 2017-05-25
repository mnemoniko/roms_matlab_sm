function [ mld ] = sta_mld_sm( run )
%STA_MLD_SM This function calculates mixed layer depth for ROMS Ross Sea
%station files based on a threshold method
%   INPUT:  run - string indicating which run number to analyze
%   OUTPUT: mld - mixed layer depth for that run at all stations.

file = romsinitialize_sm('sta',run);

%Basic data
Ipos = nc_varget(file,'Ipos');
Jpos = nc_varget(file,'Jpos');
lat = nc_varget(file,'lat_rho');
lon = nc_varget(file,'lon_rho');
staTime = nc_varget(file,'ocean_time');
staTime = staTime./3600; staTime = staTime./24; %days
staTime = staTime -2190; %set 0 to sept 15 2010
staTime = datenum(2010,9,15+staTime,0,0,0);

%To calc MLDs:
temp = nc_varget(file,'temp');
salt = nc_varget(file,'salt');
zeta = nc_varget(file,'zeta');
h = nc_varget(file, 'h'); %meters

%Remove 0 points:
good = find(Ipos~=0);
Ipos = Ipos(good); Jpos = Jpos(good); temp = temp(:,good,:);
salt = salt(:,good,:); zeta = zeta(:,good);
lat = lat(good); lon = lon(good); h = h(good);

%Set depths of layers
Vtransform = nc_varget(file,'Vtransform');
Vstretching = nc_varget(file,'Vstretching');
theta_s = nc_varget(file,'theta_s');
theta_b = nc_varget(file,'theta_b');
hc = nc_varget(file,'hc'); %critical depth
s_rho = nc_varget(file,'s_rho');
N = length(s_rho);
igrid = 1;

zeta(find(isnan(zeta)))=0; %#ok<*FNDSB>
steps = size(zeta,1); sta = size(zeta,2);
z = zeros(steps,sta,N);
for i=1:steps
    for j=1:sta
    temp2 = set_depth_sta_sm(Vtransform, Vstretching, theta_s, theta_b, hc, N, igrid, h(j), zeta(i,j),0);
    z(i,j,:)=temp2;
    end
    %display(['Finished depth i = ' num2str(i)]);
end
clear Vtransform Vstretching theta_s theta_b hc N igrid s_rho X Y steps temp2;
display(['Finished setting depth for run ' run]);

rho = zeros(size(z));
for n=1:size(z,1)
    for i=1:24
        pres = gsw_p_from_z(squeeze(z(n,:,i)),lat);
        SA = gsw_SA_from_SP(squeeze(salt(n,:,i)),pres,lon,lat);
        CT = gsw_CT_from_pt(SA,squeeze(temp(n,:,i)));
        rho(n,:,i) = gsw_rho_CT_exact(SA,CT,10.1325);
    end
end
display(['Finished setting rho for run ' run]);

%Set surface layer as reference layer
tempRef = temp(:,:,24);
rhoRef = rho(:,:,24);

%Calculate |deltaT| from surface layer and find first point greater than
%0.2 deg C
%Calculate deltaRho from surface and find first point greater than 0.03
%kg/m^3
deltaT = zeros(size(temp)); deltaRho = deltaT;
for i=1:24
    deltaT(:,:,i) = abs(temp(:,:,i)-tempRef);
    deltaRho(:,:,i) = rho(:,:,i)-rhoRef;
end
% deltaT(:,:,24)=0; deltaRho(:,:,24)=0;

mldR = zeros(length(staTime),length(Ipos)); mldT = mldR;

for n = 1:length(staTime)
    for i=1:length(Ipos)
            dT = squeeze(deltaT(n,i,:));
            layer = find(dT>.2,1,'last');
            if(isempty(layer))
                mldT(n,i)=squeeze(z(n,i,1)); %Set to bottom depth
            else
                mldT(n,i)=squeeze(z(n,i,layer));
            end
            
            dR = squeeze(deltaRho(n,i,:));
            layer = find(dR>0.03,1,'last');
            if(isempty(layer))
                mldR(n,i)=squeeze(z(n,i,1));
            else
                mldR(n,i)=squeeze(z(n,i,layer));
            end
    end
end

%Compare two depths, and choose shallower of two for fully mixed layer
mld = max(mldT,mldR); %max, because z is negative
display(['Finished calculating mld for run ' run]);
display('Saving data...');

save([run 'MLDsta.mat'],'mld', 'mldR','mldT','run','staTime');

end