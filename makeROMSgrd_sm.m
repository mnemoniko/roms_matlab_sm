function [ ] = makeROMSgrd_sm( H, lat_rho, lon_rho, mask_rho, zice, dl, filename, grdDescrip )
%makeROMSgrd_sm Creates ROMS x/y grid file (evenly spaced) based on given
%inputs
%   INPUTS:
%           H - bathymetry, positive
%           lat_rho - latitude at rho points
%           lon_rho - longitude at rho points
%           mask_rho - land mask at rho points
%           zice - ice shelf bottom bathymetry
%           dl - dx and dy in meters
%           filename - string file name for new grid file
%           grdDescrip - Brief description of grid

tmin = min(H(:));

lm = size(H,2)-2;
mm = size(H,1)-2;
f = 4.*pi./3600./24.*sind(lat_rho); %coriolis
dfdy = nanmean((f(2:end,1)-f(1:end-1,1))./dl);

Lp = lm + 2; Mp = mm + 2;
xl = lm*dl; %Length of domain in xi-direction (meters)
el = mm*dl; %Same in the eta-direction

xi_rho  = (- 0.5 : lm + 0.5)';
eta_rho = (- 0.5 : mm + 0.5)';
xi_u    = (0 : lm)';
eta_u   = eta_rho;
xi_v    = xi_rho;
eta_v   = (0 : mm)';
xi_psi  = xi_u;
eta_psi = eta_v;

pm = ones(Mp,Lp)./dl; pn = pm;

[x_rho, y_rho] = meshgrid( xi_rho * dl, eta_rho * dl );
[x_psi, y_psi] = meshgrid( xi_psi * dl, eta_psi * dl );
[x_u,   y_u  ] = meshgrid( xi_u   * dl,   eta_u * dl );
[x_v,   y_v  ] = meshgrid( xi_v   * dl,   eta_v * dl );

lon_v   = lon_rho(2 : end, :);
lat_u   = lat_rho(:, 2 : end);
lat_psi = lat_rho - 0.5 * nanmean(lat_rho(2 : end, 1)  - lat_rho(1 : end - 1, 1) );
lon_psi = lon_rho - 0.5 * nanmean(lon_rho(1, 2 : end)' - lon_rho(1, 1 : end - 1)');
lat_psi = lat_psi(2 : end, 2 : end);
lon_psi = lon_psi(2 : end, 2 : end);
lat_v   = lat_rho - 0.5 * nanmean(lat_rho(2 : end, 1)  - lat_rho(1 : end - 1, 1) );
lat_v   = lat_v(2 : end, :);
lon_u   = lon_rho - 0.5 * nanmean(lon_rho(1, 2 : end)' - lon_rho(1, 1 : end - 1)');
lon_u   = lon_u(:, 2 : end);

dndx    = nan(  Mp, Lp);
dmde    = nan(  Mp, Lp);
dndx(:, 2 : end - 1) = (pm(:, 3 : end) - pm(:, 1 : end - 2)) / 2.;
dndx(:,           1) = dndx(:,2);
dndx(:,         end) = dndx(:, end - 1);
dmde(2 : end - 1, :) = (pn(3 : end, :) - pn(1 : end - 2, :)) / 2.;
dmde(  1,         :) = dmde(2, :);
dmde(end,         :) = dmde(end - 1, :);

angler         = nan(Mp, Lp);
angler(2 : end - 1, 2 : end - 1) ...
  = atan2( (lat_rho(3 : end, 2 : end - 1) - lat_rho(1 : end - 2, 2 : end - 1)) / dl, ...
           (lat_rho(2 : end - 1, 3 : end) - lat_rho(2 : end - 1, 1 : end - 2)) / dl ) ...
  - pi / 2.;
angler(:,   1) = angler(:,       2);
angler(:, end) = angler(:, end - 1);
angler(  1, :) = angler(      2, :);
angler(end, :) = angler(end - 1, :);

[vmask,umask,pmask] = uvp_masks(mask_rho);

%% Create ROMS grid file

nc_create_empty(filename);

%Global attributes
nc_attput(filename,nc_global,'filename',filename);
nc_attput(filename,nc_global,'type','Matlab generated NETCDF4 grid file');
nc_attput(filename,nc_global,'gridid',grdDescrip);
nc_attput(filename,nc_global,'history',['Matlab version R' version('-release')]);
nc_attput(filename,nc_global,'created',datestr(now));
nc_attput(filename,nc_global,'projection','Polar Stereographic 71S on WGS84 geoid');

%Dimensions
nc_adddim(filename,'xi_psi',Lp-1);
nc_adddim(filename,'xi_rho',Lp);
nc_adddim(filename,'xi_u',Lp-1);
nc_adddim(filename,'xi_v',Lp);
nc_adddim(filename,'eta_psi',Mp-1);
nc_adddim(filename,'eta_rho',Mp);
nc_adddim(filename,'eta_u',Mp);
nc_adddim(filename,'eta_v',Mp-1);

%Variables:
v1.Name = 'xl';
v1.Datatype = 'double';
v1.Attribute(1).Name = 'long_name';
v1.Attribute(1).Value = 'domain length in the XI-direction';
v1.Attribute(2).Name = 'units';
v1.Attribute(2).Value = 'meter';
v1.Dimension = {};
v1.Chunking = [];
v1.Deflate = 0;
nc_addvar(filename,v1);
nc_varput(filename,v1.Name,xl);
clear v1

v1.Name = 'el';
v1.Datatype = 'double';
v1.Attribute(1).Name = 'long_name';
v1.Attribute(1).Value = 'domain length in the ETA-direction';
v1.Attribute(2).Name = 'units';
v1.Attribute(2).Value = 'meter';
v1.Dimension = {};
v1.Chunking = [];
v1.Deflate = 0;
nc_addvar(filename,v1);
nc_varput(filename,v1.Name,el);
clear v1

v1.Name = 'depthmin';
v1.Datatype = 'double';
v1.Attribute(1).Name = 'long_name';
v1.Attribute(1).Value = 'Shallow bathymetry clipping depth';
v1.Attribute(2).Name = 'units';
v1.Attribute(2).Value = 'meter';
v1.Dimension = {};
v1.Chunking = [];
v1.Deflate = 0;
nc_addvar(filename,v1);
nc_varput(filename,v1.Name,tmin);
clear v1

v1.Name = 'depthmax';
v1.Datatype = 'double';
v1.Attribute(1).Name = 'long_name';
v1.Attribute(1).Value = 'Deep bathymetry clipping depth';
v1.Attribute(2).Name = 'units';
v1.Attribute(2).Value = 'meter';
v1.Dimension = {};
v1.Chunking = [];
v1.Deflate = 0;
nc_addvar(filename,v1);
nc_varput(filename,v1.Name,10000);
clear v1

v1.Name = 'spherical';
v1.Datatype = 'char';
v1.Attribute(1).Name = 'long_name';
v1.Attribute(1).Value = 'Grid type logical switch';
v1.Attribute(2).Name = 'option_T';
v1.Attribute(2).Value = 'spherical';
v1.Attribute(3).Name = 'option_F';
v1.Attribute(3).Value = 'cartesian';
v1.Dimension = {};
v1.Chunking = [];
v1.Deflate = 0;
nc_addvar(filename,v1);
nc_varput(filename,v1.Name,'F');
clear v1

v1.Name = 'f0';
v1.Datatype = 'double';
v1.Attribute(1).Name = 'long_name';
v1.Attribute(1).Value = 'Coriolis parameter central value on a beta-plane';
v1.Attribute(2).Name = '_FillValue';
v1.Attribute(2).Value = 0.;
v1.Dimension = {};
v1.Chunking = [];
v1.Deflate = 0;
nc_addvar(filename,v1);
nc_varput(filename,v1.Name,nanmean(f(:)));
clear v1

v1.Name = 'dfdy';
v1.Datatype = 'double';
v1.Attribute(1).Name = 'long_name';
v1.Attribute(1).Value = 'Coriolis parameter gradient on a beta-plane';
v1.Attribute(2).Name = '_FillValue';
v1.Attribute(2).Value = 0.;
v1.Dimension = {};
v1.Chunking = [];
v1.Deflate = 0;
nc_addvar(filename,v1);
nc_varput(filename,v1.Name,dfdy);
clear v1

v1.Name = 'h';
v1.Datatype = 'double';
v1.Attribute(1).Name = 'long_name';
v1.Attribute(1).Value = 'Bathymetry at RHO-points';
v1.Attribute(2).Name = 'units';
v1.Attribute(2).Value = 'meter';
v1.Attribute(3).Name = 'field';
v1.Attribute(3).Value = 'bath, scalar';
v1.Dimension = {'eta_rho','xi_rho'};
v1.Chunking = [];
v1.Deflate = 0;
nc_addvar(filename,v1);
nc_varput(filename,v1.Name,H);
clear v1

v1.Name = 'f';
v1.Datatype = 'double';
v1.Attribute(1).Name = 'long_name';
v1.Attribute(1).Value = 'Coriolis at RHO-points';
v1.Attribute(2).Name = 'units';
v1.Attribute(2).Value = 'rad/s';
v1.Attribute(3).Name = 'field';
v1.Attribute(3).Value = 'f, scalar';
v1.Dimension = {'eta_rho','xi_rho'};
v1.Chunking = [];
v1.Deflate = 0;
nc_addvar(filename,v1);
nc_varput(filename,v1.Name,f);
clear v1

v1.Name = 'pm';
v1.Datatype = 'double';
v1.Attribute(1).Name = 'long_name';
v1.Attribute(1).Value = 'curvilinear coordinate metric in XI';
v1.Attribute(2).Name = 'units';
v1.Attribute(2).Value = 'meter-1';
v1.Attribute(3).Name = 'field';
v1.Attribute(3).Value = 'pm, scalar';
v1.Dimension = {'eta_rho','xi_rho'};
v1.Chunking = [];
v1.Deflate = 0;
nc_addvar(filename,v1);
nc_varput(filename,v1.Name,pm);
clear v1

v1.Name = 'pn';
v1.Datatype = 'double';
v1.Attribute(1).Name = 'long_name';
v1.Attribute(1).Value = 'curvilinear coordinate metric in ETA';
v1.Attribute(2).Name = 'units';
v1.Attribute(2).Value = 'meter-1';
v1.Attribute(3).Name = 'field';
v1.Attribute(3).Value = 'pn, scalar';
v1.Dimension = {'eta_rho','xi_rho'};
v1.Chunking = [];
v1.Deflate = 0;
nc_addvar(filename,v1);
nc_varput(filename,v1.Name,pn);
clear v1

v1.Name = 'dndx';
v1.Datatype = 'double';
v1.Attribute(1).Name = 'long_name';
v1.Attribute(1).Value = 'XI derivative of inverse metric factor pn';
v1.Attribute(2).Name = 'units';
v1.Attribute(2).Value = 'meter';
v1.Attribute(3).Name = 'field';
v1.Attribute(3).Value = 'dndx, scalar';
v1.Dimension = {'eta_rho','xi_rho'};
v1.Chunking = [];
v1.Deflate = 0;
nc_addvar(filename,v1);
nc_varput(filename,v1.Name,dndx);
clear v1

v1.Name = 'dmde';
v1.Datatype = 'double';
v1.Attribute(1).Name = 'long_name';
v1.Attribute(1).Value = 'ETA derivative of inverse metric factor pm';
v1.Attribute(2).Name = 'units';
v1.Attribute(2).Value = 'meter';
v1.Attribute(3).Name = 'field';
v1.Attribute(3).Value = 'dmde, scalar';
v1.Dimension = {'eta_rho','xi_rho'};
v1.Chunking = [];
v1.Deflate = 0;
nc_addvar(filename,v1);
nc_varput(filename,v1.Name,dmde);
clear v1

v1.Name = 'x_rho';
v1.Datatype = 'double';
v1.Attribute(1).Name = 'long_name';
v1.Attribute(1).Value = 'x location of RHO-points';
v1.Attribute(2).Name = 'units';
v1.Attribute(2).Value = 'meter';
v1.Dimension = {'eta_rho','xi_rho'};
v1.Chunking = [];
v1.Deflate = 0;
nc_addvar(filename,v1);
nc_varput(filename,v1.Name,x_rho);
clear v1

v1.Name = 'y_rho';
v1.Datatype = 'double';
v1.Attribute(1).Name = 'long_name';
v1.Attribute(1).Value = 'y location of RHO-points';
v1.Attribute(2).Name = 'units';
v1.Attribute(2).Value = 'meter';
v1.Dimension = {'eta_rho','xi_rho'};
v1.Chunking = [];
v1.Deflate = 0;
nc_addvar(filename,v1);
nc_varput(filename,v1.Name,y_rho);
clear v1

v1.Name = 'x_psi';
v1.Datatype = 'double';
v1.Attribute(1).Name = 'long_name';
v1.Attribute(1).Value = 'x location of PSI-points';
v1.Attribute(2).Name = 'units';
v1.Attribute(2).Value = 'meter';
v1.Dimension = {'eta_psi','xi_psi'};
v1.Chunking = [];
v1.Deflate = 0;
nc_addvar(filename,v1);
nc_varput(filename,v1.Name,x_psi);
clear v1

v1.Name = 'y_psi';
v1.Datatype = 'double';
v1.Attribute(1).Name = 'long_name';
v1.Attribute(1).Value = 'y location of PSI-points';
v1.Attribute(2).Name = 'units';
v1.Attribute(2).Value = 'meter';
v1.Dimension = {'eta_psi','xi_psi'};
v1.Chunking = [];
v1.Deflate = 0;
nc_addvar(filename,v1);
nc_varput(filename,v1.Name,y_psi);
clear v1

v1.Name = 'x_u';
v1.Datatype = 'double';
v1.Attribute(1).Name = 'long_name';
v1.Attribute(1).Value = 'x location of U-points';
v1.Attribute(2).Name = 'units';
v1.Attribute(2).Value = 'meter';
v1.Dimension = {'eta_u','xi_u'};
v1.Chunking = [];
v1.Deflate = 0;
nc_addvar(filename,v1);
nc_varput(filename,v1.Name,x_u);
clear v1

v1.Name = 'y_u';
v1.Datatype = 'double';
v1.Attribute(1).Name = 'long_name';
v1.Attribute(1).Value = 'y location of U-points';
v1.Attribute(2).Name = 'units';
v1.Attribute(2).Value = 'meter';
v1.Dimension = {'eta_u','xi_u'};
v1.Chunking = [];
v1.Deflate = 0;
nc_addvar(filename,v1);
nc_varput(filename,v1.Name,y_u);
clear v1

v1.Name = 'x_v';
v1.Datatype = 'double';
v1.Attribute(1).Name = 'long_name';
v1.Attribute(1).Value = 'x location of V-points';
v1.Attribute(2).Name = 'units';
v1.Attribute(2).Value = 'meter';
v1.Dimension = {'eta_v','xi_v'};
v1.Chunking = [];
v1.Deflate = 0;
nc_addvar(filename,v1);
nc_varput(filename,v1.Name,x_v);
clear v1

v1.Name = 'y_v';
v1.Datatype = 'double';
v1.Attribute(1).Name = 'long_name';
v1.Attribute(1).Value = 'y location of V-points';
v1.Attribute(2).Name = 'units';
v1.Attribute(2).Value = 'meter';
v1.Dimension = {'eta_v','xi_v'};
v1.Chunking = [];
v1.Deflate = 0;
nc_addvar(filename,v1);
nc_varput(filename,v1.Name,y_v);
clear v1

v1.Name = 'lat_rho';
v1.Datatype = 'double';
v1.Attribute(1).Name = 'long_name';
v1.Attribute(1).Value = 'latitude of RHO-points';
v1.Attribute(2).Name = 'units';
v1.Attribute(2).Value = 'meter';
v1.Dimension = {'eta_rho','xi_rho'};
v1.Chunking = [];
v1.Deflate = 0;
nc_addvar(filename,v1);
nc_varput(filename,v1.Name,lat_rho);
clear v1

v1.Name = 'lon_rho';
v1.Datatype = 'double';
v1.Attribute(1).Name = 'long_name';
v1.Attribute(1).Value = 'longitude of RHO-points';
v1.Attribute(2).Name = 'units';
v1.Attribute(2).Value = 'meter';
v1.Dimension = {'eta_rho','xi_rho'};
v1.Chunking = [];
v1.Deflate = 0;
nc_addvar(filename,v1);
nc_varput(filename,v1.Name,lon_rho);
clear v1

v1.Name = 'lat_psi';
v1.Datatype = 'double';
v1.Attribute(1).Name = 'long_name';
v1.Attribute(1).Value = 'latitude of PSI-points';
v1.Attribute(2).Name = 'units';
v1.Attribute(2).Value = 'meter';
v1.Dimension = {'eta_psi','xi_psi'};
v1.Chunking = [];
v1.Deflate = 0;
nc_addvar(filename,v1);
nc_varput(filename,v1.Name,lat_psi);
clear v1

v1.Name = 'lon_psi';
v1.Datatype = 'double';
v1.Attribute(1).Name = 'long_name';
v1.Attribute(1).Value = 'longitude of PSI-points';
v1.Attribute(2).Name = 'units';
v1.Attribute(2).Value = 'meter';
v1.Dimension = {'eta_psi','xi_psi'};
v1.Chunking = [];
v1.Deflate = 0;
nc_addvar(filename,v1);
nc_varput(filename,v1.Name,lon_psi);
clear v1

v1.Name = 'lat_u';
v1.Datatype = 'double';
v1.Attribute(1).Name = 'long_name';
v1.Attribute(1).Value = 'latitude of U-points';
v1.Attribute(2).Name = 'units';
v1.Attribute(2).Value = 'meter';
v1.Dimension = {'eta_u','xi_u'};
v1.Chunking = [];
v1.Deflate = 0;
nc_addvar(filename,v1);
nc_varput(filename,v1.Name,lat_u);
clear v1

v1.Name = 'lon_u';
v1.Datatype = 'double';
v1.Attribute(1).Name = 'long_name';
v1.Attribute(1).Value = 'longitude of U-points';
v1.Attribute(2).Name = 'units';
v1.Attribute(2).Value = 'meter';
v1.Dimension = {'eta_u','xi_u'};
v1.Chunking = [];
v1.Deflate = 0;
nc_addvar(filename,v1);
nc_varput(filename,v1.Name,lon_u);
clear v1

v1.Name = 'lat_v';
v1.Datatype = 'double';
v1.Attribute(1).Name = 'long_name';
v1.Attribute(1).Value = 'latitude of V-points';
v1.Attribute(2).Name = 'units';
v1.Attribute(2).Value = 'meter';
v1.Dimension = {'eta_v','xi_v'};
v1.Chunking = [];
v1.Deflate = 0;
nc_addvar(filename,v1);
nc_varput(filename,v1.Name,lat_v);
clear v1

v1.Name = 'lon_v';
v1.Datatype = 'double';
v1.Attribute(1).Name = 'long_name';
v1.Attribute(1).Value = 'longitude of V-points';
v1.Attribute(2).Name = 'units';
v1.Attribute(2).Value = 'meter';
v1.Dimension = {'eta_v','xi_v'};
v1.Chunking = [];
v1.Deflate = 0;
nc_addvar(filename,v1);
nc_varput(filename,v1.Name,lon_v);
clear v1

v1.Name = 'mask_rho';
v1.Datatype = 'double';
v1.Attribute(1).Name = 'long_name';
v1.Attribute(1).Value = 'mask on RHO-points';
v1.Attribute(2).Name = 'option_0';
v1.Attribute(2).Value = 'land';
v1.Attribute(3).Name = 'option_1';
v1.Attribute(3).Value = 'water';
v1.Attribute(4).Name = '_FillValue';
v1.Attribute(4).Value = -9999.;
v1.Dimension = {'eta_rho','xi_rho'};
v1.Chunking = [];
v1.Deflate = 0;
nc_addvar(filename,v1);
nc_varput(filename,v1.Name,mask_rho);
clear v1

v1.Name = 'mask_u';
v1.Datatype = 'double';
v1.Attribute(1).Name = 'long_name';
v1.Attribute(1).Value = 'mask on U-points';
v1.Attribute(2).Name = 'option_0';
v1.Attribute(2).Value = 'land';
v1.Attribute(3).Name = 'option_1';
v1.Attribute(3).Value = 'water';
v1.Attribute(4).Name = '_FillValue';
v1.Attribute(4).Value = -9999.;
v1.Dimension = {'eta_u','xi_u'};
v1.Chunking = [];
v1.Deflate = 0;
nc_addvar(filename,v1);
nc_varput(filename,v1.Name,umask);
clear v1

v1.Name = 'mask_v';
v1.Datatype = 'double';
v1.Attribute(1).Name = 'long_name';
v1.Attribute(1).Value = 'mask on V-points';
v1.Attribute(2).Name = 'option_0';
v1.Attribute(2).Value = 'land';
v1.Attribute(3).Name = 'option_1';
v1.Attribute(3).Value = 'water';
v1.Attribute(4).Name = '_FillValue';
v1.Attribute(4).Value = -9999.;
v1.Dimension = {'eta_v','xi_v'};
v1.Chunking = [];
v1.Deflate = 0;
nc_addvar(filename,v1);
nc_varput(filename,v1.Name,vmask);
clear v1

v1.Name = 'mask_psi';
v1.Datatype = 'double';
v1.Attribute(1).Name = 'long_name';
v1.Attribute(1).Value = 'mask on PSI-points';
v1.Attribute(2).Name = 'option_0';
v1.Attribute(2).Value = 'land';
v1.Attribute(3).Name = 'option_1';
v1.Attribute(3).Value = 'water';
v1.Attribute(4).Name = '_FillValue';
v1.Attribute(4).Value = -9999.;
v1.Dimension = {'eta_psi','xi_psi'};
v1.Chunking = [];
v1.Deflate = 0;
nc_addvar(filename,v1);
nc_varput(filename,v1.Name,pmask);
clear v1

v1.Name = 'angle';
v1.Datatype = 'double';
v1.Attribute(1).Name = 'long_name';
v1.Attribute(1).Value = 'angle between xi axis and east';
v1.Attribute(2).Name = 'units';
v1.Attribute(2).Value = 'degree';
v1.Dimension = {'eta_rho','xi_rho'};
v1.Chunking = [];
v1.Deflate = 0;
nc_addvar(filename,v1);
nc_varput(filename,v1.Name,angler);
clear v1

v1.Name = 'zice';
v1.Datatype = 'double';
v1.Attribute(1).Name = 'long_name';
v1.Attribute(1).Value = 'Ice shelf depth (below sea level)';
v1.Attribute(2).Name = 'units';
v1.Attribute(2).Value = 'meter';
v1.Dimension = {'eta_rho','xi_rho'};
v1.Chunking = [];
v1.Deflate = 0;
nc_addvar(filename,v1);
nc_varput(filename,v1.Name,zice);
clear v1

v1.Name = 'zice_Full';
v1.Datatype = 'double';
v1.Attribute(1).Name = 'long_name';
v1.Attribute(1).Value = 'Ice shelf depth (incl. land values)';
v1.Attribute(2).Name = 'units';
v1.Attribute(2).Value = 'meter';
v1.Dimension = {'eta_rho','xi_rho'};
v1.Chunking = [];
v1.Deflate = 0;
nc_addvar(filename,v1);
nc_varput(filename,v1.Name,zice);
clear v1


end

