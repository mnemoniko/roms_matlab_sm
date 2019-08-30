function [ ] = makeROMSinitial_sm( data, filename, grdfile, iniDescrip )
%makeROMSinitial Create initial file for ROMS application
%   INPUT:  data - structure with fields for zeta, ubar, vbar, u, v, temp,
%   salt
%           filename - name of initial netcdf4 file
%           grdfile - name of associated grid file
%           iniDescription - brief string description of application

%Get relevant information from grid file
h = nc_varget(grdfile,'h');
lm = size(h,2)-2;
mm = size(h,1)-2;
Lp = lm + 2; Mp = mm + 2;
N = nc_getdiminfo(grdfile,'s_rho','Length');

%Create and populate netcdf initial file
nc_create_empty(filename);

%Global attributes
nc_attput(filename,nc_global,'filename',filename);
nc_attput(filename,nc_global,'type','Matlab generated NETCDF4 initial file');
nc_attput(filename,nc_global,'description',iniDescrip);
nc_attput(filename,nc_global,'history',['Matlab version R' version('-release')]);
nc_attput(filename,nc_global,'created',datestr(now));
nc_attput(filename,nc_global,'projection','Polar Stereographic 71S on WGS84 geoid');

%Dimensions
nc_adddim(filename,'xi_rho',Lp);
nc_adddim(filename,'xi_u',Lp-1);
nc_adddim(filename,'xi_v',Lp);
nc_adddim(filename,'eta_rho',Mp);
nc_adddim(filename,'eta_u',Mp);
nc_adddim(filename,'eta_v',Mp-1);
nc_adddim(filename,'s_rho',N);
nc_adddim(filename,'N',N);
nc_adddim(filename,'s_w',N+1);
nc_adddim(filename,'time',1);

%Grid & geometry variables:
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
nc_varput(filename,v1.Name,nc_varget(grdfile,v1.Name));
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
nc_varput(filename,v1.Name,h);
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
nc_varput(filename,v1.Name,nc_varget(grdfile,v1.Name));
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
nc_varput(filename,v1.Name,nc_varget(grdfile,v1.Name));
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
nc_varput(filename,v1.Name,nc_varget(grdfile,v1.Name));
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
nc_varput(filename,v1.Name,nc_varget(grdfile,v1.Name));
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
nc_varput(filename,v1.Name,nc_varget(grdfile,v1.Name));
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
nc_varput(filename,v1.Name,nc_varget(grdfile,v1.Name));
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
nc_varput(filename,v1.Name,nc_varget(grdfile,v1.Name));
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
nc_varput(filename,v1.Name,nc_varget(grdfile,v1.Name));
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
nc_varput(filename,v1.Name,nc_varget(grdfile,v1.Name));
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
nc_varput(filename,v1.Name,nc_varget(grdfile,v1.Name));
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
nc_varput(filename,v1.Name,nc_varget(grdfile,v1.Name));
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
nc_varput(filename,v1.Name,nc_varget(grdfile,v1.Name));
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
nc_varput(filename,v1.Name,nc_varget(grdfile,v1.Name));
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
nc_varput(filename,v1.Name,nc_varget(grdfile,v1.Name));
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
nc_varput(filename,v1.Name,nc_varget(grdfile,v1.Name));
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
nc_varput(filename,v1.Name,nc_varget(grdfile,v1.Name));
clear v1

%Vertical stretching variables
v1.Name = 'theta_s';
v1.Datatype = 'double';
v1.Attribute(1).Name = 'long_name';
v1.Attribute(1).Value = 'S-coordinate surface control parameter';
v1.Dimension = {};
v1.Chunking = [];
v1.Deflate = 0;
nc_addvar(filename,v1);
nc_varput(filename,v1.Name,nc_varget(grdfile,v1.Name));
clear v1

v1.Name = 'theta_b';
v1.Datatype = 'double';
v1.Attribute(1).Name = 'long_name';
v1.Attribute(1).Value = 'S-coordinate bottom control parameter';
v1.Dimension = {};
v1.Chunking = [];
v1.Deflate = 0;
nc_addvar(filename,v1);
nc_varput(filename,v1.Name,nc_varget(grdfile,v1.Name));
clear v1

v1.Name = 'Vtransform';
v1.Datatype = 'int';
v1.Attribute(1).Name = 'long_name';
v1.Attribute(1).Value = 'vertical terrain-following transformation equation';
v1.Dimension = {};
v1.Chunking = [];
v1.Deflate = 0;
nc_addvar(filename,v1);
nc_varput(filename,v1.Name,nc_varget(grdfile,v1.Name));
clear v1

v1.Name = 'Vstretching';
v1.Datatype = 'double';
v1.Attribute(1).Name = 'long_name';
v1.Attribute(1).Value = 'Vertical terrain-following stretching function';
v1.Dimension = {};
v1.Chunking = [];
v1.Deflate = 0;
nc_addvar(filename,v1);
nc_varput(filename,v1.Name,nc_varget(grdfile,v1.Name));
clear v1

v1.Name = 'Tcline';
v1.Datatype = 'double';
v1.Attribute(1).Name = 'long_name';
v1.Attribute(1).Value = 'S-coordinate parameter surface/bottom layer width';
v1.Attribute(2).Name = 'units';
v1.Attribute(2).Value = 'meter';
v1.Dimension = {};
v1.Chunking = [];
v1.Deflate = 0;
nc_addvar(filename,v1);
nc_varput(filename,v1.Name,nc_varget(grdfile,v1.Name));
clear v1

v1.Name = 'hc';
v1.Datatype = 'double';
v1.Attribute(1).Name = 'long_name';
v1.Attribute(1).Value = 'S-coordinate parameter, critical depth';
v1.Attribute(2).Name = 'units';
v1.Attribute(2).Value = 'meter';
v1.Dimension = {};
v1.Chunking = [];
v1.Deflate = 0;
nc_addvar(filename,v1);
nc_varput(filename,v1.Name,nc_varget(grdfile,v1.Name));
clear v1

v1.Name = 's_rho';
v1.Datatype = 'double';
v1.Attribute(1).Name = 'long_name';
v1.Attribute(1).Value = 'S-coordinate at RHO-points';
v1.Attribute(2).Name = 'valid_min';
v1.Attribute(2).Value = -1.0;
v1.Attribute(3).Name = 'valid_max';
v1.Attribute(3).Value = 0.;
v1.Attribute(4).Name = 'positive';
v1.Attribute(4).Value = 'up';
v1.Attribute(5).Name = 'standard_name';
v1.Attribute(5).Value = 'ocean_s_coordinate_g2';
v1.Attribute(6).Name = 'formula_terms';
v1.Attribute(6).Value = 's: s_rho C: Cs_r eta: zeta depth: h depth_c: hc';
v1.Attribute(7).Name = 'field';
v1.Attribute(7).Value = 's_rho, scalar';
v1.Dimension = {'s_rho'};
v1.Chunking = [];
v1.Deflate = 0;
nc_addvar(filename,v1);
nc_varput(filename,v1.Name,nc_varget(grdfile,v1.Name));
clear v1

v1.Name = 's_w';
v1.Datatype = 'double';
v1.Attribute(1).Name = 'long_name';
v1.Attribute(1).Value = 'S-coordinate at W-points';
v1.Attribute(2).Name = 'valid_min';
v1.Attribute(2).Value = -1.0;
v1.Attribute(3).Name = 'valid_max';
v1.Attribute(3).Value = 0.;
v1.Attribute(4).Name = 'positive';
v1.Attribute(4).Value = 'up';
v1.Attribute(5).Name = 'standard_name';
v1.Attribute(5).Value = 'ocean_s_coordinate_g2';
v1.Attribute(6).Name = 'formula_terms';
v1.Attribute(6).Value = 's: s_w C: Cs_w eta: zeta depth: h depth_c: hc';
v1.Attribute(7).Name = 'field';
v1.Attribute(7).Value = 's_w, scalar';
v1.Dimension = {'s_w'};
v1.Chunking = [];
v1.Deflate = 0;
nc_addvar(filename,v1);
nc_varput(filename,v1.Name,nc_varget(grdfile,v1.Name));
clear v1

v1.Name = 'Cs_r';
v1.Datatype = 'double';
v1.Attribute(1).Name = 'long_name';
v1.Attribute(1).Value = 'S-coordinate stretching curves at RHO-points';
v1.Attribute(2).Name = 'valid_min';
v1.Attribute(2).Value = -1.0;
v1.Attribute(3).Name = 'valid_max';
v1.Attribute(3).Value = 0.;
v1.Attribute(4).Name = 'field';
v1.Attribute(4).Value = 'Cs_r, scalar';
v1.Dimension = {'s_rho'};
v1.Chunking = [];
v1.Deflate = 0;
nc_addvar(filename,v1);
nc_varput(filename,v1.Name,nc_varget(grdfile,v1.Name));
clear v1

v1.Name = 'Cs_w';
v1.Datatype = 'double';
v1.Attribute(1).Name = 'long_name';
v1.Attribute(1).Value = 'S-coordinate stretching curves at W-points';
v1.Attribute(2).Name = 'valid_min';
v1.Attribute(2).Value = -1.0;
v1.Attribute(3).Name = 'valid_max';
v1.Attribute(3).Value = 0.;
v1.Attribute(4).Name = 'field';
v1.Attribute(4).Value = 'Cs_w, scalar';
v1.Dimension = {'s_w'};
v1.Chunking = [];
v1.Deflate = 0;
nc_addvar(filename,v1);
nc_varput(filename,v1.Name,nc_varget(grdfile,v1.Name));
clear v1

%Initial data variables:
v1.Name = 'ocean_time';
v1.Datatype = 'double';
v1.Attribute(1).Name = 'long_name';
v1.Attribute(1).Value = 'time since initialization';
v1.Attribute(2).Name = 'units';
v1.Attribute(2).Value = 'seconds since 0001-01-01 00:00:00';
v1.Attribute(3).Name = 'calendar';
v1.Attribute(3).Value = '365.25 days in every year';
v1.Attribute(4).Name = 'field';
v1.Attribute(4).Value = 'time, scalar, series';
v1.Dimension = {'time'};
v1.Chunking = [];
v1.Deflate = 0;
nc_addvar(filename,v1);
nc_varput(filename,v1.Name,0.00);
clear v1

v1.Name = 'zeta';
v1.Datatype = 'double';
v1.Attribute(1).Name = 'long_name';
v1.Attribute(1).Value = 'free-surface';
v1.Attribute(2).Name = 'units';
v1.Attribute(2).Value = 'meter';
v1.Attribute(3).Name = 'time';
v1.Attribute(3).Value = 'ocean_time';
v1.Attribute(4).Name = 'field';
v1.Attribute(4).Value = 'free-surface, scalar, series';
v1.Dimension = {'time','eta_rho','xi_rho'};
v1.Chunking = [];
v1.Deflate = 0;
nc_addvar(filename,v1);
nc_varput(filename,v1.Name,data.zeta);
clear v1

v1.Name = 'ubar';
v1.Datatype = 'double';
v1.Attribute(1).Name = 'long_name';
v1.Attribute(1).Value = 'vertically integrated u-momentum component';
v1.Attribute(2).Name = 'units';
v1.Attribute(2).Value = 'meter second-1';
v1.Attribute(3).Name = 'time';
v1.Attribute(3).Value = 'ocean_time';
v1.Attribute(4).Name = 'field';
v1.Attribute(4).Value = 'ubar-velocity, scalar, series';
v1.Dimension = {'time','eta_u','xi_u'};
v1.Chunking = [];
v1.Deflate = 0;
nc_addvar(filename,v1);
nc_varput(filename,v1.Name,data.ubar);
clear v1

v1.Name = 'vbar';
v1.Datatype = 'double';
v1.Attribute(1).Name = 'long_name';
v1.Attribute(1).Value = 'vertically integrated v-momentum component';
v1.Attribute(2).Name = 'units';
v1.Attribute(2).Value = 'meter second-1';
v1.Attribute(3).Name = 'time';
v1.Attribute(3).Value = 'ocean_time';
v1.Attribute(4).Name = 'field';
v1.Attribute(4).Value = 'vbar-velocity, scalar, series';
v1.Dimension = {'time','eta_v','xi_v'};
v1.Chunking = [];
v1.Deflate = 0;
nc_addvar(filename,v1);
nc_varput(filename,v1.Name,data.vbar);
clear v1

v1.Name = 'u';
v1.Datatype = 'double';
v1.Attribute(1).Name = 'long_name';
v1.Attribute(1).Value = 'u-momentum component';
v1.Attribute(2).Name = 'units';
v1.Attribute(2).Value = 'meter second-1';
v1.Attribute(3).Name = 'time';
v1.Attribute(3).Value = 'ocean_time';
v1.Attribute(4).Name = 'field';
v1.Attribute(4).Value = 'u-velocity, scalar, series';
v1.Dimension = {'time','s_rho','eta_u','xi_u'};
v1.Chunking = [];
v1.Deflate = 0;
nc_addvar(filename,v1);
nc_varput(filename,v1.Name,data.u);
clear v1

v1.Name = 'v';
v1.Datatype = 'double';
v1.Attribute(1).Name = 'long_name';
v1.Attribute(1).Value = 'v-momentum component';
v1.Attribute(2).Name = 'units';
v1.Attribute(2).Value = 'meter second-1';
v1.Attribute(3).Name = 'time';
v1.Attribute(3).Value = 'ocean_time';
v1.Attribute(4).Name = 'field';
v1.Attribute(4).Value = 'v-velocity, scalar, series';
v1.Dimension = {'time','s_rho','eta_v','xi_v'};
v1.Chunking = [];
v1.Deflate = 0;
nc_addvar(filename,v1);
nc_varput(filename,v1.Name,data.v);
clear v1

v1.Name = 'temp';
v1.Datatype = 'double';
v1.Attribute(1).Name = 'long_name';
v1.Attribute(1).Value = 'potential temperature';
v1.Attribute(2).Name = 'units';
v1.Attribute(2).Value = 'Celsius';
v1.Attribute(3).Name = 'time';
v1.Attribute(3).Value = 'ocean_time';
v1.Attribute(4).Name = 'field';
v1.Attribute(4).Value = 'temperature, scalar, series';
v1.Dimension = {'time','s_rho','eta_rho','xi_rho'};
v1.Chunking = [];
v1.Deflate = 0;
nc_addvar(filename,v1);
nc_varput(filename,v1.Name,data.temp);
clear v1

v1.Name = 'salt';
v1.Datatype = 'double';
v1.Attribute(1).Name = 'long_name';
v1.Attribute(1).Value = 'salinity';
v1.Attribute(2).Name = 'units';
v1.Attribute(2).Value = 'PSU';
v1.Attribute(3).Name = 'time';
v1.Attribute(3).Value = 'ocean_time';
v1.Attribute(4).Name = 'field';
v1.Attribute(4).Value = 'salinity, scalar, series';
v1.Dimension = {'time','s_rho','eta_rho','xi_rho'};
v1.Chunking = [];
v1.Deflate = 0;
nc_addvar(filename,v1);
nc_varput(filename,v1.Name,data.salt);
clear v1

% v1.Name = 
% v1.Datatype = 
% v1.Attribute(1).Name = 
% v1.Attribute(1).Value = 
% v1.Attribute(2).Name = 
% v1.Attribute(2).Value = 
% v1.Attribute(3).Name = 
% v1.Attribute(3).Value =
% v1.Attribute(4).Name = 
% v1.Attribute(4).Value = 
% v1.Dimension = 
% v1.Chunking = [];
% v1.Deflate = 0;
% nc_addvar(filename,v1);
% nc_varput(filename,v1.Name,nc_varget(grdfile,v1.Name));
% clear v1

end

