function [  ] = addROMSvert_sm( Vtransform, Vstretching, N, theta_s, theta_b, Tcline, file )
%addROMSvert_sm Add information on vertical stretching parameters to netcdf
%file
%   INPUT:  Vtransform; Vstretching; N; theta_s; theta_b; Tcline - typical
%   parameters from ROMS grid stretching.  See ROMS wiki on appropriate
%   values.
%           file - netcdf file to write out vertical information to
%           (typically a grid file).

H = nc_varget(file,'h');

if(Vtransform==1)
    hc = nanmin([H(:); Tcline]);
else
    hc = Tcline;
end

[s_rho,Cs_r] = stretching(Vstretching, theta_s, theta_b, hc, N, 0,1);
[s_w,Cs_w] = stretching(Vstretching, theta_s, theta_b, hc, N, 1, 1);

nc_adddim(file,'s_rho',N);
nc_adddim(file,'N',N);
nc_adddim(file,'s_w',N+1);

v1.Name = 'theta_s';
v1.Datatype = 'double';
v1.Attribute(1).Name = 'long_name';
v1.Attribute(1).Value = 'S-coordinate surface control parameter';
v1.Dimension = {};
v1.Chunking = [];
v1.Deflate = 0;
nc_addvar(file,v1);
nc_varput(file,v1.Name,theta_s);
clear v1

v1.Name = 'theta_b';
v1.Datatype = 'double';
v1.Attribute(1).Name = 'long_name';
v1.Attribute(1).Value = 'S-coordinate bottom control parameter';
v1.Dimension = {};
v1.Chunking = [];
v1.Deflate = 0;
nc_addvar(file,v1);
nc_varput(file,v1.Name,theta_b);
clear v1

v1.Name = 'Vtransform';
v1.Datatype = 'int';
v1.Attribute(1).Name = 'long_name';
v1.Attribute(1).Value = 'vertical terrain-following transformation equation';
v1.Dimension = {};
v1.Chunking = [];
v1.Deflate = 0;
nc_addvar(file,v1);
nc_varput(file,v1.Name,Vtransform);
clear v1

v1.Name = 'Vstretching';
v1.Datatype = 'double';
v1.Attribute(1).Name = 'long_name';
v1.Attribute(1).Value = 'Vertical terrain-following stretching function';
v1.Dimension = {};
v1.Chunking = [];
v1.Deflate = 0;
nc_addvar(file,v1);
nc_varput(file,v1.Name,Vstretching);
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
nc_addvar(file,v1);
nc_varput(file,v1.Name,Tcline);
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
nc_addvar(file,v1);
nc_varput(file,v1.Name,hc);
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
nc_addvar(file,v1);
nc_varput(file,v1.Name,s_rho);
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
nc_addvar(file,v1);
nc_varput(file,v1.Name,s_w);
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
nc_addvar(file,v1);
nc_varput(file,v1.Name,Cs_r);
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
nc_addvar(file,v1);
nc_varput(file,v1.Name,Cs_w);
clear v1

end

