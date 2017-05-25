function [z]=set_depth_sta_sm(Vtransform, Vstretching, ...
                       theta_s, theta_b, hc, N, ...
                       igrid, h, zeta, report)
%
% SET_DEPTH:  Compute ROMS grid depth from vertical stretched variables
%
% [z]=set_depth(Vtransform, Vstretching, theta_s, theta_b, hc, N, ...
%               igrid, h, zeta)
%
% Modified from original set_depth 2015-07-16 to accommodate 1-D inputs
% from station files. igrid must be 1, always.
%
% Given a batymetry (h), free-surface (zeta) and terrain-following
% parameters, this function computes the 3D depths for the requested
% C-grid location. If the free-surface is not provided, a zero value
% is assumed resulting in unperturb depths.  This function can be
% used when generating initial conditions or climatology data for
% an application. Check the following link for details:
%
%    https://www.myroms.org/wiki/index.php/Vertical_S-coordinate
%
% On Input:
%
%    Vtransform    Vertical transformation equation:
%
%                    Vtransform = 1,   original transformation
%
%                    z(x,y,s,t)=Zo(x,y,s)+zeta(x,y,t)*[1+Zo(x,y,s)/h(x,y)]
%
%                    Zo(x,y,s)=hc*s+[h(x,y)-hc]*C(s)
%
%                    Vtransform = 2,   new transformation
%
%                    z(x,y,s,t)=zeta(x,y,t)+[zeta(x,y,t)+h(x,y)]*Zo(x,y,s)
%
%                    Zo(x,y,s)=[hc*s(k)+h(x,y)*C(k)]/[hc+h(x,y)]
%
%    Vstretching   Vertical stretching function:
%                    Vstretching = 1,  original (Song and Haidvogel, 1994)
%                    Vstretching = 2,  A. Shchepetkin (UCLA-ROMS, 2005)
%                    Vstretching = 3,  R. Geyer BBL refinement
%                    Vstretching = 4,  A. Shchepetkin (UCLA-ROMS, 2010)
%
%    theta_s       S-coordinate surface control parameter (scalar)
%
%    theta_b       S-coordinate bottom control parameter (scalar)
%
%    hc            Width (m) of surface or bottom boundary layer in which
%                    higher vertical resolution is required during
%                    stretching (scalar)
%
%    N             Number of vertical levels (scalar)
%
%    igrid         Staggered grid C-type (integer):
%                    igrid=1  => any points
%
%    h             Bottom depth, single-valued at RHO-points (m, positive),
%                    h(1)
%
%    zeta          Free-surface, single-valued at RHO-points (m), OPTIONAL,
%                    z(1)
%
%    report        Flag to report detailed information (OPTIONAL):
%                    report = 0,       do not report
%                    report = 1,       report information
%
% On Output:
%
%    z             Depths (m, negative), 1D array
%

% svn $Id: set_depth.m 647 2013-01-22 23:40:00Z arango $
%=========================================================================%
%  Copyright (c) 2002-2013 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%

z=[];

%--------------------------------------------------------------------------
%  Set several parameters.
%--------------------------------------------------------------------------

if (nargin < 8),
  disp(' ');
  disp('*** Error:  SET_DEPTH - too few arguments.');
  disp(['                     number of supplied arguments: ',          ...
       num2str(nargin)]);
  disp('                     number of required arguments: 8');
  disp(' ');
  return
end

if (Vtransform < 1 || Vtransform > 2),
  disp(' ');
  disp(['*** Error:  SET_DEPTH - Illegal parameter Vtransform = '       ...
	    num2str(Vtransfrom)]);
  disp(' ');
  return
end

if (Vstretching < 1 || Vstretching > 4),
  disp(' ');
  disp(['*** Error:  SET_DEPTH - Illegal parameter Vstretching = '      ...
        num2str(Vstretching)]);
  disp(' ');
  return
end

if (hc > min(min(h)) && Vtransform == 1),
  disp(' ');
  disp(['*** Error:  SET_DEPTH - critical depth exceeds minimum'        ...
        ' bathymetry value.']);
  disp(['                        Vtransform = ', num2str(Vtransform)]);
  disp(['                        hc         = ', num2str(hc)]);
  disp(['                        hmax       = ', num2str(min(min(h)))]);
  disp(' ');
  return
end

if (nargin < 9),
  zeta=zeros(size(h));
end

if (nargin < 10),
  report=1;
end

%--------------------------------------------------------------------------
% Compute vertical stretching function, C(k):
%--------------------------------------------------------------------------

if (report),
  disp(' ');
  if (Vtransform == 1),
    disp(['Vtransform  = ',num2str(Vtransform), '   original ROMS']);
  elseif (Vtransform == 2),
    disp(['Vtransform  = ',num2str(Vtransform), '   ROMS-UCLA']);
  end

      disp(['   igrid    = ',num2str(igrid),                            ...
            '   at horizontal RHO-points']);
end

kgrid=0;

[s,C]=stretching(Vstretching, theta_s, theta_b, hc, N, kgrid, report);

%--------------------------------------------------------------------------
%  Average bathymetry and free-surface at requested C-grid type.
%--------------------------------------------------------------------------

hr=h;
zetar=zeta;

%--------------------------------------------------------------------------
% Compute depths (m) at requested C-grid location.
%--------------------------------------------------------------------------

if (Vtransform == 1),
      for k=1:N,
        z0=(s(k)-C(k))*hc + C(k).*hr;
        z(k)=z0 + zetar.*(1.0 + z0./hr);
      end
elseif (Vtransform == 2),
      for k=1:N,
        z0=(hc.*s(k)+C(k).*hr)./(hc+hr);
        z(k)=zetar+(zeta+hr).*z0;
      end
end

return
