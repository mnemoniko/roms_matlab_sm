function [rx1]=rx1factor(z_w,rmask)
 
%
% RFACTOR:  Compute bathymetry rx1-factor
%
% [r]=rfactor(z_w,rmask)
%
% This function computes the bathymetry stiffness ratio, r-factor
%
% On Input:
%
%    z_w         depth levels at W-points (3D array)
%    rmask       Land/Sea masking at RHO-points (2D array)
%
% On Output:
%
%    rx1         R-factor (2D array)
%
%z_w = z_w; %Needs to be done in positive values

Lp = size(z_w,1); Mp = size(z_w,2); N = size(z_w,3)-1;
L=Lp-1;
M=Mp-1;

%  Land/Sea mask on U-points.

for j=1:Mp
  for i=2:Lp
    umask(i-1,j)=rmask(i,j)*rmask(i-1,j);
  end
end

%  Land/Sea mask on V-points.

for j=2:Mp
  for i=1:Lp
    vmask(i,j-1)=rmask(i,j)*rmask(i,j-1);
  end
end

umask(umask==0)=NaN;
vmask(vmask==0)=NaN; %NaNs at land points

umask=repmat(umask,1,1,N);
vmask=repmat(vmask,1,1,N);

%----------------------------------------------------------------------------
%  Compute R-factor.
%----------------------------------------------------------------------------


hx1(1:L,1:Mp,1:N)=abs(z_w(2:Lp,1:Mp,2:N+1)-z_w(1:L,1:Mp,2:N+1)+z_w(2:Lp,1:Mp,1:N)-z_w(1:L,1:Mp,1:N));
hx2(1:L,1:Mp,1:N)=z_w(2:Lp,1:Mp,2:N+1)+z_w(1:L,1:Mp,2:N+1)-z_w(2:Lp,1:Mp,1:N)-z_w(1:L,1:Mp,1:N);
hx = hx1./hx2;
    
hy1(1:Lp,1:M,1:N)=abs(z_w(1:Lp,2:Mp,2:N+1)-z_w(1:Lp,1:M,2:N+1)+z_w(1:Lp,2:Mp,1:N)-z_w(1:Lp,1:M,1:N));
hy2(1:Lp,1:M,1:N)=z_w(1:Lp,2:Mp,2:N+1)+z_w(1:Lp,1:M,2:N+1)-z_w(1:Lp,2:Mp,1:N)-z_w(1:Lp,1:M,1:N);
hy = hy1./hy2;

hx=hx.*umask;
hy=hy.*vmask;


rx1(1:L,1:M,1:N)=max(max(hx(1:L,1:M,1:N),hx(1:L,2:Mp,1:N)), ...
               max(hy(1:L,1:M,1:N),hy(2:Lp,1:M,1:N)));

rmin=min(rx1(:));
rmax=max(rx1(:));
ravg=nanmean(rx1(:));
rmed=nanmedian(rx1(:));

disp('  ');
disp(['Minimum r-value = ', num2str(rmin)]);
disp(['Maximum r-value = ', num2str(rmax)]);
disp(['Mean    r-value = ', num2str(ravg)]);
disp(['Median  r-value = ', num2str(rmed)]);

return