function [ mask_u, mask_v, mask_psi ] = uvpmask_sm( mask_rho )
%UVPMASK_SM Calculation of ROMS mask_u, mask_v, and mask_psi from a given
%mask_rho. This differs from the ROMS released version because using
%snctools to read in netcdf files switches the dimensions (I think).
%   

[Lp,Mp]=size(mask_rho);
L=Lp-1;
M=Mp-1;

%  Land/Sea mask on V-points.

mask_v(1:L,1:Mp)=mask_rho(2:Lp,1:Mp).*mask_rho(1:L,1:Mp);

%  Land/Sea mask on U-points.

mask_u(1:Lp,1:M)=mask_rho(1:Lp,2:Mp).*mask_rho(1:Lp,1:M);

%  Land/Sea mask on PSI-points.

mask_psi(1:L,1:M)=mask_rho(1:L,1:M ).*mask_rho(2:Lp,1:M ).* ...
               mask_rho(1:L,2:Mp).*mask_rho(2:Lp,2:Mp);


end

