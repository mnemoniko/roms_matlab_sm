function [ hnew, zicenew, rx1, xFix, yFix ] = rx1smoothing_sm( file, rx1max, ISswitch )
%rx1smoothing_sm Iteratively smooths a ROMS bathymetry (and ice shelf) to
%the requested rx1 value.  Outputs a new bathymetry, maximum rx1 map, and 
%list of points that still need smoothed.
%   file - ROMS grid file name
%   rx1max - Maximum value of rx1 to smooth to
%   ISswitch - How to handle IS front values: 1(ignore), 2(do not smooth),
%   3(smooth only ocean bottom).

maxIter = 500;
%deltaH = 1; %meters

%% Get appropriate parameters from ROMS file
h = nc_varget(file,'h');
zice = nc_varget(file,'zice');
Vtransform = nc_varget(file,'Vtransform');
Vstretching = nc_varget(file,'Vstretching');
theta_s = nc_varget(file,'theta_s');
theta_b = nc_varget(file,'theta_b');
mask = nc_varget(file,'mask_rho');
N = nc_getdiminfo(file,'N','Length');
Tcline = nc_varget(file,'Tcline');

ziceold = zice;

if(Vtransform==1)
    hc = min(h(:),Tcline);
else
    hc = Tcline;
end
hmin = min(h(:));
zmin = max(zice(zice<0));

X = size(h,1); Y = size(h,2);

%Find ice shelf points:
mask_ice = mask;
mask_ice(zice<0)=2; %2 at ice points, 1 at water

%Find IS Front:  Have to set this ahead of time, a 3 point check won't
%catch all front points while fixing them.

mask_ice(mask==0)=NaN;
mask_ISfront = zeros(size(mask));

for i=2:X-1
    for j=2:Y-1
        check1=nanmean(mask_ice(i-1:i+1,j));
        check2=nanmean(mask_ice(i,j-1:j+1));
        if(rem(check1,1)~=0 || rem(check2,1)~=0)
            mask_ISfront(i,j)=1;
        elseif(check1~=mask_ice(i,j) || check2~=mask_ice(i,j))
            mask_ISfront(i,j)=1;
        else
            mask_ISfront(i,j)=0;
        end
    end
end
%These IS front points cover both ice & water - need to reset to only cover
%ice points for this application
mask_ISfront(mask==0)=NaN;
mask_ISfront(mask_ice==1)=0;
%testmask = 3*mask_ISfront + mask_ice;
mask_ISfront = mask_ISfront(1:X-1,1:Y-1);

h(mask==0)=NaN;
zice(mask==0)=NaN; %Don't use land points in smoothing calc

%% Smoothing:

for q = 1:maxIter %maximum iterations
    %Calculate rx1
    z_w = set_depth(Vtransform,Vstretching,theta_s,theta_b,hc,N,5,h,zice,0);
    rx1 = rx1factor(z_w,mask);
    rx1 = max(rx1,[],3);
    
    %List of points to fix
    maskFix = zeros(size(rx1));
    maskFix(rx1>rx1max)=1;  %Modify points with high rx1
    if(ISswitch==2)
        maskFix(mask_ISfront==1)=0; %Don't modify IS front points
    end
    
    [xFix, yFix] = find(maskFix==1);
    
    disp(['Starting iteration ' num2str(q)]);
    disp(['Currently ' num2str(length(xFix)) ' points need smoothed.']);
    
    if(isempty(xFix))
        disp('There are no more points to smooth.');
        disp(['Total iterations: ' num2str(q-1)]);
        break
    end
    for p = 1:length(xFix) %For each point that needs fixed
        
        if(mask_ISfront(xFix(p),yFix(p))==1)
            isISfront=1; %Is this point on the ice shelf front?
        else
            isISfront=0;
        end
        
        hmean = nanmean([h(xFix(p),yFix(p)); h(xFix(p)+1,yFix(p)); h(xFix(p),yFix(p)+1); ...
            h(xFix(p)-1,yFix(p)); h(xFix(p),yFix(p)-1)]);
        zmean = nanmean([zice(xFix(p),yFix(p)); zice(xFix(p)+1,yFix(p)); zice(xFix(p),yFix(p)+1); ...
            zice(xFix(p)-1,yFix(p)); zice(xFix(p),yFix(p)-1)]);
        
        %If thin points, increase h, decrease zice.  If thick, adjust to
        %mean:
        thick = hmean+zmean;
        deltaH = .005.*thick; %.5% of thickness for change (?)
        
        if(thick>hmin*4)
        %Adjust that point and the two surrounding towards the mean:
            h(xFix(p),yFix(p)) = h(xFix(p),yFix(p)) + sign(hmean-h(xFix(p),yFix(p))).*deltaH;
            h(xFix(p)+1,yFix(p)) = h(xFix(p)+1,yFix(p)) + sign(hmean-h(xFix(p)+1,yFix(p))).*deltaH;
            h(xFix(p),yFix(p)+1) = h(xFix(p),yFix(p)+1) + sign(hmean-h(xFix(p),yFix(p)+1)).*deltaH;
        
            if(isISfront && ISswitch==3)
                %Do nothing
            else
            zice(xFix(p),yFix(p)) = zice(xFix(p),yFix(p)) + sign(zmean-zice(xFix(p),yFix(p))).*deltaH;
            zice(xFix(p)+1,yFix(p)) = zice(xFix(p)+1,yFix(p)) + sign(zmean-zice(xFix(p)+1,yFix(p))).*deltaH;
            zice(xFix(p),yFix(p)+1) = zice(xFix(p),yFix(p)+1) + sign(zmean-zice(xFix(p),yFix(p)+1)).*deltaH;
            end
        else
            h(xFix(p),yFix(p)) = h(xFix(p),yFix(p)) + deltaH;
            h(xFix(p)+1,yFix(p)) = h(xFix(p)+1,yFix(p)) + deltaH;
            h(xFix(p),yFix(p)+1) = h(xFix(p),yFix(p)+1) + deltaH;
            
            if(isISfront && ISswitch==3)
                %Do nothing
            else
            zice(xFix(p),yFix(p)) = zice(xFix(p),yFix(p)) + deltaH;
            zice(xFix(p)+1,yFix(p)) = zice(xFix(p)+1,yFix(p)) + deltaH;
            zice(xFix(p),yFix(p)+1) = zice(xFix(p),yFix(p)+1) + deltaH;
            end
        end
    end
    
    %Re-set minimum values:
    h(h<hmin) = hmin;
    zice(zice>zmin) = zmin;
    
    %Fix thin values:
    wc_thick = h+zice;
    wc_thick(mask==0)=hmin; %Arbitrarily set water column thickness to not modify land points
    thin = wc_thick<hmin;  %Points that are too thin 

    %Deepen thin points:
    while(~isempty(find(thin==1)))
        h(thin) = h(thin) + .5; %Deepen bathy by .5m at thin points
        zice(thin) = zice(thin) + .5; %Shallow ice by .5m at thin points
    
        %Reset h & zice mins:
        h(h<hmin) = hmin;
        zice(zice>zmin) = zmin;
    
        %Re-calculate wc thickness & thin points
        wc_thick = h+zice;
        wc_thick(mask==0)=hmin;
        thin = wc_thick<hmin;    
    end
    
    %Re-set land/open water values (just in case):
    h(mask==0)=hmin;
    zice(mask_ice==1)=0;
    
end

%%
hnew = h;
hnew(mask==0)=hmin;
zicenew = zice;
zicenew(mask_ice==1)=0.0; % Set open water to 0 zice
zicenew(mask==0)=ziceold(mask==0); %Set land points to old ice values

nc_varput(file,'h',hnew);
nc_varput(file,'zice',zicenew);

disp(['New h and zice saved in ' file]);
if(~isempty(xFix))
    disp('Please check near the following grid points and adjust by hand:');
    for i=1:length(xFix)
        if(mask(xFix(i),yFix(i))==1 && mask_ice(xFix(i),yFix(i))+mask_ice(xFix(i)+1,yFix(i)+1) ~= 3)
            disp([num2str(xFix(i)) ',' num2str(yFix(i))]);
        end
    end
    disp(['Re-save h and zice into ' file ' and re-run this function']);
end

end

