function [ h, zice ] = rx1smoothingv2_sm( file, rx1max, ISswitch )
%rx1smoothingv2_sm Smoothing rx1 values using method in Martinho & Bateen
%2006
%   INPUT:
%       file - ROMS grid file
%       rx1max - maximum desired value of rx1max (Choose something high,
%       say 10000, to not smooth and get surface and bottom rx1 maximums of
%       the current grid)
%       ISswitch - Set to 1 to smooth the ice shelf front, 0 to not.
%       Recommend to only turn on AFTER the rest of the grid is smoothed.

%% Variables & Parameters
maxIter = 3000;
tol = 0.00001;

h = nc_varget(file,'h');
zice = nc_varget(file,'zice');
Vtransform = nc_varget(file,'Vtransform');
Vstretching = nc_varget(file,'Vstretching');
theta_s = nc_varget(file,'theta_s');
theta_b = nc_varget(file,'theta_b');
mask = nc_varget(file,'mask_rho');
N = nc_getdiminfo(file,'N','Length');
Tcline = nc_varget(file,'Tcline');

Nw = N+1;
mask_zice = mask;
mask_zice(zice<0)=2; %2 at ice points, 1 at water
mask_zice = mask_zice.*mask; %Don't modify land ice points
zice_old = zice; %Preserve land values

if(Vtransform==1)
    hc = min(min(h(:)),Tcline);
else
    hc = Tcline;
end
hmin = min(h(:));
zmin = max(zice(zice<0));

X = size(h,1); Y = size(h,2);


%% Calculate initial grid volume(s)

%% Modify Bathymetry & ice shelf draft

for q = 1:maxIter %maximum iterations
    
%%%Modify bottom bathymetry%%%

    %Calculate rx1 for lowest level:
    hw = set_depth(Vtransform,Vstretching,theta_s,theta_b,hc,N,5,h,zice,0);
    
    rx1L = zeros(size(mask)); rx1R = rx1L; rx1T = rx1L; rx1B = rx1L;
    
    rx1top = abs(hw(2:X-1,2:Y-1,2) - hw(1:X-2,2:Y-1,2) + hw(2:X-1,2:Y-1,1) - hw(1:X-2,2:Y-1,1));
    rx1bot = hw(2:X-1,2:Y-1,2) + hw(1:X-2,2:Y-1,2) - hw(2:X-1,2:Y-1,1) - hw(1:X-2,2:Y-1,1);
    rx1L(2:X-1,2:Y-1) = rx1top.*mask(2:X-1,2:Y-1).*mask(1:X-2,2:Y-1)./rx1bot;
    [xC, yC] = find(rx1L>rx1max+tol); %Center points with high left rx1
    x2 = xC-1; y2 = yC; %Left points that correspond
    pointPairs = [sub2ind(size(rx1L),xC,yC) sub2ind(size(rx1L),x2,y2)]; %List of point pairs as indicies
    
    rx1top = abs(hw(2:X-1,2:Y-1,2) - hw(3:X,2:Y-1,2) + hw(2:X-1,2:Y-1,1) - hw(3:X,2:Y-1,1));
    rx1bot = hw(2:X-1,2:Y-1,2) + hw(3:X,2:Y-1,2) - hw(2:X-1,2:Y-1,1) - hw(3:X,2:Y-1,1);
    rx1R(2:X-1,2:Y-1) = rx1top.*mask(2:X-1,2:Y-1).*mask(3:X,2:Y-1)./rx1bot;
    [xC, yC] = find(rx1R>rx1max+tol); %Center points with high right rx1
    x2 = xC+1; y2 = yC; %Right points that correspond
    pointPairs = cat(1,pointPairs,[sub2ind(size(rx1L),xC,yC) sub2ind(size(rx1L),x2,y2)]); %List of point pairs as indicies
    
    rx1top = abs(hw(2:X-1,2:Y-1,2) - hw(2:X-1,1:Y-2,2) + hw(2:X-1,2:Y-1,1) - hw(2:X-1,1:Y-2,1));
    rx1bot = hw(2:X-1,2:Y-1,2) + hw(2:X-1,1:Y-2,2) - hw(2:X-1,2:Y-1,1) - hw(2:X-1,1:Y-2,1);
    rx1B(2:X-1,2:Y-1) = rx1top.*mask(2:X-1,2:Y-1).*mask(2:X-1,1:Y-2)./rx1bot;
    [xC, yC] = find(rx1B>rx1max+tol); %Center points with high bottom rx1
    x2 = xC; y2 = yC-1; %Bottom points that correspond
    pointPairs = cat(1,pointPairs,[sub2ind(size(rx1L),xC,yC) sub2ind(size(rx1L),x2,y2)]); %List of point pairs as indicies
    
    rx1top = abs(hw(2:X-1,2:Y-1,2) - hw(2:X-1,3:Y,2) + hw(2:X-1,2:Y-1,1) - hw(2:X-1,3:Y,1));
    rx1bot = hw(2:X-1,2:Y-1,2) + hw(2:X-1,3:Y,2) - hw(2:X-1,2:Y-1,1) - hw(2:X-1,3:Y,1);
    rx1T(2:X-1,2:Y-1) = rx1top.*mask(2:X-1,2:Y-1).*mask(2:X-1,3:Y)./rx1bot;
    [xC, yC] = find(rx1T>rx1max+tol); %Center points with high top rx1
    x2 = xC; y2 = yC+1; %Top points that correspond
    pointPairs = cat(1,pointPairs,[sub2ind(size(rx1L),xC,yC) sub2ind(size(rx1L),x2,y2)]); %List of point pairs as indicies

    maxrx1Bottom = max(max(max(rx1L(:)),max(rx1R(:))),max(max(rx1B(:)),max(rx1T(:))));
    
    clear x2 xC y2 yC rx1B rx1L rx1T rx1R rx1top rx1bot
    
    %Re-order point Pairs, remove doubles:
    %pointPairs = sort(pointPairs,2,'ascend');  %This is not needed. Points
    %do not actually double up at all.
    
    disp(['Starting iteration ' num2str(q)]);
    disp(['Currently ' num2str(size(pointPairs,1)) ' points have high bottom rx1.']);
    disp(['Highest rx1 bottom value: ' num2str(maxrx1Bottom)]);

    for i=1:size(pointPairs,1)
        [x1, y1] = ind2sub(size(h),pointPairs(i,1)); %e
        [x2, y2] = ind2sub(size(h),pointPairs(i,2)); %e'
        
        if(h(x1,y1)>h(x2,y2))
            h1 = abs(((rx1max+1)/(rx1max-1))*(hw(x2,y2,2) - hw(x1,y1,1)) + hw(x1,y1,2));
            h2 = abs(((rx1max-1)/(rx1max+1))*(hw(x2,y2,2) - hw(x1,y1,1)) + hw(x1,y1,2));
            h(x2,y2) = max(h1,h2);
        else
            h1 = abs(((rx1max+1)/(rx1max-1))*(hw(x1,y1,2) - hw(x2,y2,1)) + hw(x2,y2,2));
            h2 = abs(((rx1max-1)/(rx1max+1))*(hw(x1,y1,2) - hw(x2,y2,1)) + hw(x2,y2,2));
            h(x1,y1) = max(h1,h2);
        end
      
    end
    
    
    %%%Modify Ice Shelf Bathymetry%%%

    
    if(ISswitch==0)
        mask_ice = zeros(size(mask));
        mask_ice(zice<0)=1;
        mask_ice = mask_ice.*mask; %only water & ice points
    else
        mask_ice=mask; %all water points
    end
    
    rx1L = zeros(size(mask)); rx1R = rx1L; rx1T = rx1L; rx1B = rx1L;
    
    rx1top = abs(hw(2:X-1,2:Y-1,Nw) - hw(1:X-2,2:Y-1,Nw) + hw(2:X-1,2:Y-1,Nw-1) - hw(1:X-2,2:Y-1,Nw-1));
    rx1bot = hw(2:X-1,2:Y-1,Nw) + hw(1:X-2,2:Y-1,Nw) - hw(2:X-1,2:Y-1,Nw-1) - hw(1:X-2,2:Y-1,Nw-1);
    rx1L(2:X-1,2:Y-1) = rx1top.*mask_ice(2:X-1,2:Y-1).*mask_ice(1:X-2,2:Y-1)./rx1bot;
    [xC, yC] = find(rx1L>rx1max+tol); %Center points with high left rx1
    x2 = xC-1; y2 = yC; %Left points that correspond
    pointPairs2 = [sub2ind(size(rx1L),xC,yC) sub2ind(size(rx1L),x2,y2)]; %List of point pairs as indicies
    
    rx1top = abs(hw(2:X-1,2:Y-1,Nw) - hw(3:X,2:Y-1,Nw) + hw(2:X-1,2:Y-1,Nw-1) - hw(3:X,2:Y-1,Nw-1));
    rx1bot = hw(2:X-1,2:Y-1,Nw) + hw(3:X,2:Y-1,Nw) - hw(2:X-1,2:Y-1,Nw-1) - hw(3:X,2:Y-1,Nw-1);
    rx1R(2:X-1,2:Y-1) = rx1top.*mask_ice(2:X-1,2:Y-1).*mask_ice(3:X,2:Y-1)./rx1bot;
    [xC, yC] = find(rx1R>rx1max+tol); %Center points with high right rx1
    x2 = xC+1; y2 = yC; %Right points that correspond
    pointPairs2 = cat(1,pointPairs2,[sub2ind(size(rx1L),xC,yC) sub2ind(size(rx1L),x2,y2)]); %List of point pairs as indicies
    
    rx1top = abs(hw(2:X-1,2:Y-1,Nw) - hw(2:X-1,1:Y-2,Nw) + hw(2:X-1,2:Y-1,Nw-1) - hw(2:X-1,1:Y-2,Nw-1));
    rx1bot = hw(2:X-1,2:Y-1,Nw) + hw(2:X-1,1:Y-2,Nw) - hw(2:X-1,2:Y-1,Nw-1) - hw(2:X-1,1:Y-2,Nw-1);
    rx1B(2:X-1,2:Y-1) = rx1top.*mask_ice(2:X-1,2:Y-1).*mask_ice(2:X-1,1:Y-2)./rx1bot;
    [xC, yC] = find(rx1B>rx1max+tol); %Center points with high bottom rx1
    x2 = xC; y2 = yC-1; %Bottom points that correspond
    pointPairs2 = cat(1,pointPairs2,[sub2ind(size(rx1L),xC,yC) sub2ind(size(rx1L),x2,y2)]); %List of point pairs as indicies
    
    rx1top = abs(hw(2:X-1,2:Y-1,Nw) - hw(2:X-1,3:Y,Nw) + hw(2:X-1,2:Y-1,Nw-1) - hw(2:X-1,3:Y,Nw-1));
    rx1bot = hw(2:X-1,2:Y-1,Nw) + hw(2:X-1,3:Y,Nw) - hw(2:X-1,2:Y-1,Nw-1) - hw(2:X-1,3:Y,Nw-1);
    rx1T(2:X-1,2:Y-1) = rx1top.*mask_ice(2:X-1,2:Y-1).*mask_ice(2:X-1,3:Y)./rx1bot;
    [xC, yC] = find(rx1T>rx1max+tol); %Center points with high top rx1
    x2 = xC; y2 = yC+1; %Top points that correspond
    pointPairs2 = cat(1,pointPairs2,[sub2ind(size(rx1L),xC,yC) sub2ind(size(rx1L),x2,y2)]); %List of point pairs as indicies

    maxrx1Top = max(max(max(rx1L(:)),max(rx1R(:))),max(max(rx1B(:)),max(rx1T(:))));
    
    clear x2 xC y2 yC rx1B rx1L rx1T rx1R rx1top rx1bot
    
    disp(['Currently ' num2str(size(pointPairs2,1)) ' points have high surface rx1.']);
    disp(['Highest rx1 surface value: ' num2str(maxrx1Top)]);
    disp('  ');
    
    if(isempty(pointPairs2) && isempty(pointPairs))
        disp('There are no more points to smooth.');
        disp(['Total iterations: ' num2str(q-1)]);
        break
    end

    for i=1:size(pointPairs2,1)
        [x1, y1] = ind2sub(size(h),pointPairs2(i,1)); %e
        [x2, y2] = ind2sub(size(h),pointPairs2(i,2)); %e'
        
        if(zice(x2,y2) > zice(x1,y1)) %Change zice(x1,y1)
            %test1(i)=zice(x1,y1);
            z1 = ((rx1max-1)/(rx1max+1))*(hw(x1,y1,Nw-1)-hw(x2,y2,Nw)) + hw(x2,y2,Nw-1);
            z2 = ((rx1max+1)/(rx1max-1))*(hw(x1,y1,Nw-1)-hw(x2,y2,Nw)) + hw(x2,y2,Nw-1);
            zice(x1,y1) = max(z1,z2);
            %test2(i)=zice(x1,y1);
        else %Change zice(x,y2)
            %test1(i)=zice(x2,y2);
            z1 = ((rx1max-1)/(rx1max+1))*(hw(x2,y2,Nw-1)-hw(x1,y1,Nw)) + hw(x1,y1,Nw-1);
            z2 = ((rx1max+1)/(rx1max-1))*(hw(x2,y2,Nw-1)-hw(x1,y1,Nw)) + hw(x1,y1,Nw-1);
            zice(x2,y2) = max(z1,z2);
            %test2(i)=zice(x2,y2);
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
    
    %Re-set open water  & land values:
    zice(mask_zice==1)=0; %Open water has no ice shelf
    zice(mask_zice==0) = zice_old(mask_zice==0); %Restore land points
    
end

%% Re-calculate grid volume

%% Write parameters to file

 %nc_varput(file,'h',h);
 %nc_varput(file,'zice',zice);


end

