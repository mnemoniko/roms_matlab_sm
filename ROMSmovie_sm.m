function [ ] = ROMSmovie_sm( run, variable, startDay, endDay, slice, sliceCoord, sliceCoord2, varargin )
%ROMSmovie_sm Using data from ROMS simulation in folder "run", this
%function creates a .gif format animation.
% INPUT:    run - Simulation number as a string
%           variable - variable name as a string
%           startDay - Day of simulation to start (165 gives last year for
%           year and a half PRISM run)
%           endDay   - Day of simulation to end on (530 for 1.5 year Ross
%           Sea runs)
%           slice - direction for cross section: 1 = EW slice; 2 = N/S
%           slice; 3 = horizontal slice
%           sliceCoord - model coordinate for slice. If slice=3, this is
%           the model level - set to 0 to use depth instead
%           sliceCoord2 - limits to second model coordinate for slice =1&2,
%           e.g., [10 200]; model depth for slice 3
%           cmin, cmax - Option to set your own values for colorbar max &
%           min
% OUTPUT:   None.
%

% Initialize file list
file = romsinitialize_sm('avg',run);
nfile = size(file,1); steps = zeros(length(nfile),1);

%Get time variable

time = getROMSvar_sm(run,'avg','ocean_time','m');
zice = nc_varget(file(1,:),'zice');
h = nc_varget(file(1,:),'h');

switch(slice)
    case 1 %E/W
        Vtransform = nc_varget(file(1,:),'Vtransform');
        Vstretching = nc_varget(file(1,:),'Vstretching');
        theta_s = nc_varget(file(1,:),'theta_s');
        theta_b = nc_varget(file(1,:),'theta_b');
        hc = nc_varget(file(1,:),'hc'); %critical depth
        s_rho = nc_varget(file(1,:),'s_rho');
        N = length(s_rho);
        igrid = 1;
        zeta = getROMSvar_sm(run,'avg','zeta');
        zeta = squeeze(zeta(startDay:endDay,sliceCoord,sliceCoord2(1):sliceCoord2(2)));
        h = squeeze(h(sliceCoord,sliceCoord2(1):sliceCoord2(2)));
        zice = squeeze(zice(sliceCoord,sliceCoord2(1):sliceCoord2(2)));
        
        steps = size(zeta,1);
        X = size(zeta,2);
        z = zeros(steps,N,X);
        for i=1:steps
            temp2 = set_depth(Vtransform, Vstretching, theta_s, theta_b, hc, N, igrid, h, squeeze(zeta(i,:))+zice);
            z(i,:,:)=squeeze(shiftdim(temp2,2));
            display(['Finished depth i = ' num2str(i)]);
        end
        
        data = [];
        for i=1:nfile
            temp = nc_varget(file(i,:),variable);
            temp = squeeze(temp(:,:,sliceCoord,sliceCoord2(1):sliceCoord2(2)));
            data = cat(1,data,temp);
        end
        data = data(startDay:endDay,:,:); %time by depth by EW
        Y = repmat(1:X,[N,1]);
        
    case 2 %N/S
        Vtransform = nc_varget(file(1,:),'Vtransform');
        Vstretching = nc_varget(file(1,:),'Vstretching');
        theta_s = nc_varget(file(1,:),'theta_s');
        theta_b = nc_varget(file(1,:),'theta_b');
        hc = nc_varget(file(1,:),'hc'); %critical depth
        s_rho = nc_varget(file(1,:),'s_rho');
        N = length(s_rho);
        igrid = 1;
        zeta = getROMSvar_sm(run,'avg','zeta');
        zeta = squeeze(zeta(startDay:endDay,sliceCoord2(1):sliceCoord2(2),sliceCoord));
        h = squeeze(h(sliceCoord2(1):sliceCoord2(2),sliceCoord));
        zice = squeeze(zice(sliceCoord2(1):sliceCoord2(2),sliceCoord));
        
        steps = size(zeta,1);
        X = size(zeta,2);
        z = zeros(steps,N,X);
        for i=1:steps
            temp2 = set_depth(Vtransform, Vstretching, theta_s, theta_b, hc, N, igrid, h, squeeze(zeta(i,:))'+zice);
            z(i,:,:)=squeeze(shiftdim(temp2,2));
            display(['Finished depth i = ' num2str(i)]);
        end
        
        data = [];
        for i=1:nfile
            temp = nc_varget(file(i,:),variable);
            temp = squeeze(temp(:,:,sliceCoord2(1):sliceCoord2(2),sliceCoord));
            data = cat(1,data,temp);
        end
        data = data(startDay:endDay,:,:); %time by depth by EW
        Y = repmat(1:X,[N,1]);
        
    case 3 %horizontal
        if sliceCoord == 0
            temp = nc_varget(file(1,:),variable);
            data = zeros(length(time,size(temp,3),size(temp,4)));
            clear temp

            for n=1:nfile
                temp = nc_getdiminfo(file(n,:),'ocean_time');
                steps(n)=temp.Length;
            end

            for i=1:nfile
                for n=1:steps(i)
                    totalsteps = n+sum(steps(1:i))-steps(i);
                    [data(totalsteps,:,:),~,~] = roms_zslice(file(i,:),variable,n,sliceCoord2);
                    %This may be wrong - roms_zslice account for ice shelf?
                end
            end
        else
            data = [];

            for i=1:nfile
                data = cat(1,data,nc_varget(file(i,:),variable));
            end
            data = squeeze(data(:,sliceCoord,:,:));
        end
end

data = data(startDay:endDay,:,:);
time = time(startDay:endDay);

if(~isempty(varargin))
    cmin = varargin{1};
    cmax = varargin{2};
else
    cmin = min(min(min(data)));
    cmax = max(max(max(data)));
end

figure;
for i=1:1:size(data,1)
    colormap(brewermap([],'Blues'));
    if(slice==3)
        pcolor(squeeze(data(i,:,:)))
    else
        pcolor(Y,squeeze(z(i,:,:)),squeeze(data(i,:,:)));
    end
    colorbar;
    shading flat;
    if(slice==3)
        daspect([1 1 1]);
        hold on;
        contour(zice,[-10,-10],'Color',[.5 .5 .5]);
        contour(h,[100:200:1000 2000 3000],'Color','k');
    end
    caxis([cmin cmax]);
    title(['Date: ' datestr(time(i),'mmm dd, yyyy')]);
    hold off;
    if i<10
        print('-dpng', ['sm_projects/img00' num2str(i)]);
    elseif i<100
        print('-dpng', ['sm_projects/img0' num2str(i)]);
    else
        print('-dpng', ['sm_projects/img' num2str(i)]);
    end

end

% cmd = 'mogrify -format gif img*.png';
% system( cmd );
% if depthLevel ==0
%     cmd2 = ['convert -delay 30 -loop 1 img*.gif ' run variable depth 'm.gif'];
% else
%     cmd2 = ['convert -delay 30 -loop 1 img*.gif ' run variable 'Layer' depthLevel '.gif'];
% end
% system(cmd2);
% 
% cleanup = 'rm img*.gif';
% system(cleanup);
% cleanup = 'rm img*.png';
% system(cleanup);

display('Use these commands: mogrify -format gif img*.png');
display('convert -delay 30 -loop 1 img*.gif outfile.gif');
display('Then rm img*.gif and rm img*.png');


end

