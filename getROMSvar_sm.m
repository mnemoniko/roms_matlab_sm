function [ var ] = getROMSvar_sm( in_dir,type,varname,varargin )
%GETROMSVAR_SM Reads in a specific ROMS variable from a given simulation
%run and returns the array over the whole simulation time
%   INPUT:  in_dir - directory of ROMS output files
%           type - his, avg, or dia files (also a str) (might work with sta
%           too, untested)
%           varname - name of the variable to get (another string)
%           timeformat - option format for ocean_time variable, if requested. Use
%           'h' for hours; 'm' for matlab time; 's' for seconds; 'd' for
%           days
%   OUTPUT: var - Array of the appropriate variable (double)

files = romsinitialize_sm(in_dir,type);
nfiles = size(files,1);
var = [];
for i=1:nfiles
    temp = nc_varget(files(i,:),varname);
    if(i==1)
        dim = ndims(temp);
    end
    if(ndims(temp)<dim)
        temp = shiftdim(temp,-1);
    end
    var = cat(1,var,temp);
end

if(~isempty(varargin))
    format = varargin{1};
else
    format = 'none';
end

if(strcmp(varname,'ocean_time'))
    var = var./3600;
    var = var./24;
    var = var-2190; %days since model start
    switch format
        case 'h'
            var = var.*24;
        case 'm'
            var = datenum(2010,9,15+var,0,0,0);
            %var = datenum(0,0,var,0,0,0);
        case 's'
            var = var.*24.*3600;
        case 'd'
        case 'none'
            display('Please specify a time format');
        otherwise
            display('That time format was not recognized');
    end
end


end

