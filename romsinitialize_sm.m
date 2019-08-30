function [filelist] = romsinitialize_sm(in_dir, type)
%ROMSINITIALIZE_SM Sets path for a given ROMS run and reads in the list of
%netcdf files of a given type.
%   INPUT:  in_dir - directory of ROMS simulation output
%           type - type of ROMS output netcdf file to use - 'avg', 'his',
%                  or 'dia'
%   OUTPUT: filelist - list of files of the given type for the given run

%path = '/net/storage3/mack/ROMSprojects/';
%project = 'RossSea/';

%in_dir = [path project run '/'];
%in_dir = ['/media/mnemoniko/Oolong/Ross/Output/' run '/'];
%in_dir = ['/media/mnemoniko/EarlGrey/CCPO/Output/' run '/'];
cmd = ['ls ' in_dir '*_' type '*.nc > ' in_dir type 'Files.txt']; % Create a list of files
system( cmd ); clear cmd;
 
filelist = filelistread_sm(in_dir, [type 'Files.txt']);

end

