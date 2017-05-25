function [filelist] = romsinitialize_sm(type, run)
%ROMSINITIALIZE_SM Sets path for a given ROMS run and reads in the list of
%netcdf files of a given type.
%   INPUT:  type - type of ROMS output netcdf file to use - 'avg', 'his',
%                  or 'dia'
%           run - name of folder for a given ROMS run, directory is already
%                 specified in function; e.g. '016'
%   OUTPUT: filelist - list of files of the given type for the given run

%path = '/net/storage3/mack/ROMSprojects/';
%project = 'RossSea/';

%in_dir = [path project run '/'];
in_dir = ['/media/mnemoniko/EarlGrey/CCPO/Output/' run '/'];
cmd = ['ls ' in_dir 'ross_' type '*.nc > ' in_dir type 'Files.txt']; % Create a list of files
system( cmd ); clear cmd;
 
filelist = filelistread_sm(in_dir, [type 'Files.txt']);

end

