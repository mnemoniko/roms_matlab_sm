function [filelist] = filelistread_sm(directory, filename)
%FILELISTREAD Opens a ASCII file and reads a list of file names
%   INPUT: directory - full location of the file to be read
%          filename  - name of file to be read
%   OUTPUT: filelist - char type list of files
%
%Intended use: Read a list of files in that directory (avg, his, or dia)
% from ROMS output.  List of file names has been previously saved

[fid, msg] = fopen([directory filename], 'r');
 tmp  = [];  strg = 0;
 while ( strg ~= -1 )
   strg = fgetl(fid);
   if ( strg ~= -1 )
     tmp = [tmp; char( strg )];
   end
 end
 clear strg;
 filelist = tmp; clear tmp;
 fclose(fid); clear fid msg;

end

