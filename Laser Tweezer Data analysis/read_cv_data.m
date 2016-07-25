function [trap_x, trap_y, pos_x, pos_y, t, z] = read_cv_data( filename )
%[trap_x, trap_y, pos_x, pos_y, t] = READ_CV_DATA( filename )
warning('This function is depreciated! Only use with old AOD LabView Software')
%% Initialize variables.
delimiter = ',';

%% Format string for each line of text:
%   column1: double (%f) = Trap X, MHz
%	column2: double (%f) = Trap Y, MHz
%   column3: double (%f) = Bead_x, V or Nm (if conversion is used);
%	column4: double (%f) = Bead_y, V or Nm (if conversion is used);
%   column5: double (%f) = delta t, time between points, ms
%   column6: double (%f) = z-position, micron
% For more information, see the TEXTSCAN documentation.
formatSpec = '%f%f%f%f%f%f%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,  'ReturnOnError', false);

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Allocate imported array to column variable names
trap_x = dataArray{:, 1};
trap_y = dataArray{:, 2};
pos_x = dataArray{:, 3};
pos_y = dataArray{:, 4};
t = dataArray{:, 5};
z = dataArray{:, 6};

%% find all trailing zeroes due to pre-allocation
valid = find(t > 0);
trap_x = trap_x(valid);
trap_y = trap_y(valid);
pos_x = pos_x(valid);
pos_y = pos_y(valid);
t = t(valid);
z = z(valid);


