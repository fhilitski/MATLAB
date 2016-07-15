function [Fsampling,Nsamples,Velocity,MaxMove,N_traps,T_delay,Converted,Laser_Power,AODx,AODy,QPDx,QPDy] = read_acq_data(filename)
%READ_ACQ_DATA Import numeric data from a text file as column vectors.
%   [FSAMPLING,NSAMPLES,VELOCITY,MAXMOVE,N_TRAPS,T_DELAY,CONVERTED,LASER_POWER,] =
%   READ_ACQ_DATA(FILENAME) Reads data from text file FILENAME for the default
%   selection.
%   [FSAMPLING,NSAMPLES,VELOCITY,MAXMOVE,N_TRAPS,T_DELAY,CONVERTED] =
%   READ_ACQ_DATA(FILENAME, STARTROW, ENDROW) Reads data from rows STARTROW
%   through ENDROW of text file FILENAME.
%
% Example:
%   [Fsampling,Nsamples,Velocity,MaxMove,N_traps,T_delay,Converted] =
%   read_acq_data('acquisition.csv',2, 2);
%
%    See also TEXTSCAN.

% Auto-generated by MATLAB on 2014/04/23 14:32:07

%% Initialize variables.
    startRow = 2;
    endRow = 2;
    delimiter = ',';


%% Format string for each line of text:
%   columns 1-12: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
textscan(fileID, '%[^\n\r]', startRow(1)-1, 'ReturnOnError', false);
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'ReturnOnError', false);
for block=2:length(startRow)
    frewind(fileID);
    textscan(fileID, '%[^\n\r]', startRow(block)-1, 'ReturnOnError', false);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'ReturnOnError', false);
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Allocate imported array to column variable names
Fsampling = dataArray{:, 1};
Nsamples = dataArray{:, 2};
Velocity = dataArray{:, 3};
MaxMove = dataArray{:, 4};
N_traps = dataArray{:, 5};
T_delay = dataArray{:, 6};
Converted = dataArray{:, 7};
Laser_Power = dataArray{:, 8};
AODx = dataArray{:, 9};
AODy = dataArray{:, 10};
QPDx = dataArray{:, 11};
QPDy = dataArray{:, 12};