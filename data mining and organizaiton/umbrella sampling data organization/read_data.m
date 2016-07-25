function [ all_data ] = read_data(main_dir)
%Finction reads umbrella sampling data from subfolders in main folder
%

dirs = dir(main_dir);
data_dir = {};

%counts data folders
dir_count = 0;

for i = 1:length(dirs)
    name = dirs(i).name;
    if ~(strcmp(name,'.') || strcmp(name, '..'))
        dir_count = dir_count + 1;
        data_dir{dir_count} = name;
    end
end;

%%%%%%%%%%%%%%%%%%%%%%%%
%create data counter
data_counter = 0;
%create big cell for all data
all_data = {};

for i = 1 : dir_count
    
    data_path = [main_dir '\' data_dir{i}];
    
    %see what is in the directory
    dir_list = dir(data_path);
    file_list = {};
    count = 0;
    
    %only consider data files from given directory
    %create file list for this purpose
    for i = 1:length(dir_list)
        test_dir = dir_list(i).isdir;
        
        if ( test_dir ~= 1)% if the given object is not a directory
            test_ext = dir_list(i).name(end-2:end);
            if ( strcmp(test_ext, 'mat'))
                count = count + 1;
                file_list{count} = dir_list(i).name;
            end;
        end;
    end;
    
    
    %start loading data
    %use files from data list
    for i = 1:length(file_list)
        data_counter = data_counter + 1;
        fname = [data_path '\' file_list{i}];
        all_data{data_counter} = load(fname);
        
        %older version of bead tracker did not ask for overlap length input
        %once data is loaded from file fname, check if variable with name
        %Overlap_length is present
        %if not, add it with value '-1';
        %this will mean 'No overal length recorded'
        if (~isfield(all_data{data_counter},'Overlap_length'))
            all_data{data_counter}.Overlap_length = '-1';
        end;
        
        %old versions of bead tracking software recorded id in different format
        %get proper id to ensure cross-compatibility
        id = all_data{data_counter}.ID;
        if (~strcmp(id(end),'_'))
            id = [id '_'];
            all_data{data_counter}.ID = id;
        end;
        
        
    end;
    
end

end

