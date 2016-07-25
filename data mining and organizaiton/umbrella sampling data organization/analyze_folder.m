function data_list = analyze_folder(folder_name, manual_mode)
%ANALYZE_FOLDER returns data list of relevant 
%               experimental conditions (PEG and KCL pairs);
%   Data in the folder (folder_name) should be organized into subfolders
%   with naming convention as follows:
%   '(PEG_Percentage) pct PEG (KCL_concnetration)mM KCL (optional flags)'
%   where (PEG_Percentage) and (KCL_concentration) are numbers, 
%   (optional flags) are anything else in the folder name after KCL
%   currently supported flags are:
%       3-bundle - indicates that this data is for 3-MT bundle


data_list = {};
data_count = 0;

all_files = dir(folder_name);
all_dirs = {};
dir_count = 0;



%select only directories from list
for i = 1:length(all_files)
    file = all_files(i);
    if (file.isdir) && ~strcmp(file.name,'.') && ~strcmp(file.name,'..') 
        dir_count = dir_count+1;
        all_dirs{dir_count} = all_files(i);
    end;
end;
      
%now we got list of directiries in the folder
%to determine which ones are relevant, preforms scanning of all 
%sscanf cn return 2x1, 3x1 or empty vector, where the first number is
%PEG and the second is KCL, optional third corresponds to 3-bundle;

for i = 1 : dir_count
    name = all_dirs{i}.name;
    result = sscanf(name,'%f %*s %*s %d %*s %*s %d');
    
    if ( (length(result) >= 2) && (~manual_mode))
        data_count = data_count + 1;
        data_list{data_count}.PEG = result(1);
        data_list{data_count}.KCL = result(2);
        data_list{data_count}.Three_Bundle = 0;
        if (length(result) == 3)
            data_list{data_count}.Three_Bundle = 1;
        end;
        data_list{data_count}.Description = name;
    else
        data_count = data_count + 1;
        data_list{data_count}.PEG = 0;
        data_list{data_count}.KCL = 0;
        data_list{data_count}.Three_Bundle = 0;
        data_list{data_count}.Description = name;
    end;
              

end;    


end

