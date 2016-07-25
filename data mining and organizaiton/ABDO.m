%% Active Bundles Data Organization
clc;
clear all;
close all;

data_folder =...
    'D:\My Documents\G-Drive\Research\All data compiled (new)\Active MT bundles\extension velocity data.xlsx';
data_sheet = 'K356-SA';

%% Import the data
[~, ~, raw] = xlsread(data_folder, data_sheet);
raw = raw(3:end,1:8);
raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};
cellVectors = raw(:,[1,2,3,6]);
raw = raw(:,[4,5,7,8]);

% Replace non-numeric cells with NaN
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % Find non-numeric cells
raw(R) = {NaN}; % Replace non-numeric cells

% Create output variable
data1 = reshape([raw{:}],size(raw));

% Allocate imported array to column variable names
index = cellVectors(:,1);
date = cellVectors(:,2);
run = cellVectors(:,3);
speed_ums = data1(:,1);
speed_error_ums = data1(:,2);
bundle_type = cellVectors(:,4);
force_pn = data1(:,3);
k_concentration_nm = data1(:,4);

% Clear temporary variables
clearvars data1 raw cellVectors R;

%% Create a list of all concentrations and bundle types
list_concentrations = [];
list_bundle_types = {};
for i = 1:length(k_concentration_nm)
    c_tmp = k_concentration_nm(i);
    type_tmp = bundle_type{i};
    
    j = find(list_concentrations == c_tmp);
    if (~isnan(c_tmp) && isempty(j))
        list_concentrations(end+1) = c_tmp;
    end;
    
    [found, found_ind] = find_inside_cell(bundle_type, type_tmp);
    already_found = find_inside_cell(list_bundle_types, type_tmp);
    if (found && ~already_found)
        list_bundle_types{end+1} = type_tmp;
    end;
end;
clearvars('i','c_tpm','type_tmp','j','found','found_ind','already_found');
list_concentrations = sort(list_concentrations);

disp(['Available concentrations: ' num2str(list_concentrations)]);
disp(list_bundle_types);

%% make relevant plots
for j = [1,2]
    h_fig = figure;
    hold all;
    n_concentrations = length(list_concentrations);
    annotation_groups = [];
    for i = 1:n_concentrations
        annotation_groups(i) = hggroup;
    end;
    selected_bundle_type = list_bundle_types{j};
    legend_string = cell(n_concentrations,1);
    for i = 1:n_concentrations
        selected_concentration = list_concentrations(i);
        legend_string{i} = [num2str(selected_concentration) ' nM'];
        concentration_index = find(k_concentration_nm == selected_concentration);
        [type_found, type_index] = find_inside_cell(bundle_type, selected_bundle_type);
        [common_found, selected_index] = find_common_elements(type_index, concentration_index);
        
        n_points = length(selected_index);
        speed_points = speed_ums(selected_index);
        speed_errors = speed_error_ums(selected_index);
        
        % display concentration
        fprintf('Selected concentration: %i nM\n',selected_concentration);
        fprintf('Bundle type: %s\n',selected_bundle_type);
        fprintf('Measurements found: %i\n',length(selected_index));
        
        % remove NaNs from measurements
        nans_speed = find(isnan(speed_points));
        nans_errors = find(isnan(speed_errors));
        
        speed_no_nans = find(~isnan(speed_points));
        %errors_no_nans = find(~isnan(speed_errors));
        errors_no_nans = speed_no_nans;
        [common_found, relevant_points] = find_common_elements(speed_no_nans, errors_no_nans);
        
        s = speed_points(relevant_points);
        e = speed_errors(relevant_points);
        n_relevant = length(relevant_points);
        c = selected_concentration*ones(1,n_relevant);
        fprintf('After removal of NaN values: %i\n', n_relevant);
        
        
        %actually make a plot
        h_conc_plot = plot(c,s,'o');
        set(h_conc_plot,'Parent', annotation_groups(i));       
        conc_color = get(h_conc_plot, 'Color');
        h_mean_plot = errorbar(selected_concentration, mean(s), std(s),...
            'MarkerSize',30,...
            'Marker','.',...
            'LineStyle','none',...
            'LineWidth',2,...
            'Color',conc_color);
        set(h_mean_plot,'Parent', annotation_groups(i));
        set(get(get(annotation_groups(i),'Annotation'),'LegendInformation'),...
            'IconDisplayStyle','on');
        display('plot done... moving on to the next c'); 
    end;
    
    set(gca,'Xscale','log','YGrid','on','XGrid','on');
    legend(legend_string);
    
    title(['Speed: ' selected_bundle_type]);
    xlabel('Motor concentration (nM)');
    ylabel('Speed (um s^{-1})');
    %xlim([4 1200]);
    %ylim([0 0.500]);
    box(gca, 'on');
    prettify_plot(gcf);
    set(gcf,'Position',[100 100 900 400]);
end;
