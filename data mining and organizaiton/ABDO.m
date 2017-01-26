%% Active Bundles Data Organization
clc;
clear all;
close all;

data_folder =...
    'D:\My Documents\G-Drive\Research\All data\Active MT bundles\extension velocity data.xlsx';
data_sheet = 'K356-SA';
%data_sheet = 'KSA all'; %use this for k401 data

%% Import the data
[~, ~, raw_xls] = xlsread(data_folder, data_sheet);
raw = raw_xls(3:end,1:11);
raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};

%Columns are as follows:
%   1 - index
%   2 - date
%   3 - run name
%   4 - Velocity (um/s)  *
%   5 - Velocity 95% confidence interval from the fit (um/s) *
%   6 - type of bundle **
%   7 - Force max (pN) *
%   8 - Initial total length (bead-to-bead) um
%   9 - Bulk cluster concentration (nM) * (formerly #8)
%   10 - Comments (formerly #9)
%   11 - the rest has no information content 

%Save non-quantitative information in cellVectors
cellVectors = raw(:,[1,2,3,6]);
%Save velocity, conf. ints, force, initial length and concentration in raw
raw = raw(:,[4,5,7,8,9]);

% Replace non-numeric cells with NaN
raw(cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw)) = {NaN};

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
length_um = data1(:,4);
k_concentration_nm = data1(:,5);

% Clear temporary variables
%clearvars data1 raw cellVectors;

%% Create a list of all concentrations and bundle types
list_concentrations = [];
list_bundle_types = {};

[list_concentrations, ict, count_concentrations] = ...
    create_value_list(k_concentration_nm);

[list_bundle_types, ibt, count_bundle_types] = ...
    create_value_list(bundle_type);

%clearvars('i','c_tpm','type_tmp','j','found','found_ind','already_found');
%list_concentrations = sort(list_concentrations);

disp('Available concentrations: ');
disp(list_concentrations);
disp('Available bundle types: '); 
 for i = 1:length(list_bundle_types)
     disp(['   ' list_bundle_types{i} ': ' num2str(count_bundle_types(i))]);
 end;

%% make plots of speed vs concentration for all bundle types
color_index = 2;
%index j switches bundle type from list_bundle_types variable
for j = [1 2]
    h_fig = figure;
    color_order = get(gca,'ColorOrder');
    hold all;
    n_concentrations = length(list_concentrations);
    annotation_groups = [];
    for i = 1:n_concentrations
        annotation_groups(i) = hggroup;
    end;
    selected_bundle_type = list_bundle_types{j};
    legend_string = cell(n_concentrations,1);
    all_speeds = {};
    all_conc = [];
    all_means = [];
    all_errors = [];
    
    for i = 1:n_concentrations     
        selected_concentration = list_concentrations{i};
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
        %individual datapoints
        figure(h_fig);
        h_conc_plot = plot(c,s,'o','Color',color_order(color_index,:));
        %add plto to 
        set(h_conc_plot,'Parent', annotation_groups(i)); 
        set(get(get(annotation_groups(i),'Annotation'),'LegendInformation'),...
             'IconDisplayStyle','on');

        
        conc_color = color_order(color_index,:);
        %save speeds, concentrations, means and errors
        all_speeds{i} = s;
        all_conc(i) = selected_concentration;
        all_means(i) = mean(s);
        all_errors(i) = std(s);
        
        display('plot done... moving on to the next c');           
        fprintf('\n');
    end;
    
    %errorbar plot for the individual datapoints
    h_mean_plot = errorbar(all_conc,all_means,all_errors,...
            'MarkerSize',30,...
            'Marker','.',...
            'LineStyle','none',...
            'LineWidth',2,...
            'Color',conc_color);
    %set(h_mean_plot,'Parent', annotation_groups(i));
    
    set(gca,'Xscale','log','Yscale','lin','YGrid','on','XGrid','on');
    legend(legend_string);    
    title(['Speed: ' selected_bundle_type]);
    xlabel('Motor concentration (nM)');
    ylabel('Speed (um s^{-1})');
    box(gca, 'on');
    prettify_plot(gcf);
    xlim([4 1200]);
    ylim([0 0.50]);
    set(gcf,'Position',[100 100 900 400]);              
end;

%% plot force v concentration 
count = 0;
color_index = 2;
%force is only measured for buckling bundles
disp('Plotting force v concentration');
selected_bundle_type = 'Buckling';
n_concentrations = length(list_concentrations);
%create new figures and save handles
annotation_groups = [];
%h_fig is for force v concentration plot
h_fig = figure;
%save color order
color_order = get(gca,'ColorOrder');
hold all;
for i = 1:n_concentrations
    annotation_groups(i) = hggroup;
end;

legend_string = {};
all_conc = [];
all_means = [];
all_forces = {};
all_errors = [];

for i = 1:n_concentrations
    %analyze concentration 1-after-1
    selected_concentration = list_concentrations{i};
    
    concentration_index = find(k_concentration_nm == selected_concentration);
    [type_found, type_index] = find_inside_cell(bundle_type, selected_bundle_type);
    %find all buckling runs for a given concentration
    [common_found, selected_index] = find_common_elements(type_index, concentration_index);
    %before removal of NaNs record all forces and lengths
    n_points = length(selected_index);
    force_points = force_pn(selected_index);
    
    % display information about given concentration
    fprintf('Selected concentration: %i nM\n',selected_concentration);
    fprintf('Bundle type: %s\n',selected_bundle_type);
    fprintf('Measurements found: %i\n',length(selected_index));
    
    %remove NaNs from force
    relevant_points = find(~isnan(force_points));
    %number of points after NaN removal
    n_relevant = length(relevant_points);
    fprintf('After removal of NaN values: %i\n', n_relevant);
    fprintf('\n');
    
    f = force_points(relevant_points);
    c = selected_concentration*ones(1,n_relevant);
    
    if n_relevant > 0
        count = count + 1;
        legend_string{count} = [num2str(selected_concentration) ' nM'];
       
        all_means(i) = mean(f);
        all_forces{i} = f;
        all_conc(i) = selected_concentration;
        all_errors(i) = std(f);
        
        %actually make a plot
        %individual datapoints
        figure(h_fig);
        h_conc_plot = plot(c,f,'o','Color',color_order(color_index,:));
        set(h_conc_plot,'Parent', annotation_groups(i));
        conc_color = get(h_conc_plot, 'Color');
        
        set(get(get(annotation_groups(i),'Annotation'),'LegendInformation'),...
            'IconDisplayStyle','on');
         
        display('plot done... moving on to the next c');
        fprintf('\n');
    end;
end;


     
%format plot
set(gca,'Xscale','log','Yscale','lin','YGrid','on','XGrid','on');
legend(legend_string);
title(['Forces: ' selected_bundle_type]);
xlabel('Motor concentration (nM)');
ylabel('Force (pN)');
prettify_plot(gcf);
ylim([0 2]);
box(gca, 'on');
set(gcf,'Position',[100 100 900 400]);
%errorbar plot for the individual datapoints
        h_mean_plot = errorbar(all_conc, all_means, all_errors,...
            'MarkerSize',30,...
            'Marker','.',...
            'LineStyle','none',...
            'LineWidth',2,...
            'Color',conc_color);

%% plot force v length
count = 0;
%force is only measured for buckling bundles
disp('Plotting force v concentration');
selected_bundle_type = 'Buckling';
n_concentrations = length(list_concentrations);
%create new figures and save handles
annotation_groups = [];
%h_fig_fvl is for force v length plot
h_fig_fvl = figure;
%save color order
color_order = get(gca,'ColorOrder');
hold all;
for i = 1:n_concentrations
    annotation_groups(i) = hggroup;
end;

legend_string = {};
all_conc = {};
all_means = {};
all_lengths = {};
all_forces = {};

for i = 1:n_concentrations
    %analyze concentration 1-after-1
    selected_concentration = list_concentrations{i};
    
    concentration_index = find(k_concentration_nm == selected_concentration);
    [type_found, type_index] = find_inside_cell(bundle_type, selected_bundle_type);
    %find all buckling runs for a given concentration
    [common_found, selected_index] = find_common_elements(type_index, concentration_index);
    %before removal of NaNs record all forces and lengths
    n_points = length(selected_index);
    force_points = force_pn(selected_index);
    length_points = length_um(selected_index);
    
    % display information about given concentration
    fprintf('Selected concentration: %i nM\n',selected_concentration);
    fprintf('Bundle type: %s\n',selected_bundle_type);
    fprintf('Measurements found: %i\n',length(selected_index));
    
    %remove NaNs from force and length measurements
    force_no_nans = find(~isnan(force_points));
    length_no_nans = find(~isnan(length_points));
    [common_found, relevant_points] = find_common_elements(force_no_nans, length_no_nans);
    %number of points after NaN removal
    n_relevant = length(relevant_points);
    fprintf('After removal of NaN values: %i\n', n_relevant);
    fprintf('\n');
    
    f = force_points(relevant_points);
    l = length_points(relevant_points);
    c = selected_concentration*ones(1,n_relevant);
    
    if n_relevant > 0
        count = count + 1;
        all_means{i} = mean(f);
        all_forces{i} = f;
        all_conc{i} = selected_concentration;
        all_lengths{i} = l;
        
        %actually make a plot
        %individual datapoints
        
        figure(h_fig_fvl);
        hold on;
        h_fvl_plot = plot(l,f,...
            'MarkerSize',10,...
            'Marker','o',...
            'LineStyle','none',...
            'LineWidth',2,...
            'MarkerFaceColor',color_order(count,:),...
            'Color',color_order(count,:));
        
        display('plot done... moving on to the next c');
        fprintf('\n');
        legend_string{count} = [num2str(selected_concentration) ' nM'];
    end;
end;

figure(h_fig_fvl);
set(gca,'Xscale','lin','Yscale','lin','YGrid','on','XGrid','on');
legend(legend_string);
title('Forces are correlated with bundle length');
xlabel('Initial length (\mum)');
ylabel('Force (pN)');
prettify_plot(gcf);
ylim([0 2]);
box(gca, 'on');
set(gcf,'Position',[100 100 900 400]);

