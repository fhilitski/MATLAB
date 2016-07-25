
%% Collect all relevant data from files
clear all;
close all;
clc;

%define KCL and PEG concentrations
cPEG = 1; %percent
cKCL = 300; %mM


data = read_data('All data');
total = length(data);

%% Split data into groups according to ID and Day variables in data files

%read the first day
current_day = data{1}.Day;

%define colors
colors = ['r' 'g' 'b' 'c' 'm' 'k' 'r' 'g' 'b' 'c' 'm' 'k'];
current_color = 1;
legend_title = {};
legend_title{1} = current_day;

for i = 1:total
    
    if (strcmp(data{i}.Day,current_day))
        data{i}.color = current_color;
    else
        current_day = data{i}.Day;
        current_color = current_color + 1;
        data{i}.color = current_color;
        legend_title{current_color} = current_day;
    end;
end;

%one gruop for each dataset and one extra group for fitted curves
total_groups = current_color + 1;


%% Plot umbrella results on one graph, fit them and determine lambda
% and confidence intervals for each;

h_lambda_fig = figure(1);
h_lamba_vs_overlap_fig = figure(2);
h_widths_fig = figure(3);
h_potential_fig = figure(4);
h_lambda_hist = figure(5);
figure(h_potential_fig);
hold all;


%for plotting convenience, determine min and max separations
% min_x = min(bead_sep{1});
% max_x = max(bead_sep{1});


% group old and new data points into groups
% create groups for data and extra one for all fits

groups = [];
for i = 1 : (total_groups)
    groups(i) = hggroup;
end;

%create empty array for storing:
% 1. fittied values of lambda
% 2. fitting constant c
% 3. confidence intervals
% 4. overlap length from fluorescence data
% 5. widths of experimental runs
% 6. widths of calibrations

lambda = zeros(1, total);
const = zeros(1, total);
error_lambda = zeros(2, total);
overlap = zeros(1, total);
widths = zeros(1,total);
widths_cal = zeros(1,total);


%define max and mins to format plot in future
x_min = 100;
x_max = 0;
u_min = 0;
u_max = 0;

%get spearations, poetentials and weights for each run

for i = 1 : total
    
    day = data{i}.Day;
    id = data{i}.ID;
    
    %this is now implemented in read_data.
    %     %old versions of bead tracking software recorded id in different format
    %     %get proper id to ensure cross-compatibility
    %     if ~strcmp(id(end),'_')
    %         id = [id '_'];
    %         data{i}.ID = id;
    %     end;
    
    color = data{i}.color;
    
    eval(['x = data{' num2str(i) '}.Ru' id day ';']);
    eval(['u = data{' num2str(i) '}.U' id day '(:,1);']);
    eval(['weight = data{' num2str(i) '}.U' id day '(:,4);']);
    
    %run_sep contains actual recorded separations durnig experimental run
    eval(['run_sep = data{' num2str(i) '}.R' id day ';']);
    widths(i) = std(run_sep); %widths(i) is standard deviation of exp. run i
    %run_sep contains actual recorded separations durnig experimental run
    eval(['cal_sep = data{' num2str(i) '}.R' id day 'cal;']);
    widths_cal(i) = std(cal_sep); %widths(i) is standard deviation of exp. run i
    
    
    %perform weighted linear fit
    fitobject = fit(x',u,'poly1','Weight', weight);
    
    %find coefficients of the fit
    %save them into an array of values for future reference
    fit_coeffs = coeffvalues(fitobject);
    lambda(i) = fit_coeffs(1);
    const(i) = fit_coeffs(2);
    
    %get confidence interval
    conf_int = confint(fitobject);
    error_lambda(:,i) = conf_int(:,1)-lambda(i);
    
    %get overlap length
    overlap(i) = str2double(data{i}.Overlap_length);
    
    %plot potential vs separation length in current axes;
    h_plot = plot(x,u,['o' colors(color)]);
    
    %add plot to a specific group
    set(h_plot,'Parent',groups(color));
    
    hold on;
    h_fit_plot = plot(fitobject,'--r');
    set(h_fit_plot,'Parent',groups(total_groups));
    
    x_min = min(x_min,min(x));
    x_max = max(x_max,max(x));
    u_min = min(u_min,min(u));
    u_max = max(u_max,max(u));
end;

title('\fontsize{26}Experimental Data: U(d)=-\lambda(L_{total}-d)');
xlabel('\fontsize{24}Bead separation d, \mum');
ylabel('\fontsize{24}Binding enery, k_BT');
xlim([x_min x_max]);
ylim([u_min u_max]);

hold off;

%assign legend to groups
for i = 1 : total_groups-1
    set(get(get(groups(i),'Annotation'),'LegendInformation'),...
        'IconDisplayStyle','on');
    
end;

legend(legend_title);

%also, lets find mean lambda and std
lambda_mean = mean(lambda);
lambda_std = std(lambda)/sqrt(total-1);

fprintf(1,'Concentrations: PEG = %1.1f; KCL = %3.0f mM.\n',cPEG,cKCL);
fprintf(1,'Mean value of lambda +- error: %2.2f +\\- %2.2f\n',lambda_mean,lambda_std);
fprintf(1,'Total number of measurements: %u. \n',total);

hold off;



%% Plot all values of lambda on separate scatter plot

figure(h_lambda_fig);
hold on;

groups = [];
for i = 1 : (total_groups)
    groups(i) = hggroup;
end;


for i = 1 : total
    
    day = data{i}.Day;
    id = data{i}.ID;
    color = data{i}.color;
    
    
    %plot lambda;
    h_plot = errorbar(i,lambda(i),error_lambda(1,i),error_lambda(2,i),['o' colors(color)]);
    
    %add plot to a specific group
    set(h_plot,'Parent',groups(color));
    hold on;
end;

tmp_x = 1:length(lambda);
tmp_l_m = zeros(1,length(lambda));
tmp_l_m = tmp_l_m + lambda_mean;
h_fit_plot = line(tmp_x,tmp_l_m,'Color','b','LineStyle','--');
set(h_fit_plot,'Parent',groups(total_groups));

%assign legend to groups
for i = 1 : total_groups
    set(get(get(groups(i),'Annotation'),'LegendInformation'),...
        'IconDisplayStyle','on');
    
end;

legend_title{total_groups} = '<\lambda>';
legend(legend_title);
title(['\fontsize{26}PEG = ' num2str(cPEG) '%, [KCL] = ' num2str(cKCL) 'mM']);
ylabel('\fontsize{24}\lambda, k_BT/\mum');
xlabel('\fontsize{24}Measurement #');
set(gca,'Fontsize',14);

hold off;

%% Plot all values of lambda vs overlap length;

figure(h_lamba_vs_overlap_fig);
hold on;

groups = [];
for i = 1 : (total_groups)
    groups(i) = hggroup;
end;


for i = 1 : total
    
    day = data{i}.Day;
    id = data{i}.ID;
    color = data{i}.color;
    
    if (overlap(i)~= -1)
        %plot lambda;
        h_plot = errorbar(overlap(i),lambda(i),error_lambda(1,i),error_lambda(2,i),['o' colors(color)]);
        
        %add plot to a specific group
        set(h_plot,'Parent',groups(color));
        hold on;
    end;
end;


%assign legend to groups
for i = 1 : total_groups
    set(get(get(groups(i),'Annotation'),'LegendInformation'),...
        'IconDisplayStyle','on');
    
end;

legend(legend_title);
title(['\fontsize{20}\lambda vs Overlap length: \fontsize{15}PEG = ' num2str(cPEG) '%, [KCL] = ' num2str(cKCL) 'mM']);
ylabel('\fontsize{24}\lambda, k_BT/\mum');
xlabel('\fontsize{24} Bundle overlap, \mum');
set(gca,'Fontsize',14);
axis tight;

hold off;

%% Plot all values of widths on separate scatter plot
figure(h_widths_fig);
hold on;

groups = [];
for i = 1 : (total_groups)
    groups(i) = hggroup;
end;


for i = 1 : total
    
    day = data{i}.Day;
    id = data{i}.ID;
    color = data{i}.color;
    
    
    %plot widths;
    h_plot = plot(i,widths(i),['o' colors(color)],i,widths_cal(i),['*' colors(color)]);
    
    %add plot to a specific group
    set(h_plot,'Parent',groups(color));
    hold on;
end;

%average of STDs does not make much sense over separate experiments due to
%different trap strengths

tmp_x = 1:length(lambda);
tmp_l_m = zeros(1,length(lambda));
tmp_l_m = tmp_l_m + widths(1);
h_fit_plot = line(tmp_x,tmp_l_m,'Color','b','LineStyle','none');
set(h_fit_plot,'Parent',groups(total_groups));

%assign legend to groups
for i = 1 : total_groups
    set(get(get(groups(i),'Annotation'),'LegendInformation'),...
        'IconDisplayStyle','on');
    
end;

legend_title{total_groups} = 'o - experiment; * - calibration';
legend(legend_title);
title(['\fontsize{26}PEG = ' num2str(cPEG) ' [KCL] = ' num2str(cKCL) 'mM; STDs of exp. run']);
ylabel('\fontsize{24}\sigma(d), \mum');
xlabel('\fontsize{24}Measurement #');

hold off;


%% Plot histogram of lambda values to see if it resembles gaussian or has two peaks
figure(h_lambda_hist);
[n_bins, bin_size] = get_optimal_binning (lambda);
hist_start = min(lambda);
hist_end = max(lambda);

x = hist_start:bin_size:hist_end;

hist(lambda,5);


%% Save values of lambda, errors on lambda, stds of run and calibration into file

filename = ['Data ' num2str(cPEG) 'PEG ' num2str(cKCL) 'KCL' '.mat'];
save(filename, 'lambda', 'error_lambda', 'widths', 'widths_cal');
