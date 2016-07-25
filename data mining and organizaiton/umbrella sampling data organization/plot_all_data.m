function [error_code, out, figure_handles] = plot_all_data(cPEG, cKCL, Folder_Name, three_bundle)

%PLOT_ALL_DATA plots all data for given PEG and KCL concentrations from
%              provided Folder_Name;
% code of this function repeats old stitch-to-linear routine
out = {};
tb = 0;


if (nargin == 4) && (three_bundle == 1)
    tb = 1;
end;

path = [Folder_Name '\All data'];
data = read_data(path);
total = length(data);

if (total ~= 0)
    
    %Split data into groups according to ID and Day variables in data files
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
    
    h_lambda_fig = figure;
    h_lamba_vs_overlap_fig = figure(gcf+1);
    h_widths_fig = figure(gcf+1);
    h_potential_fig = figure(gcf+1);
    h_lambda_hist = figure(gcf+1);
    h_widths_difference_fig = figure(gcf+1);
    
    figure_handles = {h_lambda_fig, h_lamba_vs_overlap_fig, ...
        h_widths_fig, h_potential_fig, h_lambda_hist, h_widths_difference_fig};
    
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
        
        %get confidence interval (95% by default)
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
    
    
    %also, lets find mean lambda
    lambda_mean = mean(lambda);
    %calculate sample std of
    lambda_std = std(lambda);
    %calculate std of mean (not useful, old and probably incorrect measure
    %of error on lambda)
    lambda_std_mean = std(lambda)/sqrt(length(lambda));
    
    %weighted mean calculation
    error_squared = error_lambda.^2;
    error_squared = error_squared(1,:);
    weighted_lambda_mean = sum(lambda./error_squared) / sum(1./error_squared);
    
    
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
        h_plot = errorbar(i,lambda(i),error_lambda(1,i),error_lambda(2,i),...
            ['o' colors(color)],'MarkerSize',10,'LineWidth',2);
        
        %add plot to a specific group
        set(h_plot,'Parent',groups(color));
        hold on;
    end;
    
    %plot a line for mean lambda;
    tmp_x = 1:length(lambda);
    tmp_l_m = zeros(1,length(lambda));
    tmp_l_m = tmp_l_m + lambda_mean;
    h_fit_plot = line(tmp_x,tmp_l_m,'Color','b','LineStyle','--');
    
    %plot a line for Weighted Mean lambda
    tmp_x = 1:length(lambda);
    tmp_l_m = zeros(1,length(lambda));
    tmp_l_m = tmp_l_m + weighted_lambda_mean;
    h_fit_plot = line(tmp_x,tmp_l_m,'Color','r','LineStyle','--');
   
    
    %assign legend to groups
    for i = 1 : total_groups
        set(get(get(groups(i),'Annotation'),'LegendInformation'),...
            'IconDisplayStyle','on');
        
    end;
    
    legend_title{total_groups} = '<\lambda>';
    legend_title{total_groups + 1} = '<\lambda> weighted';
    legend(legend_title);
    titlestring = ['\fontsize{20}Polymer = ' num2str(cPEG) '%(w/w), [KCL] = ' num2str(cKCL) 'mM'];
    title({'\fontsize{30}Measured values \lambda', titlestring});
    ylabel('\fontsize{20}\lambda, k_BT/\mum');
    xlabel('\fontsize{20}Measurement #');
    set(gca,'Fontsize',20);
    axis tight;
    
    hold off;
    
    %% Plot all values of lambda vs overlap length;
    
    %proceed with plotting only if there are two or more measurements of
    %overlap length
    
    if (~isempty( find(overlap ~= -1) ))
        
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
                h_plot = errorbar(overlap(i),lambda(i),error_lambda(1,i),error_lambda(2,i),...
                    ['o' colors(color)],'MarkerSize',10,'LineWidth',2);
                
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
        
        %plot mean lambda as line on this plot, too
        tmp_x = linspace(0,max(overlap),10);
        tmp_l_m = zeros(1,10);
        tmp_l_m = tmp_l_m + lambda_mean;
        h_fit_plot = line(tmp_x,tmp_l_m,'Color','b','LineStyle','--','LineWidth',1);
        set(h_fit_plot,'Parent',groups(total_groups));
        
        
        %plot a line for Weighted Mean lambda
         tmp_x = linspace(0,max(overlap),10);
        tmp_l_m = zeros(1,10);
        tmp_l_m = tmp_l_m + weighted_lambda_mean;
        h_fit_plot = line(tmp_x,tmp_l_m,'Color','r','LineStyle','--');
       
        
        %assign legend to groups
        for i = 1 : total_groups
            set(get(get(groups(i),'Annotation'),'LegendInformation'),...
                'IconDisplayStyle','on');
            
        end;
        
        legend_title{total_groups} = '<\lambda>';
        legend(legend_title);
        
        titlestring = ['\fontsize{20}Polymer = ' num2str(cPEG) '%(w/w), [KCL] = ' num2str(cKCL) 'mM'];
        title({'\fontsize{30}Experimental \lambda vs L_{Overlap}', titlestring});
        ylabel('\fontsize{20}Measured \lambda, k_BT/\mum');
        xlabel('\fontsize{20}Bundle overlap, \mum');
        set(gca,'Fontsize',20);
        axis tight;
        hold off;
    else
        close(h_lamba_vs_overlap_fig);
    end;
    
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
    
    legend_title{total_groups} = 'o exp.; * cal.';
    legend(legend_title);
    titlestring = ['\fontsize{20}Polymer = ' num2str(cPEG) '%(w/w), [KCL] = ' num2str(cKCL) 'mM'];
    title({'\fontsize{30} Standard deviations of runs \sigma, \mum', titlestring});
    ylabel('\fontsize{20}\sigma(d), \mum');
    xlabel('\fontsize{20}Measurement #');
    set(gca,'Fontsize',20);
    axis tight;
    
    hold off;
    
    
    %% Plot histogram of lambda values to see if it resembles gaussian or has two peaks
    figure(h_lambda_hist);
    [n_bins, bin_size] = get_optimal_binning (lambda);
    hist_start = min(lambda);
    hist_end = max(lambda);
    x = hist_start:bin_size:hist_end;
    hist(lambda,n_bins);
    
    titlestring = ['\fontsize{20}Polymer = ' num2str(cPEG) '%(w/w), [KCL] = ' num2str(cKCL) 'mM'];
    title({'\fontsize{30} Histogram of measured values \lambda', titlestring});
    ylabel('\fontsize{20} Count');
    xlabel('\fontsize{20} Measured \lambda, k_BT');
    set(gca,'Fontsize',20);
    axis tight;
    
    
    %% Plot difference in standard deviations (of probability distributions)
    %  between calibration and experiemtns - in % change
    % we are going to plot sigma
    figure(h_widths_difference_fig);
    hold on;
    
    groups = [];
    for i = 1 : (total_groups)
        groups(i) = hggroup;
    end;
    
    difference = (widths_cal - widths) ./ widths_cal *100;
    for i = 1 : total
        day = data{i}.Day;
        id = data{i}.ID;
        color = data{i}.color;
        
        %plot differences widths;
        h_plot = plot(i, difference(i), ['o' colors(color)]);
        
        %add plot to a specific group
        set(h_plot,'Parent',groups(color));
        hold on;
    end;
    
    %this plots mean difference line
    tmp_l_m = tmp_l_m + mean(difference);
    h_fit_plot = line(tmp_x,tmp_l_m,'Color','b','LineStyle','--');
    
    set(h_fit_plot,'Parent',groups(total_groups));
   
     
    tmp_x = 1:length(difference);
    tmp_l_m = zeros(1,length(difference));
    %this plots 0% line
    h_fit_plot = line(tmp_x,tmp_l_m,'Color','k','LineStyle',':');
    
    %assign legend to groups
    for i = 1 : total_groups
        set(get(get(groups(i),'Annotation'),'LegendInformation'),...
            'IconDisplayStyle','on');
        
    end;
    
    legend_title{total_groups} = 'Mean';
    legend_title{total_groups + 1} = '0% line';
    legend(legend_title);
    titlestring = ['\fontsize{20}Polymer = ' num2str(cPEG) '%(w/w), [KCL] = ' num2str(cKCL) 'mM'];
    titlestring2 = '\fontsize{20}(\sigma_{CAL}-\sigma_{EXP})/\sigma_{CAL}, %';
    title({'\fontsize{30} Relative difference in St.Dev.', titlestring2, titlestring});
    ylabel('\fontsize{20}\Delta\sigma/\sigma');
    xlabel('\fontsize{20}Measurement #');
    set(gca,'Fontsize',20);
    
    hold off;
    
    
    
    out1 = sprintf('Concentrations: PEG = %1.1f; KCL = %3.0f mM.\n',cPEG,cKCL);
    out2 = sprintf('Mean value of lambda: %2.2f +\\- %2.2f\n',lambda_mean, lambda_std);
    out3 = sprintf('Std of the mean: %2.2f;\n', lambda_std_mean);
    out4 = sprintf('Std of lambda: +\\-%2.2f.\n',lambda_std);
    out5 = sprintf('Weighted mean lambda: %2.2f.\n', weighted_lambda_mean);
    out6 = sprintf('Mean of Relative difference between STDs of cal and exp: %2.2f%%',mean(difference));
    out7 = sprintf('Total number of measurements: %u. \n',total);
    if tb out0 = sprintf('3-bundle data! \n');
    else out0 = sprintf('2-bundle data! \n');
    end;
    out = {out0; out1; out2; out4; out3; out5; out6; out7};
    
    
    %Save values of lambda, errors on lambda, stds of run and calibration into file
    filename = [Folder_Name '\Data ' num2str(cPEG) 'PEG ' num2str(cKCL) 'KCL' '.mat'];
    save(filename, 'lambda', 'error_lambda', 'overlap', 'widths', 'widths_cal','difference');
    
    %save figures of histogram, lambda vs overlap and all measured lambdas
    %            h_lambda_hist  h_lamba_vs_overlap_fig  h_lambda_fig
    %fisrst, we need to generate filenames
    fname_hist = [Folder_Name '\histogram_'  num2str(cPEG) 'PEG_' num2str(cKCL) 'KCl'];
    fname_lambda = [Folder_Name '\lambdas_'  num2str(cPEG) 'PEG_' num2str(cKCL) 'KCl'];
    fname_lambda_vs_overlap = [Folder_Name '\lambda_vs_L_'  num2str(cPEG) 'PEG_' num2str(cKCL) 'KCl'];
    saveas(h_lambda_hist, [fname_hist '.tiff']);
    saveas(h_lambda_hist, [fname_hist '.fig']);
  
    saveas(h_lambda_fig, [fname_lambda '.tiff']);
    saveas(h_lambda_fig, [fname_lambda '.fig']);
    
    try
        saveas(h_lamba_vs_overlap_fig, [fname_lambda_vs_overlap '.tiff']);
        saveas(h_lamba_vs_overlap_fig, [fname_lambda_vs_overlap '.fig']);
    catch exc;
    end;
    
    %% Error code return
    error_code = 0;
else
    error_code = -1;
    out = 0;
    figure_handles = 0;
end;


