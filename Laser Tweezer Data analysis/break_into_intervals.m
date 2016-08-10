%%
clc;
clear all;
close all;
[fname,filepath,f_idx] = uigetfile('*.mat','select data file');
load([filepath fname]);

%% fix issue with many previous runs that do not have bead_dist_fd_curve variable
if (~exist('bead_dist_fd_curve', 'var'))
    bead_dist_fd_curve = trap_dist_fd_curve;
    warning('bead_dist_fd_curve does not exist');
end;

%%
%this is a temporary variable re-assignment
trap_dist_fd_curve = trap_dist_fd_curve_r;
total_force_fd_curve = total_force_fd_curve_r;
force_x_fd_curve = force_x_fd_curve_r;
force_y_fd_curve = force_y_fd_curve_r;


%% divide trap separation into intervals
%this script breaks force-displacement data into intervals and plots them\
%separately. Intervals are manually selected by the user based on the trap
%separation vs time plot

%define color for raw data
color_raw = [0.831  0.816    0.784];
%define color for smoothed data
color_smooth = [ 0.502  0.502   0.502];
%define color for block-averaged data
color_average = [0  0.447   0.741];
units_xy = 'pN';


%the neccessary variables are
%trap_dist_fd_curve = d;
total_force_fd_curve;

n_data = length(trap_dist_fd_curve); %data length

%ask the user to privide intervals for the different regions...
h_fig_ui = figure;
plot(1:1:n_data,trap_dist_fd_curve);
axis tight;
hold on;

[index_d_borders,d_borders] = ginput;
index_d_borders = round(index_d_borders);
n_intervals = length(d_borders)+1; %number of intervals provided by the user

%initialize variables for to store data for separat intervals
d_start = 0;
d_end = 0;
intervals = cell(1,n_intervals);
distances = intervals;
forces = intervals;
forces_x = intervals;
forces_y = intervals;
bead_distances = intervals;
bead_distances_bf = intervals;

%split data into the intervals according to the user input
for i = 1:n_intervals
    d_start = d_end+1;
    if i == n_intervals
        d_end = n_data;
    else
        d_end = index_d_borders(i);
    end;
    
    intervals{i} = d_start:1:d_end;
    distances{i} = trap_dist_fd_curve(intervals{i});
    forces{i} = total_force_fd_curve(intervals{i});
    forces_x{i} = force_x_fd_curve(intervals{i});
    forces_y{i} = force_y_fd_curve(intervals{i});
    bead_distances{i} = bead_dist_fd_curve(intervals{i});
    bead_distances_bf{i} = bead_dist_fd_curve_bf(intervals{i});
    
    %display regions on the plot.
    h_line = plot(intervals{i}, distances{i},'.');
    h_line.DisplayName = ['Interval: ' num2str(i)];
    h_line.Tag = [num2str(i)];
end;



%% plot separation v frame
figure;
for i = 1: n_intervals
    plot(distances{i});
    if i == 1
        hold;
    end;
end;
hold off;

% plot forces vs frame
figure;
for i = 1: n_intervals
    plot(forces{i});
    if i == 1
        hold on;
    end;
end;
hold off;

%% plot force v distance between traps for each interval
n_intervals = length(forces);
h_f = figure;
h_box = annotation(h_f,'textbox',[0.5 0.8 0.04 0.1],'FitBoxToText','on');
pause on;
for i = 1: n_intervals
    h_line = plot(distances{i},forces{i},'.-');
    % Create textbox
    set(h_box,'String',{['Interval: ' num2str(i)]});
    hold on;
    if (i < n_intervals)
        pause;
    end;
end;
hold off;
pause off;
ylabel(['Force (' units_xy ')']);
xlabel('Trap separation (nm)');
title('X-axis Force - distance curve');
axis tight

%% plot force v bead distance for each interval

h_f = figure;
h_box = annotation(h_f,'textbox',[0.5 0.8 0.04 0.1],'FitBoxToText','on');
pause on;
for i = 1: n_intervals
    h_line = plot(bead_distances{i},forces{i},'.-');
    hold on;
    h_line_bf = plot(bead_distances_bf{i},forces{i},'.-');
    % Create textbox
    set(h_box,'String',{['Interval: ' num2str(i)]});
    if (i < n_intervals)
        pause;
    end;
end;
hold off;
pause off;
ylabel(['Force (' units_xy ')']);
xlabel('Trap separation (nm)');
title('X-axis Force - distance curve');
axis tight;

%% analysis of the relevant force-displacement intervals
%define color for raw data
color_raw = [0.831  0.816    0.784];
%define color for smoothed data
color_smooth = [ 0.502  0.502   0.502];
%define color for block-averaged data
color_average = [0  0.447   0.741];
units_xy = 'pN';

d = [];
b_d = []; %bead distance
b_d_bf = []; %bead distance from brightfield
f = [];
f_y = [];
f_x = [];
if exist('relevant','var')
    choice = questdlg('Define new intervals?','Data analysis','Yes','No','Yes');
    if strcmp(choice,'Yes')    
        relevant = input('Please define relevant inervals of the exp run: ');
    else 
        disp('Keeping existing intervals...');
    end;
else
    relevant = input('Please define relevant inervals of the exp run: ');
end;

for i=relevant
    d = cat(1, d, distances{i});
    b_d = cat(1, b_d, bead_distances{i});
    b_d_bf = cat(1, b_d_bf, bead_distances_bf{i});
    f = cat(1, f, forces{i});
    f_y = cat(1, f_y, forces_y{i});
    f_x = cat(1, f_x, forces_x{i});
end;

%% plot total force v displacement
h_force_vs_displacement_corrected = figure;
plot(d, f,'.','Color',color_raw);
hold on;
%plot(d, smooth(f,5),'.','LineWidth',2,'Color',color_smooth);
ylabel(['Force (' units_xy ')']);
xlabel('Trap separation (nm)');
title('Force - distance curve');
axis tight;
grid on;

%bin force based on the fixed displacement values and plot that
[trap_dist_binned, force_binned, force_error] = binned_xy(d,f);
errorbar(trap_dist_binned, force_binned, force_error,'o','Color',color_average,'LineWidth',1);
plot(trap_dist_binned, force_binned,'.');

%find the separation that corresponds to min. force - that's the
%equilibrium separation
[f_min, i_f_min] = min(force_binned);
d_min = trap_dist_binned(i_f_min);
fprintf(1,'Equilibrium trap separation: %.0f nm\n',d_min);
hold on;
plot(d_min, f_min, 'or','LineWidth',3,'MarkerSize',10);

%max force - this is not very useful in case of total force,
%as the max force is usually reached during the extension phase
%so it will not provide us with the force plateau
% [f_max, i_f_max] = max(force_binned);
% d_max = trap_dist_binned(i_f_max);
% fprintf(1,'Max force: %.2f pN\n',f_max);
% hold on;
% plot(d_max, f_max, 'og');

%plot y-force vs displacement
h_force_y_vs_displacement_corrected = figure;
plot(d, f_y,'.','Color',color_raw);
hold on;
grid on;
%plot(d, smooth(f_y,5),'.','LineWidth',2,'Color',color_smooth);
ylabel(['Force (' units_xy ')']);
xlabel('Trap separation (nm)');
title('Force (y) - distance curve');
axis tight;
%bin force based on the fixed displacement values and plot that
[trap_dist_binned_y, force_binned_y, force_error_y] = binned_xy(d,f_y);
errorbar(trap_dist_binned_y, force_binned_y, force_error_y,'o','Color',color_average,'LineWidth',1);
plot(trap_dist_binned_y, force_binned_y,'.');
%find the max force
%since this is done in the reference frame connected with the traps,
%this will give us the force plateau
[f_max_y, i_f_max_y] = max(force_binned_y);
d_max_y = trap_dist_binned(i_f_max_y);
fprintf(1,'Max force: %.2f pN\n',f_max_y);
hold on;
plot(d_max_y, f_max_y, 'og');
plot(d_min, force_binned_y(find(trap_dist_binned_y == d_min)), 'or',...
    'LineWidth',3,'MarkerSize',10);

%plot x-force vs displacement - this is just for fun
h_force_x_vs_displacement_corrected = figure;
plot(d, f_x,'.','Color',color_raw);
hold on;
%plot(d, smooth(f_x,5),'.','LineWidth',2,'Color',color_smooth);
ylabel(['Force (' units_xy ')']);
xlabel('Trap separation (nm)');
title('Force (x) - distance curve');
axis tight;
grid on;
%bin force based on the fixed displacement values and plot that
[trap_dist_binned_x, force_binned_x, force_error_x] = binned_xy(d,f_x);
errorbar(trap_dist_binned_x, force_binned_x, force_error_x,'o','Color',color_average,'LineWidth',1);
plot(trap_dist_binned_x, force_binned_x,'.');

%% plot forces vs bead separation

%plot total force v bead separation
h_force_vs_bead_separation_corrected = figure;
plot(b_d, f,'.','Color',color_raw);
hold on;
%plot(d, smooth(f,5),'.','LineWidth',2,'Color',color_smooth);
ylabel(['Force (' units_xy ')']);
xlabel('Bead separation (nm)');
title('Force-displacement curve');
axis tight;
grid on;

%bin force based on the calculated bead separation values and plot that
[bead_dist_binned, force_binned_b, force_error_b] = binned_xy(b_d,f);
errorbar(bead_dist_binned, force_binned_b, force_error_b,'o','Color',color_average,'LineWidth',1);
plot(bead_dist_binned, force_binned_b,'.');

%plot y-force vs displacement
h_force_y_vs_bead_separation_corrected = figure;
plot(b_d, f_y,'.','Color',color_raw);
hold on;
grid on;
%plot(d, smooth(f_y,5),'.','LineWidth',2,'Color',color_smooth);
ylabel(['Force (' units_xy ')']);
xlabel('Bead separation (nm)');
title('Force (y)-displacement curve');
axis tight;
%bin force based on the fixed displacement values and plot that
[bead_dist_binned_y, force_binned_y_b, force_error_y_b] = binned_xy(b_d,f_y);
errorbar(bead_dist_binned_y, force_binned_y_b, force_error_y_b,'o','Color',color_average,'LineWidth',1);
plot(bead_dist_binned_y, force_binned_y_b,'.');
%find the max force
%this will give us the force plateau
[f_max_y_b, i_f_max_y_b] = max(force_binned_y_b);
d_max_y_b = bead_dist_binned(i_f_max_y_b);
fprintf(1,'Max force (from bead separation): %.2f pN\n',f_max_y_b);
hold on;
plot(d_max_y_b, f_max_y_b, 'og',...
    'LineWidth',3,'MarkerSize',10);
plot(d_min, force_binned_y_b(find(bead_dist_binned_y >= d_min, 1)), 'or',...
    'LineWidth',3,'MarkerSize',10);


%plot x-force vs displacement - this is just for fun
h_force_x_vs_bead_separation_corrected = figure;
plot(b_d, f_x,'.','Color',color_raw);
hold on;
%plot(b_d, smooth(f_x,5),'.','LineWidth',2,'Color',color_smooth);
ylabel(['Force (' units_xy ')']);
xlabel('Bead separation (nm)');
title('Force (x)-displacement curve');
axis tight;
grid on;
%bin force based on the fixed displacement values and plot that
[bead_dist_binned_x, force_binned_x_b, force_error_x_b] = binned_xy(b_d,f_x);
errorbar(bead_dist_binned_x, force_binned_x_b, force_error_x_b,'o','Color',color_average,'LineWidth',1);
plot(bead_dist_binned_x, force_binned_x_b,'.');
%% plot forces vs bead separation in brightfield

%plot total force v bead separation
h_force_vs_bead_separation_corrected_bf = figure;
plot(b_d_bf, f,'.','Color',color_raw);
hold on;
%plot(d, smooth(f,5),'.','LineWidth',2,'Color',color_smooth);
ylabel(['Force (' units_xy ')']);
xlabel('Brightfield bead separation (nm)');
title('Force-displacement curve');
axis tight;
grid on;

%bin force based on the calculated bead separation values and plot that
[bead_dist_binned_bf, force_binned_b_bf, force_error_b_bf] = binned_xy(b_d_bf,f);
errorbar(bead_dist_binned_bf, force_binned_b_bf, force_error_b_bf,'o','Color',color_average,'LineWidth',1);
plot(bead_dist_binned_bf, force_binned_b_bf,'.');

%plot y-force vs displacement
h_force_y_vs_bead_separation_corrected_bf = figure;
plot(b_d_bf, f_y,'.','Color',color_raw);
hold on;
grid on;
%plot(d, smooth(f_y,5),'.','LineWidth',2,'Color',color_smooth);
ylabel(['Force (' units_xy ')']);
xlabel('Brightfield bead separation (nm)');
title('Force_y-displacement curve');
axis tight;
%bin force based on the fixed displacement values and plot that
[bead_dist_binned_y_bf, force_binned_y_b_bf, force_error_y_b_bf] = binned_xy(b_d_bf,f_y);
errorbar(bead_dist_binned_y_bf, force_binned_y_b_bf, force_error_y_b_bf,'o','Color',color_average,'LineWidth',1);
plot(bead_dist_binned_y_bf, force_binned_y_b_bf,'.');
%find the max force
%this will give us the force plateau
[f_max_y_b_bf, i_f_max_y_b_bf] = max(force_binned_y_b_bf);
d_max_y_b_bf = bead_dist_binned_bf(i_f_max_y_b_bf);
fprintf(1,'Max force (from !brigtfield! bead separation): %.2f pN\n',f_max_y_b_bf);
hold on;
plot(d_max_y_b_bf, f_max_y_b_bf, 'og',...
    'LineWidth',3,'MarkerSize',10);
plot(d_min, force_binned_y_b_bf(find(bead_dist_binned_y_bf >= d_min, 1)), 'or',...
    'LineWidth',3,'MarkerSize',10);


%plot x-force vs displacement - this is just for fun
h_force_x_vs_bead_separation_corrected_bf = figure;
plot(b_d_bf, f_x,'.','Color',color_raw);
hold on;
%plot(b_d, smooth(f_x,5),'.','LineWidth',2,'Color',color_smooth);
ylabel(['Force (' units_xy ')']);
xlabel('Brightfield bead separation (nm)');
title('Force_x-displacement curve');
axis tight;
grid on;
%bin force based on the fixed displacement values and plot that
[bead_dist_binned_x_bf, force_binned_x_b_bf, force_error_x_b_bf] = binned_xy(b_d_bf,f_x);
errorbar(bead_dist_binned_x_bf, force_binned_x_b_bf, force_error_x_b_bf,'o','Color',color_average,'LineWidth',1);
plot(bead_dist_binned_x_bf, force_binned_x_b_bf,'.');

%% plot forces vs MT strain

%we now know MT equilibrium length and can re-calculate MT strain
%strain = (bead_dist - MT_length/ MT_length)
%strain based on the bead separation:
mt_strain_bead = (bead_dist_binned - d_min)./d_min;
%trap distance based strain:
mt_strain_trap = (trap_dist_binned - d_min)./d_min;
%strain based on brightfield separation:
mt_strain_bead_bf = (bead_dist_binned_bf - d_min)./d_min;

%Get strain at max force (based on bead separation):
max_strain_b = (bead_dist_binned(i_f_max_y_b) - d_min) / d_min;
fprintf(1,'Strain at max force (from bead separation): %.1f%%\n',max_strain_b*100);
%Get strain at max force (based on the trap separation):
max_strain = (trap_dist_binned(i_f_max_y) - d_min) / d_min;
fprintf(1,'Strain at max force: %.1f%%\n',max_strain*100);

%Get strain at max force (based on bead separation from brightfield):
max_strain_b_bf = (bead_dist_binned_bf(i_f_max_y_b_bf) - d_min) / d_min;
fprintf(1,'Strain at max force (from bead separation): %.1f%%\n',max_strain_b_bf*100);

%Get force at given strain +/- strain threshold
%strain_level = -0.20;
for strain_level = [-0.20 -0.10 -0.05]
    strain_threshold = 0.005; %comparing strains within interval level +/- threshold
    s1 = strain_level + strain_threshold;
    s2 = strain_level - strain_threshold;
    idx_strain = find((mt_strain_bead <= s1)&(mt_strain_bead >= s2));
    idx_strain_bf = find((mt_strain_bead_bf <= s1)&(mt_strain_bead_bf >= s2));
    force_at_strain = force_binned_b(idx_strain);
    force_at_strain_bf = force_binned_b_bf(idx_strain_bf);
    strains_level = mt_strain_bead(idx_strain);
    strains_level_bf = mt_strain_bead_bf(idx_strain_bf);
    fprintf(1,'Force at %4.1f%% strain is %.2f pN (N=%i)\n',...
        mean(strains_level) * 100, mean(force_at_strain), length(idx_strain));
    fprintf(1,'Brightfield data - Force at %4.1f%% strain is %.2f pN (N=%i)\n',...
        mean(strains_level_bf) * 100, mean(force_at_strain_bf), length(idx_strain_bf));
end;

%plot total force v straun
h_force_vs_strain = figure;

errorbar(mt_strain_trap, force_binned, force_error,'o','LineWidth',1);
hold on;
errorbar(mt_strain_bead, force_binned_b, force_error_b,'o','LineWidth',1);
errorbar(mt_strain_bead_bf, force_binned_b_bf, force_error_b_bf,'o','LineWidth',1);
plot(max_strain_b, f_max_y_b, 'og',...
    'LineWidth',3,'MarkerSize',10);
plot(max_strain, f_max_y, 'or',...
    'LineWidth',3,'MarkerSize',10);
plot(max_strain_b_bf, f_max_y_b_bf, 'ob',...
    'LineWidth',3,'MarkerSize',10);


ylabel(['Force (' units_xy ')']);
xlabel('MT strain (dimensionless)');
title('Force-strain curve');
axis tight;
grid on;
legend({'Trap separation', 'Bead Separation', 'Brightfield bead separation'});

%% plot force vs trap separation and bead displacement to compare
h_force_vs_traps_and_beads = figure;
plot(trap_dist_binned,force_binned,'.-');
hold on;
plot(bead_dist_binned,force_binned_b,'.-');
plot(bead_dist_binned_bf,force_binned_b_bf,'.-');
ylabel(['Force (' units_xy ')']);
xlabel('Bead separation (nm)');
title('Force vs bead and trap separation');
axis tight;
grid on;
legend({'Trap separation', 'Bead Separation'});


%% analysis with fraying
fraying = false;
if (fraying)
    d_fwd = [];
    d_back = [];
    f_back = [];
    f_fwd = [];
    
    f_back_y = [];
    f_back_x = [];
    f_fwd_y = [];
    f_fwd_x = [];
    
    relevant_fwd = [2];
    relevant_back = [3];
    
    for i=relevant_fwd
        d_fwd = cat(1, d_fwd, distances{i});
        f_fwd = cat(1, f_fwd, forces{i});
        f_fwd_y = cat(1, f_fwd_y, forces_y{i});
        f_fwd_x = cat(1, f_fwd_x, forces_x{i});
    end;
    
    for i=relevant_back
        d_back = cat(1, d_back, distances{i});
        f_back = cat(1, f_back, forces{i});
        f_back_y = cat(1, f_back_y, forces_y{i});
        f_back_x = cat(1, f_back_x, forces_x{i});
    end;
    
    d = d_fwd;
    f = f_fwd;
    f_y = f_fwd_y;
    f_x = f_fwd_x;
    %plot total force v displacement for forward runs
    h_force_vs_displacement_corrected_fwd = figure;
    plot(d, f,'.','Color',color_raw);
    hold on;
    plot(d, smooth(f,5),'.','LineWidth',2,'Color',color_smooth);
    ylabel(['Force (' units_xy ')']);
    xlabel('Trap separation (nm)');
    title('Forward runs: Force - distance curve');
    axis tight;
    %bin force based on the fixed displacement values and plot that
    [trap_dist_binned, force_binned, force_error] = binned_xy(d,f);
    errorbar(trap_dist_binned, force_binned, force_error,'o','Color',color_average,'LineWidth',1);
    plot(trap_dist_binned, force_binned,'.');
    
    %plot y-force vs displacement
    h_force_y_vs_displacement_corrected_fwd = figure;
    plot(d, f_y,'.','Color',color_raw);
    hold on;
    grid on;
    plot(d, smooth(f_y,5),'.','LineWidth',2,'Color',color_smooth);
    ylabel(['Force (' units_xy ')']);
    xlabel('Trap separation (nm)');
    title('Forward runs: Force (y) - distance curve');
    axis tight;
    %bin force based on the fixed displacement values and plot that
    [trap_dist_binned_y, force_binned_y, force_error_y] = binned_xy(d,f_y);
    errorbar(trap_dist_binned_y, force_binned_y, force_error_y,'o','Color',color_average,'LineWidth',1);
    plot(trap_dist_binned_y, force_binned_y,'.');
    
    f_fwd_binned = force_binned;
    d_fwd_binned = trap_dist_binned;
    e_fwd_binned = force_error;
    f_fwd_binned_y = force_binned_y;
    e_fwd_binned_y = force_error_y;
    %-------------------------------
    
    d = d_back;
    f = f_back;
    f_y = f_back_y;
    %plot total force v displacement for back runs
    h_force_vs_displacement_corrected_back = figure;
    plot(d, f,'.','Color',color_raw);
    hold on;
    plot(d, smooth(f,5),'.','LineWidth',2,'Color',color_smooth);
    ylabel(['Force (' units_xy ')']);
    xlabel('Trap separation (nm)');
    title('Back runs: Force - distance curve');
    axis tight;
    
    %bin force based on the fixed displacement values and plot that
    [trap_dist_binned, force_binned, force_error] = binned_xy(d,f);
    errorbar(trap_dist_binned, force_binned, force_error,'o','Color',color_average,'LineWidth',1);
    plot(trap_dist_binned, force_binned,'.');
    
    %plot y-force vs displacement
    h_force_y_vs_displacement_corrected_back = figure;
    plot(d, f_y,'.','Color',color_raw);
    hold on;
    grid on;
    plot(d, smooth(f_y,5),'.','LineWidth',2,'Color',color_smooth);
    ylabel(['Force (' units_xy ')']);
    xlabel('Trap separation (nm)');
    title('Back runs: Force (y) - distance curve');
    axis tight;
    %bin force based on the fixed displacement values and plot that
    [trap_dist_binned_y, force_binned_y, force_error_y] = binned_xy(d,f_y);
    errorbar(trap_dist_binned_y, force_binned_y, force_error_y,'o','Color',color_average,'LineWidth',1);
    plot(trap_dist_binned_y, force_binned_y,'.');
    
    f_back_binned = force_binned;
    d_back_binned = trap_dist_binned;
    e_back_binned = force_error;
    f_back_binned_y = force_binned_y;
    e_back_binned_y = force_error_y;
    %------------------------
    
    %plot forward abd back runs together
    h_plot_together = figure;
    errorbar(d_fwd_binned, f_fwd_binned, e_fwd_binned,'.-');
    hold on;
    errorbar(d_back_binned, f_back_binned, e_back_binned,'.-');
    ylabel(['Force (' units_xy ')']);
    xlabel('Trap separation (nm)');
    title('Total force - spearation');
    axis tight;
    legend('Forward: buckling and fraying','Back: re-bundling and extension');
    
    %plot forward abd back runs together - only y force
    h_plot_together_y = figure;
    errorbar(d_fwd_binned, f_fwd_binned_y, e_fwd_binned_y,'.-');
    hold on;
    errorbar(d_back_binned, f_back_binned_y, e_back_binned_y,'.-');
    ylabel(['Force (' units_xy ')']);
    xlabel('Trap separation (nm)');
    title('Force (y) - spearation curves: back and forth');
    axis tight;
    legend('Forward: buckling and fraying','Back: re-bundling and extension');
    
end;

%% save data
fp_initial = filepath;
filepath = [fp_initial 'intervals\'];
mkdir(filepath);
filename = [filepath 'force-displacement intervals.mat'];
if (fraying)
    %check if fowrward and back runs were separated
    save([filepath 'fraying analysis binned.mat'],...
        'f_fwd_binned','f_back_binned',...
        'd_fwd_binned','d_back_binned',...
        'e_fwd_binned','e_back_binned',...
        'f_fwd_binned_y','f_back_binned_y',...
        'e_fwd_binned_y','e_back_binned_y',...
        'intervals', 'relevant_fwd', 'relevant_back',...
        'distances', 'forces', 'forces_x','forces_y');
    
    %save fraying figures
    f_filename = [filepath 'fraying_fwd'];
    f_handle = h_force_y_vs_displacement_corrected_fwd;
    prettify_plot(f_handle);
    savefig(f_handle,f_filename);
    saveas(f_handle,f_filename,'jpeg');
    
    f_filename = [filepath 'fraying_back'];
    f_handle = h_force_y_vs_displacement_corrected_back;
    prettify_plot(f_handle);
    savefig(f_handle,f_filename);
    saveas(f_handle,f_filename,'jpeg');
    
    f_filename = [filepath 'fraying_back_and_forward'];
    f_handle = h_plot_together;
    prettify_plot(f_handle);
    savefig(f_handle,f_filename);
    saveas(f_handle,f_filename,'jpeg');
    
    f_filename = [filepath 'fraying_back_and_forward_y'];
    f_handle = h_plot_together_y;
    prettify_plot(f_handle);
    savefig(f_handle,f_filename);
    saveas(f_handle,f_filename,'jpeg');
    
else
    save(filename,...
        'intervals', 'relevant', 'distances', 'bead_distances', 'forces','forces_x','forces_y');
end;

save([filepath 'binned forces.mat'],...
    'trap_dist_binned', 'force_binned', 'force_error',...
    'trap_dist_binned_y', 'force_binned_y', 'force_error_y',...
    'trap_dist_binned_x', 'force_binned_x', 'force_error_x',...
    'bead_dist_binned', 'force_binned_b', 'force_error_b',...
    'bead_dist_binned_y', 'force_binned_y_b', 'force_error_y_b',...
    'bead_dist_binned_x', 'force_binned_x_b', 'force_error_x_b',...
    'mt_strain_bead', 'mt_strain_trap');
disp('Done saving');

%save force-y vs trap displacement
f_filename = [filepath 'force_y vs trap separation'];
prettify_plot(h_force_y_vs_displacement_corrected);
savefig(h_force_y_vs_displacement_corrected,f_filename);
saveas(h_force_y_vs_displacement_corrected,f_filename,'jpeg');

%save total force vs trap displacement
f_filename = [filepath 'force_total vs trap separation'];
prettify_plot(h_force_vs_displacement_corrected);
savefig(h_force_vs_displacement_corrected,f_filename);
saveas(h_force_vs_displacement_corrected, f_filename,'jpeg');

%save force-y vs bead separation
f_filename = [filepath 'force_y vs bead separation'];
fig_ref = h_force_y_vs_bead_separation_corrected;
prettify_plot(fig_ref);
savefig(fig_ref,f_filename);
saveas(fig_ref,f_filename,'jpeg');

%save total force vs trap displacement
f_filename = [filepath 'force_total vs bead separation'];
fig_ref = h_force_vs_bead_separation_corrected;
prettify_plot(fig_ref);
savefig(fig_ref,f_filename);
saveas(fig_ref, f_filename,'jpeg');


%save force vs mt strain
f_filename = [filepath 'force_total vs strain'];
fig_ref = h_force_vs_strain;
prettify_plot(fig_ref);
savefig(fig_ref,f_filename);
saveas(fig_ref, f_filename,'jpeg');

%% open cftool findow
cftool(trap_dist_binned_y, force_binned_y);
cftool(trap_dist_binned, force_binned);
cftool(bead_dist_binned, force_binned_b);
cftool(bead_dist_binned_bf, force_binned_b_bf);