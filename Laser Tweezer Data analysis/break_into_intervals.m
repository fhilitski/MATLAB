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

%% plot forces vs frame
figure;
for i = 1: n_intervals
    plot(forces{i}); 
    if i == 1
        hold on;
    end;
end;
hold off;

 %% plot force v distance for each interval

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


%% analysis of the relevant force-displacement intervals
 %define color for raw data
color_raw = [0.831  0.816    0.784];
%define color for smoothed data
color_smooth = [ 0.502  0.502   0.502];
%define color for block-averaged data
color_average = [0  0.447   0.741];
units_xy = 'pN';


d = [];
f = [];
f_y = [];
f_x = [];
if ~(exist('relevant','var'))
    %if variable describing relelvant intervals does not exist yet
    relevant = 7:1:11;
end;
for i=relevant
    d = cat(1, d, distances{i});
    f = cat(1, f, forces{i});
    f_y = cat(1, f_y, forces_y{i});
    f_x = cat(1, f_x, forces_x{i});
end;

%plot total force v displacement
h_force_vs_displacement_corrected = figure;
plot(d, f,'.','Color',color_raw);
hold on;
plot(d, smooth(f,5),'.','LineWidth',2,'Color',color_smooth);
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
plot(d, smooth(f_y,5),'.','LineWidth',2,'Color',color_smooth);
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
plot(d_min, force_binned_y(find(trap_dist_binned_y <= d_min ,1)), 'or','LineWidth',3,'MarkerSize',10);


%plot x-force vs displacement - this is just for fun
h_force_y_vs_displacement_corrected = figure;
plot(d, f_x,'.','Color',color_raw);
hold on;
plot(d, smooth(f_x,5),'.','LineWidth',2,'Color',color_smooth);
ylabel(['Force (' units_xy ')']);
xlabel('Trap separation (nm)');
title('Force (x) - distance curve');
axis tight;
grid on;
%bin force based on the fixed displacement values and plot that 
[trap_dist_binned_x, force_binned_x, force_error_x] = binned_xy(d,f_x);
errorbar(trap_dist_binned_x, force_binned_x, force_error_x,'o','Color',color_average,'LineWidth',1);
plot(trap_dist_binned_x, force_binned_x,'.');

%% analysis with fraying
d_fwd = [];
d_back = [];
f_back = [];
f_fwd = [];
relevant_fwd = [2,4,6,8,10];
relevant_back = [3,5,7];

for i=relevant_fwd
    d_fwd = cat(1, d_fwd, distances{i});
    f_fwd = cat(1, f_fwd, forces{i});
end;

for i=relevant_back
    d_back = cat(1, d_back, distances{i});
    f_back = cat(1, f_back, forces{i});
end;

d = d_fwd;
f = f_fwd;
%plot total force v displacement for forward runs
h_force_vs_displacement_corrected = figure;
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

f_fwd_binned = force_binned;
d_fwd_binned = trap_dist_binned;
e_fwd_binned = force_error;
%-------------------------------

d = d_back;
f = f_back;
%plot total force v displacement for back runs
h_force_vs_displacement_corrected = figure;
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
f_back_binned = force_binned;
d_back_binned = trap_dist_binned;
e_back_binned = force_error;

h_plot_together = figure;
plot(d_fwd_binned, f_fwd_binned,'.');
hold on;
plot(d_back_binned, f_back_binned,'.');
ylabel(['Force (' units_xy ')']);
xlabel('Trap separation (nm)');
title('Force - distance curve');
axis tight;

%% sort binned vars so they dots can be connected.
m = [d_fwd_binned; f_fwd_binned; e_fwd_binned];
%m = [trap_dist_binned; force_binned];
m = m';
ms = sortrows(m);
d_fwd_binned_sorted = ms(:,1);
f_fwd_binned_sorted = ms(:,2);
e_fwd_binned_sorted = ms(:,3);

m = [d_back_binned; f_back_binned; e_back_binned];
m = m';
ms = sortrows(m);
d_back_binned_sorted = ms(:,1);
f_back_binned_sorted = ms(:,2);
e_back_binned_sorted = ms(:,3);

%plot them
h_plot_together = figure;
errorbar(d_fwd_binned_sorted, f_fwd_binned_sorted,e_fwd_binned_sorted,'o-');
hold on;
errorbar(d_back_binned_sorted, f_back_binned_sorted, e_back_binned_sorted,'o-');
ylabel(['Force (' units_xy ')']);
xlabel('Trap separation (nm)');
title('Force - distance curve');
axis tight;
legend('Forward: buckling and fraying','Back: re-bundling and extension');


%%
filepath;
filename = [filepath '\force-displacement intervals.mat'];
if (exist('f_fwd_binned','var') == 1)
    %check if fowrward and back runs were separated
    save(filename,...
    'f_fwd_binned','f_back_binned','d_fwd_binned','d_back_binned',...
    'intervals', 'relevant', 'distances', 'forces','forces_x','forces_y');
else
    save(filename,...
    'intervals', 'relevant', 'distances', 'forces','forces_x','forces_y');
end;
