%% Let's begin
clc;
clear all;
close all;

%% Set-up file path and names
path ='D:\Data - MT Sliding and Friction\2016\04-26-2016\';
findcenter_path = 'findcenter 3\';
findcenter_path = [path findcenter_path];
findcenter_scan = 1;
data_path = 'force 0\';
acquisition_fname = 'acquisition.csv';
cal_fname = ['cal 1 trap' '.dat'];
cal_averages = 2^5;
number_of_traps_calibration = 1;
%calibrated stiffness is divided by the number of traps in a run
%We assume that the calibration is done with a single trap
number_of_traps_run = 2;

data_fname = 'data.csv';
img_fname = 'images.tif';

plot_trap_positions = true;
plot_findcenter_scans = true;
analyze_calibration = true;

% number of points to average to get zero-force
number_of_zero_force_points = 100;
% location of zero-force points
% accepted values are 'start', 'end' and 'custom'
zero_force_location = 'custom';
%sero force can also be specified as an interval
%this is only relevant if location is 'custom'
%otherwise, variable number_of_zero_force_points will determine the
%interval from either start or end of the subset.
zero_force_interval = 650:1:750;

analyze_subset = false;
%indicate subset indexes here
subset_start = 1;
subset_end = 3599;

overlay_scale = 20; %pre-define overlay scale; it will be calculated later
laser_power = 0; %pre-define laser power variable

%% read acquisition parameters
fname = [path data_path acquisition_fname];
try [Fsampling,Nsamples,Velocity,MaxMove,...
        N_traps,T_delay,converted,Laser_Power,AODx,AODy,QPDx_acq,QPDy_acq] = read_acq_data(fname);
catch error
    %error is encountered while reading the acquisition data
    %this usually happens when the file is not found
    %define default AOD  and QPD sensitivities instead
    disp('Error encountered while trying to get acq. options!');
    disp(error.message);
    disp('Using default AOD conversion from MHz to nm');
    AODx = -2668; %default AOD conversion of nm/MHz
    AODy = -2660; %default AOD conversion of nm/MHz
    
    QPDx_acq = NaN;
    QPDy_acq = NaN;
end;
%% read the acquired data
filename = [path data_path data_fname];
trap_y2 = [];
trap_x2 = [];

[trap_x1, trap_y1, pos_x, pos_y, t, z, trap_x2,trap_y2] = read_cv_data_twotraps(filename);
%uncomment this for compatibility with old data format
%[trap_x1, trap_y1, pos_x, pos_y, t, z] = read_cv_data(filename);

%get the total length of the data
total_data_length = length(trap_x1);
%if we need to analyze the whole run, define start and end points as such:
if (~analyze_subset)
    subset_start = 1;
    subset_end = total_data_length;
    disp('Analyzing the entire run');
end;
subset = subset_start:1:subset_end;

fprintf('Analyzing frames: %i to %i\n',subset_start, subset_end);

switch zero_force_location
    case 'end'
        %zero force points are in the end of the subset
        zero_force_points = subset_end-number_of_zero_force_points:1:subset_end;
    case 'start'
        %zero force points are in the beginning of the subset
        zero_force_points = subset_start:1:subset_start + number_of_zero_force_points-1;
    case 'custom'
        %zero force region is custom
        zero_force_points = zero_force_interval;
    otherwise
        warning('Unrecognized value for zero_force_location parameter!');
        zero_force_points = 1;
end;

%% Display information about the run
if (~converted)
    disp('Conversion NOT USED!');
else
    disp('Coversion WAS USED!');
end;

%delta_t has all time-steps in ms
delta_t = t(2:end) - t(1:end-1);
%total_run_time is the total acquisition time in ms
total_run_time = t(end) - t(1);
fprintf('\n');
fprintf('Total run time: %2.2f s\n', total_run_time/1000);
fprintf('Expected timestep: %2.2f ms\n', Nsamples/Fsampling * 1000);
fprintf('Actual mean timestep: %.2f ms\n', mean(delta_t));
%
fprintf('\nDelay, ms: %2.2f \n', T_delay);
%T_delay is not relevant for force measurements;
%
fprintf('Sampling Freq, Hz: %2.2f \n', Fsampling);
fprintf('# of samples in one data-point: %i \n', Nsamples);
fprintf('Total # of traps: %i \n', N_traps);
if (number_of_traps_run ~= N_traps)
    warning('Specified number of traps is inconsistent!');
end;

fprintf('Traps moving: %i \n', MaxMove+1);
%no traps are moving if this manual force measurement code
%ignore this parameter

fprintf('Laser Power, W: %2.2f \n', Laser_Power);

%% display findcenter results
[x_c, y_c, QPD_x_fit, QPD_y_fit, max_x_distance, max_y_distance, QPD_map_x, QPD_map_y] =...
    get_findcenter_data(findcenter_path, plot_findcenter_scans, findcenter_scan, AODx, AODy);
if ((QPD_x_fit ~= 0) && (QPD_y_fit ~= 0))
    %if get_findcenter_data returned valid fits for QPD_x and QPD_y
   
    %QPD_x_fit and QPD_y_fit are in nm/V;
  
    % it is possible to check now that QPD_x_fit
    % and QPDx (which is a fit done by
    % LabView) are close. If there is a big
    % discrepancy, it is a sign of error.
    diff_x = abs(QPDx_acq - QPD_x_fit)/min(abs(QPDx_acq),abs(QPD_x_fit));
    
    % error threshold
    diff_threshold = 0.1;
    
    if (diff_x > diff_threshold)
        disp('X QPD sensitivity error!');
    end;
    diff_y = abs(QPDy_acq - QPD_y_fit)/min(abs(QPDy_acq),abs(QPD_y_fit));
    
    if (diff_y > diff_threshold)
        disp('Y QPD sensitivity error!');
    end;
    
    disp('Defining QPD sensitivity from the fits of findcenter!');
    QPDx = QPD_x_fit;
    QPDy = QPD_y_fit;
else
    %if get_findcenter_data did not return valid fits, use the linear
    %conversion coefficient defined in the acquisition step.
    disp('Defining QPD sensitivity from the force acquisition data!');
    QPDx = QPDx_acq;
    QPDy = QPDy_acq;
end;
fprintf('QPD sensitivity X-channel: %1.2f nm/V\n', QPDx);
fprintf('QPD sensitivity Y-channel: %1.2f nm/V\n', QPDy);

%% analyze calibration file
if (analyze_calibration)
    min_points_to_average = 10; %parameter for analyze_tweezer_calibration
    high_frequency_cutoff = 10000; %parameter for analyze_tweezer_calibration
    low_frequency_cutoff = 50;
    [fc_x,fc_y, diff_x, diff_y] =...
        analyze_tweezer_calibration(findcenter_path, cal_fname, cal_averages,...
        min_points_to_average, high_frequency_cutoff, low_frequency_cutoff);
    
    D_x = diff_x*QPDx^2; %Dx in nm^2/s
    D_y = diff_y*QPDy^2; %Dy in nm^2/s
    D_x_um = D_x*10^(-6); %Dx in micron^2/s
    D_y_um = D_y*10^(-6); %Dy in micorn^2/s
    
    fprintf('Diffusion coeff in x: %1.2f micron^2/s\n', D_x_um);
    fprintf('Diffusion ceoff in y: %1.2f micron^2/s\n', D_y_um);
    
    kBT = 4; %pN*nm
    k_x = fc_x*2*pi*kBT/D_x;
    k_y = fc_y*2*pi*kBT/D_y;
    fprintf('Trap stiffness k_x: %1.4f pN/nm\n', k_x);
    fprintf('Trap stiffness k_y: %1.4f pN/nm\n', k_y);
    
    f_max_x = 0.5 * k_x * max_x_distance/(number_of_traps_run/number_of_traps_calibration);
    f_max_y = 0.5 * k_y * max_y_distance/(number_of_traps_run/number_of_traps_calibration);
    fprintf('Max reliable force in x: %1.4f pN\n', f_max_x);
    fprintf('Max reliable force in y: %1.4f pN\n', f_max_y);
    disp('-------------------------------------------------');
    
else
    %using predetermined values of k_x and k_y
    %the values are hard-coded for now
    %see the xls file with historical trap stiffness values
    %this data can be acquired directly from the Excel file in future...
    disp('! Using pre-determined trap stiffness !');
    %laser_power = input('Enter laser power (W): ');
    laser_power = Laser_Power;
    if (Laser_Power ~= laser_power) && (~analyze_calibration)
        disp('Previously specified laser power is inconsistent!');
    end;
    k_x = 0.0517*laser_power; %k_x in pN/nm;
    k_y = 0.1097*laser_power; %k_y on pN/nm;
    
    fprintf('Trap stiffness k_x: %1.4f pN/nm\n', k_x);
    fprintf('Trap stiffness k_y: %1.4f pN/nm\n', k_y);
    number_of_traps_calibration = 1; %explicitly state that when using
    %pre-recorded stiffness is based on a single trap
    
    %factor of 0.5 is because max_x_distance is the total length
    %of the linear response region. Detection trap is placed in the center
    %of this region, thus limiting bead displacement to 1/2 of the linear
    %region in each direction.
    f_max_x = 0.5* k_x * max_x_distance/(number_of_traps_run/number_of_traps_calibration);
    f_max_y = 0.5* k_y * max_y_distance/(number_of_traps_run/number_of_traps_calibration);
    fprintf('Max reliable force in x: %1.4f pN\n', f_max_x);
    fprintf('Max reliable force in y: %1.4f pN\n', f_max_y);
    disp('-------------------------------------------------');
end;


%% plot trap positions

%parameter for plotting smoothed force vs time/frame
force_smoothing = 5;

%define color for raw data
color_raw = [0.831  0.816    0.784];
%define color for smoothed data
color_smooth = [ 0.502  0.502   0.502];
%define color for block-averaged data
color_average = [0  0.447   0.741];

%determine if trap positions were recorded
traps_recorded = true;
if (numel(trap_y1) == 0) || (numel(trap_y2) == 0) || (numel(trap_x1) == 0)
    traps_recorded = false;
    warning('Two traps positions were not recorded!');
end;

if (traps_recorded)
    %convert trap positions from MHz to nm or pixels
    %get relative trap displacement
    trap_y1_MHz = trap_y1 - trap_y1(subset_start);
    trap_x1_MHz = trap_x1 - trap_x1(subset_start);
    trap_x2_MHz = trap_x2 - trap_x1(subset_start);
    trap_y2_MHz = trap_y2 - trap_y1(subset_start);
    
    %convert MHz to nm
    trap_y1_nm = trap_y1_MHz*AODy;
    trap_y2_nm = trap_y2_MHz*AODy;
    trap_x1_nm = trap_x1_MHz*AODx;
    trap_x2_nm = trap_x2_MHz*AODx;
    
    %convert nm to pixels
    %this will be used when working with graphic overlays
    camera_pixel_size = 130; %nm/pixels for Andor IKon-M on the AOD
    trap_y1_pix = trap_y1_nm/camera_pixel_size;
    trap_y2_pix = trap_y2_nm/camera_pixel_size;
    trap_x1_pix = trap_x1_nm/camera_pixel_size;
    trap_x2_pix = trap_x2_nm/camera_pixel_size;
    
    %create trap 1,2 variables for the overlay
    trap1 = [trap_x1_pix trap_y1_pix];
    trap2 = [trap_x2_pix trap_y2_pix];
    
    %calculate trap separation in nm
    trap_dist = sqrt((trap_y2_nm - trap_y1_nm).^2 + (trap_x2_nm - trap_x1_nm).^2);
    trap_separation_x = trap_x2_nm - trap_x1_nm;
    trap_dist_x = sqrt(trap_separation_x.^2);
    trap_separation_y = trap_y2_nm - trap_y1_nm;
    trap_dist_y = sqrt(trap_separation_y.^2);
    
    %rotated trap coordinates - for consistency check
    %x1,y1 should be 0 still;
    %x2 should become 0; y2 should be equal total trap distance.
    trap_x1_nm_r = trap_x1_nm;
    trap_x2_nm_r = trap_x2_nm;
    trap_y1_nm_r = trap_y1_nm;
    trap_y2_nm_r = trap_y2_nm;
    
    %calculate the angle between line connecting the traps and the
    %y-axis
    trap_line_angle_rad = atan2(trap_separation_y, trap_separation_x);
    trap_line_angle_deg = atan2d(trap_separation_y, trap_separation_x);
    for i = 1:1:length(trap_line_angle_deg)
        alpha = 90 - trap_line_angle_deg(i);
        trap_line_rotation(:,:,i) = [cosd(alpha), -sind(alpha);...
                                     sind(alpha),  cosd(alpha)];
%         rotated_trap_1 = trap_line_rotation(:,:,i) * [trap_x1_nm(i); trap_y1_nm(i)];
%         rotated_trap_2 = trap_line_rotation(:,:,i) * [trap_x2_nm(i); trap_y2_nm(i)];
%         trap_x1_nm_r(i) = rotated_trap_1(1); 
%         trap_y1_nm_r(i) = rotated_trap_1(2);
%         trap_x2_nm_r(i) = rotated_trap_2(1);
%         trap_y2_nm_r(i) = rotated_trap_2(2);
    end;
    
    %all the same quantities in the rotated coordinate frame as consistency
    %check
%     trap_dist_r = sqrt((trap_y2_nm_r - trap_y1_nm_r).^2 + (trap_x2_nm_r - trap_x1_nm_r).^2);
%     trap_separation_x_r = trap_x2_nm_r - trap_x1_nm_r;
%     trap_dist_x_r = sqrt(trap_separation_x_r.^2);
%     trap_separation_y_r = trap_y2_nm_r - trap_y1_nm_r;
%     trap_dist_y_r = sqrt(trap_separation_y_r.^2);
%     trap_line_angle_rad_r = atan2(trap_separation_y_r, trap_separation_x_r);
%     trap_line_angle_deg_r = atan2d(trap_separation_y_r, trap_separation_x_r);
    
    if plot_trap_positions
        %     h_fig_trap_traj = figure;
        %     plot(trap_y1_pix,trap_x1_pix, 'ob');
        %     hold on;
        %     plot( trap_y2_pix,trap_x2_pix, 'or');
        %   % %This plots trap positions vs time
        %     h_fig_trap_y1 = figure;
        %     plot(time_s(subset_start:subset_end), trap_y1_pix(subset_start:subset_end), '.b');
        %     title('Moving trap1 position, y axis');
        %
        %     h_fig_trap_x1 = figure;
        %     plot(time_s(subset_start:subset_end), trap_x1_pix(subset_start:subset_end), '.r');
        %     title('Moving trap1 position, x axis');
        %
        %     h_fig_trap_y2 = figure;
        %     plot(time_s(subset_start:subset_end), trap_y2_pix(subset_start:subset_end), '.b');
        %     title('Moving trap 2 position, y axis');
        %
        %     h_fig_trap_x2 = figure;
        %     plot(time_s(subset_start:subset_end), trap_x2_pix(subset_start:subset_end), '.r');
        %     title('Moving trap 2 position, x axis');
        
        %Plot trap positions vs frame index
        h_fig_trap_pos = figure;
        plotyy([subset', subset'],[trap_y1_nm(subset),trap_x1_nm(subset)],[subset', subset'],[trap_y2_nm(subset),trap_x2_nm(subset)]);
        title('Trap positions (x + y axes)');
        xlabel('Frame #');
        ylabel('Trap separation (nm)');
        legend('y_1','x_1','y_2','x_2');
        
        h_fig_trap_line_angle = figure;
        plot(subset', trap_line_angle_deg(subset));
        title('Trap line angle');
        xlabel('Frame #');
        ylabel('Angle (degrees)');
        
        h_fig_trap_dist = figure;
        plot([subset',subset'],[trap_dist_x(subset),trap_dist_y(subset)],subset',trap_dist(subset));
        title('Trap separation during run');
        ylabel('Trap separation (nm)');
        xlabel('Frame #');
        legend('Distance x', 'Distance y', 'Distance total');
        
        %plot rotated coords for consistency check
%         h_fig_trap_pos_r = figure;
%         plotyy([subset', subset'],[trap_y1_nm_r(subset),trap_x1_nm_r(subset)],[subset', subset'],[trap_y2_nm_r(subset),trap_x2_nm_r(subset)]);
%         title('Trap positions (x + y axes) ROTATED');
%         xlabel('Frame #');
%         ylabel('Trap separation (nm)');
%         legend('y_1','x_1','y_2','x_2');
%         
%         h_fig_trap_line_angle_r = figure;
%         plot(subset', trap_line_angle_deg_r(subset));
%         title('Trap line angle ROTATED');
%         xlabel('Frame #');
%         ylabel('Angle (degrees)');
%         
%         h_fig_trap_dist = figure;
%         plot([subset',subset'],[trap_dist_x_r(subset),trap_dist_y_r(subset)],subset',trap_dist_r(subset));
%         title('Trap separation during run ROTATED');
%         ylabel('Trap separation (nm)');
%         xlabel('Frame #');
%         legend('Distance x', 'Distance y', 'Distance total');
%         
    else
        disp('Trap positions not recorded');
    end;
end;

%% plot x and y bead positions in the detection beam
% if conversion used, the position is in nm,
% otherwise it is in volts;
if (converted)
    
    force_x = k_x.*pos_x;
    force_y = k_y.*pos_y;
    
    %or
    %revert back the conversion, i.e. go back to volts
    pos_x = pos_x ./ QPDx_acq;
    pos_y = pos_y ./ QPDy_acq;
    disp('conversion reverted...')
    
end;

%now all readings are in volts
%convert them into nm with the appropraite function that maps pos_x in
%V into pos_x in nm:
% x_nm = f(x_V) and y_nm = f(y_V)
%and then into force through trap_stiffness k

%the simplest possible mapping is linear:
force_x_lin = pos_x.*QPDx.*k_x/(number_of_traps_run/number_of_traps_calibration);
force_y_lin = pos_y.*QPDy.*k_y/(number_of_traps_run/number_of_traps_calibration);

%the more advanced map is from the inverse mapping function
force_x = QPD_map_x(pos_x).*k_x/(number_of_traps_run/number_of_traps_calibration);
force_y = QPD_map_y(pos_y).*k_y/(number_of_traps_run/number_of_traps_calibration);

units_xy = 'pN';
% convert time into s
time_s = t*10^-3;

%calculate zero-force in x and y directions - constant offset due to
%incomplete centering of the bead in the detection beam in the beginning of
%the experiements.
%substract zero force from the force readings (i.e shift the force by a
%constant)
f_x_0 = mean(force_x(zero_force_points));
f_y_0 = mean(force_y(zero_force_points));
force_x = force_x - f_x_0;
force_y = force_y - f_y_0;
total_force = sqrt(force_x.^2 + force_y.^2);
theta = atan2d(force_y,force_x);

force_x_r = force_x;
force_y_r = force_y;

%rotate to the trap coordinate frame
for i = 1:1:length(force_x)
 a = trap_line_rotation(:,:,i)*[force_x(i); force_y(i)];
 force_x_r(i) = a(1);
 force_y_r(i) = a(2);
end;
total_force_r = sqrt(force_x_r.^2 + force_y_r.^2);
theta_r = atan2d(force_y_r,force_x_r);

%plot F_x vs time
% h_fig_x = figure;
% plot(time_s(subset), force_x(subset), '.-');
% ylabel(['Force (' units_xy ')']);
% xlabel('Time (s)');
% title('Measured force in x-direction');
% axis tight;

%plot F_x versus frame #
% h_fig_x_frames = figure;
% plot(force_x(subset), '.-');
% ylabel(['Force (' units_xy ')']);
% xlabel('Frame #');
% title('Measured force in x-direction');
% axis tight;

%plot F_y vs time
% h_fig_y = figure;
% plot(time_s(subset),force_y(subset),'.-');
% ylabel(['Force (' units_xy ')']);
% xlabel('Time (s)');
% title('Measured force in y-direction');
% axis tight;

%plot F_y versus frame #
% h_fig_y_frames = figure;
% plot(force_y(subset),'.-');
% ylabel(['Force (' units_xy ')']);
% xlabel('Frame #');
% title('Measured force in y-direction');
% axis tight;

%plot z versus time
% h_fig_z = figure;
% plot(time_s(subset),z(subset),'.-');
% ylabel('Z-position, \mum');
% xlabel('Time (s)');
% axis tight;
% h_fig_z_frames = figure;
% plot(z(subset),'.-');
% ylabel('Z-position, \mum');
% xlabel('Frame #');
% axis tight;

%plot F_x, F_y and z vs time in the same plot using independent y-axes
%commented out on 06/07/2016 as unnecessary plot
% h_fig_all = figure;
% [h_axes, h_plot_f, h_plot_z] = plotyy([time_s(subset),time_s(subset)],...
%     [force_x(subset),force_y(subset)],time_s(subset), z(subset));
% ylabel(h_axes(2), 'Z-position (\mum)');
% ylabel(h_axes(1), 'Force (pN)');
% xlabel('Time (s)');
% title('X-force, Y-force and z-position together');
% %axis tight;

%plot F_x, F_y and z vs frame index in the same plot using independent y-axes
h_fig_all_frames = figure;
[h_axes, h_plot_f, h_plot_z] = plotyy([subset',subset'],...
    [force_x(subset),force_y(subset)],subset, z(subset));
ylabel(h_axes(2), 'Z-position (\mum)');
ylabel(h_axes(1), 'Force (pN)');
xlabel('Frame #');
%axis tight;
title('X-force, Y-force and z-position together');
legend('F_x','F_y','z');

h_fig_all_frames = figure;
[h_axes, h_plot_f, h_plot_z] = plotyy([subset',subset'],...
    [force_x_r(subset),force_y_r(subset)],subset, z(subset));
ylabel(h_axes(2), 'Z-position (\mum)');
ylabel(h_axes(1), 'Force (pN)');
xlabel('Frame #');
%axis tight;
title('X-force, Y-force and z-position ROTATED');
legend('F_x','F_y','z');


%% plot total force vs time

%plot total F vs time
%commented as deemed not essential
%
% h_force = figure;
% plot(time_s(subset),total_force(subset),'.','Color',color_raw);
% hold on;
% plot(time_s(subset),smooth(total_force(subset),force_smoothing),'-','LineWidth',2,'Color',color_smooth);
% ylabel(['Force (' units_xy ')']);
% xlabel('Time (s)');
% title('Total measured force');
% axis tight;
%
h_force_frames = figure;
plot(total_force(subset),'.','Color',color_raw);
hold on;
plot(smooth(total_force(subset),force_smoothing),'-','LineWidth',2,'Color',color_smooth);
ylabel(['Force (' units_xy ')']);
xlabel('Substack frame #');
title('Total measured force');
axis tight;

h_force_frames_r = figure;
plot(total_force_r(subset),'.','Color',color_raw);
hold on;
plot(smooth(total_force_r(subset),force_smoothing),'-','LineWidth',2,'Color',color_smooth);
ylabel(['Force (' units_xy ')']);
xlabel('Substack frame #');
title('Total measured force ROTATED');
axis tight;


%plot the angle of the force vector
h_angle = figure;
%plot(time_s(subset_start:subset_end),theta(subset_start:subset_end),'.r');
plot(theta(subset),'.');
hold on;
%add pi/2 and -pi/2 lines to the plot for reference
%plot(time_s(subset_start:subset_end),pi/2*ones(1,length(theta(subset_start:subset_end))),'-k');
%plot(time_s(subset_start:subset_end),-pi/2*ones(1,length(theta(subset_start:subset_end))),'-k');
ylabel('Angle (degrees)');
xlabel('Substack frame #');
%xlabel('Time (s)');
title('Angle of the force vector');
axis tight;

h_angle_rotated = figure;
%plot(time_s(subset_start:subset_end),theta(subset_start:subset_end),'.r');
plot(theta_r(subset),'.');
hold on;
%add pi/2 and -pi/2 lines to the plot for reference
%plot(time_s(subset_start:subset_end),pi/2*ones(1,length(theta(subset_start:subset_end))),'-k');
%plot(time_s(subset_start:subset_end),-pi/2*ones(1,length(theta(subset_start:subset_end))),'-k');
ylabel('Angle (degrees)');
xlabel('Substack frame #');
%xlabel('Time (s)');
title('Angle of the force vector ROTATED');
axis tight;

%so far we have been dealing vith vector quantities (such as force and trap
%separation) in the laboratory frame of reference connected to the camera
%pixels. 
%we need to transform these into more relevant frame of reference with one
%vector along the direction connecting the two traps and one perpendicualr
%to it. In this frame, one force is always the buckling force, and the
%other is the tangential force. 
%direction of the line that connects two traps is given by the angle
%trap_line_angle_deg or trap_line_angle_rad 
%this is the angle between the line form trap 1 to trap 2 and the x-axis


%% plot force vs displacement curve

%in the lab reference frame
total_force_subset = total_force(subset);
trap_dist_subset = trap_dist(subset);
force_y_subset = force_y(subset);
force_x_subset = force_x(subset);

%in the trap line reference frame
total_force_subset_r = total_force_r(subset);
trap_dist_subset_r = trap_dist(subset); %trap distance is the same
force_y_subset_r = force_y_r(subset);
force_x_subset_r = force_x_r(subset);

%plot a figure for user input of the force-displacement interval
h_force_vs_displacement = figure;
plot(trap_dist_subset, total_force_subset,'.','Color',color_raw);
hold on;
plot(trap_dist_subset, smooth(total_force_subset,force_smoothing),'.','LineWidth',2,'Color',color_smooth);
ylabel(['Force (' units_xy ')']);
xlabel('Trap separation (nm)');
title('Force - distance curve raw');
axis tight;
[limit_disp, limit_force] = ginput(2);
trap_dist_fd_ind = ...
    find(trap_dist_subset >= limit_disp(1) & trap_dist_subset <= limit_disp(2));

%variables in the lab reference frame
trap_dist_fd_curve = trap_dist_subset(trap_dist_fd_ind);
total_force_fd_curve = total_force_subset(trap_dist_fd_ind);
force_y_fd_curve = force_y_subset(trap_dist_fd_ind);
force_x_fd_curve = force_x_subset(trap_dist_fd_ind);

%variables in the trap line reference frame
trap_dist_fd_curve_r = trap_dist_subset_r(trap_dist_fd_ind);
total_force_fd_curve_r = total_force_subset_r(trap_dist_fd_ind);
force_y_fd_curve_r = force_y_subset_r(trap_dist_fd_ind);
force_x_fd_curve_r = force_x_subset_r(trap_dist_fd_ind);

%plot total force v displacement
h_force_vs_displacement_corrected = figure;
plot(trap_dist_fd_curve, total_force_fd_curve,'.','Color',color_raw);
hold on;
plot(trap_dist_fd_curve, smooth(total_force_fd_curve,force_smoothing),'.','LineWidth',2,'Color',color_smooth);
ylabel(['Force (' units_xy ')']);
xlabel('Trap separation (nm)');
title('Force - distance curve');
axis tight;

%bin force based on the fixed displacement values and plot that
[trap_dist_binned, force_binned, force_error] = binned_xy(trap_dist_fd_curve, total_force_fd_curve);
errorbar(trap_dist_binned, force_binned, force_error,'o','Color',color_average,'LineWidth',1);
%also bin in the rotated frame
[trap_dist_binned_r, force_binned_r, force_error_r] = binned_xy(trap_dist_fd_curve_r, total_force_fd_curve_r);

%plot y-force v displacement
h_y_force_vs_displacement_corrected = figure;
plot(trap_dist_fd_curve, force_y_fd_curve,'.','Color',color_raw);
hold on;
plot(trap_dist_fd_curve, smooth(force_y_fd_curve,force_smoothing),'.','LineWidth',2,'Color',color_smooth);
plot(trap_dist_fd_curve_r, smooth(force_y_fd_curve_r,force_smoothing),'.','LineWidth',2,'Color','blue');
ylabel(['Force (' units_xy ')']);
xlabel('Trap separation (nm)');
title('Y-axis Force - distance curve');
axis tight;
legend('raw F_y', 'smoothed F_y', 'rotated F_y');

%plot x-force v displacement
h_x_force_vs_displacement_corrected = figure;
plot(trap_dist_fd_curve, force_x_fd_curve,'.','Color',color_raw);
hold on;
plot(trap_dist_fd_curve, smooth(force_x_fd_curve,force_smoothing),'.','LineWidth',2,'Color',color_smooth);
plot(trap_dist_fd_curve_r, smooth(force_x_fd_curve_r,force_smoothing),'.','LineWidth',2,'Color','blue');
ylabel(['Force (' units_xy ')']);
xlabel('Trap separation (nm)');
title('X-axis Force - distance curve');
axis tight;
legend('raw F_x', 'smoothed F_x', 'rotated F_x');


%%
%set the force measurement point center point on the image
%do if for the first image until the user is satidfied with the location of
%the force vector and traps.
user_satisfied = false;


while (~user_satisfied)
    fname = [path data_path img_fname];
    %check if the images exist
    file_attributes = dir(fname);
    images_saved = ~isempty(file_attributes);
    if (images_saved)
        disp(['Working with video file: ' img_fname]);
        img_info = imfinfo(fname);
        total_frames = length(img_info);
        image_index = subset_start + 50;
        img = imread(fname, image_index);
        [m,n,l] = size(img); %m is y-axis, n is x-axis
        img_dimensions = min(size(img));
        
        %calculate overlay scale
        overlay_scale_x = (n/3)/max(abs(force_x(subset)));
        overlay_scale_y = (m/3)/max(abs(force_y(subset)));
        overlay_scale = min(overlay_scale_x, overlay_scale_y);
        
        disp(['Calculated overlay scale: ' num2str(overlay_scale) ' pix/pN']);
        
        
        x_0 = round(n/2);
        y_0 = round(m/2);
        if (traps_recorded)
            [overlay_img, vector_img] = overlay_force_vector_with_traps(img,x_0,y_0,force_x(image_index),force_y(image_index), overlay_scale, [1,1],trap1(image_index,:),trap2(image_index,:),10,true);
        else
            overlay_img = overlay_force_vector(img,x_0,y_0,force_x(image_index),force_y(image_index), overlay_scale, [1,1]);
        end;
        h_fig_overlay = figure;
        imshow(overlay_img);
        title('Pick correct location of the measurement bead');
        [x_0, y_0] = ginput(1);
        
        %trap1(1,:) = trap1(1,:)+ x_0;
        
        %create new overlay image
        if (traps_recorded)
            [overlay_img, vector_img] = overlay_force_vector_with_traps(img,x_0,y_0,force_x(image_index),force_y(image_index), overlay_scale, [1,1],trap1(image_index,:),trap2(image_index,:),10,true);
        else
            overlay_img = overlay_force_vector(img,x_0,y_0,force_x(image_index),force_y(image_index), overlay_scale, [1,1]);
        end;
        figure(h_fig_overlay);
        imshow(overlay_img);
        title('New location of the measurement bead');
        user_satisfied = true;
        %dbutton = questdlg('qstring','title')
    else
        disp(['Images not found: ' img_fname]);
        user_satisfied = true;
    end;
end;
%%
if (images_saved)
    %determine the size of overlay
    [m,n,l] = size(overlay_img);
    %pre-make frame array
    total_frames = length(img_info);
    total_frames = -(subset_start - subset_end + 1);
    all_frames = zeros(m,n,l,total_frames,'uint8');
    vector_frames = zeros(m,n,l,total_frames,'uint8'); %store vector images separately
    
    
    for image_index = subset_start:subset_end
        disp(['Exporting frame: ', num2str(image_index)]);
        img = imread(fname, image_index);
        if (traps_recorded)
            [overlay_img, vector_img] = ...%overlay_force_vector(img,x_0,y_0,force_x(image_index),force_y(image_index), overlay_scale, [1,1]);
                overlay_force_vector_with_traps(img,x_0,y_0,force_x(image_index),force_y(image_index), overlay_scale, [1,1],trap1(image_index,:),trap2(image_index,:),10,true);
        else
            overlay_img = overlay_force_vector(img,x_0,y_0,force_x(image_index),force_y(image_index), overlay_scale, [1,1]);
            %vector_img = overlay_img;
        end;
        all_frames(:,:,:,image_index - subset_start + 1) = overlay_img;
        if (traps_recorded)
            vector_frames(:,:,:,image_index - subset_start + 1) = vector_img;
        end;
    end;
end;

%%
%write video
if analyze_subset
    analysis_folder = ['analysis_' num2str(subset_start) '-' num2str(subset_end) '\'];
else
    analysis_folder = 'analysis\';
end;
analysis_path = [path data_path analysis_folder];
mkdir(analysis_path);
%save overlay images

if (images_saved)
    
    
    video_fname = [analysis_path img_fname '.uncompressed'];
    video_writer = VideoWriter(video_fname,'Uncompressed AVI');
    open(video_writer);
    disp('Writing video...');
    writeVideo(video_writer,all_frames);
    close(video_writer);
    msgbox('Video saved!','FilamentTracker');
    %save vector images
    video_fname = [analysis_path img_fname '.vectors.uncompressed'];
    video_writer = VideoWriter(video_fname,'Uncompressed AVI');
    open(video_writer);
    disp('Writing vectors video...');
    writeVideo(video_writer,vector_frames);
    close(video_writer);
    msgbox('Video saved!','FilamentTracker');
end;

%% save figures and all variables

% fx_filename = [analysis_path 'force_x'];
% fy_filename = [analysis_path 'force_y'];
% savefig(h_fig_x_frames,fx_filename);
% savefig(h_fig_y_frames,fy_filename);

f_filename = [analysis_path 'fx_fy_fz'];
savefig(h_fig_all_frames,f_filename)

f_filename = [analysis_path 'force'];
savefig(h_force_frames,f_filename);

% f_filename = [analysis_path 'z-position'];
% savefig(h_fig_z_frames,f_filename);

f_filename = [analysis_path 'force-displacement'];
savefig(h_force_vs_displacement_corrected,f_filename);

angle_filename = [analysis_path 'angle'];
savefig(h_angle,angle_filename);
disp('Figures saved!');

%save all variables and open figures
disp(['Saving data ...']);
var_fname = [analysis_path 'raw_data'];
save(var_fname);
disp(['Data saved in ' var_fname ' folder!']);

%save specific variables for the force-displacement curve:
%trap_dist_fd_curve
%total_force_fd_curve
%force_y_fd_curve
%force_x_fd_curve
%trap_dist_binned
%force_binned 
disp('Saving force-distance variables');
var_fname = [analysis_path 'forces-displacement_vars'];
save(var_fname,'trap_dist_fd_curve','total_force_fd_curve','force_y_fd_curve',...
    'force_x_fd_curve','trap_dist_binned','force_binned');
disp(['Data saved in ' var_fname ' folder!']);

%save variables specific to the force-displacement curve in the rotated
%frame that corresponds to line of traps
%trap_dist_fd_curve_r
%total_force_fd_curve_r
%force_y_fd_curve_r
%force_x_fd_curve_r
%trap_dist_binned_r
%force_binned_r
disp('Saving force-distance variables in the bundle reference frame');
var_fname = [analysis_path 'forces-displacement_vars_bundle_frame.mat'];
save(var_fname,'trap_dist_fd_curve_r','total_force_fd_curve_r','force_y_fd_curve_r',...
    'force_x_fd_curve_r','trap_dist_binned_r','force_binned_r');
disp(['Data saved in ' var_fname ' folder!']);