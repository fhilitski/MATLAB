function [ time, msd, error, x, t] = process_filament_data( tracking_data, sensitivity, t_scale, pixel_size)
%PROCESS_FILAMENT_DATA plots and saves bundle sliding analysis
%
%process_filament_data(tracking_data, sensitivity, t_scale, pixel_size)
%
%   tracking_data - Raw data from bundle sliding tracking has form 
%                  [frame#, Position , Total Length]
%                    where frame# is tracked frame #;
%                    position is midpoint of the bundle position along filament's
%                    contour;
%                    Total length is total contourlength (of the whole filament);
%   sensitivity   - cutoff in filament length deviation
%   t_scale       - time between frames
%   pixel_size    - pixel size in microns
% WARNING: t_scale and pixel size are not yet implemented; they are both
% effectively = 1; be careful with the implementation, as the msd
% calculation code
%   
%   
%   All measurements of tracking_data are in pixels and frames, and need to be rescaled with
%   appropriate pixel-to-um ratio and frame capture rate;

clc;


t = tracking_data(:,1);
total_length = tracking_data(:,3);
x = tracking_data(:,2);
x = x - x(1);
N = length(tracking_data);

%   It sometimes happens that skeletonized image does not capture the full
%   filament. As a result, Total Length varies. It is ok to have some
%   variability of total length due to noise in CCD, filament going out of
%   focus, etc. But frames with significant changes of length (say, >1*STD)
%   should be dropped. 
%   Following part accomplishes this noble goal
%   Sensitivity is the coefficient in front of STD of the filament, and it
%   determines how much of the length variability is acceptable.


figure(100);
plot(t,total_length,'ob');

mean_length = mean(total_length);
std_length = std(total_length);
std_mean = std_length/sqrt(N);

min_length = mean_length - sensitivity * std_length;
max_length = mean_length + sensitivity * std_length;

hold on;
plot(t,x,'or');
plot(t,mean_length*ones(N,1),'-r','LineWidth',2);
plot(t,(min_length)*ones(N,1),'-b','LineWidth',1);
plot(t,(max_length)*ones(N,1),'-b','LineWidth',1);
hold off;
title('Raw data from filament tracking');
xlabel('Time, s');
ylabel('Position, \mum');




good_points = find((total_length <= max_length)&(total_length >= min_length)); 
N_good = length(good_points);
disp(['Got ' num2str(N_good) ' good points!']);
disp(['Fraction of good points: ' num2str(N_good/N*100) '%']);

%rename vectors after removal of incorrect frames from consideration
%back to x,t (for ease of coding)
t = t(good_points);
total_length = total_length(good_points);
x = x(good_points);

figure(200);
plot(t,total_length,'ob');

mean_length = mean(total_length);

min_length = mean_length - std_length;
max_length = mean_length + std_length;

hold on;
plot(t,x,'.r');
plot(t,mean_length*ones(N_good,1),'-r','LineWidth',2);
plot(t,(min_length)*ones(N_good,1),'-b','LineWidth',1);
plot(t,(max_length)*ones(N_good,1),'-b','LineWidth',1);
hold off;
title('Tracking data - only good points');
xlabel('Time, frames');
ylabel('Position, pixels');


%With refined time and position measurements, let us now plot mean square
%displacement graph from single trajectory. This task would be elmentary
%if all recorded positions of the bundle x were separated by same time
%interval. One could have used existing script get_msd_from_trajectory.
%Since this is not the case, new script called get_msd_from_x_t is
%developed.
[time,msd,error] = get_msd_from_x_t(x,t);
time = time*t_scale;
figure(300);
h_plot = errorbar(time,msd',error,'or');
errorbar_tick(h_plot,1000);
hold on;
plot(time,.5*time,'.k');

%loglog(time,time.^2,'--k');


%In ordrer to see if removal of bad frames influences msd plot, 
% let us also extract msd without frame removal
% [msd_all_frames, error_all_frames] = get_msd_from_trajectory(tracking_data(:,2),tracking_data(:,2)*0);
% t1 = 1:1:N-1;
% 
% h_plot = errorbar(t1,msd_all_frames,'ob');
% errorbar_tick(h_plot,1000);

hold off;
set(gca,'XScale','log');
set(gca,'YScale','log');
axis tight;
set(gca,'FontSize',20);
%title('LOG-LOG MSD plot of bundle position vs time');
xlabel('\fontsize{30}Time, s');
ylabel('\fontsize{30}Bundle MSD, \mum^2');
legend('Bundle MSD','Slope 1 (diffusion)','MSD from raw data');

%plot msd with errorbars
figure(400);
h = errorbar(time,msd',error,'or');
errorbar_tick(h);
hold on;
end

