function [t, total_length, t_s, total_length_um, segment_length_um] = process_length_data(tracking_data, t_scale, pixel_size)
%PROCESS_LENGTH_DATA plots and saves total filament length analysis
%PROCESS_LENGTH_DATA(tracking_data, sensitivity, t_scale, pixel_size)
%   tracking_data - Raw data from bundle sliding tracking has form 
%                  [frame#, Position , Total Length]
%                    where frame# is tracked frame #;
%                    position is midpoint of the bundle position along filament's
%                    contour;
%                    Total length is total contourlength (of the whole filament);
%   t_scale       - time between frames
%   pixel_size    - pixel size in microns
%
%WARNING: t_scale and pixel size are not yet implemented; they are both
%effectively = 1; 
%
%   All measurements of tracking_data are in pixels and frames, and need to be rescaled with
%   appropriate pixel-to-um ratio and frame capture rate;

clc;


t = tracking_data(:,1);
total_length = tracking_data(:,3);
segment_length = tracking_data(:,4);
x = tracking_data(:,2);
x = x - x(1);
N = length(tracking_data);

t_s = t.*t_scale;
total_length_um = total_length.*pixel_size;
segment_length_um = segment_length.*pixel_size;

figure;
plot(t,total_length,'o');
p1 = 'Total length: raw data';
hold on;
smooth_total_length = smooth(total_length,15);
plot(t, smooth_total_length, '-');
p2 = 'Total length: smooth';
hold off;
title('Raw data from filament tracking');
xlabel('Time, frames');
ylabel('Position, pixels');
legend({p1,p2});

figure;
plot(t_s,total_length_um,'o');
p1 = 'Total length: raw data';
hold on;
smooth_total_length = smooth(total_length_um,15);
plot(t_s, smooth_total_length, '-');
p2 = 'Segment length: smooth';
hold off;
title('Raw data from filament tracking');
xlabel('Time (s)');
ylabel('Length (\mum)');
legend({p1,p2});

figure;
plot(t_s,segment_length_um,'o');
p1 = 'Total length: raw data';
hold on;
smooth_segment_length = smooth(segment_length_um,15);
plot(t_s, smooth_segment_length, '-');
p2 = 'Segment length: smooth';
hold off;
title('Raw data from filament tracking');
xlabel('Time (s)');
ylabel('Length (\mum)');
legend({p1,p2});



end
