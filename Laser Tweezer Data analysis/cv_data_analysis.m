%% read data from binary

%fileID = fopen('cal 0.69W.dat');
%[data,count] = fread(fileID, 'double'); % Scan the data into a vector, in this case called data
%fclose(fileID);

%% read cv data
clear all;
close all;

path = 'D:\Data - MT Sliding and Friction\2015\03-31-2015 1.0 pct PEG 200 mM K+ with KSA (1-to-8 ratio) - sliding, buckling\run3';
run_folder = '\';
fname = 'data.csv';
filename = [path run_folder fname];

[trap_x, trap_y, pos_x, pos_y, t] = read_cv_data(filename);

%% read acquisition parameters

fname = 'acquisition.csv';
filename = [path run_folder fname];

[Fsampling,Nsamples,Velocity,MaxMove,N_traps,T_delay,Converted] = read_acq_data (filename);


%% plot trap position
% only y-coordinate is changing, x

h_fig = figure;
plot(t, trap_y, 'ob');
title('Moving trap position, y axis');


% if there is delay in trap movement,
% it has to be excluded from consideration

elapsed_t = t-t(1);
t_index = find(elapsed_t >= T_delay);


linfit = fit(t(t_index),trap_y(t_index),'poly1');
hold on;
plot(linfit);

fit_coeffs = coeffvalues(linfit); 
v = fit_coeffs(1);

%v returns velocity in nm/ms 
v = v*1000;

%delta_t has all time-steps
delta_t = t(2:end) - t(1:end-1);

clc;
fprintf('Current RUN: %s \n', run_folder(1:end-1));
if (~Converted)
    disp('Conversion NOT USED!');
end;
fprintf('Actual trap velocity: %.2f nm/s \n', v);
fprintf('Set trap velocity: %2.2f nm/s \n', Velocity);

fprintf('\n');
fprintf('Expected timestep: %2.2f ms\n', Nsamples/Fsampling * 1000);
fprintf('Actual timestep: %.2f ms\n', mean(delta_t));

fprintf('\nDelay, ms: %2.2f \n', T_delay);
fprintf('Sampling Freq, Hz: %2.2f \n', Fsampling);
fprintf('# of samples: %i \n', Nsamples);
fprintf('Total # of traps: %i \n', N_traps);
fprintf('Traps moving: %i \n', MaxMove+1);

%% plot x and y bead positions in the detection beam

h_x = figure;
plot(t, pos_x, '.r');
title('Detection beam bead position, x');
axis tight;

h_y = figure;
plot(t,pos_y,'.b');
title('Detection beam bead position, y');
axis tight;

%% attempt fit of the y data

%select linear part of the plot
figure(h_y);
[times, y_coords] = ginput(2);

t1 = find(t >=times(1));
t1 = t1(1);

t2 = find( t>=times(2));
t2 = t2(1);

y1 = y_coords(1);
y1 = find(pos_y >= y1);
y1 = y1(1);

y2 = y_coords(2);
y2 = find(pos_y >= y2);
y2 = y2(2);

%plot the interval;
hold on;
plot(t(t1:t2),pos_y(t1:t2),'or');

%fit line
fit_2 = fit(t(t1:t2),pos_y(t1:t2),'poly1');
figure(h_y);
hold on;
plot(fit_2);

fit_coeffs_2 = coeffvalues(fit_2); 
v_measured = fit_coeffs_2(1);
%v_measured returns velocity in nm/ms 
%convert it to nm/s - that's the value we set in AOD GUI
v_measured = v_measured*1000;

fprintf('\nMeasured velocity: %4.2f nm/s \n', v_measured);
fprintf('Linear interval of bead movement: %4.2f nm \n', pos_y(t2)-pos_y(t1));
fprintf('Total experiment time: %4.2f s \n', 1/1000*(t(end) - t(1)));

%% plot timestep vs times
figure;
plot(delta_t,'.b')


