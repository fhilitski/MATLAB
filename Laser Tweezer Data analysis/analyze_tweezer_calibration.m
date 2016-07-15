function [corner_freq_x, corner_freq_y, D_x, D_y] = analyze_tweezer_calibration (folder_name, file_name, number_of_splits, min_average_points, high_freq_cutoff, low_freq_cutoff)

if nargin < 3
    number_of_splits = 1;
    min_average_points = 3;
    high_freq_cutoff = 10000;
    low_freq_cutoff = 100;
end;
if nargin < 4
    min_average_points = 3;
    high_freq_cutoff = 10000;
    low_freq_cutoff = 100;
end;
if nargin < 5
    high_freq_cutoff = 10000;
    low_freq_cutoff = 100;
end;
if nargin < 6
    low_freq_cutoff = 0;
end;

fname = [folder_name file_name];

[rate, d, g, z, Vx,Vy, Nx, Ny] = read_calibration_data(fname);

%plot_voltages
h_fig_vx = figure;
h_fig_vy = figure;
%plot_voltages
figure(h_fig_vx);
plot(Vx,'.','Color',[1, 0, 0]);
title('Raw data: x-asis signal');
xlabel('Datapoint #');
ylabel('Voltage (V)');

figure(h_fig_vy);
plot(Vy,'.','Color',[0, 0, 1]);
title('Raw data: y-axis signal');
xlabel('Datapoint #');
ylabel('Voltage (V)');

if (Nx ~= Ny)
    corner_freq_x = -1;
    corner_freq_y = -1;
    return;
else
    split_length = floor(Nx/number_of_splits); 
    Vx_splits = zeros([split_length, number_of_splits],'double');
    Vy_splits = Vx_splits;
    
    split_start = 1;
    for i = 1:number_of_splits
        split_end = split_start + split_length - 1;
        Vx_splits(:,i) = Vx(split_start : split_end);
        Vy_splits(:,i) = Vy(split_start : split_end);
        split_start = split_end+1;        
    end;
end;
    

% fourier transform
[corner_freq_x, D_x] = fourier_transform_routine(Vx_splits,rate,min_average_points, high_freq_cutoff, low_freq_cutoff);
[corner_freq_y, D_y] = fourier_transform_routine(Vy_splits,rate,min_average_points, high_freq_cutoff, low_freq_cutoff);
fprintf('Corner freq x: %3.0f Hz\n', corner_freq_x);
fprintf('Corner freq y: %3.0f Hz\n', corner_freq_y);
fprintf('Diffusion coeff in x: %3.0f (arb.units)^2/s\n', D_x);
fprintf('Diffusion ceoff in y: %3.0f (arb.units)^2/s\n', D_y);
