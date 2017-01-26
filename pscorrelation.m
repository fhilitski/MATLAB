%%
clc;
clear all;
close all;
cal_path = 'D:\Data - MT Sliding and Friction\2016\09-20-2016\findcenter 1\';
cal_name = 'cal 1 trap 0.79W.dat.mat';
load([cal_path cal_name]);
Vx0 = Vx;
Vy0 = Vy;
%Vx = Vx - mean(Vx);
%Vy = Vy - mean(Vy);
h_fig_rotation = figure;
colors = colormap;
i_color = 0;
sampling_frequency = rate;

do_log_bin = false;
plot_voltage_ps = true;
%%
rotation_steps = [3.5];
for alpha = rotation_steps
    color_step = floor(length(colormap)/length(rotation_steps));
    i_color = i_color + color_step;
    %create rotation matrix for a given angle alpha
    r_matrix = [cosd(alpha), -sind(alpha);...
        sind(alpha),  cosd(alpha)];
    disp(['Rotating by ' num2str(alpha) ' degrees']);
    %all N initial voltages (Vx, Vy) are packed into 2xN matrix
    v0 = [Vx0 Vy0]';
    %rotate the matrix of intitial voltages
    v = r_matrix * v0;
    %new voltages Vx and Vy are...
    Vx = v(1,:);
    Vy = v(2,:);
    %create figure for plotting the power spectrum in necessary
    if (plot_voltage_ps)
        h_ps_fig = figure;
    end;
    
    for coord_index = 1:1:2
        %perfomr FT calculation for both coordinates
        %x = 1
        %y = 2
        
        separation = v(coord_index,:);
        N = length(separation);
        data_array = separation - mean(separation);
        total_points = length(data_array);
        %fft works with 2^N datapoints best
        next2 = nextpow2(total_points);
        if ~(2^next2 == total_points)
            total_points = 2^(nextpow2(total_points)-1);
        end;
        %total time of acquisition
        total_time = total_points/sampling_frequency;
        
        %there are N discrete frequeincies
        %from 0 to f_max
        %actually, 0 frequency is never achieved, as it requires infinite
        %acquisition time; minimum frequency we achieve is 1/total_time
        %max frequency f_max is sampling frequency
        %actually, the power spectrum is only true for half of this frequency,
        % called Nyquist frequency, f_sampling/2
        
        %create vector of frequencies
        %%f = linspace(-sampling_frequency/2,sampling_frequency/2,total_points);
        f1 = linspace(1/total_time,sampling_frequency/2,total_points/2);
        f2 = fliplr(f1.*(-1));
        f = [f2 f1];
        
        %find first index after f = 0;
        f0 = find(f >= 0); %this looks like a bug since f=0 points are actually included
        f0 = f0(1);
        %freq and powerspectrum reperesent one-sided power spectrum (positive frequencies)
        freq = f(f0:end);
        
        %FT of data
        %units are arbitrary so far
        ft_sep = fft(data_array, total_points);
        abs_ft = abs(ft_sep);
        
        %normalization of the ft
        %Following Gittes and Schmidt in MCB (vol. 55)
        %units are micron^2/Hz
        abs_ft_norm = (2*total_time/(total_points^2))*(abs_ft).^2;
        
        %find mean(x^2) of the data
        mean_x2 = mean(data_array(1:total_points).^2);
        mean_x_2 = mean(data_array(1:total_points))^2;
        var_x = mean_x2 - mean_x_2;
        
        %find mean(x^2) from the power spectrum
        %sum(S_k) from 0 to N/2 /Total_time = mean(x^2)
        mean_x2_ps = sum(abs_ft_norm(1:length(freq)))/total_time;
        
        %the two should be equal
        normalization_const = mean_x2/mean_x2_ps;
        abs_ft_norm = normalization_const*(2*total_time/(total_points^2))*(abs_ft).^2;
        ps = abs_ft_norm(1:length(freq));
        
        %save FT result for further work with correlations
        R(coord_index,:) = ft_sep;
        
        if (plot_voltage_ps)
            figure(h_ps_fig);
            plot(freq(2:end), ps(2:end));
            hold on;
            axis tight;
            set(gca,'Xscale','log');
            set(gca,'Yscale','log');
            set(gca,'FontSize',12);
            xlabel('\fontsize{24} Frequency (Hz)');
            ylabel('\fontsize{24} PS density (arb.units^2/Hz)');
            title(['Voltage PS: ' num2str(alpha) '^o rotation']);
            if (do_log_bin)
                %use blogf.m to get logarithmic binning of power spectrum
                fraw = freq';
                fraw(:,2) = ps;
                result = blogf(fraw);
                
                %now plot log-binned power spectrum
                freq_binned = result(:,1);
                power_binned = result(:,2);
                power_sd_binned = result(:,3);
                
                h_binned_ps = errorbar(freq_binned(2:end),...
                    power_binned(2:end),...
                    power_sd_binned(2:end),...
                    'o','LineStyle','none','LineWidth',1.5);
            end; %if do_log_bin    
        end; %if plot_voltage_ps
    end; %for loop coord_index
    %% calculate the correlation
    %R now contains FTs of Vx and Vy such that
    %Rx = R(1,:);
    %Ry = R(2,:);
    %calculating |Rx|^2 and |Ry|^2
    P = abs(R).^2;
    Pxy = real(R(1,:) .* conj(R(2,:)));
    %finally, normalize the correlation coeff
    cc = Pxy ./ sqrt(P(1,:).*P(2,:));
    %Px = abs(Rx).^2;
    %Py = abs(Ry).^2;
    %Pxy = real(Rx .* conj(Ry));
    %cc = Pxy ./ sqrt(Px .* Py);
    
    % optionally plot Px and Py. This is not useful...
    % figure;
    % plot(freq(2:end), Px(2:length(freq)));
    % hold on;
    % plot(freq(2:end), Py(2:length(freq)));
    % set(gca,'Xscale','log');
    % set(gca,'Yscale','log');
    %%
    [fb, ccb, cce, ccsdm] = binned_xy(freq(2:end), cc(2:length(f)), true, 1000);
    figure(h_fig_rotation);
    h_series = errorbar(fb,ccb,ccsdm, '-o','Color',colors(i_color,:));
    ylabel('\fontsize{24} P_{xy}/(P_x P_y)^{-2}');
    xlabel('\fontsize{24} Frequency (Hz)');
    axis tight;
    set(gca,'FontSize',12);
    h_series.DisplayName = [num2str(alpha) ' degrees rotation'];
    hold on;
end;