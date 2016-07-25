function [ x_c, y_c, QPDx_from_fit, QPDy_from_fit, dist_x, dist_y, mapping_x, mapping_y] =...
    get_findcenter_data(data_path, plot_scans, scan_number, AODx, AODy)
%GET_FINDCENTER_DATA obtains and fits the data from the FindCenter LabView routine
% Syntax:   [ x_c, y_c, QPDx_from_fit, QPDy_from_fit, dist_x, dist_y, mapping_x, mapping_y] =
%           GET_FINDCENTER_DATA(data_path, plot_scans, scan_number, AODx, AODy)
%
% Inputs:   data_path - path of the folder that contains .csv files
%                       generated by the FindCenter LabView routine
%                       The require files are:
%                           centerXY.csv
%                           freqXY.csv
%                           voltXY.csv
%           plot_scans - boolean flag to display scan results in a separate
%                        plot; if scans are displayed, prompt user to
%                        identify two points on each scan and fit a line to
%                        determine QPD sensitivity.
%
%           scan_number - if multiple scans are saved, select the
%                         appropriate scan; set to the first scan by default
%
%           AODx,AODy  -  conversion from the Mhz to nm domain to get the
%                         displacements in nm;
%
% Outputs: x_c, y_c - center coordinates (in nm)
%          QPDx_from_fit, QPDy_from_fit - linear QPD sensitivity in nm/V
%          dist_x, dist_y - distances of linear QPD response in nm. Max allowable
%                           displacement from the center is dist_x/2.
%          mapping_x, mapping_y - function handles for mapping displacement
%                           into voltage using the entire QPD scan, not
%                           only the linear region.
%


%Revision history:  07/10/2016 Fitting of the linear part of the findcenter
%                               is supplemented by fitting the whole curve
%                               with a fourier series. Resulting function
%                               allows for a better (non-linear) conversion between
%                               coordinate in V into nm. Linear coefficiet
%                               is still kept of calibration analysis
%                               (conversion of the diffusion coeff. into
%                               nm).
%                   05/19/2015 Added display of center coordinates on scan
%                              plots; FH
%                   04/19/2015 Created; Feodor Hilitski (Dogic Lab, Brandeis)
%

mapping_x = NaN;
mapping_y = NaN;

if nargin < 3
    scan_number = 1;
end;

if nargin < 4
    AODx = 1;
    AODy = 1;
    convert_to_nm = false;
    
else
    convert_to_nm = true;
end;

dist_x = 0;
dist_y = 0;

filename_center = 'centerXY.csv';
filename_freq = 'freqXY.csv';
filename_volt = 'voltXY.csv';

%find folders with scans, which have names scan0, scan1, and so on
data_list = {};
data_count = 0;

all_files = dir(data_path);
all_dirs = {};
dir_count = 0;

for i = 1:length(all_files)
    file = all_files(i);
    if (file.isdir) && strncmp(file.name,'scan',4)
        dir_count = dir_count+1;
        all_dirs{dir_count} = all_files(i);
    end;
end;

total_scans = length(all_dirs);

if ((~isempty(all_dirs)) && (scan_number <= total_scans))
    %at least one scan is found and
    %the desired scan number is present
    i = scan_number;
    scan_path = [data_path '\' all_dirs{i}.name];
    filename_center = [scan_path '\' filename_center];
    filename_freq = [scan_path '\' filename_freq];
    filename_volt = [scan_path '\' filename_volt];
    [f_x, f_y] = import_csv_datafile(filename_freq);
    [v_x, v_y] = import_csv_datafile(filename_volt);
    [c] = import_csv_datafile(filename_center);
    
    %f_x and f_y are frequencies form AOD in Hz. Convert them in MHz
    f_x = f_x./10^6;
    f_y = f_y./10^6;
    
    %convert frequiencies into distances - um
    f_x = f_x.*AODx;
    f_y = f_y.*AODy;
    
    %remove last points as they often contain zeroes due to LabView
    %bug
    while (v_y(end) == 0)    
        f_y = f_y(1:end-1);
        v_y = v_y(1:end-1);
    end;
    
    while (v_x(end) == 0)
        f_x = f_x(1:end-1);
        v_x = v_x(1:end-1);
    end;
    
    %x_c and y_c contain center frequencies
    x_c = c(1) * AODx;
    y_c = c(2) * AODy;
    
    %shift the determined center to 0 nm purely for convenience
            f_y = f_y - y_c;
            f_x = f_x - x_c;
            x_c = 0;
            y_c = 0;
    
    %we need to find center voltages
    v_x_c_index = find(f_x >= x_c,1, 'last');
    v_y_c_index = find(f_y <= y_c,1, 'first');
    v_x_c = v_x(v_x_c_index);
    v_y_c = v_y(v_y_c_index);
    
    
    
    
    if plot_scans
        %x-scan
        h_scan_x = figure;
        plot(f_x, v_x, 'o');
        hold on;
        plot(x_c, v_x_c,'or','MarkerFaceColor','g');
        title('x-scan');
        ylabel('Voltage (V)');
        xlabel('Bead position (nm)');
        %plot(x_c, v_x_c, 'ob', 'MarkerSize', 10);
        hold off;
        axis tight;
        %select linear region of the response plot
        [f, v] = ginput(2);
        
        x1_index = find( f_x >= f(1));
        x1_index = x1_index(end);
        f_x1 = f_x(x1_index);
        v_x1 = v_x(x1_index);
        
        
        x2_index = find(f_x <= f(2));
        x2_index = x2_index(1);
        f_x2 = f_x(x2_index);
        v_x2 = v_x(x2_index);
        
        
        %plot the interval;
        hold on;
        plot(f_x(x2_index:x1_index),v_x(x2_index:x1_index),'o');
        if (convert_to_nm)
            dist_x = abs(f_x2 - f_x1);
            disp(['Max bead displacement in x (nm): ' num2str(dist_x)]);
        else
            dist_x = 0;
        end;
        
        %fit line
        fit_2 = fit(f_x(x2_index:x1_index),v_x(x2_index:x1_index),'poly1');
        figure(h_scan_x);
        x_lims = xlim;
        y_lims = ylim;
        hold on;
        h = plot(fit_2);
        legend(h, 'Linear region fit');
        fit_coeffs = coeffvalues(fit_2); %this is in V/nm
        %we are interested in the inverse mapping, i.e volts into nm
        %in case of a linear function everything is relatively simple:
        QPDx_from_fit = 1/fit_coeffs(1); %convert into nm/V
        
        %fit the fourier series
        [fData, vData] = prepareCurveData(f_x, v_x);
        % Set up fittype and options.
        ft = fittype( 'fourier4' );
        opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
        opts.Display = 'Off';
        opts.StartPoint = [0 0 0 0 0 0 0 0 0 0.001];
        % Fit model to data.
        [fit_x_fourier, gof] = fit(fData, vData, ft, opts);
        % Plot fit with data.
        h = plot(fit_x_fourier,'-b');
        set(h, 'DisplayName','Full fit - Fourier 4');
        %reset the plot limits
        xlim(x_lims);
        ylim(y_lims);
        %! HOWEVER, this function maps nm into volts. Experimentally, we
        %obtain volts and need to convert them to nm. Therefore, the
        %inverse of the mapping function is necessary.
        
        %prepare to fit the inverse function with fourier series
        %the function d(v) should be unambigously defined for each v
        %therefore, we only fit on the interval from
        %absolute max to absolute min in voltage
        [v_max_val, v_max_ind] = max(vData);
        [v_min_val, v_min_ind] = min(vData);
        start_ind = min(v_max_ind, v_min_ind);
        end_ind = max(v_max_ind, v_min_ind);
        vData = vData(start_ind:end_ind);
        fData = fData(start_ind:end_ind);
        plot(fData, vData, '*r');
        %fit the inverse data
        ft = fittype( 'fourier8' );
        opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
        opts.Display = 'Off';
        opts.StartPoint = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.4];
        % Fit model to data.
        [fit_x_fourier, gof] = fit(vData, fData, ft, opts);
        % Plot fit with data.
        %h = plot(fit_x_fourier,'-b');
        %set(h, 'DisplayName','Inverse fit - Fourier 8');
        %reset the plot limits
        xlim(x_lims);
        ylim(y_lims);
        
        %create a function handle for the fourier mapping
        mapping_x = fourier_curve_fit(fit_x_fourier);
        %test the functin handle by plotting the result
        h = plot(mapping_x(v_x),v_x,'.k');
        set(h, 'DisplayName', 'Fourier mapping function');
        
        
        
        %y-scan
        h_scan_y = figure;
        plot(f_y, v_y, 'o');
        hold on;
        plot(y_c, v_y_c,'or','MarkerFaceColor','g');
        title('y-scan');
        ylabel('Voltage (V)');
        xlabel('Bead position (nm)');
        %plot(x_c, v_x_c, 'ob', 'MarkerSize', 10);
        hold off;
        axis tight;
        %select linear region of the response plot
        [f, v] = ginput(2);
        y1_index = find( f_y <= f(1));
        y1_index = y1_index(1);
        f_y1 = f_y(y1_index);
        v_y1 = v_y(y1_index);
        
        y2_index = find(f_y >= f(2));
        y2_index = y2_index(end);
        f_y2 = f_y(y2_index);
        v_y2 = v_y(y2_index);
        
        %plot the interval;
        hold on;
        plot(f_y(y2_index:y1_index),v_y(y2_index:y1_index),'o');
        if (convert_to_nm)
            dist_y = abs(f_y2-f_y1);
            disp(['Max bead displacement in y (nm): ' num2str(dist_y)]);
        else
            dist_y = 0;
        end;
        
        %fit line
        fit_line = fit(f_y(y2_index:y1_index),v_y(y2_index:y1_index),'poly1');
        figure(h_scan_y);
        x_lims = xlim;
        y_lims = ylim;
        hold on;
        h = plot(fit_line);
        legend(h, 'Linear region fit');
        fit_coeffs = coeffvalues(fit_line);
        QPDy_from_fit = 1/fit_coeffs(1);
        
        %fit the fourier series
        [fData, vData] = prepareCurveData(f_y, v_y);
        % Set up fittype and options.
        ft = fittype( 'fourier4' );
        opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
        opts.Display = 'Off';
        opts.StartPoint = [0 0 0 0 0 0 0 0 0 0.001];
        % Fit model to data.
        [fit_y_fourier, gof] = fit(fData, vData, ft, opts);
        % Plot fit with data.
        h = plot(fit_y_fourier,'-b');
        set(h, 'DisplayName','Full fit - Fourier 4');
        %reset the plot limits
        xlim(x_lims);
        ylim(y_lims);
        %! HOWEVER, this function maps nm into volts. Experimentally, we
        %obtain volts and need to convert them to nm. Therefore, the
        %inverse of the mapping function is necessary.
        
        %prepare to fit the inverse function with fourier series
        %the function d(v) should be unambigously defined for each v
        %therefore, we only fit on the interval from
        %absolute max to absolute min in voltage
        [v_max_val, v_max_ind] = max(vData);
        [v_min_val, v_min_ind] = min(vData);
        start_ind = min(v_max_ind, v_min_ind);
        end_ind = max(v_max_ind, v_min_ind);
        vData = vData(start_ind:end_ind);
        fData = fData(start_ind:end_ind);
        plot(fData, vData, '*r');
        %fit the inverse data
        ft = fittype( 'fourier8' );
        opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
        opts.Display = 'Off';
        opts.StartPoint = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.2];
        % Fit model to data.
        [fit_y_fourier, gof] = fit(vData, fData, ft, opts);
        % Plot fit with data.
        %h = plot(fit_y_fourier,'-b');
        %reset the plot limits
        xlim(x_lims);
        ylim(y_lims);
        
        mapping_y = fourier_curve_fit(fit_y_fourier);
        h = plot(mapping_y(v_y),v_y,'.k');
        set(h, 'DisplayName', 'Fourier mapping function');
        
        %         close(h_scan_x);
        %         close(h_scan_y);
    else
        %if plotting scans was not necessary, QPDx and QPDy are not
        %determined
        QPDx_from_fit = 0;
        QPDy_from_fit = 0;
    end;
    
    
    
else
    %no scans found
    x_c = 0;
    y_c = 0;
    QPDx_from_fit = 0;
    QPDy_from_fit = 0;
end;

end

