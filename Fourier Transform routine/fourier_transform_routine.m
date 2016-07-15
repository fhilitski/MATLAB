function [fc, D] = fourier_transform_routine(separation, sampling_freq, n_min, f_max, f_min)
%FOURIER_TRANSFORM_ROUTINE performs fourier transform of given time-series,
%log binning and lorentzian fit;
%Refer to K. Berg-Sørensen and H. Flyvbjerg, Review of Scientific Instruments 75, 594 (2004)
%for information on lornetzian coeff. calculation
%FOURIER_TRANSFORM_ROUTINE(separation)
%displays user interfce dialog toenter sampling frequency;
%FOURIER_TRANSFORM_ROUTINE(separation,sampling_freq) 
%does not display the dialog, perfomrs FT using provided sampling
%frequency;
%FOURIER_TRANSFORM_ROUTINE(separation,sampling_freq, n_min)
%does not display the dialog, and uses only data resulted from averaging at least n_min
%power-spectrum points during the log-binning process for the Lorentzian fit;
%by default, n_min = 1;
%FOURIER_TRANSFORM_ROUTINE(separation,sampling_freq,  n_min, f_max)
%in additon to minimum averaging creterion, there is upper frequency cutoff f_max
%used in the Lornetzina fit. By default, the cut-off is half the sampling
%frequency (aka Nyquist frequency).
%INPUTS:
%       separation - time series of a given variable (first was designed for
%                    separation between beads, hence the name);
%                    If separation is Nx1 vector, performs FT, log-binning and fit;
%                    if separation is NxM matrix, it treats the matrix as M independent 
%                    runs, performs M Fourier transforms, finds RMS-average of them and 
%                    then performs log-binning and fits.
%RETURNS:
%       fc - corner frequency from lorentzian fit;
%       D - the diffusion coefficient

%Last revision: 04/16/2015 Feodor Hilitski (Dogic lab, Brandeis)

warning('off','MATLAB:Axes:NegativeDataInLogAxis');
 
min_points_to_average = 1;

%we need time between images, i.e. frame rate of the aquisition
if nargin == 1
    
    prompt = 'Enter number of frames per second: ';
    dialog_title = 'FT options...';
    user_input = inputdlg(prompt, dialog_title);
    if (length(user_input) > 0)
        frames_per_sec = str2num(user_input{1});
    else
        frames_per_sec = 0;
    end;
    n_min = min_points_to_average;
    sampling_freq = frames_per_sec;
end;

if nargin >= 2
    frames_per_sec = sampling_freq;
    
end;

if nargin >= 3
    min_points_to_average = n_min;
end;

if nargin >= 4
    freq_cutoff = f_max;
else
    
    freq_cutoff = sampling_freq/2;
    f_max = freq_cutoff;
end;

if nargin < 5
    f_min = 0;
else
   %disp(['Min freq: ' num2str(f_min)]); - for debug only
end;

if (frames_per_sec ~= 0)
    
    sampling_freq = frames_per_sec;
    
    %if the separation is comprised of N separate runs, average the
    %power spectra
    figure;
    hold on;
    
    [M,N] = size(separation);
    for i = 1:N
        
        single_run = separation(:,i) - mean(separation(:,i));
        [f(:,i),ps(:,i), n_tot(:,i), norm(:,i)] = power_spectrum_FT(single_run, sampling_freq);
        %disp(['Power spectrum normalization: ' num2str(norm)]);
        %first, plot raw data - the original power spectrum
        
        plot(f(2:end,i), ps(2:end,i), '.','Color',[(135*i)/(N*255) (135*i)/(N*255) (135*i)/(255*N)]);
        
    end;
    axis tight;
    set(gca,'Xscale','log');
    set(gca,'Yscale','log');
    set(gca,'FontSize',12);
    xlabel('\fontsize{24} Frequency (Hz)');
    ylabel('\fontsize{24} PS density (arb.units^2/Hz)');
    
    f = mean(f,2);
    %average the power spectrum
    ps_reg = mean(ps,2);
    %use RMS average:
    ps = sqrt(sum(ps.^2,2)/N);
    
    plot(f(2:end), ps(2:end),'.','Color',[1 0 0]);
    plot(f(2:end), ps_reg(2:end),'.','Color',[1 1 0]);
    
    ps = ps_reg;
    %use blogf.m to get logarithmic binning of power spectrum
    result = blogf(cat(2,f,ps));
    
    %plot one-sided power spectrum
    h_power_spectrum = figure;
    
    %first, plot raw data - the original power spectrum
    plot(f(2:end), ps(2:end), '-','Color',[1 75/255 75/255]);
    hold on;
    
    %now plot log-binned power spectrum
    freq = result(:,1);
    power = result(:,2);
    power_sd = result(:,3);
      
        h_binned_ps = errorbar(freq(2:end),power(2:end),power_sd(2:end),'ok',...
            'LineStyle','none','LineWidth',1.5);
        axis tight;
        set(gca,'Xscale','log');
        set(gca,'Yscale','log');
        set(gca,'FontSize',12);
        xlabel('\fontsize{24} Frequency (Hz)');
        ylabel('\fontsize{24} PS (V^2/Hz)');
        
                
        %fit Lorentzian of the form  P(f) = a / (1 + (f/fc)^2),
        %where fc is the corner frequency.
        %drop points in the log binning that result from a single raw data,
        %these are noisy and have no statistics.
       
        good_for_fit = find(result(:,4) > min_points_to_average);
        result_for_fit = result(good_for_fit,:);
        fit_freq = result_for_fit(:,1);
        fit_power= result_for_fit(:,2);
        result_sd = result_for_fit(:,3);
        
        %introduce high frequency cut-off, above which the ps is inaccuarte  
        %freq_cutoff;
        good_for_fit = find(fit_freq <= freq_cutoff);
        result_for_fit = result_for_fit(good_for_fit,:);
        fit_freq = result_for_fit(:,1);
        fit_power= result_for_fit(:,2);
        result_sd = result_for_fit(:,3);
        
        %plot good_for_fit datapoints on the same plot;
        h_fit_ps = errorbar(fit_freq,fit_power,result_sd,'oy',...
            'LineStyle','none','LineWidth',2);
        
        %introduce low frequency cut-off, below which the ps is inaccuarte
        %this is on top of the high_freq cut-off
        %f_min;
        good_for_fit = find(fit_freq >= f_min);
        result_for_fit = result_for_fit(good_for_fit,:);
        fit_freq = result_for_fit(:,1);
        fit_power= result_for_fit(:,2);
        result_sd = result_for_fit(:,3);
        
        %plot good_for_fit datapoints on the same plot;
        h_fit_ps = errorbar(fit_freq,fit_power,result_sd,'og',...
            'LineStyle','none','LineWidth',2);
        
        
        %check if there is enough data points
        if (length(fit_power) >=2 )
           
            %This part is trying to fit Lorentzian into the power spectrum
            %using standard Matlab fitting methods (mean square)
            %it has not been very sucessfull, and has been commented out on
            % 06/10/2016
%             
%             %power is usually very large number, and fit is not doing good
%             %job with such numbers; 
%             %we normalize with by its max value, and fit normalized one.
%             max_power = max(fit_power);
%             normalized_fit_power = fit_power./max_power;
%             
%             %%[fit_result, how_bad] = LorentzianFitWeighted(fit_freq, normalized_fit_power,1./result_sd);
%             [fit_result, how_bad] = LorentzianFit(fit_freq, normalized_fit_power);
%             %[fit_result, how_bad] = LorentzianFit(fit_freq, fit_power);
%             
%             coeffs = coeffvalues(fit_result);
%             confidence = confint(fit_result,0.69);
%             a = coeffs(1)* max_power; %recover real a from normalization;
%             %a = coeffs(1);
%             b = coeffs(2);
%             conf_int_b = confidence(:,2);
%             
%             fc_fit = 1/b;
%             dfc = mean(abs(1./conf_int_b - fc_fit));
%             
%             %create variable to plot lorentzian fit together with power spectrum;
%             ps_fit = a ./(1 + b^2 *freq.^2);
%             figure(h_power_spectrum);
%             loglog(freq,ps_fit,'--k','LineWidth',2);
%             
            %calcualtion of coefficients - 
            %It is possible to directly calculate the coefficients and 
            %thus fit the Lorentzian
            %see Lorentzian_coeffs function for details
            [f1,D1] = Lorentzian_coeffs(fit_freq,fit_power);
            ps_fit2 = D1 ./(f1^2 + freq.^2);
            figure(h_power_spectrum);
            loglog(freq,ps_fit2,'-b','LineWidth',2);
                       
            fc = f1; % = k/(2*pi*gamma)
            D = D1*2*pi^2; % = k_B T / (gamma)
            %see the paper by Flyvberg,Born-Sorensen for the explanation;
            %Interestingly, another popular paper (Gittes, Shmidt MCB, vol
            %55) which I refer to when the fourier transfor of the raw data
            %is calculated, implies (see page 138) that 
            D = D1*pi^2; 
            % i.e. without the factor of 2.
            %This value seems to produce physically correct diffusion
            %coefficient values and it is used in the subsequent analysis
            
                        
            title_string = {['N=' num2str(n_tot(1))];...
                ['Corner frequency calculated: ' num2str(f1,3) ' Hz'];...
                ['Diffusion coefficient: ' num2str(D,3) '(arb.units)^2/s']};
            title(title_string);
            axis tight;
        else 
            fc = 0;
            D = 0;
            title('Power spectrum not Lorentzian, no corner frequency!');
        end;        
else
    disp('Invalid sampling freq. input!!');
    fc = 0;
    D = 0;
end;
warning('on','MATLAB:Axes:NegativeDataInLogAxis');


