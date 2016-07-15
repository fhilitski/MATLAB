function [freq, powerspectrum, total_points, normalization_const] = power_spectrum_FT(data_array, sampling_frequency)
%POWER_SPECTRUM_FT Returns one-sided power spectrum of data
%
%[freq, powerspectrum, total_points, normalization_const] = power_spectrum_FT(data_array, sampling_frequency)
%  Returns:
%           freq - a vector of positive frequencies;
%           powerspectrum - vector of corresponding values of FT power spectrum
%                           of data;
%           total_points - total number of data-points included in the spectrum;
%           normalization_const - ratio between variance of the <x^2> of
%                                 the data and sum of fourier coeffs. 
%                                 If the data has 0 mean, this coefficient
%                                 should be 1. Otherwise, 
%                                 See Gittes and Schmidt in MCB (vol. 55)
%        
%  Inputs: 
%           data_array - vector of values of time-dependent quantity;
%           sampling_frequency - frequency of acquisition for values in the
%           data_array;
%


%total number of datapoints
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
f0 = find(f >= 0);
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


powerspectrum = abs_ft_norm(1:length(freq));


