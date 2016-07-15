%% Example from Matlab Help 
clc;
clear all;

Fs = 1000;                    % Sampling frequency
T = 1/Fs;                     % Sample time
L = 1000;                     % Length of signal
t = (0:L-1)*T;                % Time vector
% Sum of a 50 Hz sinusoid and a 120 Hz sinusoid
x = 0.7*sin(2*pi*50*t) + sin(2*pi*120*t); 
y = x + 2*randn(size(t));     % Sinusoids plus noise
figure(1);
plot(Fs*t(1:50),y(1:50))
title('Signal Corrupted with Zero-Mean Random Noise')
xlabel('time (milliseconds)')
 

%It is difficult to identify the frequency components by looking at the original signal. Converting to the frequency domain, the discrete Fourier transform of the noisy signal y is found by taking the fast Fourier transform (FFT):
NFFT = 2^nextpow2(L); % Next power of 2 from length of y
Y = fft(y,NFFT)/L;
f = Fs/2*linspace(0,1,NFFT/2+1);

% Plot single-sided amplitude spectrum.
figure(2);
plot(f,2*abs(Y(1:NFFT/2+1))) 
title('Single-Sided Amplitude Spectrum of y(t)')
xlabel('Frequency (Hz)')
ylabel('|Y(f)|')

%% my FFT applied to same signal
f = linspace(-Fs/2, Fs/2, length(t));
y_w = fft(y);

y_w_abs = abs(fft(y))/length(t);
y_w_abs_sq = y_w_abs.^2;

figure(3);
f0 = find(f == 0);
plot(f,2*y_w_abs,'b');
axis tight;
%%plot(f,2*y_w_abs,'b');

std1 = std(y);
std2 = std(y_w_abs);