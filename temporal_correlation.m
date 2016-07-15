function [ time, correlations, errors ] = temporal_correlation(x)
%TEMPORAL_CORRELATION Produces autocorrelation function for x(t)
%	TEMPORAL_CORRELATION Takes the input x_i, i = 1..N and calculates time 
%                           autocorrelation function R(dt); 
%                           time intervals are assumed to be 0,1,...,N
%[ time, correlations, errors ] = temporal_correlation(x)


N = length(x);
mju = mean(x);
sigma = std(x);

x = x - mju;

%time vector
time = 0:1:N-1;
%initialize msd vector
corr = zeros(N,N);
%initialize errorbar vector
errors = zeros(1,N);


for i = 1:N
    
   
    counter = 0;
    
    for k = i:N
        
        dt = k-i;   
        corr(k-i+1,i) = x(i)*x(k);
        counter = counter + 1;
        
    end;
    if (mod(i,100) == 0) 
        disp(['Frame ' num2str(i) ' out of ' num2str(N) ' done...']);
    end;
end;

correlations = zeros(1,N);
disp('Averaging... Almost done...');

for i = 1:N
    all_corrs = corr(i,:);
    non_zero_corrs = find(all_corrs ~= 0);
    real_corrs = all_corrs(non_zero_corrs);
    correlations(i) = mean(real_corrs)/sigma^2;
    errors(i) = std(real_corrs)/sqrt(length(real_corrs));
end;



