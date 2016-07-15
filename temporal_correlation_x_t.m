function [time, correlation, errors] = temporal_correlation_x_t(x, t)
%TEMPORAL_CORRELATION_X_T Computes autocorrelation function for 
%the signal x(t); 
%
% [time, correlation, errors] = temporal_correlation_x_t(x, t)


%total number of timesteps
N_x = length(x);
N_t = length(t);

if ~(N_x == N_t)
    error('X and T must have unequal lengths!');
else
    N = N_x;
end;


%initialize correlation and time vector
time = zeros(1,N);
correlation = zeros(1,N);
%initialize errorbar vector
errors = zeros(1,N);

%initialize temporry corr vector, that will turn into matrix and hold all
%possible correlation values x(t)*x(t+dt)
corr = zeros(1,N);

%start with the first elements of the vectors, x1 and t1
%calculate all possible dt(1i) and x(1)x(i), i = 2:N
%here, x(1)x(1) will correspond to time = 0;
corr = x * x(1);
time = t - t(1);

%this makes vectors dt and corr
%Now calculate all possible dt(2,i), x(2)*x(i) for elements i=2..N
%Here, time dt(2,2)=0 corresponds to x(2)x(2)
for i = 2:(N)
    
    dx2 = x * x(i);
    dt = t - t(i);
    %this calculates everything including x2*x2, but this was already
    %done, so we need to disregard it;
    dx2 = dx2(i:end);
    dt = dt(i:end);
    
    %now we need to loop over all timesteps in dt and find same timesteps
    %in time vector
    %once we found corresponding time vector, record appropriate sd in
    %table of sds;
    for j = 1:length(dt)
        index = find( dt(j) == time);
        if ~isempty(index)
            
            total_displacements = length(find(corr(index,:) ~= 0));
            corr(index,total_displacements+1) = dx2(j);
        else
            %we have encountered dt not previously registered!
            %it needs to be added to the end of the time and correlation matrices;
            last = length(time);
            index = last+1;
            time(index) = dt(j);
            total_displacements = 0;
            corr(index,total_displacements+1) = dx2(j);
            disp('New time!');
        end;
    end;
end;

%at this point, we have vector of all time intervals, called time
%and matrix of all square displacements.
%for given time time(i),
%all square displacements are put into row sd(i,:)
%and empty spaces are zeroes.
%
%one needs to calculate mean square displacement for each time,
%disregarding zeroes;

total_time_points = length(time);
for i = 1:total_time_points
    all_nonzero_sds = find(corr(i,:) ~= 0);
    real_sds = corr(i,all_nonzero_sds);
    if isempty(real_sds)
        real_sds = 0;
    end;
    correlation(i) = mean(real_sds);
    errors(i) = std(real_sds)/sqrt(length(real_sds));
end;






