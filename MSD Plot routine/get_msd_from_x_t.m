function [time, msd, errors] = get_msd_from_x_t( x,t )
% GET_MSD_FROM_X_T Given trajectory x(t) in one dimension, returns...
%      msd as a function of time-step dt <x(dt)^2>
% [time, msd, errors] = get_msd_from_x_t( x,t )
%   This code is especially useful is time-steps t are of unequal size;
%   If timesteps are equal, it is faster to use get_msd_from trajectory


%total number of timesteps
N_x = length(x);
N_t = length(t);

if ~(N_x == N_t)
    error('X and T must have unequal lengths!');
else
    N = N_x;
end;


%initialize msd and time vector
time = zeros(1,N);
msd = zeros(1,N);
%initialize error vector
errors = zeros(1,N);

%initialize temporry sd vector, that will turn into matrix and hold all
%possible square displacements (sd)
sd = zeros(1,N);

%start with the first elements of the vectors, x1 and t1
%calculate all possible dt(1i) and dx(1i), i = 2:N
sd = (x - x(1)).^2;
time = t - t(1);
%this makes vectors dt and dx
%Now calculate all possible dt(2,i), dx2(2,i) for elements i=2..N
%obviously, dt(2,2)=0 and dx(2,2) = 0

for i = 2:(N-1)
    
    dx2= (x - x(i)).^2;
    dt = t - t(i);
    %this calculates everything including x1 - x2, but this was already
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
            
            total_displacements = length(find(sd(index,:) ~= 0));
            sd(index,total_displacements+1) = dx2(j);
        else
            %we have encountered dt not previously registered!
            %it needs to be added to the end of the time and msd matrices;
            last = length(time);
            index = last+1;
            time(index) = dt(j);
            total_displacements = 0;
            sd(index,total_displacements+1) = dx2(j);
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
    all_nonzero_sds = find(sd(i,:) ~= 0);
    real_sds = sd(i,all_nonzero_sds);
    if isempty(real_sds)
        real_sds = 0;
    end;
    msd(i) = mean(real_sds);
    errors(i) = std(real_sds)/sqrt(length(real_sds));
end;






