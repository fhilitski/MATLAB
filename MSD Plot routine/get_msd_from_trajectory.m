function [ msd, errors ] = get_msd_from_trajectory( x, y )
%GET_MSD_FROM_TRAJECTORY Returns MSD from particle trajectory
%   Returns mean square displacement vector from trajectory of a diffusing
%   particle given by x and y, where each element is displacement of a particle
%   from intial position r0^2 at time t.
%   All points x(i) and y(i) must be separated by same value of delta*t;

%total number of timesteps
N_x = length(x);
N_y = length(y);

if ~(N_x == N_y)
    error('Trajectories have unequal lengths!');
else
    N = N_x;
end;


%initialize msd vector
msd = zeros(1,N-1);
%initialize errorbar vector
errors = zeros(1,N-1);


for i = 1:(N-1)
    %the smallest possible time-step is 1*t
    %and the largest possible time-step in a given trajectory is N-1
    %calculate msd for each possible step length
    
    
    %this variable holds displacement for given time-step i (i = 1*t_step,
    % 2*t_step, ...,N*t_step;
    displacement = [];
    counter = 0;
    
    for k = 1:(N - i);
        
        x0 = x(k);
        y0 = y(k);
        
        j = k + i;
        r2 = (x(j) - x0)^2 + (y(j) - y0)^2;
        displacement(k) = r2;
        counter = counter + 1;
        
    end;
    
    msd(i) = mean(displacement);
    errors(i) = std(displacement)/sqrt(length(displacement));
end;


