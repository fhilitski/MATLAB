function [ x,y ] = generate_random_walk( n_steps, step_size, v_x, v_y, random_v, constrain_r)
%GENERATE_RANDOM_WALK generates 2-dimensional random walk
%RETURNS coordinates x,y of the particle performing random walk
%Usage:
%GENERATE_RANDOM_WALK(n_steps, step_size)
%generates random walk in 2D
%GENERATE_RANDOM_WALK( n_steps, step_size, v_x, v_y, random_v)
%generates 2D walk with ballistic components, which can be random
%Returns:
%        x,y -  coordinates for 2D random walk with given Step_size
%               and Number_of_steps
%               optionally, ballistic components V_x and V_y (in terms of
%               distance traveled per step)
%               if flag random_V is set to true, code adds random velocity
%               (from standard normal distribution) to both x and y
%
if nargin == 1
    %if only n_steps is defined
    step_size = 1;
elseif nargin == 2
    %only n_steps and step_size are defined
    %no balistic movement and no random velocity
    v_x = 0;
    v_y = 0;
    random_v = false;
    constrain_r = 0;
elseif nargin == 5
    constrain_r = 0;
end;

%initialize x and y arrays
x = zeros(1,n_steps);
y = x;
%the initial x and y are always 0;
%generate x and y for elements 2:n_steps
for i = 2 : n_steps
    satisfactory = false;
    while ~satisfactory
        angle = random('unif',0,2*pi);
        if random_v
            random_V_x = random('normal',0,1);
            random_V_y = random('normal',0,1);
        else
            random_V_x = 0;
            random_V_y = 0;
        end;
        x(i) = x(i-1) + step_size * cos(angle)+ v_x + random_V_x;
        y(i) = y(i-1) + step_size * sin(angle) + v_y + random_V_y;
        %check if constrained radius is non-zero
        if constrain_r > 0
            %check if displacement fits within constrained region
            r = sqrt(x(i)^2 + y(i)^2);
            if r>constrain_r
                %repeat displacement generation until constrain condition is
                %satisfied
                satisfactory = false;
            else
                satisfactory = true;
            end;
        else
            %if constrain radius is 0, do not constrain the particle;
            satisfactory = true;
        end;
    end;
end;

