%% Define parameters
%solving differential equations
%
%For the first order the equation dy/dt = f(y, t)
%initial conditions y(t0) = t(0)

%Second order ODE (Newton's equation, beam equation)
%d2y/dt2 = f(y, y', t)
%with initial conditions
% y(t0) = y0;
% y'(t0) = v0;


clear all;
%close all;
clc;

%initial conditions
v0 = 0; %initial velocity v(t0) - for Newton's equations
        %initial tangent t(x0) - for beam bending
        %initial curvature (at s = 0) - for beam equations, c = 1/R 
        
y0 = 1; %initial position - for Newton's equations
        %initial deflection y(x0) - for beam bending
        %initial tangent (at s = 0) - for beam equation T

%timespan or filament contour length
%(x and t are used interchangibly)
t = 0:0.01:30;

funct = @Newtons_equations;
%funct = @beam_equation;
%funct = @beam_equation_small_defl;
%funct = @beam_equation_no_approx;

[time,solution] = ode45(funct,t,[v0,y0]);
%time is the solvers t (or x) - it is the variable in the y(t)
v = solution(:,1); % solved velocity y'(t)
y = solution(:,2); % solved position y(t)

h_fig_v = figure;
h_v = plot(time,v);
h_fig_y = figure;
h_y = plot(time,y);

set(h_v,'Marker','.');
set(h_v,'DisplayName','y\prime = y\prime(t)');

set(h_y,'Marker','.');
set(h_y,'DisplayName','y = y(t)');

prettify_plot(h_fig_v);
prettify_plot(h_fig_y);




