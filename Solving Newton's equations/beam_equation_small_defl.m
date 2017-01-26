function [result] = beam_equation_small_defl(x, values, parameters)
%BEAM_EQUATION_SMALL_DEFL
%Beam equation for small deflections
%RETURNS [d(dy)/dx, dy/dx]
%
%Define following variables for a beam:
% t - coordinate along the beam
% y - deflection from equilibrium
%
%Beam equation has form:
% d2y/dt2 = M(t)
%
%Rewriting this as two first order ODEs
%dy/dt = theta (tan of angle theta between the tangent and x-axis)
%dtheta/dt = M(t)

%these are solutions of the immediately preceeding step
theta = values(1);
y = values(2);

w = 1; %force on the end of the beam
EI = 1; %flexural rigidiy
l = 1; %total length of the beam
m = 1; %mass of the beam
g = 1; %gravity acceleration

%this is the first derivative
dydt = theta;

%equation for a beam with weight w at the end
d2ydt2 = w/EI*(l - x);

%equation for a of mass m
%d2ydt2 = m*g/(2*EI*l)*(x^2 + l^2 - 2*x*l);

result = [d2ydt2; dydt];
end

