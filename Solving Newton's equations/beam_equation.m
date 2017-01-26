function [result] = beam_equation(x, values, parameters)
%BEAM_EQUATION
%Beam equation for all deflections; See BEAM_EQUATION_SMALL_DEFL for small
%deflection approximation. 
%RETURNS [d(dy)/dx, dy/dx]
%For given x and values ([dT/ds,s]) (1/R and contour length)
%Parameters are stiffnes and daping [k,gamma] 
%
%Define following variables for a beam:
% s - coordinate along the beam  - former t
% theta - tangent at a point s  - former y
%
%Beam equation has form:
% d2theta/ds2 = -f(theta) 
%
%Rewriting this as two first order ODEs
%dtheta/ds = c (curvature)      theta == y1, c == y2, y1' = y2
%dc/ds = d2theta/ds = f(theta)  dcds = y2' = f(y1)

%parameters that might go into the equation f(y,t)
w = 1; %force on the end of the beam
EI = 1; %flexural rigidiy
l = 1; %total length of the beam 
m = 1; %mass of the beam
g = 1; %gravity acceleration
F = 1.02*pi^2*EI/l^2; %buckling force


%these are solutions of the immediately preceeding step
c = values(1);
theta = values(2);

%this is the first derivative
dthetads = c;

%this is the second derivative
%dcds = -F/EI*sin(theta); %buckling equation
dcds = -F/EI * sin(theta);
%dcds = -abs(theta); - test case from the example

result = [dcds; dthetads];
end

