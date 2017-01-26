function result = beam_equation_bc(ya, yb)
%BEAM_EQUATION_BC
%Boundary conditions for the buckling beam equation BEAM_EQUATION
%Beam equation has form:
%d2theta/ds2 = f(theta) 
%
%Rewriting this as two first order ODEs
%dtheta/ds = c (curvature)      theta == y1, c == y2, y1' = y2
%dc/ds = d2theta/ds = f(theta)  dcds = y2' = f(y1)

%Example 1.
%Boundary conditions (in a familiar form) set initial angle theta:
%theta(s_i) = theta_i;
%theta(s_f) = theta_f;
%In our y1/y2 notation this becomes:
% ya1 = y1_i;
% yb1 = y1_f;
%Note: a,b are just labels for s_i and s_f, respectively.

%For the MATLAB solver, boundary conditions should have form:
%ya(1) -y1_i = 0;
%yb(1) -y1_f = 0;

%define boundary values y1_i, y1_f:
y1_i = pi/400;
y1_f = -pi/400;
%-----------------------------

%Example 2.
%Pinned boundary conditions assume that the angle is not constrained
%Instead, there is no curvature at the boundary:
%theta'(s_i) = 0;
%theta'(s_f) = 0;
%In our y1/y2 notation this becomes:
% ya2 = 0;
% yb2 = 0;
%Note: a,b are just labels for s_i and s_f, respectively.

%For the MATLAB solver, boundary conditions should have form:
%ya(2) -y2_i = 0;
%yb(2) -y2_f = 0;
r = 0.001;
F = pi^2*1/1^2;
y2_i = -0*F*r*cos(ya(1));
y2_f = -0*F*r*cos(yb(1));


%return boundary conditions:
%result = [ya(1)-y1_i yb(1)-y1_f]; %for clamped boundary (fixed angle)
result = [ya(2)-y2_i; yb(2)-y2_f];  %for pinned boundary (no torque)
end

