function [result] = Newtons_equations(t, values, parameters)
%NEWTONS_EQUATIONS 
%Returns [dv/dt, dx/dt] (acceleration and velocity)
%For given t and values ([v,x]) (velocity and position)
%Parameters are stiffnes and daping [k,gamma] 
%
%Newtons equation for bead in the optical trap is 
% dv/dt =  -\gamma*v+ k*x^2
% dx/dt = v
% where 
% \gamma = 6 \pi \eta a / m
% k = stifness / m 

v = values(1);
x = values(2);

k = 1;
g = 1;

a = -k*x - g*v;

result = [a; v];

end

