%% Define parameters
%newtons equation for bead in the optical trap is 
% a = v' =  - 2*\gamma*v+ k*x^2
% v = dx/dt

% where 
% \gamma = 6 \pi \eta a / m
% k = stifness / m 
clear all;
clc;

gamma = 1;  
k = 1;
params = [k, gamma];

%initial conditions
x0 = 0;
v0 = 1;

%timespan
t = 0:0.01:10;

funct = @Newtons_equations;


[T,solution] = ode45(funct,t,[v0,x0]);

v = solution(:,1);
x = solution(:,2)
figure(1);
hold on;
plot(T,v,'or');
plot(T,x,'ob');
hold off;


