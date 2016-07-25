function [ lb ] = bjerrum_length(epsilon, T)
%screening_length Returns Bjeerum length
%   (lengt at which energy between two charges is equal to thermal energy) 
%   based on relative permittivity epsilon,
%   and temperatue T in Kelvin;
%   returns Bjerrum length lb in nm;

epsilon_0 = 8.85*10^(-12); %F/m
D = epsilon; %water relative permittivity
kB = 1.38*10^(-23); %J/K at room T

e = 1.61*10^(-19); %coulomb;

lb = e^2 / (4*pi*epsilon_0*D*kB*T);

%convert into nm
lb = lb*10^9;
end

