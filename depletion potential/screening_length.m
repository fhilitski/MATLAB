function [ ls ] = screening_length(c_ions, epsilon)
%screening_length Returns Debye length (electrostatic screening length) 
%   based on molar concentration of monovalent ions, and
%   relative permittivity epsilon;
%   returns electrostatic screening length (debye length) ls in nm;

epsilon_0 = 8.85*10^(-12); %F/m
D = epsilon; %water relative permittivity
kB = 1.38*10^(-23); %J/K at room T

e = 1.61*10^(-19); %coulomb;

c = c_ions * 6.02*10^23/(10^(-3)); %concentration in particles per m^3

ls = sqrt(D*epsilon_0*kB*300/(2*e^2*c));

%convert into nm
ls = ls*10^9;
end

