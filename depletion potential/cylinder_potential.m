function [ phi, phi_V ] = cylinder_potential(r, R_colloid, Sigma, L_screening, Epsilon)
%CYLINDER_POTENTIAL Electrostatic potential of charged cylinder in ionic
%solution
%   Colloidal particle with cylindrical symmetry, radius R_colloid (in nm), 
%   surface charge Sigma (in e/nm^2), suspended in ionic solution with Debye length 
%   L_screeing (in nm), and dielectric constant Epsilon (dimensionless). Distance from axis to
%   the point where potential is calculated, r.
%   Returns potential in k_B*T/e and in volts

epsilon_0 = 8.85*10^(-12); %F/m
D = Epsilon; %water relative permittivity
kB = 1.38*10^(-23); %J/K
T = 300; %K

e = 1.61*10^(-19); %coulomb;

s = Sigma*e*10^18; %convert sigma from e/nm^2 to C/m^2;

phi_V = s*(L_screening*10^-9) /(D*epsilon_0) * besselk(0,r./L_screening) /(besselk(1,R_colloid./L_screening));
%this is potential in volts

%convert volts into k_BT/e
phi = phi_V /(kB * T/ e);


end

