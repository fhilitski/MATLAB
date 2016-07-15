function [ U_electrostatic ] = electrostatic_energy(d, R_c, sigma,ld, D)
%electrostatic_energy  Electrostatic energy of two charged cylinders in ionic
%solution
%   Colloidal particle with cylindrical symmetry, radius R_c (in nm), 
%   surface charge sigma (in e/nm^2), suspended in ionic solution with Debye length 
%   ld(in nm), and dielectric constant D (dimensionless). Distance between axises of
%   the cylinders is d.
%   Returns potential energy per micron overlap in k_B*T/micron 
%   Approximation of very long cylinders neglects edge effects, 
U_electrostatic = zeros(1,length(d));

N = 1000;
theta = linspace(0,2*pi,N);

for i = 1:length(d)
    r_n = sqrt(R_c^2*sin(theta).^2+(d(i) - R_c*cos(theta)).^2);
    [phi_n, phi_n_V] = cylinder_potential(r_n,R_c,sigma,ld,D);
    phi_tot = sum(phi_n);
    
    %multiplication by 1000 is necessary to convert kT/nm into kT/micron
    U_electrostatic(i) = 2*pi*R_c*sigma*phi_tot/length(theta)*1000;
end;

end

