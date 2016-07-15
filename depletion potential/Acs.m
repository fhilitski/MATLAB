function area = Acs(d, R_colloid, Rg)
% Returns cross-sectional area of excluded volume overlap between two
% cylindrical objects of radius R_colloid suspended in polymers with radius
% of gyration Rg.
% optput has same units as input R_c and R_g, so if radii are in nm, output
% is in nm^2.

r = R_colloid + Rg;

atn = atan(sqrt(4*r^2-d.^2)./d);
sroot = sqrt(4*r^2 - (d).^2);
c =  atn - (d .* sroot)/(4*r^2);
area = 2*r^2 * c;

end