function[Vdep]=depletion(c,R_g,R_c, Polymer_MW, d)
%Gives the Depletion Energy/micron in units of kT/micron, between two
%colloidal cylinders (microtubules)of radius R_c; suspended in solution of
%micromolecules of with gyration radius R_g; and with center-to center
%distance d; polymer concentration cp is in mg/ml.

%PegConc returns particles/nm^3 and Acs returns nm^2, so we have to
%multiply by 1000 to convert resfogirult from kT/nm to kT/mum
Vdep = - PegCon(Polymer_MW, c)* Acs(d,R_c,R_g) * 1000;


end
