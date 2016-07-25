function[Vrod]=rodpot(lb,b,r,ls)
%calculate the electrostatic repulsion between two cylinders. taken
%from Brenner and Parsegian 1974 paper. is valid as long as the seperation
%is between the surfaces is larger than the screening length
%lb is bjjerem length
%b is the charge spacing
%r is the center to center distance
% ls is the screening length

Vrod=2*lb/(b^2)* besselk(0,r./ls)*1000;% Units are in kT/micron (?)
end
