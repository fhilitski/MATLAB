function[c]=PegCon(Mw,cp)
%Outputs concentration of polymer based on Molecular Wegith (Mw)[gm/mole]
%and concentration cp in [mg/mL] 
%Output units are particles/nm^3

%convert  cp from mg/ml into g/m^3
conc = cp*1000;
c = conc*(6.02*10^23/Mw)/10^(27);

end