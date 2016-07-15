function [fc,D] = Lorentzian_coeffs(freq,power)
%LORENTZIAN_COEFFS Calculates coefficients for lorentzian fit of binned
%power spectrum power with frequencies freq;
%
%[fc,D] = LORENTZIAN_COEFFS(f,p)
%
%Based on:
%K. Berg-Sørensen and H. Flyvbjerg, Review of Scientific Instruments 75, 594 (2004).
%specifically, p.595-596
%
% 

function s = Summation(f,p,l,q)

    s = 0;
   N = length(f);
   for i=1:N
       s = s + (f(i)^(2*l))*(p(i)^q);
   end;       
end


s01 = Summation(freq,power,0,1);
s11 = Summation(freq,power,1,1);
s12 = Summation(freq,power,1,2);
s02 = Summation(freq,power,0,2);
s22 = Summation(freq,power,2,2);

fc = sqrt((s01*s22 - s11*s12)/(s11*s02 - s01*s12));
D = (s02*s22 - s12^2)/(s11*s02 - s01*s12);    


end

