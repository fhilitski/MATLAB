function [mapping_function] = fourier_curve_fit(fit_result)
%[MAPPING_FUNCTION] = FOURIER_CURVE_FIT(FIT_RESULT);
%RETURNS mapping_function -  function handle for the fourier mapping function:
%fit_y_fourier(x) =  a0 + a1*cos(x*w) + b1*sin(x*w) +
%                    a2*cos(2*x*w) + b2*sin(2*x*w) + a3*cos(3*x*w) + b3*sin(3*x*w) +
%                    a4*cos(4*x*w) + b4*sin(4*x*w)
% This function is designed to map displacement onto voltage. 

coeffs = coeffvalues(fit_result);
a0 = coeffs(1);
a1 = coeffs(2);
b1 = coeffs(3);
a2 = coeffs(4);
b2 = coeffs(5);
a3 = coeffs(6);
b3 = coeffs(7);
a4 = coeffs(8);
b4 = coeffs(9);
a5 = coeffs(10);
b5 = coeffs(11);
a6 = coeffs(12);
b6 = coeffs(13);
a7 = coeffs(14);
b7 = coeffs(15);
a8 = coeffs(16);
b8 = coeffs(17);
w = coeffs(18);
mapping_function =  @(x) a0 + a1*cos(x*w) + b1*sin(x*w) + ...
               a2*cos(2*x*w) + b2*sin(2*x*w) + a3*cos(3*x*w) + b3*sin(3*x*w) + ...
               a4*cos(4*x*w) + b4*sin(4*x*w) + a5*cos(5*x*w) + b5*sin(5*x*w) +  ...
               a6*cos(6*x*w) + b6*sin(6*x*w) + a7*cos(7*x*w) + b7*sin(7*x*w) + ...
               a8*cos(8*x*w) + b8*sin(8*x*w);

