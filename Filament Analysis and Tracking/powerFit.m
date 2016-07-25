function [fitresult, gof] = powerFit(t, l)
%CREATEFIT(T_GOOD,L_GOOD)
%  Create a power-law fit l = a*t^b
%
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 29-Apr-2014 11:21:27

[xData, yData] = prepareCurveData( t, l );

% Set up fittype and options.
ft = fittype( 'a*x^b', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( ft );
opts.Display = 'Off';
opts.Lower = [-Inf -Inf];
opts.Upper = [Inf Inf];
opts.StartPoint = [0.795199901137063 0.186872604554379];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
figure( 'Name', 'untitled fit 1' );
h = plot( fitresult, xData, yData );
legend( h, 'L_good vs. t_good', 'untitled fit 1', 'Location', 'NorthEast' );
% Label axes
xlabel( 't_good' );
ylabel( 'L_good' );
grid on


