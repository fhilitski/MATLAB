clc;
clear all;
close all;
%%
c = [5; 10; 25; 38; 50; 500; 1000]; %motor concentrations
buckling_speed = [0.28; 0.32; 0.21; 0.09; 0.08; 0.07; 0.07]; %for KSA in um/s
extension_speed = [0.14; 0.21; 0.17; 0.18; 0.12; 0.04; 0.07]; %for KSA in um/s
buckling_error = [NaN; 0.12; 0.18; 0.01; 0.04; 0.04; 0.02]; %for KSA in um/s
extension_error = [0.04; 0.09; 0.02; 0.03; 0.06; 0.01; 0.02]; %for KSA in um/s
%buckling_measurements = [1; 7; 1; 3; 17; 2; 11]; %for KSA - number of measurements
% extension_measurements = [3; 6; 5; 3; 11; 2; 11]; %for KSA - number of measurements
% force_measurement = [1.5; 0.8; NaN; NaN; 1.1; NaN; 0.52]; %for KSA in pN
% force_error = [0; 0.321455025; NaN; NaN; 0.593295879; NaN; 0.083666003]; %for KSA in pN
import_bundle_speeds;
plot_k365 = false;
buckling_speed_k365 = [NaN; NaN; 0.29; NaN; 0.24; 0.18; 0.12]; %for K365-SA in um/s
buckling_error_k365 = [NaN; NaN; 0.08; NaN; 0.05; 0.05; 0.03]; %for K365-SA in um/s
extension_speed_k365 = [NaN; NaN; 0.19; NaN; 0.30; 0.16; 0.11]; %for K365-SA in um/s
extension_error_k365 = [NaN; NaN; 0.11; NaN; 0.07; 0.07; 0.01]; %for K365-SA in um/s
%buckling_measurements_k365 = [0; 0; 2; 0; 3; 2; 2]; %for KSA - number of measurements
%extension_measurements_k365 = [0; 0; 0; 0; 4; 3; 2]; %for KSA - number of measurements
%%
%create_speed_plot([c c],[extension_speed extension_speed_k365]...
%    ,[extension_error extension_error_k365],[extension_error extension_error_k365]);

create_speed_plot([c c],[buckling_speed buckling_speed_k365]...
    ,[buckling_error buckling_error_k365],[buckling_error buckling_error_k365]);


%%
h_fig = figure;
h_buckling_plot = errorbar(c,buckling_speed,buckling_error,'o');
hold on;
h_extension_plot = errorbar(c, extension_speed, extension_error, 'o');
if plot_k365
    h_buckling_plot_k365 = errorbar(c, buckling_speed_k365, buckling_error_k365, 'o');
    h_extension_plot_k365 = errorbar(c, extension_speed_k365, extension_error_k365, 'o');
end;
set(gca,'XScale','log');
prettify_plot(h_fig);
if plot_k365
    legend({'KSA Buckling'; 'KSA extension'; 'K365 buckling'; 'K365 extension'});
else
    legend({'KSA Buckling'; 'KSA extension';});
end;
xlabel ('K365 concentration (nM)');
ylabel ('Speed (\mum s^{-1})');
xlim([9,1200]);


%%
h_force_fig = figure;
h_force_plot = errorbar(c, force_measurement, force_error,'o');
set(gca,'XScale','log');
prettify_plot(h_force_fig);