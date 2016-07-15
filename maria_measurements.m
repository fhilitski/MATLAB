%% initialize variables
%data has to be added manually
birthday = datenum('13-Jun-2014');
measurement_dates = {'15-Aug-2014','26-Jun-2014','18-Jun-2014','18-Sep-2014','22-Sep-2014','17-Oct-2014'...
    '15-Dec-2014','27-Feb-2015','24-Mar-2015','18-Jun-2015','20-Jul-2015','25-Aug-2015','14-Dec-2015','9-Jun-2016'};
measurement_weight = [5.04, 3.54, 3.35, 5.78, 5.87, 6.73, 7.90, 8.82, 9.28, 9.52, 9.69,9.996,10.796,12.61]; %in kg
measurement_height = [0.58, 0.53, 0.527, NaN, NaN, 0.635, 0.673, NaN, 0.743, 0.775, NaN,0.80,0.864,0.889]; %in m
measurement_head_circ = [0.38, 0.35, 0.345, NaN, NaN, 0.40, 0.42, NaN, 0.435, 0.45, NaN,0.44,0.46, NaN]; %in m %head circumference

today = datenum(date);
days_old = datenum(measurement_dates)-birthday;


%%
close all;
h_fig_w = figure;
h_weight = plot(days_old, measurement_weight, 'o');
title('Maria''s weight gain');

h_fig_h = figure;
h_height = plot(days_old, measurement_height, 'o');
title('Maria''s height');

h_fig_head = figure;
h_head = plot(days_old, measurement_head_circ, 'o');
title('Maria''s head circumference');


%%
h_fig = figure;
normalized_w = measurement_weight / max(measurement_weight);
normalized_h = measurement_height / max(measurement_height);
h_plot_n = plot(days_old, normalized_w, 'o', days_old, normalized_h,'o');
legend('Normalized weight','Normalized height');

