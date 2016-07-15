function [ output_args ] = prettify_plot(figure_handle)
%PRETTIFY_PLOT makes the plot look nice! 
%   It makes title font HUGE, axes labels LARGE, ticks acceptable for
%   Zvonimir's critical eye, axes TIGHT and etc. 
%   Make sure your plot is active figure (click on it before running!) or
%   pass the figure handle

if (nargin == 0)
    figure_handle = gcf;
end;
figure(figure_handle);

%get figures axes;
h_axes = get(figure_handle, 'Children');
if (length(h_axes) >= 2)
    h_axes = h_axes(2);
end;

%Set title font;
h_title = get(h_axes, 'Title');
set(h_title, 'FontName', 'Calibri');
set(h_title, 'FontSize', 35);

%Set axes labels' font to something large
h_label = get(h_axes, 'YLabel');
set(h_label, 'FontName', 'Calibri');
set(h_label, 'FontSize', 30);

h_label = get(h_axes, 'XLabel');
set(h_label, 'FontName', 'Calibri');
set(h_label, 'FontSize', 30);

%Ticks labels should also be visible, right. Zvonimir like-y...
set(h_axes, 'FontName', 'Calibri');
set(h_axes, 'FontSize', 20);

end

