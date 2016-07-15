function [ output_args ] = plot_msd( x, t, errors )
%PLOT_MSD Plots MSD of particle
%   PLOT_MSD (x,t,e)


h_fig = figure;
h_msd = errorbar(t, x, errors, '.r');
errorbar_tick(h_msd);
hold on;
plot(t,t,'-.k');
%plot(t,t.^2,'--k');
hold off;
set(gca,'XScale','log');
set(gca,'YScale','log');
axis tight;

title('LOG-LOG MSD plot of bundle position vs time');
xlabel('Time, frames');
ylabel('Bundle MSD, pixel^2');
legend('Bundle MSD','Slope 1 (diffusion)','Slope 2 (ballistic)');


