h_fig = figure;
h = animatedline('LineStyle','-','Marker','.');
x = trap_dist_fd_curve;
y = smooth(total_force_fd_curve,force_smoothing);

for k = 1:length(x)
    addpoints(h,x(k),y(k));
    drawnow
end