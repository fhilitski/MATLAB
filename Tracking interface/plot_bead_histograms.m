function [h1,h2,h3,h4] = plot_bead_histograms( fcoord )
%plot_bead_histograms Plots histograms of tracked beads
%   Input coordinates fcoord

bead_1 = find(fcoord(:,4) == 1);
bead_1_coords = fcoord(bead_1,:);

bead_2 = find(fcoord(:,4) == 2);
bead_2_coords = fcoord(bead_2,:);

%% Plot histograms of bead positions and find stds

h1 = figure(gcf+1);
hist(bead_1_coords(:,1),100);
bead_1_x_std = std(bead_1_coords(:,1));
h = findobj(gca,'Type','patch');
set(h,'FaceColor','r','EdgeColor','w');
title(['X_{1} histogram; \delta(x_1) = ' num2str(bead_1_x_std)]);

h2 = figure(gcf+1);
bead_2_x_std = std(bead_2_coords(:,1));
title(['X_{2} histogram; \delta(x_2) = ' num2str(bead_2_x_std)]);
hist(bead_2_coords(:,1),100);
title(['X_{2} histogram; \delta(x_2) = ' num2str(bead_2_x_std)]);


h3 = figure(gcf+1);
hist(bead_1_coords(:,2),100);
bead_1_y_std = std(bead_1_coords(:,2));
h = findobj(gca,'Type','patch');
set(h,'FaceColor','r','EdgeColor','w');
title(['Y_{1} histogram; \delta(y_1) = ' num2str(bead_1_y_std)]);

h4 = figure(gcf+1);
hist(bead_2_coords(:,2),100);
bead_2_y_std = std(bead_2_coords(:,2));
title(['y_{2} histogram; \delta(y_2) = ' num2str(bead_2_y_std)]);
end;

