%Please import manual tracking data into variable data!
clc;
clearvars -except data dt time;
close all;


length_per_pixel = 0.13;
h_figure = figure;
hold on;
colors = ['r','g','b','c','m','k'];
color_index = 0;
v = [];
for i = 1:max(data(:,1))
    track1 = data(find(data(:,1) == i),:);
    if (size(track1) > 0)
        x = track1(:,3)-track1(1,3);
        y = track1(:,4)-track1(1,4);
        d = length_per_pixel  * sqrt(x.^2 + y.^2);
        slice_number = track1(:,2);
        %t = 0.118*(track1(:,2)-track1(1,2)); %this is calculation of time
        %with known time step
        %in case time steps are saved in separate variable dt:
        t = dt(slice_number)./1000; %this is time in seconds
        t = t';
        
        color_index = color_index + 1;
        if (color_index > size(colors))
            color_index = 1;
        end;
        h_plot = plot(t,d,['.-' colors(color_index)]);
        hAnnotation = get(h_plot,'Annotation');
        hLegendEntry = get(hAnnotation','LegendInformation');
        set(hLegendEntry,'IconDisplayStyle','on');
        %set(h_plot,'DisplayName','Tracked Dot position');
        
        [fit,gof] = linearFit(t,d);
        h_fit = plot(fit,colors(color_index));
        coeffs = coeffvalues(fit); 
        v = [v, coeffs(1)];
        %turn off annotation for the fit;
        hAnnotation = get(h_fit,'Annotation');
        hLegendEntry = get(hAnnotation','LegendInformation');
        set(hLegendEntry,'IconDisplayStyle','on');
    else
        disp(['Track ' num2str(i) ' has zero size!']);
    end;
end;
axis tight;
legend('Tracked dot position', 'Line fit');
xlabel('Time (s)');
ylabel('Displacement (um)');
title('Tracking TIRF videos');
prettify_plot(h_figure);

fprintf('Fitted velocity: %4.3f um/s\n', v);
