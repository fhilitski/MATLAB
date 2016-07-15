function [ fc ] = msd_plot_routine(separation)
%MSD_PLOT_ROUTINE plots mean square displacement from given time-series,
%on log-log scale

% separation - time series of a given variable (first was designed for
% separation between beads, hence the name);

% returns 0 if succesfull, -1 otherwise;


%we need time between images, i.e. frame rate of the aquisition
%and number of images to analyze
N = length(separation);
prompt = {'Enter number of frames per second: '; 'Number of images to analyze:'};
dialog_title = 'MSD options...';
user_input = inputdlg(prompt, dialog_title);

if (length(user_input) > 0)
    
    frames_per_sec = str2num(user_input{1});
    N_frames = str2num(user_input{2});
    
    if (~isempty(frames_per_sec) && ~isempty(N_frames) )
        
        if (N_frames > N) 
            N_frames = N;
        end;
        
        separation = separation(1:N_frames);
        %dt = time between the images
        dt = 1/frames_per_sec;
        
        %let us create time vector
        %first images correspond to t1, t2, t3, so
        %delta_t = 0, t2-t1, t3-t1,....,tN-t1
        %        = (0,1,2,3,...N-1)*dt
        N = length(separation);
        time_vector = (1:1:(N-1))*dt;
        
        
        [msd, errors] = get_msd_from_trajectory (separation, 0*separation);
                
        %first, plot MSD
        h_fig = figure;
        h_msd = errorbar(time_vector, msd, errors, 'or');
        errorbar_tick(h_msd);
        hold on;
        
        %plot lines of slope 1
        N_points = min(20,N_frames-1);
        plot(time_vector(1:N_points),(msd(1)/time_vector(1))*time_vector(1:N_points),'--k');
        hold off;
        set(gca,'XScale','log');
        set(gca,'YScale','log');
        axis tight;
        fc = 0;     
    end;
    
else
    fc = -1;
end;


