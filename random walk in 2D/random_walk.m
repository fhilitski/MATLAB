%% Random Walk Simulation

clc;
clear all;
close all;

N_experiments = 10;

theoretical_msds = zeros(N_experiments, 1);
msd_values = zeros(N_experiments, 1);
number_steps = zeros(N_experiments, 1);

for t = 1:N_experiments
    
    N_steps = 1*t;
    N_iterations = 100;
    step_size = 1;
    
    square_displacements = zeros(N_iterations,1);
    single_track = cell(N_iterations,1);
    
    fprintf(1, 'Starting parallel simultaion: %d iterations \n', N_iterations);
    fprintf(1, 'Steps: %d \n', N_steps);
    
    parfor j = 1 : (N_iterations)
        
        if (mod(j,1) == 0 )
            fprintf(1,'Currently doing %i step out of %d\n', ...
                j, N_iterations);
        end;
        
        
        
        
        x_0 = 0;
        y_0 = 0;
        % start = [x_0, y_0];
        
        traj = zeros(N_steps + 1,2);
        traj(1,:) = [x_0,y_0];
        
        x = x_0;
        y = y_0;
        
        
        for i = 1 : N_steps
            
            angle = random('unif',0,2*pi);
            x = x + step_size * cos(angle);
            y = y + step_size * sin(angle);
            
            traj(i+1,:) = [x, y];
            
        end;
        
        single_track{j} = traj;
        
        square_displacements(j) = ((x - x_0)^2 + (y - y_0)^2);
        
    end;
    
    if (gcf == 1)
        figure(gcf);
    else
        figure();
    end;
    
    traj = single_track{1};
    plot(traj(:,1),traj(:,2),'-b');
    hold on;
    
    start = traj(1,:);
    endpoint = traj(end,:);
    
    plot(traj(1,1),traj(1,2),'sg');
    plot(traj(end,1),traj(end,2),'sr');
    
    line_x = [start(1), endpoint(1)];
    line_y = [start(2), endpoint(2)];
    line(line_x,line_y,'Color','red','LineWidth',2);
    axis auto;
    hold off;
    
    figure();
    plot(square_displacements,'or');
    axis auto;
    
    msd = mean(square_displacements);
    fprintf(1,'Exp. MSD: %f \n',msd);
    fprintf(1,'Theoretical. MSD: %f \n',step_size^2*N_steps);
    
    theoretical_msds(t) = step_size^2*N_steps;
    msd_values(t) = msd;
    number_steps(t) = N_steps;
    
end;

figure(gcf+1);
loglog(number_steps,msd_values,'ob',number_steps,theoretical_msds,'sr');
grid on;
axis tight;
legend('Simulation','1D diffusion');



