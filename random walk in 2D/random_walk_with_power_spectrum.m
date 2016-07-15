%% Generates random walk and plots the power spectrum of the trajectory
clear all;
close all;
clc;

N_steps = 10000;
N_particles = 1;
msds = zeros(N_particles,N_steps);
sigmas = zeros(N_particles,N_steps);

for i = 1 : N_particles
    
[x,y] = generate_random_walk(N_steps,1,0,0,false,100);

figure (1);
plot(x,y,'-b');

[m,s] = get_msd_from_trajectory(x,y);

msds(i,2:end)= m;
sigmas(i,2:end) = s;

time = 1:1:N_steps;

%[t_2x,m_x,sx_2] = get_msd_from_x_t(x',time');
%[t_2y,m_y,sy_2] = get_msd_from_x_t(y',time');
%m_2 = m_x + m_y;

figure(2);
subplot(2,1,1);
h = errorbar(time,msds(i,:),sigmas(i,:),'o','LineStyle','none');
%errorbar_tick(h,300);
hold on;

subplot(2,1,2);
h = errorbar(time - time(1),msds(i,:),sigmas(i,:),'.');
%errorbar_tick(h, 3000);
hold on;

%plot a line of slope 1 for reference
plot(time,time,'-b');

set(gca,'Xscale','log');
set(gca,'Yscale','log');
axis tight;
end;

msd = mean(msds,1);
t = 0:1:N_steps-1;


figure(3);
loglog(t,msd,'og');
figure(4);
lt = log(t(2:end));
lm = log(msd(2:end));
plot(lt,lm,'og');
fitobj = fit(lt',lm','poly1');
coeffs = coeffvalues(fitobj);

fprintf('Slope of the MSD graph is %2.1f.\n',coeffs(1));

hold on;
plot(fitobj,'-r');
hold off;
%% Use fourier transoform routine
sampling_f = 1000;
[fc_x, d_x] = fourier_transform_routine(x',sampling_f);
[fc_y, d_y] = fourier_transform_routine(y',sampling_f);
