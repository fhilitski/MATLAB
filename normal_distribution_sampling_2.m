% Code for The Data Incubator challenge - exercise 1
%% initialize workspace
clc;
close all;
clear all;

%starting position
start_key = 0;

%sum of keys
s = 0;

%intialize random generator
rng('shuffle');

%initialize parallel computation
par_pool = gcp;

%% Part 1. After T=10 moves, what is the expected value of the quantity Smod10 ?
%let's run multipe trials n_moves, histogram the quantity and obtain
%statistics about it.
%the number of trials to run:
n_trials = 10000;
h_fig_s = figure;
title('S');
h_fig_q = figure;
title('S mod 10');


for j = 1:1000:n_trials
    disp(['Starting simulation with ' num2str(j) ' trials']);
    %initialize variable to store s and quantity of interest(s mod 10)
    all_s = zeros(1, j);
    all_q = zeros(1, j);
    
    parfor i = 1:1:j
        n_moves = 10;
        s = randn;
        quantity = mod(s,10); %that's the quantity we are looking for
        
        all_s(i) = s;
        all_q(i) = quantity;
        if mod(i,1000) == 0
            %disp('.');
        end;
        %fprintf('Trial %i finished\n', i);
    end;
    
    %fprintf('\n');
    %disp('Done with trials!');
    %figure; histogram(all_q); title('Histogram of s mod 10');
    %figure; histogram(all_s); title('Histogram of s');
    
    mean_q = mean(all_q);
    mean_s = mean(all_s);
    
    std_q = std(all_q);
    std_s = std(all_s);
    
    
    fprintf('<s>: %14.12f \n', mean_s);
    fprintf('<(s mod 10)>: %14.12f \n', mean_q);
    fprintf('std(s): %14.12f \n', std_s);
    fprintf('std((s mod 10)): %14.12f \n', std_q);
    
    figure(h_fig_s);
    errorbar(j,mean_s,std_s,'ob');
    hold on;
    drawnow;
    figure(h_fig_q);
    errorbar(j,mean_q,std_q,'ob');
    hold on;
    drawnow;

end;
