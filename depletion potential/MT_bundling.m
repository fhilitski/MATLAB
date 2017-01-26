%% This script calculates enery of MT bundles and longitudinal depletion force
clc;
clear all;
close all;

r = 12.5; %radius of mts in nm;

mw = 20000; %molecular weight of PEG in Kda
f = 2/sqrt(pi); 
r_g = f*8; %nm, PEG radius of gyration;


c = 1.0; % %mg/ml, PEG concentration
cK = 0.100; % Molar, concentration of K+ ions

epsilon = 81;%dielectic constant of water
T = 300; % K room tem

%sigma = 0.48; %electron charges/nm^2 of surface - from protein sequence
sigma = 0.23; %electron charges/nm^2 of surface - measured by electrophoresis experiments

%% Calculate depletion potential
figure(1);

min_d = 2*r;
max_d = 2*(r+r_g);

d = linspace(min_d, max_d, 1000);
U_depletion = depletion(c,r_g,r,mw,d);

h_depletion = plot(d,U_depletion,'og');

%% Calculate electrostatic repulsion

U_es = electrostatic_energy(d,r,sigma,screening_length(cK,epsilon),epsilon);

hold on;
h_es = plot(d,U_es,'or');

%% Total potential
U_total = U_es + U_depletion;
h = plot(d,U_total,'ok');

[U_min, index] = min(U_total);
d_min = d(index);

ylim([U_min - 20, U_min+100]);
xlim([d_min-3, d_min+2]);

%% Label axis in figure 1
xlabel('\fontsize{24} Center-to-center spearation, \mum');
ylabel({'\fontsize{24} Bundle energy '; 'per unit overlap, k_BT/\mum'});
legend_string = {'Depletion','Electrostatic','Total'};
legend(legend_string,'FontSize',20);
title('\fontsize{28} [PEG] = 1%, [K^+] = 100 mM ');


%% Get depletion strength;
lambda = depletion(c,r_g,r,mw,d_min);
disp(['lambda = ' num2str(lambda)]);

%% Fix [PEG], change [K+]
h_fig = figure(2);
%color = ['b','r','g'];
figure(h_fig);
hold all;


counter = 0;
polymer = [5, 10, 15, 20, 30, 100];
cK = [0.1, 0.5,1:1:50, 100, 200, 300, 400, 500, 650] .* 0.001;

force = zeros(length(polymer),length(cK));
%gr = [];

for c = polymer
    
    counter = counter + 1;
    U_depletion = depletion(c,r_g,r,mw,d);
    %gr(counter) = hggroup;
        
    for i = 1:length(cK)
        conc_K = cK(i);
        U_es = electrostatic_energy(d,r,sigma,screening_length(conc_K,epsilon),epsilon);
        U_total = U_es + U_depletion;
        [U_min, index] = min(U_total);
        d_min = d(index);
        lambda = -(depletion(c,r_g,r,mw,d_min));
        force(counter,i) = lambda;
       
    end;
        %h_plot = plot(cK,force(counter,:),['-o', color(counter)],'MarkerSize',8,'LineWidth',2);
        h_plot = plot(cK,force(counter,:),'-o','MarkerSize',8,'LineWidth',2);
        %set(h_plot,'Parent',gr(counter));
        hold on;
    
end;

%% Add Labels in figure 2
xlabel('\fontsize{24} [K+], Molar');
ylabel({'\fontsize{24} Bundle energy '; 'per unit overlap, k_BT/\mum'});

title({'\fontsize{28} Calculation: ','Depletion Froce vs [K+]'});

% for i = 1 : length(gr)
%     set(get(get(gr(i),'Annotation'),'LegendInformation'),...
%         'IconDisplayStyle','on');    
% end;
legend_string = {'0.5% PEG','1.0% PEG ','1.5% PEG'};
legend(legend_string,'FontSize',20);


%% Fix [K+], change [PEG]
figure(3);
color = ['k','r','g','b'];

counter = 0;
polymer = [0, 5, 10, 15];
cK = [0.1, 100, 300, 500] .* 0.001;

for i = cK
   counter = counter + 1;
   lambda = [];
   for c = polymer
       
        U_depletion = depletion(c,r_g,r,mw,d);
        U_es = electrostatic_energy(d,r,sigma,screening_length(i,epsilon),epsilon);
        U_total = U_es + U_depletion;
        [U_min, index] = min(U_total);
        d_min = d(index);
        lambda = [lambda, -(depletion(c,r_g,r,mw,d_min))];
       
        hold on;
    end;
     plot(polymer,lambda,['o', color(counter)]);
end;




