%%  Files to track

clc;
clear all;
close all;

Name = 'cropped\run05.tif';
Namec = 'cropped\run05cal.tif';

%make sure to change ID and day variables for different runs
ID = 'run05_';
Day = '12102012';

%% Start tracking
disp(['Tracking ' Name '...']);

fcoord = BeadTrackPar(Name);

total_tracked = fcoord(end,4);
disp(['Tracked ' num2str(total_tracked) ' beads.']);

%% Remove unwanted tracks
%here we might need to add code to remove unwanted tracks
%
% fc = fcoord;
%fcoord = fcoord(5323:end,:);
%fcoord(:,4) = fcoord(:,4) - 1;
% fcoord1 = fc(1:5415,:);
% fcoord2 = fc(10831:end,:);
% fcoord2(:,4)=fcoord2(:,4) - 1;
% fcoord = cat(1,fcoord1,fcoord2);
% fcoord = fc(1:9732,:);
%

%% Proceed with tracking

%plot histogram of modulus of coordinates.
%We strive for sub-pixel resolution, thus coordintes should better be
%fractional numbers, NOT integers.
%Thus, histogram on the next figure should be approximately flat.
%

h_f1 = figure;
hist(mod(fcoord(:,2),1));

%plot_bead_histograms(fcoord);
%% Sorting tracks for convenience

%save beads coordinates in a separate variable
eval([ID Day '_beads_xy = fcoord;']);

traj = SortBeads(fcoord);
disp('X1        Y1      X2      Y2');
disp(num2str(traj(1,:)));

%% Introduce scale ...


%t1=input('Enter Trap1 index: ');

t1 = 1; %Trap 1 X coord is in the column #1
%t2=input('Enter Trap2 index: ');
t2 = 3; %Trap 2 X coord is in the column #3

%Introduce magnification scale to convert pixels into nanometers
mag = 0.115/2.5; % microns per pixel

%% Evaluating separation between the beads



eval([ID,Day,'=traj;']);
eval(['R',ID,Day,'=magR(traj,t1,t2) * mag;']); %Separation of two beads in
disp('Separation Measured Succesfully!');

%plot separation vs time
current_fig = get(0,'CurrentFigure');
h_fig = figure(current_fig + 1);
eval(['h_p1 = plot (R',ID,Day,');']);
axis tight;
set(h_p1,'LineStyle','none','Marker','.','Color','red');

%% Track calibration beads

disp(['Tracking ' Namec '...']);

fcoord = BeadTrackPar(Namec);

total_tracked = fcoord(end,4);
disp(['Tracked ' num2str(total_tracked) ' beads.']);

%% Remove unwanted tracks
%here we might need to add code to remove unwanted tracks

% fc = fcoord;
% fcoord1 = fc(1:5580,:);
% fcoord2 = fc(8652:end,:);
% fcoord2(:,4)=fcoord2(:,4) - 1;
% fcoord = cat(1,fcoord1,fcoord2);

%% Proceed with tracking

%plot histogram of modulus of coordinates.
%We strive for sub-pixel resolution, thus coordintes should better be
%fractional numbers, NOT integers.
%Thus, histogram on the next figure should be approximately flat.
%
fh = figure(gcf+1);
hist(mod(fcoord(:,2),1));

%% Sorting tracks for convenience

%save calibration coordinates in a separate variable
eval([ID Day '_beads_xy_cal = fcoord;']);

% plot_bead_histograms(fcoord);

traj_cal = SortBeads(fcoord);
disp('X1        Y1      X2      Y2');
disp(num2str(traj_cal(1,:)));

%% next...


%t1=input('Enter Trap1 index: ');

t1 = 1; %Trap 1 X coord is in the column #1
%t2=input('Enter Trap2 index: ');
t2 = 3; %Trap 2 X coord is in the column #3

%Introduce magnification scale to convert pixels into nanometers
mag = 0.115/2.5; % microns per pixel

%% Evaluating separation between the beads


eval([ID,Day,'calibration','=traj_cal;']);
eval(['R',ID,Day,'cal = magR(traj_cal,t1,t2) * mag;']); %Separation of two beads in
disp('Separation Measured Succesfully!');

figure(h_fig);
hold on;
eval(['h_p2 = plot (R',ID,Day,'cal);']);
axis tight;
set(h_p2,'LineStyle','none','Marker','.','Color','blue');
hold off;

%% Perform comparison of distributions to get potential?

eval(['[U',ID,Day,',Ru',ID,Day,'] = Umbrella(R',ID,Day,',R',ID,Day,'cal,0.02);']);

%% Save all data variables

save_filename = ['data\' ID(1:end-1) 'data'];
if (exist([save_filename '.mat'], 'file')==0)
    
    save(save_filename);
    disp([save_filename '.mat saved successfully.']);
else
    
    disp('File already exists! Check ID...');
end;

close all;



