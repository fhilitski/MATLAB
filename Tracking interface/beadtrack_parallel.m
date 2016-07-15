function tracked_coords = beadtrack_parallel(name,bit_depth, Bsmall, Blarge)
%Bead Tracking with parallel analysis of images
%INPUTS
% Name - filename of tif file with images to track
% bit_depth - bit depth of the image
% Bsmall - small bandpass filter value
% Blarge - large bandpass setting

%Added Parallel Functionality 7-27-2012 Andy Ward

info = imfinfo(name);
num_images = numel(info);

param.mem = 100;
param.quiet = 1;
param.good = round(num_images/2);
param.dim = 2;
coordstruct = {};
proceed = 0;


coord = [];
fcoord = [];
j = 1;
cnt = [];

%This part makes use of possibility to filter particles by axes
%it is desabled for now

% FilA = 0;
% FilA = input('Filter Particles by Axes: ');
% if FilA==1
%     Ax = input('Enter Experimental Coordinate which is ~fixed (1=X  2=Y) : ');
%     Center = input('Enter Center coord of beads: ');
% else
%     Ax = 0;
%     Center = 0;
% end


%Max intenisty threshold M
%All pixels with intensity below 30% of image max intensity will not be
%tracked
M=0.3;

%Bc=B+5;
%Bsmall is always pixel size, i.e. 1, filtering out pixel noise
%Bsmall=input('Enter Band Pass Small Setting:  ');


%h_bar = waitbar(0,'Tracking Beads');
%progress_counter = 0;
disp('Start tracking...');

parfor i=1 : num_images
    % waitbar(i/num_images)
    
    if (mod(i,100)) == 0
        disp(['Tracking frame: ' num2str(i) '... ']);
    end;
    
    a = double(imread(name,i,'Info',info));
    
    % invert the image   
    a = 2^bit_depth-a;
       
    b = bpass(a,Bsmall,Blarge);
    maxIntensity = max(max(b));
    
    thresh = maxIntensity*M;
    
    
    pk = pkfnd(b,thresh,Blarge);
    %pkfind finds peaks in filtered image b
    
    %if we know beads to be along one axis,
    %we only need to consider these peaks,
    %so we are going to look at peaks +/- beads_range pixels
    %from axis
    
    beads_range = 15;
    
    %         if FilA==1
    %             f=find(abs(pk(:,Ax)-Center) < beads_range );
    %             pk=pk(f,:);
    %         end
    
    cnt = cntrd(b,pk,Blarge,0); %this gets subpixel resolution from
    %preliminary results of pkfind
    
    s = size(cnt);
    
    
    if s(1)<1
        %only continue if at least one bead is found
        coordstruct{i}=0;
        continue
    end
    coordstruct{i} = cnt(:,1:2);
    %coord(j:j+length(cnt(:,1))-1,1:2)=cnt(:,1:2);
    %coord(j:j+length(cnt(:,1))-1,3)=i;
    %j=j+length(cnt(:,1));
end;

%Sort through the structure coordstruct to organize coordinates in  the
%format that track likes them.
coord=[];
for j=1:length(coordstruct)
    tempcoord=coordstruct{j};
    tempcoord(:,3)=j;
    coord=[coord;tempcoord];
    
end;

disp('Tracking Is Complete!');

tracked_coords = track(coord,30,param);
end

