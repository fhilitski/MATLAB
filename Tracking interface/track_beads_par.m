function tracked_coords = BeadTrackPar(Name)
%BeadTrack with parallel computation mode
%INPUTS 
% Name - string filename of multipage tiff;

%   Detailed explanation goes here
%Added Parallel Functionality 7-27-2012 Andy Ward 
info = imfinfo(Name);
num_images = numel(info);

param.mem=100;
param.quiet=1;
param.good=round(num_images/2);
param.dim=2;
coordstruct={};
proceed = 0;


%get figure handles
h_f1 = figure;
h_f2 = h_f1+1;
h_f3 = h_f2+1;

while (proceed == 0)
    a=double(imread(Name,1,'Info',info));
    %bq=input('Enter Bit Level (default is 8 bit): ');
    bq = 8;
    
    
    addbuff = 0;
    % addbuff=input('add buffer pixels? (1=Y,0=N) ');
    % adding buffer pixels is optional.
    if addbuff==1
        dimms=size(a);
        buff=zeros(300,dimms(2))+mean(mean(a));
        a=[buff;a];
    else 
        buff = 0;
    end;
    
    %plot first frame of the original image
    figure(h_f1); imagesc(a);colormap('gray');axis image; title ('Original Image');
    
    Inv = 1;
    disp('Inverting image');
    if Inv==1
        a=2^bq-a;
    end
    
    %plot inverted image
    figure(h_f2); imagesc(a);colormap('gray');axis image; title ('Inverted Image');
    
    B = input('Enter Max Bandpass setting (usually twice the size of bead in pixels):');
    min_bp = input('Enter Max Bandpass setting (usually 1):');
    b = bpass(a,min_bp,B);
    
    %show image after bandpass
    figure(h_f3); imagesc(b);colormap('gray');axis image; title ('Bandpass Image');
    
    %if image looks ok, proceed = 1 and we go to line 48, else go back to
    %the beginning of the loop to optimize bandpass settings.
    
    proceed = input('Check your image? How is the filtering? \n If it works, then proceed (1), otherwise repeat (0): ');
end;

    
    coord=[];
    fcoord=[];
    j=1;
    cnt=[];
    
    FilA = 0;
    FilA = input('Filter Particles by Axes: ');
    if FilA==1
        Ax = input('Enter Experimental Coordinate which is ~fixed (1=X  2=Y) : ');
        Center = input('Enter Center coord of beads: ');
    else 
        Ax = 0;
        Center = 0;
    end
    
    maxIntensity = max(max(b));
    disp(['Max Intensity of the frame: ' num2str(maxIntensity)]);
    M=input('Enter Max Intensity Multiplier (recommended value 0.3): ');%what?
    
    %Bc=B+5;
    
    %Bsmall is always pixel size, i.e. 1, filtering out pixel noise
    %Bsmall=input('Enter Band Pass Small Setting:  ');
    Bsmall = min_bp;
    
    %h_bar = waitbar(0,'Tracking Beads');
    %progress_counter = 0;
    disp('Start tracking...');
    
    parfor i=1:num_images
       % waitbar(i/num_images)
      
      if (mod(i,100)) == 0
            disp(['Tracking frame: ' num2str(i) '... ']);
      end;
      
        a = double(imread(Name,i,'Info',info));
        if addbuff == 1
            a = [buff;a];
        end
        
        if Inv == 1
            a = 2^bq-a;
        end
        
        b=bpass(a,Bsmall,B);
        maxIntensity = max(max(b));
        
        thresh=maxIntensity*M;
        
        
        pk=pkfnd(b,thresh,B);
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
        
        cnt=cntrd(b,pk,B,0); %this gets subpixel resolution from 
                             %preliminary results of pkfind
        
        s=size(cnt);
        
        
        if s(1)<2 
            %only continue if two beads are found
            coordstruct{i}=0;
            continue
        end
        coordstruct{i}=cnt(:,1:2); 
        %coord(j:j+length(cnt(:,1))-1,1:2)=cnt(:,1:2);
        %coord(j:j+length(cnt(:,1))-1,3)=i;
        %j=j+length(cnt(:,1));
    end
    
    %Sort through the structure coordstruct to organize coordinates in  the
    %format that track likes them.
    coord=[];
    for j=1:length(coordstruct)
        tempcoord=coordstruct{j};
        tempcoord(:,3)=j;
    coord=[coord;tempcoord];

    end
    
    disp('Tracking Is Complete!');
    %close(h_bar);
    close(h_f1);
    close(h_f2);
    close(h_f3);
    tracked_coords = track(coord,20,param);
end

