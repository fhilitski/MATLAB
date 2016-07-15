function[data]=blogf(rawPSD)
%BLOGF performs logarithmic binning of power spectrum rawPSD
%
%data is an array of size n x 4, where n us the number of resulting
%frequency bins;
%data(i,1) - mean frequency of the data in the bin; 
%data(i,2) - mean value of power spectrum for the bin;
%data(i,3) - standard deviation of power spectrum for the bin;
%data(i,4) - starting frequency of bin (left edge of the bin in the
%            histogram);
%data(i,5) - total number of data points averaged for the bin;
%

fraw=rawPSD(:,1);
praw=rawPSD(:,2);

%sometimes initial frequency returned by power spectrum function is 0, and
%this have to be avoided.

if (fraw(1) == 0 )
    bs = floor(log10(fraw(2)));
else
    bs = floor(log10(fraw(1)));
end;
be = ceil(log10(fraw(end)));
fb =[0];

for i=bs:be
    f1=10^(i);
    df=f1;
    for j=0:8
        fb(end+1)=f1+j*df;
    end
end;

[c,x]=histc(fraw,fb);
fd=find(c==0);
fb(fd)=[];

% %temp. histogram
% figure(100);
% hist(fraw,fb);

j=1;
fc=0;
k=0;
PS=[];
data=[];

for i=1:length(fraw)
    
    f_current = fraw(i);
    
    bin_start = fb(j);
    
    if j < length(fb)
        bin_end = fb(j+1);
    else bin_end = max(fraw); %last bin is from the edge of the bins to infinity
    end;
    
    if ( (f_current >= bin_start) && (f_current < bin_end) )
        %if current frequency value falls inside current bin
        fc = fc + f_current;
        PS(end+1) = praw(i);
        k=k+1;
    else
        %if we need to switch bins,
        %first, record data for previous frequency bin
        data(j,1)= fc/k;
        data(j,2)= mean(PS);
        data(j,3)= std(PS)/sqrt(length(PS));
        data(j,4)= bin_start;
        data(j,5)= k;
        
        %then switch bin index
        j = j+1;
        
        %and start analyzing new frequency
        fc = f_current;
        PS = [];
        PS(1) = praw(i);
        k = 1;
        
    end
end
end

