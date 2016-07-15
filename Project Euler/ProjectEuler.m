%%Project Euler Problems.
clear all;
clc;

%% Problem 1. Sum of multiples
disp('Running Problem 1.');
result = 0;
  for j=0:1:999
    if ( rem(j,3) == 0  || rem(j,5) == 0)
        result = result + j;
    end;
  end 
  
%% Problem 2. Fibonacci sequence

disp('Running Fibonacci sequence...');
current = 2;
previous = 1;
new = current + previous;
result = 2;

while (new < 4)
    if (rem(new,2) == 0) result = result + new; 
    end;
    previous = current;
    current = new;
    new = previous + current;
end;
disp('Done!');

%% Problem 11. Reading from file and finding sequence with the largest sum.
disp('Running Problem 11...');
clear all;
data = importdata('problem11.txt');

 %%analyze strings
 maxval = 1;
 length = 4;

for i = 1:20
    for j = 1:20-(length-1)
      currentval = 1;  
      for k = j : j+(length-1)
          readval = data(i,k);
          currentval = currentval*readval;   
      end    
     maxval = max(maxval, currentval);   
        
    end
end

 %%analyze columns
 
for j = 1:20
    for i = 1:20-(length-1)
      currentval = 1;  
      for k = i : i+(length-1)
          readval = data(i,k);
          currentval = currentval*data(i,k);   
      end    
     maxval = max(maxval, currentval);   
        
    end
end

 %%analyze adjacent diagonal elements
 %%left to right
 for i = 1: 20-(length-1)
    for j = 1 : 20-(length-1)
      currentval = 1;  
      for k = 0 : (length-1)
          readval = data(i+k,j+k);
          currentval = currentval*data(i+k,j+k);   
      end    
     maxval = max(maxval, currentval);   
        
    end
 end
 
 %%analyze adjacent diagonal elements
 %%rigth to left
 for i = 1:20 - (length-1)
    for j = 20 :-1:length
      currentval = 1;  
      for k = 0 : (length-1)
          readval = data(i+k,j-k);
          currentval = currentval*data(i+k,j-k);   
      end    
     maxval = max(maxval, currentval);   
        
    end
 end
 
 disp('Max product of four adjacent numbers: ');
disp(maxval);

%% Problem 14
% Finds the number that produces longest chain in the sequence
% n -> n/2 (n is even)
% n -> 3n + 1 (n is odd)

ms = 0;
num = 0;
clear all;
disp('Running Problem 14...');
[ms,num] = Problem14(101);

fprintf(1,'Maximum number of steps: %d \nProduced by number: %d \n', ms,num);

