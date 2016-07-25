%% Problem 14
% Finds the number that produces longest chain in the sequence
% and the number of steps in the sequence
% n -> n/2 (n is even)
% n -> 3n + 1 (n is odd)
% Returns:  [maxSteps,number]


function [maxSteps,number] = Problem14(limit)
maxCount = 1;
numberWithMC = 1;

for i = limit:-1:1

j = i;
count = 1;
while (j >1)
 if (rem(j,2) == 0) 
   j = j/2;    
 else j = 3*j + 1;
 end;
 count = count + 1;
end
maxCount = max(maxCount,count);
if (count >= maxCount) numberWithMC = i; end;

end;

maxSteps = maxCount;
number = numberWithMC;
end



