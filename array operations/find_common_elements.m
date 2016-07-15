function [ found, common ] = find_common_elements( vector1, vector2 )
%FIND_COMMON_ELEMENTS Summary of this function goes here
%   Detailed explanation goes here
 common = [];
 
 n = length(vector1);
 m = length(vector2);
 
 if (m*n) > 0
     for i = 1:n
         for j = 1:m
             if vector1(i) == vector2(j)
                 common(end+1) = vector1(i);
             end;
         end;
     end;   
 end;
 found = ~isempty(common);
end

