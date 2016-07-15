function [ found, index ] = find_inside_cell( target, expression)
%FIND_INSIDE_CELL Summary of this function goes here
%   Detailed explanation goes here
 found = false;
 index = [];
 
 n = length(target);
 if n > 0
     for i = 1:n
      if (strcmp(target{i}, expression))
          index = [index i];
      end;
     end; 
 end;
 
 found = ~isempty(index);
end

