function [ x_binned, y_binned, error_binned ] = binned_xy( x, y )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
  
[vl,il,cl] = create_value_list(x);

n = length(vl);
x_binned = vl;
y_binned = NaN(1,n);
error_binned = NaN(1,n);

for i = 1:length(il)
    y_group = y(il{i}); 
    y_binned(i) = mean(y_group);
    error_binned(i) = std(y_group);
end;

end

