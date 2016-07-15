function [ y ] = linear_map(x,x_min,x_max,y_min,y_max)
%LINEAR_MAP Summary of this function goes here
%   linearly maps orgument x from the range [x_min; x_max] into y on the interval
%   [y_min; y_max]
 
k = (y_min - y_max)/(x_min-x_max);
b = 1/2 * ((y_min+y_max) - k*(x_min+x_max));
y = (k*x+b);

end

