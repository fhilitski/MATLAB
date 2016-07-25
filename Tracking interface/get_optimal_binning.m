function [number_of_bins, bin_size] = get_optimal_binning (array)
%[number_of_bins, bin_size] = get_optimal_binning (array)
%array - vector of data;
%number_of_bins - number of bins in normally distributed
%                 array using  square root  rule 
%bin_size - returns optimal bin size using Scott's rule 

%get std of the array
%Scott's rule

st_dev = std(array);
n = length(array);

bin_size =  3.5*st_dev/(n^(1/3));

min_value = min(array);
max_value = max(array);


%number_of_bins = round((max_value - min_value)/bin_size);

%squareroot rule
number_of_bins = ceil(sqrt(n));





