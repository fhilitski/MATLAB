function [ value_list, index_list, count_list ] = create_value_list( input_array )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

value_list = [];
index_list = {};
count_list = [];

n = length(input_array);
last_value = false;

value_counter = 1;
value_index = 1;

while (~last_value)
    value = input_array(value_index);
    if ~isnan(value)
        indexes = find(input_array == value);
        value_list = [value_list value];
        index_list{value_counter} = indexes;
        count_list = [count_list length(indexes)];
        input_array(indexes) = NaN;
        value_counter = value_counter + 1;
    end;
    value_index = value_index + 1;
    if value_index > n
        last_value = true;
    end;
end;
end

