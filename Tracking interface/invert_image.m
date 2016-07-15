function [ inverted ] = invert_image (img_array, bit_depth)
%invert_image  inverts grayscale image
%img_array  double-precision image to be inverted
%bit_depth  bit depth of the image

inverted = 2^bit_depth - img_array;
    

end

