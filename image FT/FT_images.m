clear all;
close all;
clc;

%% read image

image = imread('bundled MTs and traps.jpg');
%image = image';

%% check if image is CMYK and convert to sRGB
dimensions = size(image);

if (dimensions(3) == 4)
    form = makecform('cmyk2srgb');
    srgb_image = applycform(image,form);
    image = rgb2gray(srgb_image);
else
    image = rgb2gray(image);
end;


%% plot grayscale image
figure(1);
imagesc(image);
colormap('gray');
title('Initial Image');
axis off;

%% proceed with FT
ft_no_shift = fft2(double(image));
ft_image = fftshift(ft_no_shift);
amp_ft = abs(ft_image);
phase_ft = angle(ft_image);

%log_ft = log(phase_ft);
log_ft = uint8(amp_ft);
%log_ft = amp_ft;

figure(2);
colormap(gray);
imagesc(log_ft);
title('FT of the initial image');
axis off;

%% make the filter
dim = size(image);
filter_1 = makeFilter(100,dim(1),dim(2),0);

figure(3);
colormap('gray');
imagesc(filter_1);
axis off;
title('Filter for FT');
%% display filtered Fourier images
filtered_ft = filter_1 .* ft_image;

amp_ft = abs(filtered_ft);
phase_ft = angle(filtered_ft);

log_ft = uint8(amp_ft);

figure(4);
colormap(gray);
imagesc(log_ft);
title('Filtered Fourier plane image');
axis off;
%% perform the inverse Fourier transform

result = ifft2(ifftshift(filtered_ft));

figure(5);
colormap(gray);
imagesc(abs(result));
title('Resulting Image with Blocked Spatial Frequencies');
axis off;


