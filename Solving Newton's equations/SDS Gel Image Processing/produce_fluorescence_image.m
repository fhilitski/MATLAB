clc;
clear all;
%% Import Typhoon scanner tables here manually

%% Concatenated all the imported sheets

gel_raw = cat(2,a647_1,a647_2);
imagesc(gel_raw);
colormap(gray);
rgb = uint16(round(gel_raw));
imagesc(rgb);
imwrite(rgb,'gel_a647.tiff');