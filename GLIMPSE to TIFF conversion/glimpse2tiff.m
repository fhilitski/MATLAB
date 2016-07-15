clc;
clear all;
close all;


folder_path = 'D:\Data - MT Sliding and Friction\2016\03-23-2016 TIRF active bundles\jchung00012\';
tiff_folder = [folder_path 'TIFF images'];
tiff_folder_top = [folder_path 'TIFF images top'];
tiff_folder_bottom = [folder_path 'TIFF images bottom'];
split_channels = true; %split top and bottom views;
fix_channel_offset = true;%two channels appear to be offset spatially 

%% obtain acquisition information
header = open([folder_path 'header.mat']);
gheader = header.vid;
number_of_images = gheader.nframes; 
img_width = gheader.width;
img_height = gheader.width;
time = gheader.ttb; %absolute times
dt = time - time(1); %relative time from the first frame



%% extract first image
i = 1;
extracted_img = glimpse_image(folder_path,gheader,i);
max_pixel = max(max(extracted_img));
min_pixel = min(min(extracted_img));
h_fig = figure;
h_img = imshow(extracted_img,[min_pixel max_pixel]);

%split fields of view if necessary
if (split_channels)
    img_fov_height = img_height / 2;
    %the top half of the image
    top_img_fov = extracted_img(1:img_fov_height,:);
    %bottom half of the image
    bottom_img_fov = extracted_img(img_fov_height+1:end,:);
   
    %display the top and bottom separately
    h_fig_top = figure;
    max_pixel = max(max(top_img_fov));
    min_pixel = min(min(top_img_fov));
    h_img_top = imshow(top_img_fov,[min_pixel max_pixel]);
   
    h_fig_bottom = figure;
    max_pixel = max(max(bottom_img_fov));
    min_pixel = min(min(bottom_img_fov));
    h_img_bottom = imshow(bottom_img_fov,[min_pixel max_pixel]);
    
    merged_img = imfuse(bottom_img_fov,top_img_fov, 'ColorChannels',[1,2,0]);
    h_fig_merged = figure;
    h_merged_img = imshow(merged_img);
end;

%% prepare to save the images
mkdir_status = mkdir(tiff_folder);
mkdir_status = mkdir(tiff_folder_top);
mkdir_status = mkdir(tiff_folder_bottom);

%% extract and save the remaining images without showing
for i = 1 : number_of_images
    extracted_img = glimpse_image(folder_path,gheader,i);
    tiff_filename = [tiff_folder '\' num2str(i) '.tiff'];
    imwrite(extracted_img, tiff_filename,  'TIFF');
    if (split_channels)
        img_fov_height = img_height / 2;
        %the top half of the image
        top_img_fov = extracted_img(1:img_fov_height,:);
        %bottom half of the image
        bottom_img_fov = extracted_img(img_fov_height+1:end,:);
        tiff_filename_top = [tiff_folder_top '\' num2str(i) '.tiff'];
        tiff_filename_bottom = [tiff_folder_bottom '\' num2str(i) '.tiff'];
        imwrite(top_img_fov, tiff_filename_top,  'TIFF');
        imwrite(bottom_img_fov, tiff_filename_bottom,  'TIFF');
    end;
    disp(['Working on image ' num2str(i) ' out of ' num2str(number_of_images)]);
end;
disp('Done!');

%% save time as a separate variable
save([folder_path 'time.mat'], 'time','dt');

