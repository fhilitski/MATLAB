%% Create Tiff object

filename = 'run10.tif';
mode = 'r+' %open file for reading and writing

images = Tiff(filename, mode);

%% Read image information
filename = 'run11.tif';
info = imfinfo(filename);

number_of_images = size(info);

for i=1:number_of_images
    image = imread(filename,'Index',i);
    
    figure(1);
    imagesc(image);
    colormap(gray);
    
    noise_height = 29 + (380-350);
    noise_region = image(150:150+noise_height,:);
    
    figure(2);
    imagesc(noise_region);
    colormap(gray);
    
    new_image = image;
    
    new_image(1:29,:) = noise_region(1:29,:);
   
    new_image(350:380,:) = noise_region(30:end,:);
    
    figure(3);
    imagesc(new_image);
    colormap(gray);
    
    imwrite(new_image,'e:\run11_new.tif','WriteMode','append');
    
    %if ( mod(i,100) == 0)
     %   disp(i);
    %end;
end;



