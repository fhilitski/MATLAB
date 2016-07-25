function [ merged_img ] = addContour( img, contour, style, debug_mode)
%ADDCONTOUR adds a contour to the image
%   [ merged_img ] = ADDCONTOUR( img, contour, style, debug_mode)
%   img - grayscale uint16 image, i.e. fluorescent image of a filament
%   contour - n-by-2 matrix that has x-y coordinates of the tracked
%              filament.
%   style - 1-by-2 matrix that contains width of the contour overlay.
%           Example: [1,1]; [0,1];
%   debug_mode - boolean expression to enable temporary output;
%
%   merged_image - m-by-n-by-3 uint8 RGB image of the overlay

%the image os m-by-n-by-l
[m,n,l] = size(img);

%style defines width and height of the overlay line of contour
w = style(1);
h = style(2);

%for proper scaling of images, difeine max and min intensities
max_intensity = max(max(max(img)));
min_intensity = min(min(min(img)));

%create template for contour overlay image
contour_img = uint16(zeros(m,n,l));

%coordinates of the filament
x_coords = round(contour(:,1));
y_coords = round(contour(:,2));
bad_x = find(x_coords <= 1);
bad_y = find(y_coords <= 1);
x_coords(bad_x) = 2;
y_coords(bad_y) = 2;

%create an overlay image
 for i = 1:length(x_coords)
     contour_img(y_coords(i),x_coords(i)) = max_intensity;
     contour_img(y_coords(i)+ w,x_coords(i)) = max_intensity;
     contour_img(y_coords(i)- w,x_coords(i)) = max_intensity;
     contour_img(y_coords(i),x_coords(i)+ h) = max_intensity;
     contour_img(y_coords(i),x_coords(i)- h) = max_intensity;
 end;

%merge image with the tracked contour
merged_img = imfuse(contour_img, img, 'ColorChannels',[2 1 2]);
 
%plot out images if in debug mode
if (debug_mode)
    h_fig_img = figure;
    imshow(img,[min_intensity max_intensity]);
    hold on;
    %plot(x_coords,y_coords, style);
    %frame_image = frame2im(getframe(h_fig));
    
    h_fig_contour = figure;
    imshow(contour_img,[min_intensity max_intensity]);
    
    h_fig_merge = figure;
    imshow(merged_img);
    
    close(h_fig_img);
    close(h_fig_contour);
    close(h_fig_merge);
end;

end

