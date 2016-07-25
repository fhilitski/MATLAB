function [ merged_img, vector_img ] = overlay_force_vector_with_traps( fluor_image, point_x, point_y, force_x, force_y, scaling, style, trap1, trap2, radius_pix, overlay_force_labels)
%OVERLAY_FORCE_VECTOR Overlays force vector onto an image + displays
%optical traps
%INPUTS:
%       fluor_image - grayscale uint16 image, i.e. fluorescent image of a filament
%       point_x = x-coordinate (in the image reference frame) of the bead
%       point_y = y-coordinate (in the image reference frame) of the bead
%       force_x = x-component of the force vector
%       force_y = y-component of the force vector
%       scaling = scaling in pixels/pN to make the vectors
%       style = [width, height] 
%               where width, height define line thichness in pixels; use [0,0] for
%               the thinnest line
%       trap1 = [x,y] vector of trap1 coordinates in pixels
%       trap2 = [x,y] vector of trap2 coordinates in pixels
%       radius = radius of the trap in pixels    
%       overlay_force_labels = boolean flag to enable the overlay of text
%                              for the force magnitude
%RETURNS:
%       merged_img
%       vector_img
 
%the image is m-by-n-by-l
[m,n,l] = size(fluor_image);

% %style defines width and height of the overlay line of contour
w = style(1);
h = style(2);

%for proper scaling of images, difeine max and min intensities
max_intensity = max(max(max(fluor_image)));
min_intensity = min(min(min(fluor_image)));
%max allowed intensity in the uint8 based RGB image
max_allowed_intensity = 255; 

%create template for contour overlay image
vector_img = uint8(zeros(m,n,3));
%length of the force in pN
force_length = sqrt(force_x.^2 + force_y.^2); 
%insert text bar with instantaneous force value in pN
if (overlay_force_labels)
    vector_img = insertText(vector_img, [10,10], [num2str(force_length,'%3.1f') ' pN']);
end;

%convert force from pN to pixels;
force_x = round(force_x*scaling);
force_y = round(force_y*scaling);
%length of the force in pixels
force_length = round(sqrt(force_x.^2 + force_y.^2)); 
%force_length = scaling*force_total;
%angle of the line
alpha = asin(force_y/force_length);

x1 = point_x;
y1 = point_y;
x2 = x1 + force_x;
y2 = y1 + force_y;

%length of the force in pixels
force_length = round(sqrt(force_x.^2 + force_y.^2)); 

%increments in x and y;
dx = force_x / force_length;
dy = force_y / force_length;


%create a vector at (point_x, point_y) with components (force_x, force_y)
%create a line of a given length and orientation
for i = 1:force_length;
     x = round(x1 + (i-1)* dx);
     y = round(y1 + (i-1)* dy);    
     for j = -w:1:w
      vector_img(y+j,x,:) = max_allowed_intensity;
     end;
     for j = -h:1:h
      vector_img(y,x+j,:) = max_allowed_intensity;
     end;
end;

%create an arrowhead
ah_length = round(0.2 * force_length);
%endpoints of a line are starting points for the arrowhead
ah_x1 = x2;
ah_y1 = y2; 

%angle of the arrowhead with respect to the line
beta = pi/6;
%draw one side of the arrowhead
% ah_angle = (alpha - pi + beta);
% if ((force_x) < 0) 
%     ah_angle = (alpha - pi + beta)+pi/2+beta; 
% end;
% if (force_y < 0)
%     ah_angle = (alpha - pi + beta)-pi/2+beta; 
% end;

ah_angle1 = ((alpha + beta)-pi);
ah_angle2 = ((alpha - beta)-pi);
ah_angle = [ah_angle1, ah_angle2];
for k = 1:2
    ah_index = ah_angle(k);
    %disp(alpha * 180/pi);
    %disp(ah_index * 180/pi);
    
    if (force_x > 0)
        dx = cos(ah_index);
    else
        dx = -cos(ah_index);
    end;
    
    dy = sin(ah_index);
    for i = 1 : ah_length
        x = round(ah_x1 + (i-1)* dx);
        y = round(ah_y1 + (i-1)* dy);
        vector_img(y,x,:) = max_allowed_intensity;
        for j = -w:1:w
           vector_img(y+j,x,:) = max_allowed_intensity;
        end;
        for j = -h:1:h
           vector_img(y,x+j,:) = max_allowed_intensity;
        end;
    end;
end;

% %draw another side
% beta = -beta;
% ah_angle = alpha - pi + beta;
% if ((force_x) < 0)
%     ah_angle = (alpha - pi + beta)+pi/2+beta;
% end;
% if (force_y < 0)
%     ah_angle = (alpha - pi + beta)-pi/2+beta;
% end;
% 
% dx = cos(ah_angle);
% dy = sin(ah_angle);
% 
% for i = 1 : ah_length
%     x = round(ah_x1 + (i-1)* dx);
%     y = round(ah_y1 + (i-1)* dy);
%     for j = -w:1:w
%         vector_img(y+j,x) = max_intensity;
%     end;
%     for j = -h:1:h
%         vector_img(y,x+j) = max_intensity;
%     end;
% end;

%merged_img = imfuse(fluor_image, vector_img, 'ColorChannels',[1,2,0]);

%DRAW trap1 circle;
%theta = linspace(0,2*pi,100);
radius = radius_pix;
%trap radius in pixels;

for theta = 0:0.1:2*pi
    x_circle = round(point_x + trap1(1) + radius*cos(theta));
    y_circle = round(point_y + trap1(2) + radius*sin(theta));
    vector_img(y_circle, x_circle,:) = max_allowed_intensity;
    
    x_circle = round(point_x + trap2(1) + radius*cos(theta));
    y_circle = round(point_y + trap2(2) + radius*sin(theta));
    %traps can be located outside the image (or the selected ROI)
    %make sure these coordinates are handled properly
    if (x_circle <= 0)
        x_circle = 1;
    end;
    if (x_circle >= n)
        x_circle = n;
    end;
    if (y_circle <= 0)
        y_circle = 1;
    end;
    
    if (y_circle >= m)
        y_circle = m;
    end;
    vector_img(y_circle, x_circle,:) = max_allowed_intensity;
end;

merged_img = imfuse(vector_img,fluor_image, 'ColorChannels',[1,2,0]);


%  h_fluor_fig = figure;
%  imshow(fluor_image, [min_intensity, max_intensity]);
%  colormap(gray);
%  
%  h_vector_fig = figure;
%  imshow(vector_img, [min_intensity, max_intensity]);
 
 
%  h_merged_fig = figure;
%  imshow(merged_img);

end

