function tform = register_DMD(camera_image, ctrl_points)

%Send DMD image through same process as camera image, just in case there is
%any shift at all that could occur
DMD_pts = ctrl_points; 

% Smooth the image with Gaussian filter before making mask
npix = 15;       % smoothing function is an npix by npix square Gaussian weighting with a stdev of sigma
sigma = npix/3;     % standard deviation of gaussian for smoothing
GaussFilt = fspecial('gaussian', npix, sigma);      % generates the gaussian filter
cam_image = imfilter(camera_image,GaussFilt,'replicate','same','conv');    % applies the filter to generate the smoothed background

% Convert to grayscale 
camera_image = mat2gray(cam_image); 

% Automatically threshhold the image
level_cam = graythresh(camera_image); 

BW_cam = im2bw(camera_image, level_cam);

% Get rid of noise
se2 = strel('square', 5); 
BW_cam = imerode(BW_cam, se2);

% Amplify signal even greater
se2 = strel('square', 20); 
BW_cam = imopen(BW_cam, se2); 

% Show image for debugging purposes
figure(1)
imshow(BW_cam, [], 'InitialMagnification', 'fit')

% Find the connected regions 
CC_cam = bwconncomp(BW_cam); 

num_obj = CC_cam.NumObjects;

% Calculate the center of mass of each detected object
cam_points = regionprops(CC_cam, BW_cam, {'WeightedCentroid'});

% Acess to 'Centroid' field of return structure of regionprops is only
% possible with the following syntax
cam_roi = [cam_points.WeightedCentroid];

% Make the format of each row [x y] 
cam_pts = zeros(num_obj, 2);

for i = 1:num_obj
    cam_pts(i, 1:2) = cam_roi(1, (2*i-1):(2*i)); 
end
 
cam_points = uint16(cam_pts); 

% Display the calculated center points by pointing to them

% This code creates an object that looks like a pointer: forgive the
% hard-coding of values
pointer_obj = {1.0*10^3*[1.2901 0.5914; 1.3255 0.6799; 1.3476 0.6533; 1.4140 0.7108; 1.4449 0.6577; 1.3786 0.6312; 1.3786 0.6091; 1.2901 0.5914]}; 
pointer = roi(pointer_obj);

% Plot the image with the pointers showing you where the centers of mass
% were calculated
% camera_image = int16(camera_image); 
figure(2)
for i = 1:num_obj
    pointer.top_leftx = cam_points(i, 1);
    pointer.top_lefty = cam_points(i, 2);
    
    % Prvent collision between pointer and the DMD camera window (copied
    % from fov.m
    object = pointer.data;
    obj = pointer;
    
    if ((obj.top_lefty + obj.height - 1) > size(camera_image, 1))
        yindices = obj.top_lefty:size(camera_image, 1);
        object = object(1:numel(yindices), :);
    elseif (obj.top_lefty < 1)
        yindices = 1:(obj.height - abs(obj.top_lefty) - 1);
        object = object((abs(obj.top_lefty) + 2):end,:);
    else
        yindices = obj.top_lefty:(obj.top_lefty + obj.height - 1);
    end
    
    
    if ((obj.top_leftx + obj.width - 1) > size(camera_image, 2))
        xindices = obj.top_leftx:size(camera_image, 2);
        object = object(:, 1:numel(xindices));
    elseif (obj.top_leftx < 1)
        xindices = 1:(obj.width - abs(obj.top_leftx) - 1);
        object = object(:,(abs(obj.top_leftx + 2):end));
    else
        xindices = obj.top_leftx:(obj.top_leftx + obj.width - 1);
    end
    
    % Add the pointer to the camera image
    camera_image(yindices, xindices) = camera_image(yindices, xindices) + (max(camera_image(:)))*double(object);
end 
imshow(camera_image, [0 1], 'InitialMagnification', 'fit')

%Calculate transformation 
tstruct = maketform('affine',cam_pts,DMD_pts);

% return transformation
tform = tstruct.tdata.T';

