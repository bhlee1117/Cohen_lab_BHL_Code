function tform = register_DMD_man(camera_image, ctrl_points)

%Send DMD image through same process as camera image, just in case there is
%any shift at all that could occur
DMD_pts = ctrl_points; 
camera_image = mat2gray(camera_image); 

% Use clicky to get the camera points
temp = clicky_stripped(camera_image); 

% Access the camera points in an array
cam_pts = temp{1};

% Acess only the first three rows
cam_pts = cam_pts(1:3,:); 

% sort the element along x in ascending order
cam_pts = sortrows(cam_pts, 1); 
DMD_pts = sortrows(DMD_pts, 1); 

%Calculate transformation 
tstruct = maketform('affine',cam_pts,DMD_pts);

% return transformation
tform = tstruct.tdata.T';

