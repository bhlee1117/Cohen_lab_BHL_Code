function tform = adjust_tform2(im, DMD, ctrl_pts, tform)

% map the DMD coordinates onto camera coordinates
[y, x] = find(DMD);  

% put those points into a format on which the transform can act
temp = [x'; y'; ones(1, numel(x))]; 

% map those points onto a camera field of view 
cam = zeros(2048); 
dmd_tformpts = tform\temp; 

temp = uint16(dmd_tformpts); 

for i = 1:numel(dmd_tformpts(1,:))
    cam(temp(2,i), temp(1,i)) = 1; 
end

% close the boxes
se = strel('disk', 25); 
BW_cam = imclose(cam, se);

% Find the connected regions 
CC_cam = bwconncomp(BW_cam); 
L = labelmatrix(CC_cam);

% find the centers of 3 connected regions
cam_pts(1,:) = find_center(im, L, 1); 
cam_pts(2,:) = find_center(im, L, 2); 
cam_pts(3,:) = find_center(im, L, 3); 


% sort the element along x in ascending order
cam_pts = sortrows(cam_pts, 1); 
cam_pts = double(cam_pts);
ctrl_pts = sortrows(ctrl_pts, 1); 

%Calculate transformation 
tstruct = maketform('affine',cam_pts,ctrl_pts);

% return transformation
tform = tstruct.tdata.T';
end