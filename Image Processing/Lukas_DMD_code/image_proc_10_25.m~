% Create square and move it across field of view of some size

% Create square

width = 4; 
height = 4; 
for i = 1:width
    for j = 1:height
        square(i,j) = 1; 
    end
end

%Create field of view 

fovwidth = 10; 
fovheight = 10; 
for i = fovwidth
    for j = fovheight
        fov(i,j) = 0; 
    end
end

% Add the square to the field of view with a center at some point

tleft = [2; 3]; 

fov(tleft(1,:):(tleft(1,:) + height - 1), tleft(2,:):(tleft(2,:) + width -1)) = square; 

% Use clicky to specify an arbitrary shape of field of view
refimg = ones(20, 20); 

figure(1)
imshow(refimg, [], 'InitialMagnification', 'fit') %scales the reference image to fit the clicky image window
hold on;

[ysize, xsize] = size(refimg(:,:,1)); %because it will return the number of rows and then columns
npts = 1;
colorindex = 0;
order = get(gca,'ColorOrder');
nroi = 1;
intens = [];
  [x, y] = meshgrid(1:xsize, 1:ysize);
 while(npts > 0) %where does npts go to <= 0?

  subplot(1,2,1)
  [xv, yv] = (getline(gca, 'closed')); %getline(fig) lets you select a polyline in the current axes of figure fig using the mouse. 
  %Coordinates of the polyline are returned in X and Y. Use normal button clicks to add points to the polyline. A shift-, right-, or double-click adds a final point and ends the polyline selection. 
  %Pressing Return or Enter ends the polyline selection without adding a final point. Pressing Backspace or Delete removes the previously selected point from the polyline.
  if size(xv,1) < 3  % exit loop if only a line is drawn
  break
  end
  inpoly(:,:,nroi) = inpolygon(x,y,xv,yv);
  
 %draw the bounding polygons and label them
  currcolor = order(1+mod(colorindex,size(order,1)),:); %wrap the colors so that if there are more plots than "order" holds, 
  %Then the current color will wrap around to a color already used
  plot(xv, yv, 'Linewidth', 1,'Color',currcolor);
  text(mean(xv),mean(yv),num2str(colorindex+1),'Color',currcolor,'FontSize',12);
  
  
  colorindex = colorindex+1;
  
  roi_points{nroi} = [xv, yv];
  nroi = nroi + 1;
 end

 
% find a bounding rectangle that surrounds the polygon - a bounding box
sum_rows = sum(inpoly, 1); 
sum_column = sum(inpoly, 2); 
coords_columns = find(sum_rows > 0); 
coords_rows = find(sum_column > 0); 
box_locale = [coords_rows(1); coords_columns(1)]; 
box_width = size(coords_columns, 2);  
box_height = size(coords_rows, 1);  

% Create this box
box = inpoly(box_locale(1,:):(box_locale(1,:) + box_height - 1), box_locale(2,:):(box_locale(2,:) + box_width -1));

% Insert the box into an image, fov, at an abitrary location
fov2 = zeros(40); 
tleft = [2; 3]; 

fov2(tleft(1,:):(tleft(1,:) + box_height - 1), tleft(2,:):(tleft(2,:) + box_width -1)) = box; 

% Now, think about what needs to happen if there is a linear transformation
% of the box before it is re-inserted into the field of view...actually kit
% already has the code on this so it's not as important

%Use the classes roi and fov to add polygons to an image

roi1 = roi(inpoly);
fov3 = fov(40, 40);
fov3.add_obj(roi1, tleft); 

tleft = [20; 3]; 
fov3.add_obj(roi1);
fov3.set_loc(roi1, tleft); 

image1 = zeros(1024, 1024); 
image1 = fov5.data;
[roi_points inpoly] = clicky_stripped(image1); 

for i = 1:size(roi_points, 2)
    roi_array(i) = roi(roi_points(i), inpoly(i));
end

fov5 = fov(1024, 1024); 
for i = 1:size(roi_array, 2)
    fov5.add_obj(roi_array(i));
end

% Maybe I can specify the time data of the objects in long arrays? 
% Problem 1: alternate two objects back and forth
clf
clear fov1
fov1 = fov(1024, 1024, 8); 

t1 = [0 1 0 1 0 1 0 1]; 
t2 = [1 0 1 0 1 0 1 0]; 

image1 = zeros(1024, 1024);
[roi_points inpoly] = clicky_stripped(image1); 

for i = 1:size(roi_points, 2)
    roi_array(i) = roi(roi_points(i), inpoly(i));
end

% Find the frames that a given object will be active
obj1 = find(t1 > 0); 
obj2 = find(t2 > 0); 

% Add the object to those frames
for i = 1:size(obj1, 2)
    fov1.add_obj(roi_array(1), obj1(i)); 
end

for i = 1:size(obj2, 2)
    fov1.add_obj(roi_array(2), obj2(i)); 
end

% Visualize if the object has actually been added
for i = 1:fov1.frames
    figure(i)
    imshow(fov1.data(:,:,i))
end


% Success! 

% Now try to use the function makewave_DMD: 

% Try using a transform matrix
tformMat = [cos(pi/8) sin(pi/8) 0; -sin(pi/8) cos(pi/8) 0; 0 0 1]; 

tformMat = [0 1 0; 1 0 0; 0 0 1]; 

clf
clear fov1
fov1 = fov(1024, 1024, 8); 

t1 = [0 1 0 1 0 1 0 1]; 
t2 = [1 0 1 0 1 0 1 0]; 
t3 = [0 1 0 1 0 1 0 1]; 

image1 = zeros(1024, 1024);
roi_points = clicky_stripped(image1); 

for i = 1:size(roi_points, 2)
    roi_array(i) = roi(roi_points(i), tformMat);
end

makewave_DMD(fov1, roi_array(1), t1, roi_array(2), t2); 

for i = 1:fov1.frames
    figure(i+1)
    imshow(fov1.data(:,:,i))
end

% Play around with transformation matrices
x = [1 1; 2 4; 3 2]; 
y = [1 1; 3 -1; 2 -3]; 

tstrcut = maketform('affine', x, y);

x = [x' ; ones(1,3)]; 
y = [y' ; ones(1,3)]; 

trans = tstruct.tdata.T'; 

% Now it's time to make some time vectors! 

peak1 = peak_obj(1, 10000, 10, 100); 

figure(1)
plot(peak1.data)

time_vect1 = time_vect(50000, peak1); 

figure(2) 
plot(time_vect1.data)

peak2 = peak_obj(40000, 7500, 1, 7500);
figure(3)
plot(peak2.data)

time_vect1.add_peak(peak2)

figure(4)
plot(time_vect1.data)



peak1 = peaks(1, 10000, 10, 100); 
time_vect1 = time_vect(50000, peak1); 
peak2 = peaks(40000, 7500, 1, 7500);
time_vect1.add_peak(peak2)
t1 = time_vect1.data; 

peak1 = peaks(1, 10000, 10, 100); 
time_vect2 = time_vect(50000, peak1); 
peak2 = peaks(40000, 7500, 1, 7500);
time_vect2.add_peak(peak2)
t2 = time_vect2.data; 

peak1 = peaks(1, 10000, 10, 100); 
time_vect3 = time_vect(50000, peak1); 
peak2 = peaks(40000, 7500, 1, 7500);
time_vect3.add_peak(peak2)
t2 = time_vect3.data; 

%Make a waveform to go to the DMD
fov1 = fov(1080, 1920, 1); 

red_region = [512 844; 512 1204; 1536 1204; 1536 844];
all_on = true(1080, 1920); 
movie_to_bin(all_on, 'C:\data\Lukas_data\all_on.bin')

fid = fopen('C:\data\Lukas_data\all_on.bin', 'w');
fwrite(fid, permute(all_on,[2 1 3]), 'ubit1','b');     % DMD reads by rows, Matlab saves by columns, big endian ordering
fclose(fid);

dlmwrite('C:\data\Lukas_data\red_roi.txt', red_region, 'delimiter', '\t'); 

square = {[100 200; 200 200; 200 100; 100 100]}; 
roi1 = roi(square); 
roi1.top_leftx = 600; 
roi1.top_lefty = 300; 

roi4 = roi(square); 
roi4.top_leftx = 500; 
roi4.top_lefty = 200; 

roi2 = roi(square); 
roi2.top_leftx = 900; 
roi2.top_lefty = 700; 

roi5 = roi(square); 
roi5.top_leftx = 800; 
roi5.top_lefty = 600; 

roi3 = roi(square); 
roi3.top_leftx = 1200; 
roi3.top_lefty = 300; 

roi6 = roi(square); 
roi6.top_leftx = 1100; 
roi6.top_lefty = 200; 

fov1.add_obj(roi1, 1); 
fov1.add_obj(roi2, 1); 
fov1.add_obj(roi3, 1);
fov1.add_obj(roi4, 1); 
fov1.add_obj(roi5, 1); 
fov1.add_obj(roi6, 1); 

figure(1)
imshow(fov1.data(:,:,1), [], 'InitialMagnification', 'fit')

ctl_points = [600 300; 900 700; 1200 300]; 

% Write that waveform to binary
movie_to_bin(fov1.data, 'C:\data\Lukas_data\image_reg.bin')
im = imread('C:\data\Lukas_data\image_reg.tif'); 
tform = register_DMD(im, ctl_points); 

mpath = 'C:\data\Lukas_data\movie.bin';
% Write that movie to binary
movie_to_bin(fov1.data, mpath); 

fileID = fopen('C:\data\Lukas_data\image_reg.txt','w');
fprintf(fileID,'%.2f %.2f\n',ctl_points);
fclose(fileID);

dlmwrite('C:\data\Lukas_data\image_reg.txt', ctl_points, 'delimiter', '\t'); 

x = dlmread('C:\data\Lukas_data\image_reg.txt'); 

% Work through some image registration codee
im = imread('C:\data\Lukas_data\im_test1_at101545.tiff'); 
ctl_points = dlmread('C:\data\Lukas_data\image_reg.txt'); 

register_DMD(im, ctl_points);

% Debug labview code by running it here
% clear previous roi_points cell array and roi objects

data = imread('C:\data\Lukas_data\im_test4_at092119.tiff'); 
data2 = imread('C:\data\Lukas_data\im_test3_at091749.tiff'); 
tpath = 'C:\data\Lukas_data\transform.txt';
xoffset = 0;
yoffset = 0; 

figure(1)
subplot(1,2,1)
imshow(data, [], 'InitialMagnification', 'fit');
subplot(1,2,2)
imshow(data2, [], 'InitialMagnification', 'fit');
clear roi_points; 
clear roi_objects; 

% define roi points
roi_points = clicky_stripped(data); 

% read in transformation matrix
tform = dlmread(tpath); 

% Add xoffset and yoffset and create roi objects
for i = 1:numel(roi_points)
    roi_points{i}(:,1) = roi_points{i}(:,1) + xoffset; 
    roi_points{i}(:,2) = roi_points{i}(:,2) + yoffset; 
    roi_object(i) = roi(roi_points(i), tform); 
end

% create fov object
fov1 = fov(1080, 1920, 1); 

% Add the roi objects to the field of view
for i = 1:numel(roi_object)
fov1.add_obj(roi_object(i), 1); 
end

% Write that movie to binary
movie_to_bin(fov1.data, mpath); 

% test the manual image registration
im = imread('C:\data\Lukas_data\im_test1_at101545.tiff'); 
ctl_points = dlmread('C:\data\Lukas_data\image_reg.txt'); 
tform = register_DMD_man(im, ctl_points); 

tpath = 'C:\data\Lukas_data\transform.txt';
tform_check = dlmread(tpath);

data = imread('C:\data\Lukas_data\_at133612.tiff');

% define roi points
roi_points = clicky_stripped(data); 

% read in transformation matrix

% Add xoffset and yoffset and create roi objects
for i = 1:numel(roi_points)
    roi_points{i}(:,1) = roi_points{i}(:,1); 
    roi_points{i}(:,2) = roi_points{i}(:,2); 
    roi_objects(i) = roi(roi_points(i)); 
end

% Make the roi objects from the transformed points
for i = 1:numel(roi_points)
    roi_points{i}(:,1) = roi_points{i}(:,1); 
    roi_points{i}(:,2) = roi_points{i}(:,2); 
    roi_objects(i) = roi(roi_points(i)); 
end

% create fov object
fov1 = fov(1080, 1920, 1); 

% Add the roi objects to the field of view
for i = 1:numel(roi_objects)
fov1.add_obj(roi_objects(i), 1); 
end

% create fov object
fov2 = fov(2048, 2048, 1); 

% Add the roi objects to the field of view
for i = 1:numel(roi_objects)
fov2.add_obj(roi_objects(i), 1); 
end

figure(3)
imshow(fov1.data, [], 'InitialMagnification', 'fit'); 

figure(4)
imshow(fov2.data, [], 'InitialMagnification', 'fit'); 

% Write that movie to binary
movie_to_bin(fov1.data, mpath); 

% Saving to a text file
exposure = 100; 
frames = 10; 
bin = 4; 
xoffset = 200; 
yoffset = 300; 
xwidth = 500; 
ywidth = 700; 

c = cell([7 2]); 
c{1,1} = 'exposure'; 
c{1,2} = exposure;  
c{2,1} = 'frames'; 
c{2,2} = frames;
c{3,1} = 'bin'; 
c{3,2} = bin;

fid = fopen('C:\data\Lukas_data\test_file_io.txt', 'w');
for i = 1:7
    fprintf(fid, '%s\t%g\n', c{i,:}); 
end
fclose(fid);

x = imread('C:\data\Lukas_data\all_on_at215224.tiff'); 
x = double(x);
roi_points = clicky_stripped(x); 

% Debug Labview code 

%load in image registration movie and control points
ctrl_points = [600 300; 900 700; 1200 300]; 
im = imread('/Users/lukasgemar/Documents/MATLAB/Lukas_DMD_code/_at223259.tiff'); 

%Calculate the transform; tform in terms of 1x1 binning
%for i = 1:10
tform(:,:,i) = register_DMD_man(im, ctrl_points);
%end
bin = 4;
tform(1, 1, :) = (1/bin)*tform(1,1,:);
tform(2,2, :) = (1/bin)*tform(2,2,:);

tform2 = adjust_tform(

roi_points = clicky_stripped(im); 
num_obj = numel(roi_points); 
xoffset = 0; 
yoffset = 0; 

clear roi_objects; 
% Add xoffset and yoffset and create roi objects (in terms of 1x1 bin)
for j = 1:size(tform, 3)
    for i = 1:num_obj
        if (j == 1)
            roi_points{i}(:,1) = (roi_points{i}(:,1) + xoffset) * bin;
            roi_points{i}(:,2) = (roi_points{i}(:,2) + yoffset) * bin;
        end
        roi_objects(i) = roi(roi_points(i), tform(:,:,j));
        corner(j,:, i) = [roi_objects(i).top_leftx, roi_objects(i).top_lefty];
    end
end

% for each object calculate the standard deviation of the mapping onto DMD
for i = 1:num_obj
    std_devs_dmd(i,:) = std(corner(:,:,i)); 
end

%Calculate the mean value of the control points
for i = 1:num_obj
    means(i,:) = mean(corner(:,:,i)); 
end

% Calculate how many pixels off you would have to be to produce that type
% of variation in the DMD coords
tform_inv = inv(tform); 
test1 = tform_inv(:,:,1)*[600; 300; 1];
test2 = tform_inv(:,:,1)*[630; 330; 1];
test3 = tform_inv(:,:,1)*[570; 270; 1];

% Check to make sure that the locations are actually where the objects were
% placed: [600 300; 900 700; 1200 300]

% create fov object
fov1 = fov(1080, 1920, 1); 

% Add the roi objects to the field of view
for i = 1:numel(roi_objects)
fov1.add_obj(roi_objects(i), 1); 
end

figure(3)
imshow(fov1.data)

% Map the DMD boxes back onto the camera using the inverse transform
% Create a fov object representing the camera 
fov2 = fov(2048, 2048, 1); 

% Create roi objects for each of the boxes on the DMD
square = {[100 200; 200 200; 200 100; 100 100]}; 
roi1 = roi(square); 
roi1.top_leftx = 600; 
roi1.top_lefty = 300; 
roi1.apply_tform(inv(tform)); 
fov2.add_obj(roi1, 1); 

roi4 = roi(square); 
roi4.top_leftx = 500; 
roi4.top_lefty = 200; 
roi4.apply_tform(inv(tform)); 
fov2.add_obj(roi4, 1); 

roi2 = roi(square); 
roi2.top_leftx = 900; 
roi2.top_lefty = 700; 
roi2.apply_tform(inv(tform)); 
fov2.add_obj(roi2, 1); 

roi5 = roi(square); 
roi5.top_leftx = 800; 
roi5.top_lefty = 600; 
roi5.apply_tform(inv(tform)); 
fov2.add_obj(roi5, 1); 

roi3 = roi(square); 
roi3.top_leftx = 1200; 
roi3.top_lefty = 300; 
roi3.apply_tform(inv(tform)); 
fov2.add_obj(roi3, 1); 

roi6 = roi(square); 
roi6.top_leftx = 1100; 
roi6.top_lefty = 200; 
roi6.apply_tform(inv(tform)); 
fov2.add_obj(roi6, 1); 

figure(1)
imshow(fov2.data(:,:,1), [], 'InitialMagnification', 'fit')

% Apply the convolution of a selected part of the inverse matrix 
bin = 4; 
im = imresize(im, bin); 

% call the image matrix A and the convolution matrix B
top_left = [roi4.top_leftx, roi4.top_lefty]; 
btm_rgt = [(roi1.top_leftx + roi1.width), (roi1.top_lefty + roi1.height)]; 
A = fov2.data(top_left(2):btm_rgt(2), top_left(1):btm_rgt(1)); 
ycoords_vect = uint16((top_left(2) - .5*(btm_rgt(2) - top_left(2))):(btm_rgt(2) + .5*(btm_rgt(2) - top_left(2)))); 
xcoords_vect = uint16((top_left(1) - .5*(btm_rgt(1) - top_left(1))):(btm_rgt(1) + .5*(btm_rgt(1) - top_left(1)))); 
B = im(ycoords_vect, xcoords_vect); 
A = double(A); 
B = double(B); 
C = conv2(B, A, 'same'); 

% find the indices of the maximum value of the convolution
[xmax, ymax] = find(C == max(C(:)));

% Try a more general method of finding the mapping of DMD boxes onto camera
DMD = fov1.data; 
[y, x] = find(DMD); 

% put those points into a format on which the transform can act
temp = [x'; y'; ones(1, numel(x))]; 
dmd_tformpts = tform\temp; 

% map those points onto a camera field of view 
cam = zeros(2048); 

temp = uint16(dmd_tformpts); 

for i = 1:numel(dmd_tformpts(1,:))
    cam(temp(2,i), temp(1,i)) = 1; 
end

% close the boxes
se = strel('square', 5); 
BW_cam = imclose(cam, se);

% Find the connected regions 
CC_cam = bwconncomp(BW_cam); 
L = labelmatrix(CC_cam);

% make images out of the connected regions
[y1, x1] = find(L == 1); 
[y2, x2] = find(L == 2); 
[y3, x3] = find(L == 3); 

obj1 = zeros(2048); 
for i = 1:numel(y1)
    obj1(y1(i), x1(i)) = 1; 
end

obj1 = obj1(min(y1):max(y1), min(x1):max(x1)); 

obj2 = zeros(2048); 
for i = 1:numel(y2)
    obj2(y2(i), x2(i)) = 1; 
end

obj2 = obj2(min(y2):max(y2), min(x2):max(x2));

obj3 = zeros(2048); 
for i = 1:numel(y3)
    obj3(y3(i), x3(i)) = 1; 
end

obj3 = obj3(min(y3):max(y3), min(x3):max(x3));

% reshape the image registration image
im = imresize(im, 4); 

% do the same thing with the function ctrl_points
pts1 = find_center(im, L, 1); 
pts2 = find_center(im, L, 2); 
pts3 = find_center(im, L, 3); 

% perform the convolutions to find the new image registration control
% points

% call the image matrix A and the convolution matrix B
top_left = [min(x1), min(y1)]; 
btm_rgt = [max(x1), max(y1)]; 
A = obj1; 
ycoords_vect = uint16((top_left(2) - .5*(btm_rgt(2) - top_left(2))):(btm_rgt(2) + .5*(btm_rgt(2) - top_left(2)))); 
xcoords_vect = uint16((top_left(1) - .5*(btm_rgt(1) - top_left(1))):(btm_rgt(1) + .5*(btm_rgt(1) - top_left(1)))); 
B = im(ycoords_vect, xcoords_vect); 
A = double(A); 
B = double(B); 
C = conv2(B, A, 'same'); 

% find the indices of the maximum value of the convolution
[xmax, ymax] = find(C == max(C(:)));

% convert those indices to real camera coordinates
xcam = xmax + min(xcoords_vect); 
ycam = ymax + min(ycoords_vect); 


%see if it's working! 
ctrl_points = [600 300; 900 700; 1200 300]; 
im = imread('/Users/lukasgemar/Documents/MATLAB/Lukas_DMD_code/_at223259.tiff'); 

%Calculate the transform; tform in terms of 1x1 binning
tform(:,:,i) = register_DMD_man(im, ctrl_points);
bin = 4;
tform(1, 1, :) = (1/bin)*tform(1,1,:);
tform(2,2, :) = (1/bin)*tform(2,2,:);

DMD = fov1.data;
im = imresize(im, 4); 
tform2 = adjust_tform(im, DMD, ctrl_points, tform); 






















