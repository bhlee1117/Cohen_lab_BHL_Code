app.registration_data = [0 0;1 0;0 1];
figure
imshow(img,[])
hold on
title('Please clicky ROIs')
[ysize, xsize] = size(img);
npts = 1;
nroi = 1;

while(npts > 0)
[xv, yv] = (getline(gca, 'closed'));
if size(xv,1) < 3  % exit loop if only a line is drawn
break
end
%draw the bounding polygons and label them
plot(xv, yv, 'Linewidth', 1);
roi_points{nroi} = [xv, yv];
nroi = nroi + 1;
end
if nroi>1
L = 1024;
W = 768;
pat = false(L,W);
dmd_cam_trans = app.registration_data(1,:);
dmd_rot_dil_mat = app.registration_data(2:3,:)'-[dmd_cam_trans;dmd_cam_trans]';
[x, y] = meshgrid(1:W, 1:L);
dmd_pixel_pos = cellfun(@(xy) round(dmd_rot_dil_mat*(xy'-repmat(dmd_cam_trans',[1 length(xy)]))),roi_points,'uniformoutput',false);

in = cellfun(@(pixel_pos) inpolygon(x,y,pixel_pos(1,:),pixel_pos(2,:)),dmd_pixel_pos,'uniformoutput',false);
for i=1:length(in), pat = pat|in{i}; end
end
figure;imshow(pat)