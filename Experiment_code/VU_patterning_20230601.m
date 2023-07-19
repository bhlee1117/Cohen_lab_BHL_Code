addpath('C:\Users\Labmember\Data\ByungHun')
fpath='C:\Users\Labmember\Data\ByungHun\Optopatch\20230603\Snaps';
files = dir(fullfile(fpath, '*.mat'));
files = files(~[files.isdir]);
[~, idx] = max([files.datenum]);
recent_file = files(idx).name;
%recent_file = '184857BHLm026_10x_FOV.mat';
load(fullfile(fpath,recent_file));
%load(fullfile(fpath,'143506BHLm025_optopatch.mat'));
%load(['C:\Users\Labmember\Data\ByungHun\dmd_cal_230601.mat'])
ref_im=snap.img; sz_ref=size(ref_im);
radius=10;
centers=Cell_segment_circle_10x_VU(ref_im,0.85);
centers=cell_detection_manual(ref_im,centers,[]);
figure;
imshow2(ref_im,[]);
hold all
plot(centers(:,1),centers(:,2),'ro')

dmd_now=xx.GetDevice('DMD_Device');
mask_im = zeros(sz_ref);
for i = 1:size(centers, 1)
    [columnsInImage, rowsInImage] = meshgrid(1:sz_ref(2), 1:sz_ref(1));
    circlePixels = (rowsInImage - centers(i,2)).^2 ...
        + (columnsInImage - centers(i,1)).^2 <= radius.^2;
    mask_im(circlePixels) = 1; % change pixel intensity to 255 (white) inside the circle
end

%Orange
transform=dmd_now(1,1).tform;
%transform=snap.tform;
Rfixed = imref2d([dmd_now(1,1).Dimensions(2) dmd_now(1,1).Dimensions(1)]);

newSize = [2304 2304];
rmask = imresize(mask_im, newSize);
% apply transformation from image space to patterning device space
tmask=imwarp(rmask,transform,'OutputView',Rfixed);

dmd_now(1,1).Target=tmask;
dmd_now(1,1).Write_Static();

%Blue
radius_blue=15;
mask_im_blue = zeros(sz_ref);
transform_blue=dmd_now(1,2).tform;
Rfixed_blue = imref2d([dmd_now(1,2).Dimensions(2) dmd_now(1,2).Dimensions(1)]);
imshow2(ref_im,[])
centers_blue=ginput();
for i = 1:size(centers_blue, 1)
    [columnsInImage, rowsInImage] = meshgrid(1:sz_ref(2), 1:sz_ref(1));
    circlePixels = (rowsInImage - centers_blue(i,2)).^2 ...
        + (columnsInImage - centers_blue(i,1)).^2 <= radius_blue.^2;
    mask_im_blue(circlePixels) = 1; % change pixel intensity to 255 (white) inside the circle
end

rmask_blue = imresize(mask_im_blue, newSize);
%rmask_blue = imresize(mask_im, newSize);
% apply transformation from image space to patterning device space
tmask_blue=imwarp(rmask_blue,transform_blue,'OutputView',Rfixed_blue);

dmd_now(1,2).Target=tmask_blue;
dmd_now(1,2).Write_Static();