

fn1= '202734spot_z_3.7536.tif';
fn2= '203103spot_z_3.8075.tif';
clim = [1e+3 20e+3];

im = imread(fn1);
figure(19);clf
subplot(1,2,1)
imshow2(double(im),[]);title(fn1(end-9:end-3))
max(im(:))
colormap(jet);colorbar;set(gca,'clim',clim)
im = imread(fn2);
max(im(:))
subplot(1,2,2)
imshow2(double(im),[]);title(fn2(end-9:end-3))
colormap(jet);colorbar;set(gca,'clim',clim)

%%
fn1= '160105im0.tif';
fn2= '160116im1.tif';

im1 = imread(fn1);
im2 = imread(fn2);

% figure;imshowpair(ffzpad(im1,size(im2)),im2)
figure;imshow2(ffzpad(im1,size(im2))-im2,[])