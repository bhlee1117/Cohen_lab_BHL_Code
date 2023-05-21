function dmd_register_correct(app)

im_live = mat2gray(app.ref_im);
im_test = zeros(size(im_live));
[X,Y] = meshgrid(1:size(im_test,2),1:size(im_test,1));
im_test(inpolygon(X,Y,app.current_rois{:}(:,1),app.current_rois{:}(:,2)))=1;
im_test = mat2gray(imgaussfilt(im_test,1));
[optimizer, metric] = imregconfig('multimodal');
% tform = imregtform(im_test,im_live, 'affine', optimizer, metric);
tform = imregcorr(im_live,im_test);
app.correction_tform = tform.T;

dlmwrite('dmd_registration_correction_data.txt',app.correction_tform,'precision',16);
% im_corr = imwarp(im_test,tform);
% figure;imshowpair(im_corr, im_test,'diff')
% figure;imshowpair(im_test, im_live)
% roi_corr = bwboundaries(im_corr);
