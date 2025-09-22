function [aligned_avgImg tform]=imReg_faster(avgImg,avgImg2)

% Normalize input images
fixed  = mat2gray(avgImg2);  % target image (reference)
moving = mat2gray(avgImg);   % image to be aligned

% Automatically estimate rotation + translation
tform = imregcorr(moving, fixed, 'rigid');

% Apply transformation to avgImg
Rfixed = imref2d(size(fixed));
aligned_avgImg = imwarp(moving, tform, 'OutputView', Rfixed);

% Display results
nexttile([1 1]);
imshowpair(aligned_avgImg, fixed, 'falsecolor');
title('Aligned avgImg with avgImg2');