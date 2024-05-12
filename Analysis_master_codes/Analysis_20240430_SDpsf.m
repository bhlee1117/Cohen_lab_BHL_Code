clear;
fpath='/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/20240429/F8887_20nm_10E4dilution_500msexp_0p5umpersec_30um_2106.tif';
beadImg=double(readtiff(fpath));
pixelSz=6.5/25; %um
zstep=0.25; %um

%%
[roi, intens]=clicky(beadImg);
%%
[~, focus_z]=max(intens);
cent=round(mean(roi{1},1));
crop_im = beadImg(cent(2)+[-20:20],cent(1)+[-20:20],focus_z);

%%

gaussianFunc = @(p, data) ...
    p(1) * exp(-((cos(p(6))^2/(2*p(4)^2) + sin(p(6))^2/(2*p(5)^2)) * (data{1}-p(2)).^2 + ...
    (sin(2*p(6))/(4*p(5)^2) - sin(2*p(6))/(4*p(4)^2)) * (data{1}-p(2)).*(data{2}-p(3)) + ...
    (sin(p(6))^2/(2*p(4)^2) + cos(p(6))^2/(2*p(5)^2)) * (data{2}-p(3)).^2)) + p(7);

gaussianFunc1d = @(p,data) exp(-(data - p(1)).^2 / (2 * p(2)^2));

[xData, yData] = meshgrid(1:size(crop_im, 2), 1:size(crop_im, 1));
zData = double(crop_im);

% Initial guesses for A, x0, y0, sigmaX, sigmaY, theta, and offset
initialGuess = [max(zData(:)), size(image, 2)/2, size(image, 1)/2, 20, 20, 0, min(zData(:))];

% Lower and upper bounds for the parameters
lb = [0, 1, 1, 1, 1, -pi, 0];
ub = [Inf, size(crop_im, 2), size(crop_im, 1), size(crop_im, 2), size(crop_im, 1), pi, max(zData(:))];

% Perform the fitting
params = lsqcurvefit(gaussianFunc, initialGuess, {xData, yData}, zData, lb, ub);

fittedImage = gaussianFunc(params, {xData, yData});

[zsigma, zmu] = gaussfit( [1:size(beadImg,3)], rescale(intens),std(intens),focus_z);
fittedZ = gaussianFunc1d([zmu, zsigma],[1:size(beadImg,3)]);

figure;
subplot(1,3,1);
imshow2(crop_im,[])
subplot(1,3,2);
imshow2(fittedImage,[])
title(['sigmaX : ' num2str(params(5)*pixelSz,2) ' \mum ,sigmaY : ' num2str(params(4)*pixelSz,2) ' \mum'])
subplot(1,3,3);
plot([1:size(beadImg,3)]*zstep,rescale(intens));
hold all
plot([1:size(beadImg,3)]*zstep,fittedZ)
title(['sigmaZ : ' num2str(zsigma*zstep) ' \mum'])

