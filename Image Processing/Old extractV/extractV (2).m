function [Vout, corrimg, weightimg, offsetimg] = extractV(imgs, Vin)
% function [Vout, corrimg, weightimg] = extractV(imgs, Vin)
% Calculates an estimated voltage, from a movie of input images and a
% control voltage (or estimate thereof)
% Uses a linear model at each pixel, i:
% S(i,t) = V(t)*m(i) + b(i) + e(i,t)
% where S(i,t) is the movie, m(i) is a position-dependent weight factor,
% b(i) is a position-dependent offset and e(i,t) is the noise.
% AEC 22 Jan 2011

avgimg = mean(imgs, 3);
avgV = mean(Vin);
dV = Vin - avgV;
L = length(Vin);

% subtract off background
for j = 1:L
    imgs(:,:,j) = imgs(:,:,j) - avgimg;
end;

% correlate the changes in intensity with the applied voltage
[ysize, xsize] = size(avgimg);
corrimg = zeros(ysize, xsize);
if ysize*xsize < L;
    for y = 1:ysize;
        for x = 1:xsize;
            corrimg(y,x) = mean(dV.*squeeze(imgs(y,x,:)));
        end;
    end;
else
    for j = 1:L;
        corrimg = corrimg + dV(j)*imgs(:,:,j);
    end;
    corrimg = corrimg/L;
end;
corrimg = corrimg/mean(dV.^2);

% calculate a dV estimate at each pixel, based on the linear regression.
imgs2 = zeros(size(imgs));
if ysize*xsize < L;
    for y = 1:ysize;
        for x = 1:xsize;
            imgs2(y,x,:) = imgs(y,x,:)./corrimg(y,x);
        end;
    end;
else
    for j = 1:L;
        imgs2(:,:,j) = imgs(:,:,j)./corrimg;
    end;
end;

% Look at the residuals to get a noise at each pixel
sigmaimg = zeros(ysize, xsize);
if ysize*xsize < L;
    for y = 1:ysize;
        for x = 1:xsize;
            sigmaimg(y,x) = mean((squeeze(imgs2(y,x,:))-dV).^2);
        end;
    end;
else
    for j = 1:L
        sigmaimg = sigmaimg + (imgs2(:,:,j) - dV(j)).^2;
    end;
    sigmaimg = sigmaimg/L;
end;

weightimg = 1./sigmaimg;
weightimg = weightimg/mean(weightimg(:));

dVout = zeros(L,1);
for j = 1:L
    dVout(j) = mean(mean(imgs2(:,:,j).*weightimg));
end

Vout = dVout + avgV;
offsetimg = avgimg - avgV*corrimg;

    
    