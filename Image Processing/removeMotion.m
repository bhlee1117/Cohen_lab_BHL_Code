function [out, dX, dY] = removeMotion(in);
% function [out, dX, dY] = removeMotion(in);
% Remove motion artifacts from a movie by projecting against spatial
% gradients
%
% 1 June 2018 AEC and HT

[ySize, xSize, nFrames] = size(in);

avgImg = mean(in, 3);
gradX = avgImg(:,2:end) - avgImg(:,1:end-1);
gradX = [gradX, gradX(:,end)];  % keep the original dimensions
gradY = avgImg(2:end,:) - avgImg(1:end-1,:);
gradY = [gradY; gradY(end,:)];  % keep the original dimensions

xNorm = gradX(:)'*gradX(:);
yNorm = gradY(:)'*gradY(:);

dMov = in - repmat(avgImg, [1 1 nFrames]);

dMovV = tovec(dMov);

dX = gradX(:)'*dMovV/xNorm;
dY = gradY(:)'*dMovV/yNorm;

dMovVCorr = dMovV - gradX(:)*dX - gradY(:)*dY;
out = toimg(dMovVCorr, ySize, xSize) + avgImg;


