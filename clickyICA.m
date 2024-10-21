function [ROI, PCA, ICA] = clickyICA(mov, refImg, nPCs,hideOutput)
% function [ROI, intens, eigVals, pcImgs] = clickyPCA(mov, refimg, nPCs);
% mov: input movie.
% refImg: optional image (black and white or color) on which to click to
% select ROIs.
% nPCs: # of principal components to compute
% Click to define a single ROI, right click to complete.
% Show the principal components in space and time domain over the ROI.
% Returns coordinates of ROI, s and the intensity traces.
%
% ROI: x and y coordinates of the boundary of the ROI
% intens: intensity traces of the principal components
% eigVals: eigenvalues of the principal components
% pcImgs: spatial components of the principal components
%
% AEC 23 March 2022
% Byung Hun Lee, 19 October 2024, Added ICA analysis

arguments
    mov
    refImg = []
    nPCs = 2
    hideOutput = 0
end
if isempty(refImg)
    refImg = mean(mov, 3);
end

nFrames = size(mov, 3);
[ySize, xSize] = size(refImg(:,:,1)); 
intens = [];
[x, y] = meshgrid(1:xSize, 1:ySize);


figure
imshow(refImg, [], 'InitialMagnification', 'fit')
title('Click one ROI')
hold on;
[xv, yv] = (getline(gca, 'closed'));
ROI = [xv, yv];
if size(xv,1) < 3  % exit loop if only a line is drawn
	return
end
inpoly = inpolygon(x,y,xv,yv);
plot(xv, yv, 'Linewidth', 1);

movVec = tovec(mov);
movVec = movVec(inpoly(:), :);
nPix = size(movVec, 1);

if nFrames >= nPix
    movVec_centered = movVec-mean(movVec,2);
    covMat = (movVec_centered*movVec_centered')./numel(movVec(1,:));
    [V, D] = eig(covMat);
    V = V(:,end:-1:1);
    D = diag(D);
    D = D(end:-1:1);
    eigVals = D(1:nPCs);
    
    eigImgsVec = zeros(ySize*xSize, nPCs);
    eigImgsVec(inpoly(:),:) = V(:, 1:nPCs);
    pcImgs = toimg(eigImgsVec, ySize, xSize);
    
    intens = eigImgsVec(inpoly(:),:)'*movVec;
else
    movVec_centered = movVec-mean(movVec,2);
    covMat = (movVec_centered'*movVec_centered)./numel(movVec(1,:));
    [V, D] = eig(covMat);
    V = V(:,end:-1:1);
    D = diag(D);
    D = D(end:-1:1);
    eigVals = D(1:nPCs);
    
    intens = V(:,1:nPCs)';
    
    pcSpace = movVec*V(:,1:nPCs);
    eigImgsVec = zeros(ySize*xSize, nPCs);
    eigImgsVec(inpoly(:),:) = pcSpace;
    pcImgs = toimg(eigImgsVec, ySize, xSize);
end

PCA.intens=intens;
PCA.eigVals=eigVals;
PCA.pcImgs=pcImgs;

[icsTrace, ~, sepmat]=sorted_ica(intens(1:nPCs,:)',nPCs);
icsSpace = movVec*icsTrace;
icaImgsVec = zeros(ySize*xSize, size(icsTrace,2));
icaImgsVec(inpoly(:),:) = icsSpace;
icaImgs=toimg(icaImgsVec,ySize,xSize);
icsTrace=icsTrace';

ICA.intens=icsTrace;
ICA.icaImgs=icaImgs;

if ~hideOutput
figure;
tiledlayout(nPCs,nPCs);

    for j = 1:nPCs
        nexttile(1+(j-1)*nPCs,[1 1])
        imshow2(pcImgs(:,:,j), []);
        nexttile(2+(j-1)*nPCs,[1 floor(nPCs/2)-1])
        plot(mat2gray(intens(j,:))-(j-1)); hold all
    end

    for j = 1:size(icaImgs,3)
        nexttile(floor(nPCs/2)+1+(j-1)*nPCs,[1 1])
        imshow2(icaImgs(:,:,j), []);
        nexttile(floor(nPCs/2)+2+(j-1)*nPCs,[1 floor(nPCs/2)-1])
        plot(mat2gray(icsTrace(j,:))-(j-1)); hold all
    end
    hold off
end
end
