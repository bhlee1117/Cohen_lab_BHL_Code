function [intens, eigVals, pcImgs] = maskPCA(mask, mov, nPCs)
% function [intens, eigVals, pcImgs] = maskPCA(mask, mov, nPCs);
% mask: region of movie to calculate PCA
% mov: input movie.
% nPCs: # of principal components to compute
% 
% Show the principal components in space and time domain over the ROI.
%
% intens: intensity traces of the principal components
% eigVals: eigenvalues of the principal components
% pcImgs: spatial components of the principal components
%
% AEC 23 March 2022

if isempty(nPCs);
    nPCs = 1;
end

nFrames = size(mov, 3);
[ySize, xSize] = size(mask); 
intens = [];
[x, y] = meshgrid(1:xSize, 1:ySize);

movVec = tovec(mov);
movVec = movVec(mask(:), :);
nPix = size(movVec, 1);

if nFrames >= nPix;
    covMat = movVec*movVec';
    [V, D] = eig(covMat);
    V = V(:,end:-1:1);
    D = diag(D);
    D = D(end:-1:1);
    eigVals = D(1:nPCs);
    
    eigImgsVec = zeros(ySize*xSize, nPCs);
    eigImgsVec(mask(:),:) = V(:, 1:nPCs);
    pcImgs = toimg(eigImgsVec, ySize, xSize);
    
    intens = eigImgsVec'*movVec;
else
    covMat = movVec'*movVec;
    [V, D] = eig(covMat);
    V = V(:,end:-1:1);
    D = diag(D);
    D = D(end:-1:1);
    eigVals = D(1:nPCs);
    
    intens = V(:,1:nPCs)';
    
    pcSpace = movVec*V(:,1:nPCs);
    eigImgsVec = zeros(ySize*xSize, nPCs);
    eigImgsVec(mask(:),:) = pcSpace;
    pcImgs = toimg(eigImgsVec, ySize, xSize);
end;

figure
imgIdx = ((1:nPCs)-1)*nPCs + 1;
plotIdx = setdiff(1:nPCs^2, imgIdx);
for j = 1:nPCs
    subplot(nPCs, nPCs, imgIdx(j));
    imshow2(pcImgs(:,:,j), []);
    subplot(nPCs, nPCs, plotIdx);
    plot(mat2gray(intens(j,:))-(j-1)); hold all
end;
hold off
