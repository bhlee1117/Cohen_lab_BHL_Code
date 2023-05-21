function [intens, eigVals, pcImgs] = maskPCA(mask, mov, nPCs,hideOutput)
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
% AEC 23 March 2022 (FPB cleaned up slightly and added hideOutput option
% since graphics display takes ~forever for large nPCs 13 Apr 2022)
arguments
    mask
    mov
    nPCs = 2
    hideOutput = 0
end

nFrames = size(mov, 3);
[ySize, xSize] = size(mask); 
%intens = [];
%[x, y] = meshgrid(1:xSize, 1:ySize);

movVec = tovec(mov);
movVec = movVec(mask(:), :);
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
    eigImgsVec(mask(:),:) = V(:, 1:nPCs);
    pcImgs = toimg(eigImgsVec, ySize, xSize);
    
    intens = eigImgsVec(mask(:),:)'*movVec;
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
    eigImgsVec(mask(:),:) = pcSpace;
    pcImgs = toimg(eigImgsVec, ySize, xSize);
end
if ~hideOutput
    figure
    imgIdx = ((1:nPCs)-1)*nPCs + 1;
    plotIdx = setdiff(1:nPCs^2, imgIdx);
    for j = 1:nPCs
        subplot(nPCs, nPCs, imgIdx(j));
        imshow2(pcImgs(:,:,j), []);
        subplot(nPCs, nPCs, plotIdx);
        plot(mat2gray(intens(j,:))-(j-1)); hold all
    end
    hold off
end
end
