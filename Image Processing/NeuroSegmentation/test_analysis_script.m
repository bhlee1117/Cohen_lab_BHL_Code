% path = 'X:\Lab\Data and Analysis\Electrochromic Protein\coke_can_microscope\2013-10-07_TTX_incubated_neurons\Plate 2_T1_r1';
path = 'C:\Users\Adam\Desktop\500Hz_640_488_1x1 bin3__stack.tif';
mov = loadTIFFStack(path, 360, 1024, 549);
[ySize, xSize, nFrames] = size(mov);
avgImg = mean(mov, 3);
imshow2(avgImg, [])


% Define the frames when the red and blue lasers are on
startFrames = 1:600:10800;
offsetIndx = 1:548;
nreps = 18; 

goodFrames = repmat(offsetIndx, nreps, 1);
goodFrames = bsxfun(@plus, offsetIndx, startFrames');

% Only load the highest three intensities in each rep
toLoad = goodFrames([7:9 16:18], :);
toLoad = sort(toLoad(:));
[dat, ySize, xSize] = load_movie(path, toLoad(:));

% convert to 3-D array
mov = toimg(dat, ysize, xSize);

% look at the average
avgImg = mean(mov, 3);
imshow2(avgImg, [])

% look at the whole-field intensity trace
intens = squeeze(mean(mean(mov)));
plot(intens)

% filter the movie spatially, then downsample
sigma = 3;
filt = fspecial('Gaussian', [2*sigma 2*sigma], sigma);
movFilt = imfilter(mov, filt, 'replicate');
avgFiltImg = mean(movFilt, 3);
imshow2(avgFiltImg, [])

% down-sample movie
movDownSamp = movFilt(1:sigma:end,1:sigma:end,:);
[ySizeSmall, xSizeSmall, nFrames] = size(movDownSamp);
clicky(movDownSamp)


%%
% Want to remove spurious signals from the movie.  Run EITHER the first or
% second block of code below
mov2Dds = tovec(movDownSamp);  % ds = Down Sampled in space.
dmov2D = bsxfun(@minus, mov2Dds, mean(mov2Dds, 2));  % Subtract each pixel's mean.

% No filtering
mov2DFilt = dmov2D;

% *** OR ***

% project out average intensity from intensity trace at each pixel
% intens2D = mean(mov2Dds, 1);  % whole-field intensity
% dIntens = intens2D - mean(intens2D); % mean-subtracted whole-field intensity
% slope = (dmov2D*dIntens')/(dIntens*dIntens');  % correlate each pixel with whole-field intensity
% imshow2(toimg(slope, ySizeSmall, xSizeSmall), [])
% mov2DFilt = dmov2D -slope*dIntens;  % project out whole-field intensity
[dF, scaleImgs] = SeeResiduals(movDownSamp, intens, 1);
slope = scaleImgs(:,:,2);
mov2DFilt = tovec(dF);

clicky(dF, slope);

% *** OR ***

% to run both segments:
mov2Dds = mov2DFilt;

% high-pass filter in time
filtT = 30; % time window for median filtering
% mov2DFilt = mov2Dds - medfilt2(mov2Dds, [1 filtT], 'symmetric');  % median filtering is super slow--high pass is quicker.
mov2DFilt = mov2Dds - imfilter(mov2Dds, ones(1, filtT)/filtT, 'replicate');
clicky(toimg(mov2DFilt, ySizeSmall, xSizeSmall), slope);

%%

nEig = 100; % # of eigenvalues to calculate.
[u, s, v] = pca_eig(mov2DFilt, nEig);
stackplot(v(:,1:20))
% Look at the spatial eigenfunctions
figure(4); clf
for j = 1:20;
    subplot(5,4,j);
    imshow2(toimg(u(:,j), ySizeSmall, xSizeSmall), [])
end;
figure(5); clf
imshow2(v', [])

%%

% Insert a function to autocorrelate each PC and to estimate # of cells
Lags = 0:60;
nLags = length(Lags);
autoCorrs = zeros(nLags, nEig);
for j = 1:nLags;
    autoCorrs(j,:) = mean(v((1+Lags(j)):end,:).*v(1:(end-Lags(j)),:), 1);
end;
figure(6); clf
stackplot(autoCorrs(:,1:50));
SNR = autoCorrs(2,:)./autoCorrs(1,:);
plot(SNR)
nIcs = find(SNR < 0.05, 1, 'first');
nEigUse = nIcs+4;


% Try spatial ICA.  This works very well!
[icsSpace, mixmat, sepmat] = sorted_ica(u(:,1:nEigUse),nIcs);
icsTime = v(:,1:nEigUse)*sepmat';
figure(1); clf
stackplot(icsTime(:,1:min(30, nIcs)));  % plot up to 30 independent components

figure(2); clf
for j = 1:30;
    subplot(6,5,j);
    imshow2(toimg(icsSpace(:,j), ySizeSmall, xSizeSmall), []);
end;


%% Try temporal ICA.  This does not work so well.
[icsTime, mixmat, sepmat] = sorted_ica(v(:,1:nEigUse),nIcs);
if nIcs ~= size(icsTime, 2)
    'number of ics is less than requested'
    nIcs = size(icsTime, 2);
end;
figure(1); clf
stackplot(icsTime(:,1:min(30, nIcs)));  % plot up to 30 independent components
icsSpace = u(:,1:nEigUse)*sepmat';

figure(6); clf
for j = 1:20;
    subplot(5,4,j);
    imshow2(toimg(icsSpace(:,j), ySizeSmall, xSizeSmall), []);
end;

figure(8)  % Check my understanding of mixmat and sepmat
subplot(1,2,1); plot(v(:,1:nEigUse) - icsTime*mixmat')
subplot(1,2,2); plot(icsTime - v(:,1:nEigUse)*sepmat')



%% Try mixed ICA.  This works well for small mu.
mu = .2;
sigUse = [mu*v(:,1:nEigUse); (1-mu)*u(:,1:nEigUse)];   % combine the spatial and temporal parts
sigUse = bsxfun(@times, sigUse, 1./mean(sigUse.^2, 1).^.5);  % make sure each column has standard deviation 1. 
[ics, mixmat, sepmat] = sorted_ica(sigUse, nIcs);
icsTime = ics(1:nFrames,:);
icsSpace = ics(nFrames+1:end,:);

figure(5); clf
stackplot(icsTime(:,1:min(20, nIcs)));  % plot up to 20 independent components

figure(6); clf
for j = 1:20;
    subplot(5,4,j);
    imshow2(toimg(icsSpace(:,j), ySizeSmall, xSizeSmall), []);
end;


% Later use a spike-finding algorithm to determine # of spiking cells.
nCells = nIcs;


%% Move to high spatial resolution domain
cellBlobSmall = zeros(ySizeSmall, xSizeSmall, nCells);
cellImgsSmall = toimg(icsSpace(:,1:nCells), ySizeSmall, xSizeSmall);
for j = 1:nCells;
    oneImg = cellImgsSmall(:,:,j);
    % Watershed around region of maximum intensity
    L = watershed(-oneImg);
    props = regionprops(L, oneImg, 'maxIntensity'); 
    [~, indx] = max([props(:).MaxIntensity]);
    cellBlobSmall(:,:,j) = oneImg.*(L == indx).*(oneImg > 0);
%     imshow2(cellBlovSmall(:,:,j), [])
%     pause;
end;
overlaps = sign(tovec(cellBlobSmall)'*tovec(cellBlobSmall));  % Overlap matrix of cell images
imshow2(overlaps)

% project out components correlated with whole-field mean
dIntens = intens - mean(intens);
mov2D = tovec(mov);
dmov2D = bsxfun(@minus, mov2D, avgImg(:));
corrImg2D = (dmov2D*dIntens)/(dIntens'*dIntens);
fieldImg = toimg(corrImg2D, ySize, xSize); 
imshow2(fieldImg, [])
dmov2DnoBkg = dmov2D - corrImg2D*dIntens';
clear dmov2D

clicky(toimg(dmov2DnoBkg, ySize, xSize), toimg(corrImg2D, ySize, xSize));

% remove median-filtered baseline
% dmov2DFilt = dmov2DnoBkg - medfilt2(dmov2DnoBkg, [1 filtT], 'symmetric'); % too slow
tic
dmov2DFilt = dmov2DnoBkg - imfilter(dmov2DnoBkg, ones(1, filtT)/filtT, 'replicate');
toc
clear dmov2DnoBkg
clicky(toimg(dmov2DFilt, ySize, xSize), toimg(corrImg2D, ySize, xSize));

% normalize the independent components
icsNorm = bsxfun(@times, icsTime, 1./sum(icsTime.^2).^.5);
vTraces = zeros(nFrames, nCells);
flatAvgTraces = zeros(nFrames, nCells);
corrImgs = zeros(ySize, xSize, nCells);
weightImgs = zeros(ySize, xSize, nCells);
offsetImgs = zeros(ySize, xSize, nCells);
cellImgs = zeros(ySize, xSize, nCells);
for j = 1:nCells;
    j
    cellBlob = imresize(cellBlobSmall(:,:,j), [ySize xSize], 'nearest');
    cellImgs(:,:,j) = fieldImg.*sign(cellBlob);
    pixList = find(tovec(cellBlob));
    pixTraces = dmov2DFilt(pixList,:);
    % For each cell, project out all the other cells and then look at the residuals
    % I want to project into an orthogonal subspace of all but one IC.  I
    % think there is something better than the iterative approach, but I
    % need to look it up.

    for k = setdiff(find(overlaps(:,j)), j)  % for every cell that partially overlaps cell j
       pixTraces = pixTraces - (pixTraces*icsNorm(:,k))*icsNorm(:,k)';  % project out the component of the movie along that cell
    end;
    flatAvgTraces(:,j) = mean(pixTraces,1);
    [vTraces(:,j), corrimg, weightimg, offsetimg] = extractV2d(pixTraces, icsNorm(:,j));
    tmp = zeros(ySize, xSize);    tmp(pixList) = corrimg;
    corrImgs(:,:,j) = tmp;
    tmp = zeros(ySize, xSize);    tmp(pixList) = weightimg;
    weightImgs(:,:,j) = tmp;
    tmp = zeros(ySize, xSize);    tmp(pixList) = offsetimg;
    offsetImgs(:,:,j) = tmp;
end;

figure(10); clf
for j = 1:min(30, nCells);
    subplot(5, 6, j);
    imshow2(corrImgs(:,:,j), []);
end;
figure(11); clf
for j = 1:min(30, nCells);
    subplot(5, 6, j);
    imshow2(weightImgs(:,:,j), []);
end;
figure(12); clf
for j = 1:min(30, nCells);
    subplot(5, 6, j);
    imshow2(offsetImgs(:,:,j), []);
end;
figure(13); clf
for j = 1:min(30, nCells);
    subplot(5, 6, j);
    imshow2(cellImgs(:,:,j), []);
end;
figure(14); clf
stackplot(vTraces)
figure(15); clf
stackplot(icsTime)
figure(16); clf
stackplot(flatAvgTraces)











