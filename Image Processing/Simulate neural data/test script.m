%% Start with a single cell, and try making spatial and temporal filters to extract the intensity trace with maximum SNR
cellCoords = [50, 50; 55, 55];  % x- and y-coordinates of the center of the cell
cellSize = [5, 5];  % Cell radius
spikeExponent = [20, 1];  % intensity made from rand^spikeExponent to generate spiky data.  Intensity is in [0, 1]
cellOffset = [1, 1];      % cell baseline brightness in A.U.
cellAmp = [1, 0];         % maximum amplitude of cell spikes in A.U.
bkg = 1;               % background brightness
scale = 100;           % multiplicative scale (important for determining size of shot noise
addNoise = 0;          % whether to add noise
ySize = 100;           % dimensions of the field of view
xSize = 100;
nFrames = 500;         % number of frames to simulate
[mov, trueImg, trueIntens] = sim_camera(cellCoords, cellSize, cellOffset, cellAmp, spikeExponent, bkg, scale, addNoise, ySize, xSize, nFrames);
nPix = ySize*xSize;

% Look at the data
avgImg = mean(mov, 3);
imshow2(avgImg, [])
clicky(mov);
intens = squeeze(mean(mean(mov)));
plot(intens)

% only keep the intensity of the first cell
trueIntens = trueIntens(1,:);
% only keep the first cell image
trueImg = trueImg(:,:,1);

% Calculate the SNR in the raw means
[spaceNoise, timeNoise, spaceErrMap, timeErrMap] = calc_noise(avgImg, trueImg, intens, trueIntens);
spaceNoise  % .43
timeNoise  % .66;
plot(timeErrMap)
imshow2(spaceErrMap, [])

%% Compare to the SNR from PCA
covmat = tovec(mov)'*tovec(mov);
[V1, D1] = eig(covmat);
V1 = fliplr(V1);
D1 = diag(D1); D1 = flipud(D1);
semilogy(D1)
stackplot(V1(:,1:4))
plot(V1(:,2), intens, '.')
R1 = tovec(mov)*V1;
for j = 1:4;
    subplot(2,2,j);
    imshow2(toimg(R1(:,j), ySize, xSize), []);
end;
[spaceNoise, timeNoise, spaceErrMap, timeErrMap] = calc_noise(avgImg, trueImg, intens, trueIntens);
spaceNoise
timeNoise
[spaceNoise, timeNoise, spaceErrMap, timeErrMap] = calc_noise(toimg(R1(:,2), ySize, xSize), trueImg, V1(:,2), trueIntens);
spaceNoise
timeNoise
% You do much better starting with PCA than starting with raw means.

% Calculate the noise.
nEig = 2;
D = R1(:,1:nEig); A = V1(:,1:nEig);
for j = 1:nEig;  % Make the eigenvectors have positive peaks.
    [~, maxIdx] = max(abs(D(:,j)));
    D(:,j) = D(:,j) * sign(D(maxIdx, j));
    A(:,j) = A(:,j) * sign(D(maxIdx, j));
end;
movRecon = D*A';
movNoise = tovec(mov) - movRecon;
SigmaD = var(movNoise, [], 2);
imshow2(toimg(SigmaD, ySize, xSize), [])  % noise
imshow2(toimg(D(:,2), ySize, xSize), [])  % signal

spaceNoise = zeros(21,21);
for j = 0:.25:5;
    for k = 0:.25:5;
        SNR = (abs(D(:,2)).^j)./(SigmaD.^k);
%     imshow2(toimg(SNR, ySize, xSize), [])  % SNR
%     imshow2(toimg(D(:,2).*SNR, ySize, xSize), [])  % weight image
       [spaceNoise(1 + 4*j, 1+4*k), ~, ~, ~] = calc_noise(D(:,2).*SNR, tovec(trueImg), [], []);
%     pause
    end;
end;
imshow2(spaceNoise, [])

j = 1.75; k = 1.25;
SNR = (abs(D(:,2)).^j)./(SigmaD.^k);
 imshow2(toimg(SNR, ySize, xSize), [])  % SNR
 imshow2(toimg(D(:,2).*SNR, ySize, xSize), [])  % weight image
[spaceNoise, ~, ~, ~] = calc_noise(D(:,2).*SNR, tovec(trueImg), [], [])


[~, timeNoise, ~, ~] = calc_noise([], [], A(:,2), trueIntens)
timeNoise = zeros(9,9);
for j = 0:.25:2;
    for k = 0:.25:2;
        DWeight = D.*((abs(D).^j)./((SigmaD*ones(1,nEig)).^k));
        Ahat = (inv(DWeight'*DWeight)*DWeight'*tovec(mov))';
        [~, timeNoise(1+4*j, 1+4*k), ~, ~] = calc_noise([], [], Ahat(:,2), trueIntens);
%         [j, k, timeNoise]
    end;
end;

% timeNoise == min(timeNoise(:))
[r c] = ind2sub(size(timeNoise),find(timeNoise == min(timeNoise(:))))
jmin = (r - 1)/4
kmin = (c - 1)/4

j = jmin; k = kmin; 
DWeight = D.*((abs(D).^j)./((SigmaD*ones(1,nEig)).^k));
Ahat = (inv(DWeight'*DWeight)*DWeight'*tovec(mov))';
plot(Ahat)
[~, timeNoise, ~, ~] = calc_noise([], [], Ahat(:,2), trueIntens)
% Use Ahat to get a better estimate of Dhat
% Assume that shot-noise is not time varying.  Then shot noise is a
% constant, which we can set to 1.
[spaceNoise, ~, ~, ~] = calc_noise(D(:,2), tovec(trueImg), [], [])
for j = 0:.25:2;
    AWeight = Ahat.*(abs(Ahat).^j);
    Dhat = (inv(AWeight'*AWeight)*AWeight'*tovec(mov)')';
    [spaceNoise(1 + 4*j), ~, ~, ~] = calc_noise(Dhat(:,2), tovec(trueImg), [], []);
%     [j spaceNoise]
end;
spaceNoise

imshow2(-toimg(Dhat(:,2), ySize, xSize), []);

movRecon = DWeight*Ahat';
movNoise = tovec(mov) - movRecon;

SigmaD = var(movNoise, [], 2);




Sigma = movNoise'*movNoise;
% SigmaInv = diag(1/mean(diag(Sigma))*ones(nFrames,1));
SigmaInv = diag(ones(nFrames,1));
imshow2(Sigma, [])
imshow2(SigmaInv, [])
plot(diag(Sigma))
% Use Scott's formula:
D = R1(:,1:2); A = V1(:,1:2);



% Dhat = inv(A'*inv(Sigma)*A)*A'*inv(Sigma)*tovec(mov)';
Dhat = inv(A'*SigmaInv*A)*A'*SigmaInv*tovec(mov)';
plot(Dhat(2,:)-D(:,2)', D(:,2), '.')
subplot(1,2,1); imshow2(toimg(-D(:,2), ySize, xSize), []);
subplot(1,2,2); imshow2(toimg(-Dhat(2,:), ySize, xSize), []);
[spaceNoise, ~, ~, ~] = calc_noise(toimg(D(:,2), ySize, xSize), trueImg, [], [])
[spaceNoise, ~, ~, ~] = calc_noise(toimg(Dhat(2,:), ySize, xSize), trueImg, [], [])



noiseImg = var(toimg(movNoise, ySize, xSize), [], 3);
imshow2(noiseImg, [])


%%
% Try using **SPARSE PCA**
% Start with a very small movie

cellCoords = [7, 7; 10, 10];  % x- and y-coordinates of the center of the cell
cellSize = [3, 3];  % Cell radius
spikeExponent = [20, 20];  % intensity made from rand^spikeExponent to generate spiky data.  Intensity is in [0, 1]
cellOffset = [1, 1];      % cell baseline brightness in A.U.
cellAmp = [1, 1];         % maximum amplitude of cell spikes in A.U.
bkg = 1;               % background brightness
scale = 100;           % multiplicative scale (important for determining size of shot noise
addNoise = 1;          % whether to add noise
ySize = 15;           % dimensions of the field of view
xSize = 15;
nFrames = 500;         % number of frames to simulate
[mov, trueImg, trueIntens] = sim_camera(cellCoords, cellSize, cellOffset, cellAmp, spikeExponent, bkg, scale, addNoise, ySize, xSize, nFrames);
nPix = ySize*xSize;
% Look at the data
avgImg = mean(mov, 3);
imshow2(avgImg, [])
clicky(mov);
intens = squeeze(mean(mean(mov)));
plot(intens)

movVec = tovec(mov);
dmov = movVec - mean(movVec, 2)*ones(1, nFrames);
imshow2(dmov, [])
dmovN = dmov./(std(dmov, [], 2)*ones(1, nFrames));
imshow2(dmovN, [])
std(dmovN, [], 2)
[B SV L D PATHS] = spca(dmovN', [], 4, .5, 4) ;
for j = 1:4;
    subplot(2,2,j)
    imshow2(toimg(B(:,j), 15, 15), [])
end;




%% try iterating extractV
nIter = 10;
Vouts = zeros(nIter + 1, nFrames);
Vouts(1,:) = intens;
tNoises = zeros(nIter + 1,1); % time-domain noise
tNoises(1) = timeNoise;
sNoises = zeros(nIter + 1,1);  % space-domain noise
sNoises(1) = spaceNoise;
corrImgs = zeros(ySize, xSize, nIter + 1); 
corrImgs(:,:,1) = avgImg;
for j = 2:nIter;
    [Vouts(j,:), corrImg,~,~] = extractV2d(tovec(mov), Vouts(j-1,:));
    corrImgs(:,:,j) = toimg(corrImg, ySize, xSize);
    [sNoises(j), tNoises(j), spaceErrMap, timeErrMap] = calc_noise(corrImgs(:,:,j), trueImg, Vouts(j,:), trueIntens);
end;
% The best we can hope to do is by training with the trueIntens
j = nIter + 1;
    [Vouts(j,:), corrImg2,~,~] = extractV2d(tovec(mov), trueIntens);
    corrImgs(:,:,j) = toimg(corrImg2, ySize, xSize);
    [sNoises(j), tNoises(j), spaceErrMap, timeErrMap] = calc_noise(corrImgs(:,:,j), trueImg, Vouts(j,:), trueIntens);

figure(2); clf
plot(1:nIter+1, tNoises, 1:nIter+1, sNoises)  % Iteration helps!
legend('Time noise', 'Space noise')

% Now try to get a better spatial image by preferentially weighting the
% correlation at times when the SNR is good
dmov = tovec(mov) - repmat(avgImg(:), [1 nFrames]);
imshow2(corrImgs(:,:,nIter), [])
[imgOut, corrTrace, weightTrace, offsetTrace] = extractV2d(dmov', corrImg);
figure(3); clf
plot(offsetTrace)
plot(weightTrace)
plot(1:nFrames, mat2gray(corrTrace), 1:nFrames, mat2gray(Vouts(nIter,:)))
plot(1:nPix, corrImg, 1:nPix, imgOut, 1:nPix, avgImg(:))
imshow2(toimg(imgOut, ySize, xSize), [])

[spaceNoise, timeNoise, spaceErrMap, timeErrMap] = calc_noise(corrImgs(:,:,j), trueImg, Vouts(nIter,:), trueIntens);
timeNoise
[spaceNoise, timeNoise, spaceErrMap, timeErrMap] = calc_noise(corrImgs(:,:,j), trueImg, corrTrace, trueIntens);
timeNoise

%% Try a ping-pong approach, alternating between spatial and temporal extractV:
nIter = 6;
Vouts = zeros(nIter + 1, nFrames);
Vouts(1,:) = intens;
corrImgs = zeros(nIter + 1, nPix); 
corrImgs(1,:) = tovec(avgImg);
imgOuts = zeros(nIter+1, nPix);
imgOuts(1,:) = tovec(avgImg);
corrTraces = zeros(nIter+1, nFrames);
corrTraces(1,:) = intens;
dmov = tovec(mov) - repmat(avgImg(:), [1 nFrames]);
% dmovT = tovec(mov)' - repmat(intens, [1, nPix]);
for j = 2:nIter
    [Vouts(j,:),corrimg,weightimg,offsetimg] = extractV2d(dmov, corrTraces(j-1,:));
    corrImgs(j,:) = (dmov*Vouts(j,:)')'/nFrames;
    [imgOuts(j,:),corrtrace,weighttrace,offsettrace] = extractV2d(dmov', corrImgs(j,:));
    corrTraces(j,:) = (imgOuts(j,:)*dmov)/nPix;
end;

imshow2(toimg(offsetimg, ySize, xSize), [])

plot(corrTraces')
plot(Vouts')
plot(imgOuts' - corrImgs')

% Calculate the noise from each round of ping-pong
sNoises = zeros(nIter,1); 
tNoises = zeros(nIter,1);
sNoises2 = zeros(nIter,1); 
tNoises2 = zeros(nIter,1);
for j = 1:nIter;
    [sNoises(j), tNoises(j), ~, ~] = calc_noise(corrImgs(j,:), trueImg(:)', Vouts(j,:), trueIntens);
    [sNoises2(j), tNoises2(j), ~, ~] = calc_noise(imgOuts(j,:), trueImg(:)', corrTraces(j,:), trueIntens);
end;
figure;
plot(1:nIter, sNoises, 1:nIter, tNoises, 1:nIter, sNoises2, 1:nIter, tNoises2)  % winners are: corrTrace, imgOuts


[Vout, corrimg, weightimg, offsetimg] = extractV2d(tovec(mov), intens);
corrimg = toimg(corrimg, ySize, xSize);
weightimg = toimg(weightimg, ySize, xSize);
offsetimg = toimg(offsetimg, ySize, xSize);
plot(Vout)
subplot(2,3,1); imshow2(corrimg, []); 
subplot(2,3,2); imshow2(weightimg, []);
subplot(2,3,3); imshow2(offsetimg, []);
subplot(2,3,4:6);
plot(1:nFrames, intens, 1:nFrames, Vout)

[spaceNoise, timeNoise, spaceErrMap, timeErrMap] = calc_noise(weightimg, trueImg, Vout, trueIntens);
spaceNoise
timeNoise





%% Now consider two fluctuating cells
cellCoords = [50, 50; 65 65];
cellSize = [20, 10];
spikeExponent = [20, 20];
cellOffset = [1, .4];
bkg = 1;
scale = 100;
addNoise = 1;
ySize = 100;
xSize = 100;
nFrames = 500;
[mov, trueImg, trueIntens] = sim_camera(cellCoords, cellSize, cellOffset, spikeExponent, bkg, scale, addNoise, ySize, xSize, nFrames);

avgImg = mean(mov, 3);
imshow2(avgImg, [])
clicky(mov);
intens = squeeze(mean(mean(mov)));
plot(1:nFrames, trueIntens)


%%
% Now test my cell finding algorithm!
% filter the movie spatially, then downsample
sigma = 5;
filt = fspecial('Gaussian', [2*sigma 2*sigma], sigma);
imgBdry = (ySize*mean(mean(mov(:,1,:))) + xSize*mean(mean(mov(1,:,:))) + ySize*mean(mean(mov(:,end,:))) + xSize*mean(mean(mov(end,:,:))))/(2*xSize + 2*ySize);
movFilt = imfilter(mov, filt, imgBdry);  % replace values outside of the array by the mean around all edges.
avgFiltImg = mean(movFilt, 3);
imshow2(avgFiltImg, [])
imshow2(var(movFilt(:,1:25,:), [], 3), [])

% down-sample movie
movDownSamp = movFilt(1:sigma:end,1:sigma:end,:);
[ySizeSmall, xSizeSmall, nFrames] = size(movDownSamp);
clicky(movDownSamp)
imshow2(var(movDownSamp(:,1:5,:), [], 3), [])


%%
% Want to remove spurious signals from the movie.  Run EITHER the first or
% second block of code below
mov2Dds = tovec(movDownSamp);  % ds = Down Sampled in space.
dmov2D = bsxfun(@minus, mov2Dds, mean(mov2Dds, 2));  % Subtract each pixel's mean.

% No filtering
mov2DFilt = dmov2D;

% *** OR ***

% project out average intensity from intensity trace at each pixel
intens2D = mean(mov2Dds, 1);  % whole-field intensity
dIntens = intens2D - mean(intens2D); % mean-subtracted whole-field intensity
slope = (dmov2D*dIntens')/(dIntens*dIntens');  % correlate each pixel with whole-field intensity
imshow2(toimg(slope, ySizeSmall, xSizeSmall), [])
mov2DFilt = dmov2D -slope*dIntens;  % project out whole-field intensity

dF = toimg(mov2DFilt, ySizeSmall, xSizeSmall);
clicky(dF, toimg(slope, ySizeSmall, xSizeSmall));

% *** OR ***

% to run both segments:
mov2Dds = mov2DFilt;

% high-pass filter in time
filtT = 30; % time window for median filtering
% mov2DFilt = mov2Dds - medfilt2(mov2Dds, [1 filtT], 'symmetric');  % median filtering is super slow--high pass is quicker.
mov2DFilt = mov2Dds - imfilter(mov2Dds, ones(1, filtT)/filtT, 'replicate');
clicky(toimg(mov2DFilt, ySizeSmall, xSizeSmall), toimg(slope,ySizeSmall, xSizeSmall));

%%

nEig = 30; % # of eigenvalues to calculate.
[u, s, v] = pca_eig(mov2DFilt, nEig);
figure(6); clf; stackplot(v(:,1:20))
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


nIcs = 6;
nEigUse = 8;



% Try spatial ICA.
[icsSpace, mixmat, sepmat] = sorted_ica(u(:,1:nEigUse),nIcs);
icsTime = v(:,1:nEigUse)*s(1:nEigUse,1:nEigUse)*mixmat; % Note: mixmat == sepmat'*inv(sepmat*sepmat')'
figure(1); clf
stackplot(icsTime(:,1:min(30, nIcs)));  % plot up to 30 independent components

figure(2); clf
for j = 1:min(30, nIcs);
    subplot(6,5,j);
    imshow2(toimg(icsSpace(:,j), ySizeSmall, xSizeSmall), []);
end;

figure(8); clf  % Check my understanding of mixmat and sepmat
subplot(1,2,1); plot(u(:,1:nEigUse) - icsSpace*mixmat')  % This one is approximate
subplot(1,2,2); plot(icsSpace - u(:,1:nEigUse)*sepmat')  % This one is exact


%% Try temporal ICA.
[icsTime, mixmat, sepmat] = sorted_ica(v(:,1:nEigUse),nIcs);
if nIcs ~= size(icsTime, 2)
    'number of ics is less than requested'
    nIcs = size(icsTime, 2);
end;
figure(1); clf
stackplot(icsTime(:,1:min(30, nIcs)));  % plot up to 30 independent components
icsSpace = u(:,1:nEigUse)*s(1:nEigUse,1:nEigUse)*mixmat; %  Note: mixmat == sepmat'*inv(sepmat*sepmat');

figure(6); clf
for j = 1:min(30, nIcs);
    subplot(6,5,j);
    imshow2(toimg(icsSpace(:,j), ySizeSmall, xSizeSmall), [-1 3]);
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
stackplot(icsTime(:,1:min(30, nIcs)));  % plot up to 20 independent components

figure(6); clf
for j = 1:30;
    subplot(6,5,j);
    imshow2(toimg(icsSpace(:,j), ySizeSmall, xSizeSmall), []);
end;

