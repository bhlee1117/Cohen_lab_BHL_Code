function [superlocColormov dtimg dFF]=generate_SNAPTmov(STAmovie,mask,StrImg,tformReg)

if nargin<4
    tformReg=[];
end

if size(mask,3)>1
    error('Check the mask')
end
spikeMov = STAmovie;
dSpikeMov = spikeMov - mean(spikeMov,3);  % look at the deviations from the mean

bkgImg=mean(STAmovie(:,:,end-5:end),3);
dFMov=STAmovie-bkgImg;
APamp=imgaussfilt(max(dFMov,[],3),5);
dFF=max(APamp./imgaussfilt(bkgImg,5),[],3);

% perform spatial and PCA filtering to clean up movie
filtsize = 7; sigma = 3;
spikemovfilt1 = spatialfilt(dSpikeMov, filtsize, sigma);  % Spatially filter to remove speckle noise
spikemovfilt1(isnan(spikemovfilt1))=0;
[spikemovfilt2, eigvecs, eigvals] = pcafilt(spikemovfilt1, 5);  % PCA filter to remove temporal noise

%%
disp('fit a spline to the data (computationally demanding)')
tic;
thresh = 0.1; dir = 1; %it was 0.3
dtimg = splinezeros(spikemovfilt2, thresh, dir);
dtimg = dtimg - nanmean(dtimg(:));
dtimg = set_edge(dtimg,5,NaN);
validImg=imfill(~isnan(dtimg),'holes');
toc;

if ~isempty(tformReg)
    dtimg_Reg=imwarp(dtimg, tformReg, 'OutputView', imref2d(size(StrImg)));
    spikeMov_Reg=imwarp(dSpikeMov, tformReg, 'OutputView', imref2d(size(StrImg)));
    dFF_Reg=imwarp(dFF, tformReg, 'OutputView', imref2d(size(StrImg)));
    validImg_Reg=double(imwarp(validImg, tformReg, 'OutputView', imref2d(size(StrImg))));
    validImg_Reg(validImg_Reg==0)=NaN;
    dtimg_Reg=dtimg_Reg.*validImg_Reg;

    dFF_Reg=dFF_Reg.*mask;
    dFF_Reg(dFF_Reg==0)=NaN;
else
    dtimg_Reg=dtimg;
    spikeMov_Reg=dSpikeMov;
    dFF_Reg=dFF;
end
%%
figure; clf;
nexttile([1 1])
dtimg_Reg=dtimg_Reg.*mask;
dtimg_Reg(dtimg_Reg==0)=NaN;
imshow2(dtimg_Reg-prctile(dtimg_Reg(:),5), [0 3]); title('Timing')
hold all
Maskboundary = cell2mat(bwboundaries(mask));
plot(Maskboundary(:,2),Maskboundary(:,1),'r.')
colormap('turbo')
nexttile([1 1])
imshow2(dFF_Reg,[]); hold all
plot(Maskboundary(:,2),Maskboundary(:,1),'r.')
colormap('turbo')



subframeT = 0.025; % ms
initialT = -2; % ms
finalT = 2; % ms
sigma = 0.05; % ms.  This is how long the flash lasts at each pixel.  Adjust to get a reasonable-looking movie.
%stdimg = std(spikeMov, [], 3);
stdimg = mean(spikeMov_Reg, 3);
[ysize, xsize, ~] = size(spikeMov_Reg);
times = initialT:subframeT:finalT;
nSNAPT = length(times);

GaussPeaksmov = zeros(ysize,xsize,nSNAPT);
for q = 1:nSNAPT
 GaussPeaksmov(:,:,q) = exp(-(dtimg_Reg-times(q)*ones(ysize,xsize)).^2/(2*sigma^2)).*double(mask); %mult'iplying mask to keep it at the neuron
end

superlocmov = GaussPeaksmov.*repmat(mat2gray(stdimg), [1 1 nSNAPT]).*mat2gray(dFF_Reg);  
superlocmov(isnan(superlocmov))=0;

superlocColormov = zeros(ysize, xsize, 3, nSNAPT);
Colormov = zeros(ysize,xsize,3);

for j = 1:length(times)
    %Colormov(:,:,1) = superlocmov(:,:,j)*5; %Tune this to enhance this color in RGB
    Colormov = grs2rgb(double(superlocmov(:,:,j)*4),colormap("hot"));
    superlocColormov(:,:,:,j) = grs2rgb(double(StrImg),colormap(gray))+Colormov;
end
end
