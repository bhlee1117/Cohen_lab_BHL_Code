
clear
clc;
cd '/Volumes/BHL_WD18TB/PP72_PlaceCellResults';
[~, fpath, raw] = xlsread(['/Volumes/BHL_WD18TB/' ...
    'PrismPCdata_Arrangement.xlsx'], 'Sheet1', 'C5:C7');

% [~, ~, NeuronsToUse]=xlsread(['/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/' ...
%     'PlaceCellData_Arrangement.xlsx'], 'Sheet1', 'L8:M42');
% 
% NeuronsToUse=cellfun(@(x) (str2num(num2str(x))),NeuronsToUse,'UniformOutput',false);

save_figto='/Volumes/BHL_WD18TB/PP72_PlaceCellResults';/Volumes/cohen_lab-1/Lab/Labmembers/Pojeong Park/Data/220525_PP071_P15_cre_10x_CA1/Cell3/195204PP071_P15_cell3_opto_soma/SNAPT

place_bin=150; time_segment=15000;


load(fullfile(save_figto,'Result_PC_Prism_20240106.mat'))
%%
 f=3;%1:3
staMov = rescale2(tovec(PC_Result{f}.STA_movie_align),1);
staMov = double(toimg(mean(staMov,2),size(PC_Result{f}.ref_im,1),size(PC_Result{f}.ref_im,2)));
spikeMov = -staMov(:,:,10:35);
dSpikeMov = spikeMov - mean(spikeMov,3);  % look at the deviations from the mean

% perform spatial and PCA filtering to clean up movie
filtsize = 7; sigma = 4;
spikemovfilt1 = spatialfilt(spikeMov, filtsize, sigma);  % Spatially filter to remove speckle noise
[spikemovfilt2, eigvecs, eigvals] = pcafilt(spikemovfilt1, 5);  % PCA filter to remove temporal noise

%%
% fit a spline to the data (computationally demanding)
thresh = 0.1; dir = 1; %it was 0.3
dtimg = splinezeros(spikemovfilt2, thresh, dir);
dtimg = dtimg - nanmean(dtimg(:));
%%

% stdimg = std(spikemovfilt2, [], 3);
% hist(stdimg(:), 50)  % Based on this standard deviation image, define a threshold for pixels that lie on the neuron
% 
% thresh = 0.003;
% mask = stdimg > thresh;
mask = mask_dendrite(avgimg,6,'Skeleton',0.2);
figure(12); clf
imshow(mask, []);
avgimg = mean(spikeMov,3);
%%
figure(14); clf; colormap('gray');
subplot(1,3,1); imshow2(avgimg, []); title('Average'); freezeColors;
subplot(1,3,2); imshow2(std(spikeMov, [], 3), []); title('Standard deviation'); freezeColors;
subplot(1,3,3); imshow2(dtimg.*mask, [-1 2]); title('Timing')
colormap('jet')

% saveas(gcf, 'Spike_timing_map.fig')
% saveas(gcf, 'Spike_timing_map.png')

%%
subframeT = 0.025; % ms
initialT = -2; % ms
finalT = 2; % ms
sigma = 0.05; % ms.  This is how long the flash lasts at each pixel.  Adjust to get a reasonable-looking movie.
%stdimg = std(spikeMov, [], 3);
stdimg = mean(spikeMov, 3);
[ysize, xsize, ~] = size(spikeMov);
times = initialT:subframeT:finalT;
nSNAPT = length(times);



GaussPeaksmov = zeros(ysize,xsize,nSNAPT);
for q = 1:nSNAPT
 GaussPeaksmov(:,:,q) = exp(-(dtimg-times(q)*ones(ysize,xsize)).^2/(2*sigma^2)).*double(mask); %mult'iplying mask to keep it at the neuron
end

superlocmov = GaussPeaksmov.*repmat(mat2gray(stdimg), [1 1 nSNAPT]);  



superlocColormov = zeros(ysize, xsize, 3, nSNAPT);
Colormov = zeros(ysize,xsize,3);
img_z=mat2gray(std(staMov,0,3)); img_z=img_z-medfilt2(img_z,[35 35]);
for j = 1:length(times)
    Colormov(:,:,1) = superlocmov(:,:,j)*5; %Tune this to enhance this color in RGB
    superlocColormov(:,:,:,j) = grs2rgb(double(img_z),colormap(gray))+Colormov;
end


close all
figure(20)
v = VideoWriter('SNAPT_movie','MPEG-4');

open(v);
ww = 1; %Frame counter

for j = 1:length(times)
    imshow(superlocColormov(:,:,:,j),[])
    pbaspect([size(double(superlocColormov(:,:,:,j)),2) size(double(superlocColormov(:,:,:,j)),1) 1]),colormap(gray)
    
    
    axis off
    text(2,20,[num2str(times(j)+1.7) 'ms'], 'FontSize', 20, 'color', [0.99 0.99 0.99])% the value 1. is to adjust timing by eyes       

    %annotation(gcf,'line',[0.7125 0.8281],...
    %[0.2 0.2],'Color',[1 1 1],'LineWidth',2);%length of bar is 46 pixels = 46*6.5/10 = 30 um
   

     pause(0.1)
    
    set(gcf,'color','w')    % Sets background to white
    frame = getframe(gcf);
    writeVideo(v,frame);
     ww = ww + 1;
    
end;
close(v);
