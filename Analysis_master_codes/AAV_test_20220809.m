% Test AAV expression, in house YQ201
% 2022/08/09 
clear
%filename='/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/20220808_YQ201_IH_bistability_first_explore/231910_YQ201_IH_M2_expression_p3';
[fname,fpath] = uigetfile('*.*','MultiSelect','on'); if fname == 0, return;end
mov=vm(fpath);
%mov=double(readBinMov_times(filename,304,512,[1:10000]));
moviefixsc(mov)
load(fullfile(fpath,'settings.mat'))
%%
mov_test=mov(:,:,150:250);
mov_test = single(mov_test)./single(max(mov_test.data(:)));
mov_test = movmean(mov_test,10,3);
mov_ref = squeeze(median(mov_test,3));
[mov_mc,xyField]=optical_flow_motion_correction_LBH(mov,mov_ref);

%%
[cell_mask]=clicky(mov_mc);
%Blue=interp1(DAQ_waves.time*1.25/1.27,DAQ_waves.amplitude(4,:),[1:size(mov,3)]*1.25*1e-3);
[m Blue]=clicky(mov_mc); Blue=lowpass(Blue,5,800);
mov_rm=SeeResiduals(double(mov_mc),Blue);
%%
figure;
volt=-apply_clicky(cell_mask,mov_rm,'no')';
scale=70;
for i=1:size(volt,1)
plot([1:size(volt,2)]*1.25*1e-3,volt(i,:)-median(volt(i,:))+i*scale)
hold all
end
yyaxis right
plot(DAQ_waves.time*1.25/1.27,DAQ_waves.amplitude(4,:),'b')
%%

%mov_mc = mov_mc - min(mov_mc(:));
%options_rigid = NoRMCorreSetParms('d1',size(mov_mc,1),'d2',size(mov_mc,2),'bin_width',200,'max_shift',30,'us_fac',50,'init_batch',200);
%tic; [mov_mc2,shifts1,template1,options_rigid] = normcorre(mov_mc,options_rigid); toc
mov([1 end],:,:) = 0;
mov(:,[1 end],:) = 0;

%% High-pass filter
mov_res=SeeResiduals(mov_mc,squeeze(mean(movmean(mov_mc,200,3),[1 2])));
im_G=imgaussfilt(mean(mov_mc,3),2);
[centers radii]=Cell_segment_circle_080222(im_G);
centers=cell_detection_manual(mean(mov_mc,3),centers);
%%
rad=5;
img_cent=zeros(size(mov_mc,1),size(mov_mc,2),size(centers,1));
for i=1:size(centers,1)
 [x y]=meshgrid([1:size(mov,2)],[1:size(mov,1)]);
 dist=sqrt(sum(([x(:) y(:)]-centers(i,[1 2])).^2,2));
 m=reshape(dist<rad,size(mov,1),[]);
mask_cent(i,:)=m(:);
end
D_cent=-bwdist(~(toimg(sum(mask_cent,1),size(mov,1),size(mov,2))>0));
L=watershed(D_cent); L(~(toimg(sum(mask_cent,1),size(mov,1),size(mov,2))>0))=0;
nROI = max(L(:))
imshow2(L,[])
%%
stats2 = regionprops(L, D_cent, 'Area', 'MaxIntensity', 'MinIntensity', 'PixelIdxList');
for j=1:nROI2; area(j)=stats2(j).Area; end; stats2(find(area<50))=[];
nROI2 = length(stats2);
mov_r_vec = tovec(mov_mc);
mov_res_vec = tovec(mov_res(:,:,1:3000));
clear traces c_ftprnt
n_comp = 3; g=1; nFrames=size(mov_mc,3);
traces = zeros(1, nFrames, nROI2); % take 3 PCs per ROI 
c_ftprnt = zeros(numel(im_G),nROI2);
for j = 1:nROI2
%     subMov = butterworth_filt(double(mov_r_vec(stats2(j).PixelIdxList,idx_seg)'),4,[10 inf],800)';  % take the pixels in the ROI

    subMov = mov_res_vec(stats2(j).PixelIdxList,idx_seg);

    covMat = subMov*subMov';  % PCA within each region
    [V, D] = eig(covMat);
    D = diag(D); 
    D = D(end:-1:1);
    V = V(:,end:-1:1);
    vSign = sign(max(V) - max(-V));  % make the largest value always positive
    V = V.*vSign;

%     mask = mat2gray(sum(V(:,1:n_comp).*D(1:n_comp)'.*(V(:,1:n_comp)>0),2));
    mask = mat2gray(mean(abs(V(:,1:n_comp)).*D(1:n_comp)',2));

    traces(:,:,j) = -mask' * mov_r_vec(stats2(j).PixelIdxList,:);
    c_ftprnt(stats2(j).PixelIdxList,j) = mask;
    
    j
%     pause

end;

c_ftprnt = reshape(c_ftprnt,size(ref_im,1),size(ref_im,2),[]);

% save('results.mat','c_ftprnt','traces')

%%
n=squeeze(sum(sum(c_ftprnt>0,1),2));
idx_t = 100:size(traces,2)-100;
traces_plot = squeeze(traces) - movmean(squeeze(traces),250,1);
idx_range = 10:7000;
traces_plot = traces_plot ./ range(traces_plot(idx_range,:));
t_sk = 10:7000;
sk = squeeze(skewness(butterworth_filt(traces_plot(t_sk,:)-median(traces_plot(t_sk,:),2),4,[20,inf],800), [], 1));
[~, idx] = sort(sk, 'descend');
idx(isnan(sk(idx))) = [];
idx(n(idx)<15)=[];
%idx100 = idx(1:150);
 idx100 = 1:size(traces_plot,2);
dt = 1.27e-3;

colr = flip(max(colormap(jet(length(idx100)))-.3,0),1);
% colr = flip(max(colormap(jet(size(traces_plot,2)))-.3,0),1);
figure;

traces_plot(:,idx100) = traces_plot(:,idx100)-median(traces_plot,2);
lines = plot(idx_t*dt,traces_plot(idx_t,idx100)+cumsum(range(traces_plot(idx_t,idx100))));
% lines = plot(idx_t*dt,traces_plot+cumsum(range(traces_plot)));
arrayfun(@(l,c) set(l,'Color',c{:}),lines,num2cell(colr,2))
% ylim([-1 102])
xlabel('Time (s)')
saveas(gcf,['skew' num2str(length(idx100)) '_pc'])
figure;imshow2(squeeze(sum(c_ftprnt(:,:,idx100).*reshape(colr,1,1,[],3),3)),[]);
% figure;imshow2(squeeze(sum(c_ftprnt.*reshape(colr,1,1,[],3),3)),[]);
% figure;imshow2(mat2gray(mov.mean)+squeeze(sum(c_ftprnt(:,:,idx100).*reshape(colr-max(mat2gray(mov.mean),[],'all'),1,1,[],3),3)),[]);
saveas(gcf,['skew' num2str(length(idx100)) '_map' '_pc'])


%%

load('STA_waveform.mat','sta_w'); 
%sta_w=(sta_w-min(sta_w))/(max(sta_w)-min(sta_w));
mov_gauss=imgaussfilt3(mov_mc,[1 1 0.1]);
a=normxcorr2(double(-sta_w)',tovec(mov_gauss(:,:,1:3000))); 
%m=movsum(tovec(mov_gauss),length(sta_w),2);
%a=a(:,length(sta_w)+1:end-length(sta_w))./m(:,floor(length(sta_w)/2)+2:end-floor(length(sta_w)/2)-1);
mov_bin_s=reshape(a,304,512,[]);%./max(mov_gauss,[],3);
%%

Fs = 800;            % Sampling frequency                    
T = 1/Fs;             % Sampling period       
L = 2899;             % Length of signal
t = (0:L-1)*T;
Y = fft(squeeze(sum(max(mov_bin_s2(:,:,100:3000),[],2),1)));
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;
plot(f,P1) 
xlabel('f (Hz)')
ylabel('Power')
xlim([0 50]); ylim([0 10])
%%

mcTrace = squeeze(mean(xyField,[1 2]));
mcTrace_hi = mcTrace-movmean(mcTrace,50,1);
mov_res = SeeResiduals(mov_mc,squeeze(mean(mov_mc,[1 2])));
%mov_res = SeeResiduals(mov_res,mcTrace);
%mov_res = SeeResiduals(mov_res,mcTrace_hi);
mov_r_vec = tovec(mov_res);

%%
mov=double(mov);
frames = 1:size(mov,3);
mov_bin = imresize(mov(:,:,frames),1/2,'box');
mov_bin_s = single(mov_bin)-movmean(single(mov_bin),200,3);
mov_bin_s = movmean(mov_bin_s,10,3);

tic
n_pc = 16;
%mov_bin_s=mov_bin_s(15:end-15,15:end-15,1:3000);
%mov_xcorr=toimg(tovec(mov_bin_s(:,:,:)).*tovec(mean(mov,3)),124,152);
%covMat = tovec(mov_bin_s(:,:,1:3000))*tovec(mov_bin_s(:,:,1:3000))';
 %[V_pc, D_pc,U_pc] = eigs(double(covMat),n_pc);
%in = gpuArray(double(tovec(mov_bin_s)));
[V_pc, D_pc,U_pc] = svds(double(tovec(mov_bin_s)),n_pc);
%[V_pc, D_pc,U_pc] = svds(in,n_pc);
clear in
V_pc = gather(V_pc);
D_pc = gather(D_pc);
U_pc = gather(U_pc);
 toc
D_pc = diag(D_pc); 
% D = D(end:-1:1);
% V = V(:,end:-1:1);
vSign = sign(max(V_pc) - max(-V_pc));   % make the largest value always positive
V_pc = V_pc.*vSign;
toc

idx_pc = 1:16;
% pc_sub = (V_pc(:,idx_pc)'*tovec(mov_bin_s))'; %PCA
% pc_sub = pc_sub-mean(pc_sub); %PCA
pc_sub = U_pc; %svds
dt = 1.27e-3;
figure(9);clf
plot((1:2900)*dt,pc_sub./range(pc_sub(501:end-500,:))+(1:length(idx_pc)))
ylim([0 length(idx_pc)+1])

figure(10);clf
for ii = 1:16
    subplot(4,4,ii)
    imshow2(toimg(V_pc(:,idx_pc(ii)),size(mov_bin_s,1),size(mov_bin_s,2)),[prctile(V_pc(:,idx_pc(ii)),1) prctile(V_pc(:,idx_pc(ii)),99)])
    title(num2str(ii))
end

%%
idx_pc = [2:3];
[icasig, A,W] = fastica((V_pc(:,idx_pc)'*tovec(mov_bin_s)));
n_ic = size(W,1);
V_ic = (W * V_pc(:,idx_pc)')';
ic_sub = (V_ic(:,1:n_ic)'*tovec(mov_bin_s))';
mov_bin_s = SeeResiduals(mov_bin_s,ic_sub);
%%
% mov_bin = imresize(mov.data,1/2,'box');
for ii = 1:10
frames = 1:size(mov,3);
mov_bin = imresize(mov_res(:,:,frames),1/2,'box');

% mov_bin_s = single(mov_bin);%-movmean(single(mov_bin),200,3);
mov_bin_s = single(mov_bin)-movmean(single(mov_bin),200,3);
mov_bin_s = movmean(mov_bin_s,10,3);
% mov_bin_s = squeeze(mean(reshape(mov_bin_s(:,:,1:floor(mov.frames/8)*8),size(mov_bin_s,1),size(mov_bin_s,2),[],8),4));

% covMat = tovec(mov_bin_s)*tovec(mov_bin_s)';  % PCA within each region

tic
n_pc = 16;
% [V_pc, D_pc,U_pc] = eigs(double(covMat),n_pc);
%in = gpuArray(double(tovec(mov_bin_s)));
 [V_pc, D_pc,U_pc] = svds(double(tovec(mov_bin_s)),n_pc);
%[V_pc, D_pc,U_pc] = svds(in,n_pc);
clear in
V_pc = gather(V_pc);
D_pc = gather(D_pc);
U_pc = gather(U_pc);
% toc
D_pc = diag(D_pc); 
% D = D(end:-1:1);
% V = V(:,end:-1:1);
vSign = sign(max(V_pc) - max(-V_pc));   % make the largest value always positive
V_pc = V_pc.*vSign;
toc

% plot PCs and footprint
idx_pc = 1:16;
% pc_sub = (V_pc(:,idx_pc)'*tovec(mov_bin_s))';
% pc_sub = pc_sub-mean(pc_sub);
pc_sub = U_pc;
dt = 1.27e-3;
figure(9);clf
plot((frames)*dt,pc_sub./range(pc_sub(501:end-500,:))+(1:length(idx_pc)))
ylim([0 length(idx_pc)+1])
% saveas(gcf,'mov_pc')
% figure;plot((1:mov.frames)*dt,pc_sub + cumsum(range(pc_sub(501:end,:))))

% plot PC footprint
figure(10);clf
for ii = 1:16
    subplot(4,4,ii)
    imshow2(toimg(V_pc(:,idx_pc(ii)),size(mov_bin_s,1),size(mov_bin_s,2)),[prctile(V_pc(:,idx_pc(ii)),1) prctile(V_pc(:,idx_pc(ii)),99)])
    title(num2str(ii))
end
% saveas(gcf,'mov_pc_img')
% Do ICA on first 16 PCs
idx_pc = [1:3];
[icasig, A,W] = fastica((V_pc(:,idx_pc)'*tovec(mov_bin_s)));
n_ic = size(W,1);
V_ic = (W * V_pc(:,idx_pc)')';
ic_sub = (V_ic(:,1:n_ic)'*tovec(mov_bin_s))';
% 
% figure;plot((frames)*dt,ic_sub./range(ic_sub(501:end,:))+(1:n_ic))
% 
% figure
% % V_ic = (ic_sub' * tovec(mov_bin_s(:,:,frames))')';
% for ii = 1:n_ic
%     subplot(4,4,ii)
% %     imshow2(toimg(V_ic(:,ii),size(mov_bin_s,1),size(mov_bin_s,2)),[]);%[prctile(V(:,ii),1) prctile(V(:,ii),99)])
%     imshow2(toimg(V_ic(:,ii),size(mov_bin_s,1),size(mov_bin_s,2)),[prctile(V_ic(:,ii),1) prctile(V_ic(:,ii),99)])
%     title(num2str(ii))
% end

mov_res = SeeResiduals(mov_res,ic_sub);
%pause
end


%%

[ySize, xSize, nFrames] = size(mov_mc);

gSig = 12; 
gSiz = 15; 
psf = fspecial('gaussian', round(gSiz), gSig);
ind_nonzero = (psf(:)>=max(psf(:,1)));
psf = psf-mean(psf(ind_nonzero));
psf(~ind_nonzero) = 0;   % only use pixels within the center disk
imgFilt = imfilter(squeeze(mean(std(mov_res,0,3),3)),psf,'same');
% imgFilt = imfilter(squeeze(mean(mov.data,3)),psf,'same');

ref_im = imgFilt;
% fluctImg(fluctImg<prctile(fluctImg(:),90)) = 0;
% subplot(2,1,2);
% imgFilt = imgaussfilt(ref_im,1);  % Smooth the reference image a bit
% imshow2(imgFilt,[])%[prctile(fluctImg(:),0.1), prctile(fluctImg(:),99.9)])

% Try simple intensity-based segmentation
% imgFilt = imfilter(fluctImg, fspecial('gaussian', [20 20], 4), 'replicate');  % Smooth the reference image a bit
% imshow2(imgFilt, [])
% imgFilt = mov.mean;
L2 = watershed(-imgFilt);  % find the watershed regions

 figure(7);clf
 imshow2(L2, [])

nROI = max(L2(:))
stats = regionprops(L2, imgFilt, 'Area', 'MaxIntensity', 'MinIntensity','weightedcentroid');  % get the properties of each region
areas = [stats(:).Area];
iMin = [stats(:).MinIntensity]';
iMax = [stats(:).MaxIntensity]';
ctr = cell2mat({stats(:).WeightedCentroid}');

%% PCA demix
 %mov_r_vec = SeeResiduals(single(mov.data),[squeeze(mean(mov.data,[1 2])) mcTrace]);
 %mov_r_vec = single(mov_r_vec)-movmean(mov_r_vec,200,3);
 %mov_r_vec = tovec(mov_r_vec);
mov_res2 = SeeResiduals(mov_mc(:,:,1:2900),squeeze(mean(mov_mc(:,:,1:2900),[1 2])));
%mov_res = SeeResiduals(mov_res,mcTrace);
%mov_res = SeeResiduals(mov_res,mcTrace_hi);
mov_r_vec = tovec(mov_res2(:,:,1:2900));


[ySize, xSize, nFrames] = size(mov_res);
ref_im = mean(mov,3);

idx_seg = 201:2900;
% Segment each ROI
stats2 = regionprops(L2, imgFilt, 'Area', 'MaxIntensity', 'MinIntensity', 'PixelIdxList');
nROI2 = length(stats2);

mov_res_vec = tovec(mov_res);
clear traces c_ftprnt
n_comp = 1; g=1;
traces = zeros(1, nFrames, nROI2); % take 3 PCs per ROI 
c_ftprnt = zeros(numel(imgFilt),nROI2);
for j = 1:nROI2
%     subMov = butterworth_filt(double(mov_r_vec(stats2(j).PixelIdxList,idx_seg)'),4,[10 inf],800)';  % take the pixels in the ROI

    subMov = mov_res_vec(stats2(j).PixelIdxList,idx_seg);

    covMat = subMov*subMov';  % PCA within each region
    [V, D] = eig(covMat);
    D = diag(D); 
    D = D(end:-1:1);
    V = V(:,end:-1:1);
    vSign = sign(max(V) - max(-V));  % make the largest value always positive
    V = V.*vSign;

%     mask = mat2gray(sum(V(:,1:n_comp).*D(1:n_comp)'.*(V(:,1:n_comp)>0),2));
    mask = mat2gray(mean(abs(V(:,1:n_comp)).*D(1:n_comp)',2));

    traces(:,:,j) = -mask' * mov_r_vec(stats2(j).PixelIdxList,:);
    c_ftprnt(stats2(j).PixelIdxList,j) = mask;
    
    j
%     pause

end;

c_ftprnt = reshape(c_ftprnt,size(ref_im,1),size(ref_im,2),[]);

% save('results.mat','c_ftprnt','traces')

%%
n=squeeze(sum(sum(c_ftprnt>0,1),2));
idx_t = 100:size(traces,2)-100;
traces_plot = squeeze(traces) - movmean(squeeze(traces),250,1);
idx_range = 10:2900;
traces_plot = traces_plot ./ range(traces_plot(idx_range,:));
t_sk = 10:2900;
sk = squeeze(skewness(butterworth_filt(traces_plot(t_sk,:)-median(traces_plot(t_sk,:),2),4,[20,inf],800), [], 1));
[~, idx] = sort(sk, 'descend');
idx(isnan(sk(idx))) = [];
idx(n(idx)<15)=[];
idx100 = idx(1:150);
% idx100 = 1:size(traces_plot,2);
dt = 1.27e-3;

colr = flip(max(colormap(jet(length(idx100)))-.3,0),1);
% colr = flip(max(colormap(jet(size(traces_plot,2)))-.3,0),1);
figure;

traces_plot(:,idx100) = traces_plot(:,idx100)-median(traces_plot,2);
lines = plot(idx_t*dt,traces_plot(idx_t,idx100)+cumsum(range(traces_plot(idx_t,idx100))));
% lines = plot(idx_t*dt,traces_plot+cumsum(range(traces_plot)));
arrayfun(@(l,c) set(l,'Color',c{:}),lines,num2cell(colr,2))
% ylim([-1 102])
xlabel('Time (s)')
saveas(gcf,['skew' num2str(length(idx100)) '_pc'])
figure;imshow2(squeeze(sum(c_ftprnt(:,:,idx100).*reshape(colr,1,1,[],3),3)),[]);
% figure;imshow2(squeeze(sum(c_ftprnt.*reshape(colr,1,1,[],3),3)),[]);
% figure;imshow2(mat2gray(mov.mean)+squeeze(sum(c_ftprnt(:,:,idx100).*reshape(colr-max(mat2gray(mov.mean),[],'all'),1,1,[],3),3)),[]);
saveas(gcf,['skew' num2str(length(idx100)) '_map' '_pc'])
