[~,dataDir] = uigetfile('*.*');
cd(dataDir)
time = 0;
mov = vm(dataDir,[round(time/1.27e-3)+1 round((time+10)/1.27e-3)]);
    % mov = vm(dataDir);
dt = 1.27e-3;
%%
[mov(1:mov.frames),flow_xy]=optical_flow_motion_correction_gpu(mov);
mov([1 end],:,:) = 0;
mov(:,[1 end],:) = 0;
%% 
[mov_mc,flow_xy]= optical_flow_motion_correction_gpu(mov(1:3000));
% [mov_mc,flow_xy]= optical_flow_motion_correction_gpu_subsampl(mov(401:1000),1);
mov_mc([1 end],:,:) = 0;
mov_mc(:,[1 end],:) = 0;
%%
mcTrace = squeeze(mean(flow_xy,[1 2]));
mcTrace_hi = mcTrace-movmean(mcTrace,50,1);
mov_res = SeeResiduals(single(mov.data),squeeze(mean(mov.data,[1 2])));
mov_res = SeeResiduals(mov_res,mcTrace);
mov_res = SeeResiduals(mov_res,mcTrace_hi);
mov_r_vec = tovec(mov_res);
%%
%frames = 1:mov.frames;
%mov_bin = imresize(mov_res(:,:,frames),1/2,'box');

% mov_bin_s = single(mov_bin);%-movmean(single(mov_bin),200,3);
mov_bin_s = single(mov_mc_crop)-movmean(single(mov_mc_crop),200,3);
mov_bin_s = movmean(mov_bin_s,10,3);
frames=1:size(mov_mc_crop,3);

tic
n_pc = 16;
 %[V_pc, D_pc,U_pc] = eigs(double(covMat),n_pc);
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

figure(10);clf
for ii = 1:16
    subplot(4,4,ii)
    imshow2(toimg(V_pc(:,idx_pc(ii)),size(mov_bin_s,1),size(mov_bin_s,2)),[prctile(V_pc(:,idx_pc(ii)),1) prctile(V_pc(:,idx_pc(ii)),99)])
    title(num2str(ii))
end
%%
idx_pc = [1:4];
[icasig, A,W] = fastica((V_pc(:,idx_pc)'*tovec(mov_bin_s)));
n_ic = size(W,1);
V_ic = (W * V_pc(:,idx_pc)')';
ic_sub = (V_ic(:,1:n_ic)'*tovec(mov_bin_s))';
mov_res = SeeResiduals(mov_mc_crop,ic_sub);
%%
% mov_bin = imresize(mov.data,1/2,'box');
for ii = 1:10
frames = 1:mov.frames;
mov_bin = imresize(mov_res(:,:,frames),1/2,'box');

% mov_bin_s = single(mov_bin);%-movmean(single(mov_bin),200,3);
mov_bin_s = single(mov_bin)-movmean(single(mov_bin),200,3);
mov_bin_s = movmean(mov_bin_s,10,3);
% mov_bin_s = squeeze(mean(reshape(mov_bin_s(:,:,1:floor(mov.frames/8)*8),size(mov_bin_s,1),size(mov_bin_s,2),[],8),4));

% covMat = tovec(mov_bin_s)*tovec(mov_bin_s)';  % PCA within each region

tic
n_pc = 16;
% [V_pc, D_pc,U_pc] = eigs(double(covMat),n_pc);
in = gpuArray(double(tovec(mov_bin_s)));
% [V_pc, D_pc,U_pc] = svds(double(tovec(mov_bin_s)),n_pc);
[V_pc, D_pc,U_pc] = svds(in,n_pc);
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
pause
end
%% PCA demix
% mov_r_vec = SeeResiduals(single(mov.data),[squeeze(mean(mov.data,[1 2])) mcTrace]);
% mov_r_vec = single(mov_r_vec)-movmean(mov_r_vec,200,3);
% mov_r_vec = tovec(mov_r_vec);
[ySize, xSize, nFrames] = size(mov_res);
ref_im = mov.mean;

idx_seg = 201:6000;
% Segment each ROI
stats2 = regionprops(L2, imgFilt, 'Area', 'MaxIntensity', 'MinIntensity', 'PixelIdxList');
nROI2 = length(stats2);

mov_res_vec = tovec(mov_res);

n_comp = 1;
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
idx_t = 100:size(traces,2)-100;
traces_plot = squeeze(traces) - movmean(squeeze(traces),250,1);
idx_range = 10:6000;
traces_plot = traces_plot ./ range(traces_plot(idx_range,:));
t_sk = 10:6000;
sk = squeeze(skewness(butterworth_filt(traces_plot(t_sk,:)-median(traces_plot(t_sk,:),2),4,[20,inf],800), [], 1));
[~, idx] = sort(sk, 'descend');
idx(isnan(sk(idx))) = [];
idx100 = idx(1:size(traces_plot,2));
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

%% use mask only

mov_r = single(mov_res)-movmean(mov_res,200,3);
[ySize, xSize, nFrames] = size(mov_r);
ref_im = mov.mean;

idx_seg = 201:7800;
% Segment each ROI
stats2 = regionprops(L2, imgFilt, 'Area', 'MaxIntensity', 'MinIntensity', 'PixelIdxList');
nROI2 = length(stats2);

mov_r_vec = tovec(mov_r);

n_comp = 1;
traces = zeros(1, nFrames, nROI2); % take 3 PCs per ROI 
c_ftprnt = zeros(numel(imgFilt),nROI2);
for j = 1:nROI2

    subMov = mov_r_vec(stats2(j).PixelIdxList,idx_seg);

    mask = c_ftprint(stats2(j).PixelIdxList,j);
    traces(:,:,j) = -mask' * mov_r_vec(stats2(j).PixelIdxList,:);
%     c_ftprnt(stats2(j).PixelIdxList,j) = mask;
    j
%     pause
end
% c_ftprnt = reshape(c_ftprnt,size(ref_im,1),size(ref_im,2),[]);

% save('results.mat','c_ftprnt','traces','s_idx','e_idx')

idx_t = 100:size(traces,2)-100;
traces_plot = squeeze(traces) - movmean(squeeze(traces),250,1);
idx_range = 2000:6000;
traces_plot = traces_plot ./ range(traces_plot(idx_range,:));
t_sk = 101:7900;
sk = squeeze(skewness(butterworth_filt(traces_plot(t_sk,:)-median(traces_plot(t_sk,:),2),4,[20,inf],800), [], 1));
[~, idx] = sort(sk, 'descend');
idx(isnan(sk(idx))) = [];
idx100 = idx(1:200);
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