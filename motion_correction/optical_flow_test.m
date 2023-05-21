[~,dataDir] = uigetfile('*.*');
cd(dataDir)
mov = vm(dataDir);
% %%
% mov = single(mov./max(mov.data(:)));
%%

% rows = 101:200;
% cols = 151:250;
% ts = 1501:2500;

rows = 1:size(mov,1);
cols = 1:size(mov,2);
ts = 1:size(mov,3);

xField = zeros(length(rows),length(cols),length(ts),'single');
yField = zeros(length(rows),length(cols),length(ts),'single');
% figure(1);clf
mov_mc = zeros(length(rows),length(cols),length(ts),'single');
cnt = 1;
mov_test = mov.data(rows,cols,ts);
mov_test = single(mov_test)./single(max(mov_test(:)));
mov_test = movmean(mov_test,10,3);
mov_temp = squeeze(median(mov_test(:,:,51:150),3));
tic

% figure(1);clf
opticFlow = opticalFlowFarneback;
parfor ii = ts
%     clc
    frame_ii = squeeze(mov_test(:,:,ii));
    flow = estimateFlow(opticFlow,mov_temp);
    flow = estimateFlow(opticFlow,frame_ii);
    xField(:,:,ii) = flow.Vx;
    yField(:,:,ii) = flow.Vy;
    mov_mc(:,:,ii) = imwarp(single(squeeze(mov.data(rows,cols,ts(ii)))),cat(3,flow.Vx,flow.Vy),'cubic');
%     mov_mc(:,:,ii) = apply_shift_freq(single(squeeze(mov.data(rows,cols,ts(ii)))),cat(3,flow.Vx,flow.Vy));
%     clf
%     imshow2(mov_test(:,:,ii),[])
%     hold on
%     plot(flow,'DecimationFactor',[2 2],'ScaleFactor',10);
%     hold off
%     pause(10^-3)
%     pause
%     tic
%     if mod(ii,100) == 0
% %         timeNow = toc;
%         clc
%         fprintf('%d/%d, time elapsed %.1d seconds\n',ii,ts(end)); 
%     end
%     toc
end
toc
% mov_mc(:,:,1) = mov.data(rows,cols,ts(1));
%%
figure;moviefixsc(mov_mc,[],100)
% figure;moviefixsc(mov(rows,cols,ts))
%%
figure;moviefixsc(movmean(mov_mc-movmean(mov_mc,100,3),10,3),[-300 300])
%%
mov_res = regress_motion(mov_mc,xField,yField);
%%
figure;moviefixsc(mov_res,[],100)
%%
figure;moviefixsc(movmean(mov_res-movmean(mov_res,100,3),10,3),[-300 300])

%%
function mov_res = regress_motion(mov,shift_x,shift_y)
% do some simple segmentation    
    gSig = 10; 
    gSiz = 15; 
    psf = fspecial('gaussian', round(gSiz), gSig);
    ind_nonzero = (psf(:)>=max(psf(:,1)));
    psf = psf-mean(psf(ind_nonzero));
    psf(~ind_nonzero) = 0;   % only use pixels within the center disk
    Y = imfilter(single(squeeze(mean(mov,3))),psf,'same');
    Y2 = zeros(size(Y));
    Y2(11:end-10,11:end-10) = Y(11:end-10,11:end-10);
    [~,L] = bwboundaries(imbinarize(mat2gray(Y2)));
    pix = regionprops(L, Y2,'PixelIdxList');
    
% do motion regression on segmented regions    
mov_res = tovec(mov);
mov_res_cell = arrayfun(@(x) mov_res(x.PixelIdxList,:),pix,'uniformoutput',false);

shift_x = tovec(shift_x);
shift_y = tovec(shift_y);
shift_x_cell = arrayfun(@(x) mean(shift_x(x.PixelIdxList,:)),pix,'uniformoutput',false);
shift_y_cell = arrayfun(@(x) mean(shift_y(x.PixelIdxList,:)),pix,'uniformoutput',false);

% mov_res_cell = cellfun(@(mov,x,y,xy,xx,yy) SeeResiduals_vec(mov,[x;y],1,0) + mean(mov,2),...
%     mov_res_cell,shift_x_cell, shift_y_cell,...
%     'uniformoutput',false);
reg_coeffs_cell = cellfun(@(mov,x,y) SeeResiduals_vec(movmean(mov-mean(mov,2),50,2),[x;y],0,1),...
    mov_res_cell,shift_x_cell, shift_y_cell,...
    'uniformoutput',false);
mov_res_cell = cellfun(@(mov,x,y,c) mov-c*[x;y],...
    mov_res_cell,shift_x_cell, shift_y_cell, reg_coeffs_cell,...
    'uniformoutput',false);
    for ii = 1:length(pix)
        mov_res(pix(ii).PixelIdxList,:) = mov_res_cell{ii};
    end
    mov_res = toimg(mov_res,size(mov,1),size(mov,2));
end