function [mov_mc,flow_xy_out]= optical_flow_motion_correction_gpu(mov)
% for it = 1:3

rows = mov.rows;
cols = mov.cols;
ts = mov.frames;

n_colrep = 16;

ts_p = ceil(mov.frames/n_colrep)*n_colrep;

mov = padarray(mov.data,[0 0 ts_p-ts],0,'post');


% xField = zeros(rows,cols,ts,'single');
% yField = zeros(rows,cols,ts,'single');
% figure(1);clf
mov_mc = zeros(rows,cols,ts,'single');
mov_test = toimg(sgolayfilt(tovec(double(mov)),3,9,ones(9,1),2),rows,cols);
% mov_test = single(mov_test)./single(max(mov_test(:)));

flow_xy_out = zeros(rows,cols*n_colrep,ts_p/n_colrep,2);
n_avg = 1;
mov_test = reshape(movmean(mov_test,n_avg,3),rows,cols*n_colrep,[]);
if size(mov_test,3)>150
    mov_template = squeeze(median(mov_test(:,:,51:150),3));
else
    mov_template = squeeze(median(mov_test(:,:,1:10),3));
end
% mov = mov.data;

maxGPUBuffer = 4e+9/4;
maxT = floor(maxGPUBuffer/rows/cols/4/n_colrep);
tic

% opticFlow = opticalFlowFarneback;
mov_template = gpuArray(single(mov_template/2^8));
for ii = 1:maxT:(ts_p/n_colrep)
    clear movGpu flow_xy
    movGpu = gpuArray(single(mov_test(:,:,ii:min((ii+maxT-1),ts_p/n_colrep))/2^8));
    nframes = size(movGpu,3);
    flow_xy = gpuArray(zeros(rows,cols*n_colrep,nframes,2));
%     tic
    for jj = 1:nframes
        flow = optical_flow_motion_correction_OCV_GPU(mov_template,squeeze(movGpu(:,:,jj)));
        flow_xy(:,:,jj,1) = flow(:,1:2:end);
        flow_xy(:,:,jj,2) = flow(:,2:2:end);
    end
%     toc
    flow_xy_out(:,:,ii:min((ii+maxT-1),ts_p/n_colrep),:) = gather(flow_xy);
end

flow_xy_out = reshape(flow_xy_out,rows,cols,[],2);
flow_xy_out = flow_xy_out(:,:,1:ts,:);
%     toc
for ii = 1:ts
    mov_mc(:,:,ii) = imwarp(single(squeeze(mov(:,:,ii))),-squeeze(flow_xy_out(:,:,ii,:)),'cubic');
%     mov(:,:,ii) = imwarp(single(squeeze(mov(:,:,ii))),-squeeze(flow_xy_out(:,:,ii,:)),'cubic');
end
% mov = vm(mov);
% end
toc
% mov_mc = mov.data;
clear movGpu flow_xy flow nframes mov_template
%%
% for ii = 1:13
%     mov_mc(:,:,ii) = imwarp(single(squeeze(mov(:,:,ii+107))),-squeeze(flow_xy_out(:,:,ii+107,:)),'cubic');
% end
% 
% 
% %%
% mov_mc = zeros(rows,cols,ts,'single');
% movGpu = gpuArray(single(mov(:,:,108:120))/2^8);
% flow_xy = gpuArray(zeros(rows,cols,13,2));
% for jj = 1:13-1
%         flow = optical_flow_motion_correction_OCV_GPU(squeeze(movGpu(:,:,jj)),squeeze(movGpu(:,:,jj+1)));
%         flow_xy(:,:,jj+1,1) = flow(:,1:2:end);
%         flow_xy(:,:,jj+1,2) = flow(:,2:2:end);
% end
% 
% flow_xy = gather(flow_xy);
% %%
% flow_xy(:,:,:,1) = toimg((triu(ones(13)) * tovec(flow_xy(:,:,:,1))')' ,rows,cols);
% flow_xy(:,:,:,2) = toimg((triu(ones(13)) * tovec(flow_xy(:,:,:,2))')' ,rows,cols);
% %%
% for ii = 1:13
%     mov_mc(:,:,ii) = imwarp(single(squeeze(mov(:,:,ii+107))),squeeze(flow_xy(:,:,ii,:)),'cubic');
% %     mov_mc(:,:,ii) = mov(:,:,ii+107);
% end
% %%
% 
% for ii = 1:13
%     mov_mc(:,:,ii) = imwarp(single(squeeze(mov(:,:,ii+107))),-squeeze(flow_xy(:,:,ii,:)),'cubic');
%     for jj = 2:ii
%         mov_mc(:,:,ii) = imwarp(single(mov_mc(:,:,ii)),-squeeze(flow_xy(:,:,jj,:)),'cubic');
%     end
% %     mov_mc(:,:,ii) = mov(:,:,ii+107);
% end
% %%
% 
% mov_template = gpuArray(single(mov(:,:,108)/2^8));
% movGpu = gpuArray(single(mov_mc(:,:,1:13))/2^8);
% flow_xy = gpuArray(zeros(rows,cols,13,2));
% for jj = 1:13-1
%         flow = optical_flow_motion_correction_OCV_GPU(mov_template,squeeze(movGpu(:,:,jj+1)));
%         flow_xy(:,:,jj+1,1) = flow(:,1:2:end);
%         flow_xy(:,:,jj+1,2) = flow(:,2:2:end);
% end
% 
% flow_xy = gather(flow_xy);
% %%
% for ii = 1:13
%     mov_mc(:,:,ii) = imwarp(mov_mc(:,:,ii),squeeze(flow_xy(:,:,ii,:)),'cubic');
% %     mov_mc(:,:,ii) = mov(:,:,ii+107);
% end