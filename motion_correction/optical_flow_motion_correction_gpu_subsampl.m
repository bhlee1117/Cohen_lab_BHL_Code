function [mov_mc,flow_xy_out]= optical_flow_motion_correction_gpu_subsampl(mov,nsampl)

rows = mov.rows;
cols = mov.cols;
ts_total = mov.frames;
ts = mov.frames/nsampl;

if floor(ts)<ts
    error('nframes should be divisible by nsampl')
end


mov_mc = zeros(rows,cols,ts,'single');
mov_test = mov.data;
% mov_test = single(mov_test)./single(max(mov_test(:)));

% n_avg = 1;
mov_test = squeeze(mean(reshape(mov_test,rows,cols,nsampl,[]),3));
flow_xy_out = zeros(rows,cols,ts_total,2);
tic

    movGpu = gpuArray(single(mov_test/2^8));
%     nframes = size(movGpu,3);
    flow_xy = gpuArray(zeros(rows,cols,ts,2));

    for jj = 1:ts-1
        flow = optical_flow_motion_correction_OCV_GPU(squeeze(movGpu(:,:,jj)),squeeze(movGpu(:,:,jj+1)));
        flow_xy(:,:,jj+1,1) = flow(:,1:2:end);
        flow_xy(:,:,jj+1,2) = flow(:,2:2:end);
    end
    clear movGpu flow movGpu 
    flow_xy= gather(flow_xy);
flow_xy_out(:,:,:,1) = toimg((triu(ones(ts_total)) * interp1(1:nsampl:ts_total,tovec(flow_xy(:,:,:,1))',1:ts_total))',...
                        rows,cols);
flow_xy_out(:,:,:,2) = toimg((triu(ones(ts_total)) * interp1(1:nsampl:ts_total,tovec(flow_xy(:,:,:,2))',1:ts_total))',...
                        rows,cols);
    toc
for ii = 1:ts_total
%     if ii>1
%         flow_xy_out(:,:,ii,:) = sum(flow_xy_out(:,:,(1:ii),:),3);
%     end
    mov_mc(:,:,ii) = imwarp(single(squeeze(mov(:,:,ii))),squeeze(flow_xy_out(:,:,ii,:)),'cubic');
%     ii
end
toc

