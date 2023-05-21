function [mov_mc,flow_xy]= optical_flow_motion_correction(mov)

rows = mov.rows;
cols = mov.cols;
ts = mov.frames;

xField = zeros(rows,cols,ts,'single');
yField = zeros(rows,cols,ts,'single');
% figure(1);clf
mov_mc = zeros(rows,cols,ts,'single');
mov_test = mov.data;
mov_test = single(mov_test)./single(max(mov_test(:)));
mov_test = movmean(mov_test,10,3);
mov_temp = squeeze(median(mov_test(:,:,51:150),3));
mov = mov.data;
tic

opticFlow = opticalFlowFarneback;
% parfor ii = 1:ts
for ii = 1:ts
    frame_ii = squeeze(mov_test(:,:,ii));
    flow = estimateFlow(opticFlow,mov_temp);
    flow = estimateFlow(opticFlow,frame_ii);
    xField(:,:,ii) = flow.Vx;
    yField(:,:,ii) = flow.Vy;
    
    mov_mc(:,:,ii) = imwarp(single(squeeze(mov(:,:,ii))),cat(3,flow.Vx,flow.Vy),'cubic');

end
flow_xy = cat(4,xField,yField);
toc