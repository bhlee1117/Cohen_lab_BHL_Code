function [mov_mc,flow_xy]= optical_flow_motion_correction_LBH(mov,mov_temp,method)

switch method
    case 'optic_flow'
rows = mov.rows;
cols = mov.cols;
ts = mov.frames;
disp('optical_flow in progress')
xField = zeros(rows,cols,ts,'single');
yField = zeros(rows,cols,ts,'single');
mov_mc = zeros(rows,cols,ts,'single');
  mov_test = mov.data;
  mov_test = single(mov_test)./single(max(mov_test(:)));
  mov_test = movmean(mov_test,10,3);
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
    case 'normcorre'
 mov=single(mov);
[d1,d2,T] = size(mov);                                % dimensions of dataset
d = d1*d2;          

gSig = 7; 
gSiz = 17; 
psf = fspecial('gaussian', round(gSiz), gSig);
ind_nonzero = (psf(:)>=max(psf(:,1)));
psf = psf-mean(psf(ind_nonzero));
psf(~ind_nonzero) = 0;   % only use pixels within the center disk
Y = imfilter(single(mov),psf,'same');
%Ypc = Yf - Y;
bound = 2*ceil(gSiz/2);


%total number of pixels
options_rigid = NoRMCorreSetParms('d1',d1-bound,'d2',d2-bound,'bin_width',200, ...
    'grid_size',[128,128],'mot_uf',4,'correct_bidir',false, ...
    'overlap_pre',32,'overlap_post',32,'max_shift',20);
    
tic; 
[M2,shifts2,template2] = normcorre(Y(bound/2+1:end-bound/2,bound/2+1:end-bound/2,:),options_rigid); 
toc 

tic; mov_mc = apply_shifts(mov,shifts2,options_rigid,bound/2,bound/2); toc 

shifts_nr = cat(ndims(shifts2(1).shifts)+1,shifts2(:).shifts);
shifts_nr = reshape(shifts_nr,[],ndims(Y)-1,T);
shifts_x = squeeze(shifts_nr(:,2,:))';
shifts_y = squeeze(shifts_nr(:,1,:))';

flow_xy=[shifts_x shifts_y];
% apply the shifts to the removed percentile


% register filtered data


%mov=vm(mov_mc);
% 
% rows = mov.rows;
% cols = mov.cols;
% ts = mov.frames;
% disp('normcorre in progress')
% xField = zeros(rows,cols,ts,'single');
% yField = zeros(rows,cols,ts,'single');
%  %figure(1);clf
% mov_mc = zeros(rows,cols,ts,'single');
%   mov_test = mov.data;
%   mov_test = single(mov_test)./single(max(mov_test(:)));
%   mov_test = movmean(mov_test,10,3);
% % mov_temp = squeeze(median(mov_test(:,:,51:150),3));
% mov = mov.data;
% tic
% 
% opticFlow = opticalFlowFarneback;
% % parfor ii = 1:ts
% for ii = 1:ts
%     frame_ii = squeeze(mov_test(:,:,ii));
%     flow = estimateFlow(opticFlow,mov_temp);
%     flow = estimateFlow(opticFlow,frame_ii);
%     xField(:,:,ii) = flow.Vx;
%     yField(:,:,ii) = flow.Vy;
%     
%     mov_mc(:,:,ii) = imwarp(single(squeeze(mov(:,:,ii))),cat(3,flow.Vx,flow.Vy),'cubic');
% 
% end
% flow_xy = cat(4,xField,yField);
% toc
end
end