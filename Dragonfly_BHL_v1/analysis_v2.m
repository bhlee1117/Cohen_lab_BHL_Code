
[fname,fpath] = uigetfile('*.*'); if fname == 0, return;end
% im = imread(fullfile(fpath,fname));
% figure;imshow(im,[])
% mov1 = imread(fullfile(fpath,fname),1);
% mov = zeros([size(mov1) 100],'uint8');
% mov(:,:,1) = mov1;
% for i=2:100
%  
% mov(:,:,i) = imread(fullfile(fpath,fname),i);
% 

% end
% clicky_faster(mov)

% % mov(:,:,420:end)=[];
% mov_m = double(mean(mov,3));
% mov_std = double(std(single(mov),0,3));
% 
% [N,X,Y]=histcounts2(mov_m,mov_std,[500,500]);
% 
% figure;imagesc(N','XData',X,'ydata',Y);set(gca,'ydir','normal')



mov = vm(fpath);
%%
% mov = readBinMov(fullfile(fpath,fname),664,664);

% mov_hp = toimg(butterworth_filt(double(mov.tovec()),4,[20,inf],1/2.5e-3),mov.rows,mov.cols);
% figure(11);clf;moviesc(mov)
% mov = mov_125mm;

% ref_im = imread('145446p1_div18_FOV1_YQ0104_GFP.tif');
% figure;imshow(mov.mean,[],'initialmagnification','fit')

% mov = mov.tovec;
% [~,ind]=max(mov.data(:,1));
% max_p_trace = mov.data(ind,:);
% std(double(max_p_trace))/mean(max_p_trace)

% figure;plot(double(max_p_trace)/mean(max_p_trace))
% figure;imshow(mov.mean,[])
% ref = mean(mov(:,:,1:end-1).*mov(:,:,2:end),3);
% [row,col,~] = mov.size;
% mov = fft_clean(tovec(single(mov.data))');
% mov = toimg(mov',[row col]);
[~,intens]=clicky_faster(mov);
% 
% mov = mov-mov.mean;
% mov = SeeResiduals(mov,squeeze(mean(mov.data,[1 2])));
%  figure;moviefixsc(mov)

% mov_std = mov(101:200);
% mov_std = mov_std-mov_std.mean;
% mov_std = mov(401:end)-mov(401:end).mean;
% mov_std = mov_std(1:end-1).*mov_std(2:end);
% % % cd(fpath)
% mov2 = mov-mov.mean;
% mov2 = mov2(1:end-1).*mov2(2:end);

% figure;plot(squeeze(mean(mov(:,:,250:end).data,[1 2])))
% [~,intens]=clicky_faster(mov,ref_im);
% [~,intens]=clicky_faster(mov,mov_std.mean);
%%
% figure;moviefixsc(mov-mov(1:800).mean)
% figure;plot(squeeze(mean(mov.data,[1 2])))
% clear mov
% figure;plot(fft_clean(intens(101:end,:)))
ylim([100 max(intens(:))*1.1])
% dt = 1e-3;
% L = length(intens);
% intens = intens(L*5e-2+1:L*95e-2);
% figure(98);clf;plot(intens)
% % calculate fft
%  L = length(intens);
% X = fft(intens); 
% P = abs(X/L); P = P(1:L/2+1); P(2:end-1) = 2*P(2:end-1);
% F = 1/dt*(0:L/2)/L;
% [amp,loc] = findpeaks(P(F>400 & F<450),F(F>400 & F<450),'NPeaks',1,'sortstr','descend');
% figure(99);clf;plot(F,P);
% clear mov
%%
[b,a] = butter(4,10/(1000/2),'high');
intens_filt = filtfilt(b,a,intens);
% noise = [];
% for i=1:10
%     noise = [noise std(intens(i*1000-100:i*1000))];
% end
% mov_m = double(mean(mov.data,3));
% mov_std = double(std(single(mov.data),0,3));
% 
% [N,X,Y]=histcounts2(mov_m,mov_std,[500,500]);
% 
% figure;imagesc(N','XData',X,'ydata',Y);set(gca,'ydir','normal')
% % figure
% clear mov
%%

frame_size = 200;
stitch_overlap = 0;
std_lim=40;
all_rois_mask = zeros(size(mov(1)));
stitch_std = [];
n_sub_mov = ceil(size(mov,1)/frame_size)*ceil(size(mov,2)/frame_size);
id11 = [1:frame_size-stitch_overlap:size(mov,1)];
id12 = id11+200-1; id12(id12>size(mov,1))=size(mov,1);
id12 = unique(id12);id11 = id11(1:length(id12));
id21 = [1:frame_size-stitch_overlap:size(mov,2)];
id22 = id21+200-1; id22(id22>size(mov,2))=size(mov,2);
id22 = unique(id22);id21 = id21(1:length(id22));

for i=1:length(id11)
    for j=1:length(id21)
%         Y = vm(zeros(size(mov)));
%     Y(id11(i):id12(i),id21(j):id22(j),:) = vm(mov.data(id11(i):id12(i),id21(j):id22(j),:));
    Y = vm(mov.data(id11(i):id12(i),id21(j):id22(j),1:min(100,mov.frames)));
    Y2 = Y-100-Y.mean;
%     figure;moviesc(Y2,1,'fixed')
%     Y3 = mean(Y2(1:std_lim-1).*Y2(2:std_lim));
%     stitch_std{i,j} = mat2gray(Y3);
%     stitch_std(id11(i):id12(i),id21(j):id22(j)) = mat2gray(mean(Y2(1:std_lim-1).*Y2(2:std_lim)));
    
    im0 = mat2gray(mean(Y2(1:std_lim-1).*Y2(2:std_lim)));
    im_seg = zeros(size(mov(1)));
    
    params.equalization_cliplim.value = 0.01;
    params.equalization_cliplim.on = false;

    params.background_size.value = 30;
    params.background_size.on = false;

    params.median_size.value = 21;
    params.median_size.on = false;

    params.gaussian_sigma.value = 2;
    params.gaussian_sigma.on = true;

    params.minimum_area.value = 7^2;
    params.minimum_area.on = true;

    params.maximum_area.value = 30^2;
    params.maximum_area.on = true;

    params.minimum_signal.on = true;

    % background subtraction
    figure(203);clf;imshow(im0)
%     figure(201);clf;histogram(im0);hold on;

    medfilt2_opts = 'symmetric';
    im_bg = medfilt2(im0, [1 1]*params.background_size.value,medfilt2_opts);
    im0 = imsubtract(im0, im_bg);
    figure(204);clf;imshow(im0,[])
    
%     im0 = imfilter(im0, ...
%         fspecial('gaussian', round(params.gaussian_sigma.value*3), params.gaussian_sigma.value));
    im0 = cutoff_filter(im0,'bandpass',[sqrt(numel(im0))/50,sqrt(numel(im0))/30]);
%     im0 = cutoff_filter(im0,'lowpass',[sqrt(numel(im0))/30]);
    
    im_seg(id11(i):id12(i),id21(j):id22(j))=im0;
    figure(199);clf;imshow(im_seg,[])
    
    figure(202);clf;imshow(im_bg,[])
    
    im0 = mat2gray(im0);
    figure(201);clf; [f,x]=ecdf(im0(:)); cfit = fit(x(1:50:end),f(1:50:end),'cubicinterp');x_f = linspace(x(1),x(end),100);
    cfit_d = fit(x_f(1:end-1)',diff(feval(cfit,x_f)),'smoothingspline');x_f = linspace(x(1),x(end),1000);f_d=feval(cfit_d,x_f);
    c1_idx = find(f_d==max(f_d))+find(f_d(find(f_d==max(f_d)):end)<max(f_d)*.8 & f_d(find(f_d==max(f_d)):end)>max(f_d)*.2);
    cfit_lin = fit(x_f(c1_idx)',f_d(c1_idx),'poly1'); a = coeffvalues(cfit_lin);
    plot(x_f,f_d);hold on;plot((f_d(end)-a(2))/a(1),feval(cfit_d,(f_d(end)-a(2))/a(1)),'o')
    params.minimum_signal.value = (f_d(end)-a(2))/a(1);
    
    
    
%     dlg_opt.WindowStyle = 'normal';
%     params.minimum_signal.value = str2num(cell2mat(inputdlg('','',1,{''},dlg_opt)));
    
    
    [label_matrix, ~] = find_regions(im0, params);
    label_matrix_m = zeros(size(label_matrix));
    for ii=1:max(label_matrix,[],'all')
        roi = im0;
        roi(label_matrix~=ii)=0;
        thres = multithresh(roi,2);
        seg_I = imquantize(roi,thres);
        RGB = label2rgb(seg_I); 	 
        if length(find(roi>thres(2)))> params.minimum_area.value
            label_matrix_m(roi>thres(2)) = ii;
        end
%         pause
    end
    seg_rois_mask = zeros(size(mov(1)));
    seg_rois_mask(id11(i):id12(i),id21(j):id22(j))=label_matrix_m;
    figure(200);clf;imshow(seg_rois_mask,[])
    all_rois_mask = all_rois_mask|seg_rois_mask; 
%     pause
    end
end
%%
all_rois = bwboundaries(all_rois_mask);
figure;moviesc(mov);hold on
for k = 1:length(all_rois)
   boundary_k = all_rois{k};
   plot(boundary_k(:,2), boundary_k(:,1), 'r', 'LineWidth', 2)
end

%% multiple file mean
fpath = uigetdir();
flist = dir([fpath '\*YQ0303L_JFx608_DIV16_FOV20_cell13_1V_rep*']);

intens = [];
for ii=1:length(flist)
    mov = vm([fpath '\' flist(ii).name]);
    if ii==1
        [roi,~] = clicky_faster(mov);
    end
    intens = [intens apply_clicky_faster(roi,mov)];
    
end

figure; l=plot(1:length(intens),intens);
c = jet(length(l));
for ii=1:length(l);l(ii).Color = c(ii,:);end
%%
figure;plot(mean(intens,2))
%%
