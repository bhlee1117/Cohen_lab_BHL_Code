function mask = get_roi_mask(mov,start_frame,frame_duration,seg_size)
if isempty(start_frame),start_frame=1;end
if isempty(frame_duration),frame_duration=40;end
if isempty(seg_size),seg_size=100;end

frame_size = seg_size;
std_lim = frame_duration;
figure;moviesc(mov)
%%
all_rois_mask = zeros(size(mov(1)));


id11 = [1:frame_size:size(mov,1)];
id12 = id11+seg_size-1; id12(id12>size(mov,1))=size(mov,1);
id12 = unique(id12);id11 = id11(1:length(id12));
id21 = [1:frame_size:size(mov,2)];
id22 = id21+seg_size-1; id22(id22>size(mov,2))=size(mov,2);
id22 = unique(id22);id21 = id21(1:length(id22));

for i=1:length(id11)
    for j=1:length(id21)

        Y = vm(mov.data(id11(i):id12(i),id21(j):id22(j),start_frame:min(100,mov.frames)));
        if Y.mean==0,continue;end
        Y2 = Y-100-Y.mean;


        im0 = mat2gray(mean(Y2(start_frame:std_lim-1).*Y2(start_frame+1:std_lim)));


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
    %     figure(203);clf;imshow(im0)
    %     figure(201);clf;histogram(im0);hold on;

        medfilt2_opts = 'symmetric';
        im_bg = medfilt2(im0, [1 1]*params.background_size.value,medfilt2_opts);
        im0 = imsubtract(im0, im_bg);
        figure(204);clf;imshow(im0,[])

    %     im0 = imfilter(im0, ...
    %         fspecial('gaussian', round(params.gaussian_sigma.value*3), params.gaussian_sigma.value));
        im0 = cutoff_filter(im0,'bandpass',[sqrt(numel(im0))/50,sqrt(numel(im0))/30]);
    %     im0 = cutoff_filter(im0,'lowpass',[sqrt(numel(im0))/30]);
    
        im_seg = zeros(size(mov(1)));
        im_seg(id11(i):id12(i),id21(j):id22(j))=im0;
        figure(199);clf;imshow(im_seg,[])

    %     figure(202);clf;imshow(im_bg,[])

        im0 = mat2gray(im0);
    %     
        [f,x]=ecdf(im0(:)); cfit = fit(x(1:50:end),f(1:50:end),'cubicinterp');x_f = linspace(x(1),x(end),100);
        cfit_d = fit(x_f(1:end-1)',diff(feval(cfit,x_f)),'smoothingspline');x_f = linspace(x(1),x(end),1000);f_d=feval(cfit_d,x_f);
        c1_idx = find(f_d==max(f_d))-1+find(f_d(find(f_d==max(f_d)):end)<max(f_d)*.8 & f_d(find(f_d==max(f_d)):end)>max(f_d)*.2);
        cfit_lin = fit(x_f(c1_idx)',f_d(c1_idx),'poly1'); a = coeffvalues(cfit_lin);
%         figure(201);clf;plot(x_f,f_d);hold on;plot((f_d(end)-a(2))/a(1),feval(cfit_d,(f_d(end)-a(2))/a(1)),'o')
        params.minimum_signal.value = (f_d(end)-a(2))/a(1);



    %     dlg_opt.WindowStyle = 'normal';
    %     params.minimum_signal.value = str2num(cell2mat(inputdlg('','',1,{''},dlg_opt)));


        [label_matrix, ~] = find_regions(im0, params);
        label_matrix_m = zeros(size(label_matrix));
        for ii=1:max(label_matrix,[],'all')
            roi = im0;
            roi(label_matrix~=ii)=0;
            thres = multithresh(roi,2); 
            if length(find(roi>thres(2)))> params.minimum_area.value
                label_matrix_m(roi>thres(2)) = ii;
            end
    %         pause
        end
        seg_rois_mask = zeros(size(mov(1)));
        seg_rois_mask(id11(i):id12(i),id21(j):id22(j))=label_matrix_m;
        figure(200);clf;imshow(seg_rois_mask,[])
        all_rois_mask = all_rois_mask|seg_rois_mask; 
%         pause
    end
    mask = all_rois_mask;
end