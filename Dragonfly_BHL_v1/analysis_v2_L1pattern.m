clear;close all;
[fname,fpath] = uigetfile('*.*'); if fname == 0, return;end
mov = vm(fpath);
cd(fpath)
% mov=setSaturationLimits(mov,[0 700]);
figure;moviesc(mov,size(mov,3),'fixed')%imshow(mov.mean,[])
clicky_autobg(mov.data,mov(1:end/2).mean);
% vm2avi(mov(400:600),'400-600.avi',25)
%%
% [d1,d2,T] = mov.size; bound = 0;
% options_r = NoRMCorreSetParms('d1',d1-bound,'d2',d2-bound,'bin_width',200,'max_shift',20,'iter',1,'correct_bidir',false);
% 
% % register using the high pass filtered data and apply shifts to original data
% tic; [~,shifts1,template1] = normcorre_batch(single(mov.toimg.data),options_r); toc % register filtered data
% % exclude boundaries due to high pass filtering effects
% tic; mov_mc = apply_shifts(single(mov.toimg.data),shifts1,options_r,0,0); toc % apply shifts to full dataset
% figure;moviesc(mov_mc,1,'fixed')


fig=uifigure;
stim_type = uiconfirm(fig,'','','option',{'wide field','center surround','invert center surround','ramp center surround'},'closeFcn',@(h,e) close(fig));
switch stim_type
    case 'wide field'
        mov_2 = mov-mov.mean-100;
        rois = clicky(double(mov));
        cm_noise = apply_clicky(rois,double(mov_2.data),'no');

        mov_2 = vm(SeeResiduals(double(mov_2.data),cm_noise));

        rois = clicky(double(mov_2),mov.mean);

        wf_intensities = apply_clicky(rois,double(mov_2.data),'no');
    case 'center surround'
        rois = clicky(double(mov),mean(mov(1:500)));

%         loc_noise = mean(apply_clicky(rois,double(mov.data),'no'),2);

%         mov_2 = vm(SeeResiduals(double(mov.data),cm_noise));

%         rois = clicky(double(mov_2),mov(1:end).mean);
        
%         cs_intensities = apply_clicky(rois,double(mov.data),'no');
        
        cs_intensities_and_noise = apply_clicky(rois,double(mov.data),'no');
        cs_intensities=cs_intensities_and_noise(:,1)-cs_intensities_and_noise(:,2);
        cs_intensities(1:499) = cs_intensities(1:499)*cs_intensities(501)/cs_intensities(499);
        t = 0:.03:(length(cs_intensities)-1)*.03;
        figure;plot(t,cs_intensities)
    case 'invert center surround'
%         rois = clicky(double(mov(1:500)));
% 
%         cm_noise = mean(apply_clicky(rois,double(mov.dat a),'no'),2);
% 
%         mov_2 = vm(SeeResiduals(double(mov.data),cm_noise));
% 
%         rois = clicky(double(mov_2),mov(500:end).mean);
% 
%         ics_intensities = apply_clicky(rois,double(mov_2.data),'no');
        
        rois = clicky(double(mov),mean(mov(1:500)));
        
        ics_intensities_and_noise = apply_clicky(rois,double(mov.data),'no');
        ics_intensities=ics_intensities_and_noise(:,1)-ics_intensities_and_noise(:,2);
        ics_intensities(1:499) = ics_intensities(1:499)*ics_intensities(501)/ics_intensities(499);
        t = 0:.03:(length(ics_intensities)-1)*.03;
        figure;plot(t,ics_intensities)
    case 'ramp center surround'
        rois = clicky(double(mov(1:500)));

        cm_noise = mean(apply_clicky(rois,double(mov.data),'no'),2);

        mov_1 = vm(SeeResiduals(double(mov(50:333).data),cm_noise(50:333)));
        mov_2 = vm(SeeResiduals(double(mov(334:666).data),cm_noise(334:666)));
        mov_3 = vm(SeeResiduals(double(mov(700:end).data),cm_noise(700:end)));

        rois = clicky(double(mov_2),mov(1:end).mean);

        rcs_intensities = [apply_clicky(rois,double(mov_1.data),'no'); ...
                            apply_clicky(rois,double(mov_2.data),'no'); ...
                            apply_clicky(rois,double(mov_3.data),'no')];
        figure;plot(rcs_intensities)
end
%%
figure('Position',[0 0 2560 1440]); 
subplot(2,2,1)
imshow(mov.mean,[],'initialmagnification','fit')

subplot(2,2,2)
imshow(mov.mean,[],'initialmagnification','fit')

% subplot(2,3,3)
% imshow(mov,[],'initialmagnification','fit')
hold on
cellfun(@(x) plot(x(:,1),x(:,2)),rois);

subplot(2,2,3:4)
t = 0:.03:(length(intensities)-1)*.03;
% delF_F = (intensities-intensities(1,:))./intensities(1,:);
% offsets = range(delF_F,1);offsets(end) = 0;offsets = circshift(offsets,1);
% plot(t,delF_F+offsets*triu(ones(length(rois))))

offsets = range(intensities,1);offsets(end) = 0;offsets = circshift(offsets,1);
plot(t,intensities+offsets*triu(ones(length(rois))))
xlabel('Time (s)')
ylabel('Residual Counts')
% saveas(gca,'residual_clicky_1.png')
% saveas(gca,'delF_F_clicky_1.png')
figure;imagesc([0 30],[1:size(intensities,2)],intensities')
colorbar
%%
mov_50 = mov(1:50);
vm2avi(mov_50,'50_frames.avi',10)
% saveas(gca,'original.png')
% mov_2 = mov-mov.mean-100;
% figure;moviesc(mov_2)
%%
mov_2 = vm(SeeResiduals(mov_2,squeeze(mean(mov_2.data,[1 2])),0));
clear mov
% cd(fpath)
% mov_r.saveavi('residual.avi')

% figure;moviesc(mov_2)

% 
% [d1,d2,T] = mov.size; bound = 0;
% options_r = NoRMCorreSetParms('d1',d1-bound,'d2',d2-bound,'bin_width',200,'max_shift',20,'iter',1,'correct_bidir',false);
% 
% % register using the high pass filtered data and apply shifts to original data
% tic; [~,shifts1,template1] = normcorre_batch(single(mov.toimg.data),options_r); toc % register filtered data
% % exclude boundaries due to high pass filtering effects
% tic; mov_mc = apply_shifts(single(mov.toimg.data),shifts1,options_r,0,0); toc % apply shifts to full dataset
% 
% figure(2);clicky(double(mov.toimg.data)); figure(3);clicky(double(mov_mc))
%%
mov_2 = vm(uint16(mov_2.data));
mov_2 = mov_2(1:end-1).*mov_2(2:end);
% d_mov_bar_filt = d_mov.blur(1);
mov_2 = mov_2.mean;

% figure
% subplot(2,3,2)
% imshow(mov_2,[],'initialmagnification','fit')
% axis equal
% colormap(gray(256))
% set(gca,'clim',[mean(d_mov_bar_filt.mean,'all') mean(d_mov_bar_filt.mean,'all')+std(d_mov_bar_filt.mean,0,'all')])
%%
mov = vm(fpath);
mov = mov-mov.mean-100;
% figure;moviesc(mov_2)

mov = vm(SeeResiduals(mov,squeeze(mean(mov.data,[1 2])),0));
rois = clicky(double(mov),mov_2);

intensities = apply_clicky(rois,double(mov.data),'no');
% figure;
% t = 0:.03:(length(intensities)-1)*.03;
% offsets = range(intensities,1);offsets(end) = 0;offsets = circshift(offsets,1);
% plot(t,intensities+offsets*triu(ones(length(rois))))
% xlabel('Time (s)')
% ylabel('Residual Counts')
%%

figure('Position',[0 0 2560 1440]); 
subplot(2,3,1)
imshow(mov.mean,[],'initialmagnification','fit')

subplot(2,3,2)
imshow(mov_2,[],'initialmagnification','fit')

subplot(2,3,3)
imshow(mov_2,[],'initialmagnification','fit')
hold on
cellfun(@(x) plot(x(:,1),x(:,2)),rois);

subplot(2,3,4:6)
t = 0:.03:(length(intensities)-1)*.03;
offsets = range(intensities,1);offsets(end) = 0;offsets = circshift(offsets,1);
plot(t,intensities+offsets*triu(ones(length(rois))))
xlabel('Time (s)')
ylabel('Residual Counts')
saveas(gca,'residual_std_ref_clicky_1.png')