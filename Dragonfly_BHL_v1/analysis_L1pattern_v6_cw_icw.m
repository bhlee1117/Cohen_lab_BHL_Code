intensities=[];
cw_idx = 1:2:20;
icw_idx = 2:2:20;

% clear;close all;
cur_dir = pwd;
[fname,fpath_cs] = uigetfile('*.*'); if fname == 0, return;end
% cd(fpath_cs)
mov = vm(fpath_cs);
%%
n_stim_frames=50;
n_frames =301;
%% motion correction
% mov_mc = motion_cor(mov,n_frames,n_stim_frames);
% mov=vm(mov_mc);
% figure;moviesc(mov_mc)
%% CW
% [fname_cal,fpath_cal] = uigetfile('*.*'); if fname == 0, return;end
% mov_cal = vm(fpath_cal);
% cd(fpath)
% mov=setSaturationLimits(mov,[0 700]);
% figure(1);clf;moviesc(mov)%imshow(mov.mean,[])
% [mask,roi] = get_roi_mask(mov(1:40),1,40,max(size(mov,[1 2])),.4);
% [intens_cs,roi]=apply_clicky_autobg(mov,roi,'no');
% if length(roi)>3
    
    
    idx=[];
    stim_type = 'cs';
    switch stim_type
        case 'cs'
            F0_idx = 1;
            for i=1:length(cw_idx)-1,idx = [idx (i)*n_frames*2+1:(i)*n_frames*2+1+9];end
        case 'ics'
            F0_idx = n_stim_frames+1;
            for i=1:length(icw_idx)-1,idx = [idx (i)*n_frames*2-n_frames+1:(i)*n_frames*2-n_frames+1+20];end
    end
    [roi,intens_all,mask]=clicky_autobg(mov,mov(idx).mean);
%     intens_cs=apply_clicky_autobg(mov,roi);
% end


sub_fig_size = 100;se = strel('square',100);win_idx = imdilate(mask,se);[rows,cols]=ind2sub(size(win_idx),find(win_idx));
figure(100);clf;
subplot(2,5,2);imshow(mov(1:n_stim_frames).mean,[])
subplot(2,5,[1 6])
imshow(mov.mean,[]);hold on;cellfun(@(x) plot(x(:,1),x(:,2),'r'),roi)
subplot(2,5,7)
im1 = mov.mean;
imshow(im1(min(rows):max(rows),min(cols):max(cols)),[],'initialmagnification','fit');hold on;cellfun(@(x) plot(x(:,1)-min(cols),x(:,2)-min(rows),'r'),roi)
title('center->wide field illum (zoomed-in view)')


intens_all = reshape([intens_all; nan],n_frames,[]);
intens_cw = intens_all(:,cw_idx);
t = 0:.03:(length(intens_all(:,1))-1)*.03;
% figure(10);clf;
subplot(2,5,[4 9]);plot(t,intens_cw)
subplot(2,5,[5 10]);plot(t,(intens_cw-intens_cw(F0_idx,1))./intens_cw(F0_idx,1))

% ICW

intens_icw = intens_all(:,icw_idx);

figure(100)
subplot(2,5,3);imshow(mov(n_frames+1:n_frames+n_stim_frames).mean,[])
subplot(2,5,8)
im1 = mov.mean;
imshow(im1(min(rows):max(rows),min(cols):max(cols)),[],'initialmagnification','fit');hold on;cellfun(@(x) plot(x(:,1)-min(cols),x(:,2)-min(rows),'r'),roi)
title('hollow center->wide field illum (zoomed-in view)')

% offsets = range(intensities,1);offsets(end) = 0;offsets = circshift(offsets,1);
% figure;plot(intensities+offsets*triu(ones(size(intensities,2))))
t = 0:.03:(length(intens_all(:,1))-1)*.03;


subplot(2,5,[4 9]);hold on;plot(t,intens_icw)
set(groot,'defaultAxesColorOrder',[winter(size(intens_cw,2));autumn(size(intens_icw,2))])
xlabel('Time (s)')
ylabel('Counts')
c = colorbar(gca,'Ticks',[0 1],'ticklabels',{'1','10'});
c.Location='northoutside';
c.Label.String = 'center->wide field trials';colormap(gca,'winter')

subplot(2,5,[5 10]);hold on;plot(t,(intens_icw-intens_cw(F0_idx,1))./intens_cw(F0_idx,1))
xlabel('Time (s)')
ylabel('{\Delta}F/F')
set(groot,'defaultAxesColorOrder',[winter(size(intens_cw,2));autumn(size(intens_icw,2))])
c = colorbar(gca,'Ticks',[0 1],'ticklabels',{'1','10'});
c.Location='northoutside';
c.Label.String = 'hollow center->wide field trials';colormap(gca,'autumn')
% figure(2);clf;moviesc(mov)%imshow(mov.mean,[])

% cd(fpath_cs)
% saveas(gca,'cs_wide')
% cd(cur_dir)
%%
function mov_mc = motion_cor(mov,n_frames,n_stim_frames)
mov_fake=mov;
for i=1:19
cs_idx = (i-1)*n_frames+1:(i-1)*n_frames+n_stim_frames;
mov_fake(cs_idx)=mov_fake(cs_idx+n_stim_frames);
mov_fake(i*n_frames)=mov_fake(i*n_frames-1);
end

[d1,d2,T] = mov.size; bound = 0;
options_r = NoRMCorreSetParms('d1',d1-bound,'d2',d2-bound,'bin_width',200,'max_shift',20,'iter',1,'correct_bidir',false);

% register using the high pass filtered data and apply shifts to original data
tic; [~,shifts1,template1] = normcorre_batch(single(mov_fake.toimg.data),options_r); toc % register filtered data
clear mov_fake
% exclude boundaries due to high pass filtering effects
tic; mov_mc = apply_shifts(single(mov.toimg.data),shifts1,options_r,0,0); toc % apply shifts to full dataset
end