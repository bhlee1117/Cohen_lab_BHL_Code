% Analysis on AAV expression sample and plot, in house YQ201
% 2022/11/09, Byung Hun Lee

%clear
[fpath] = uigetfile_n_dir;
%%
for i=11:length(fpath)
load([fpath{i} '/output_data.mat'])
sz=double(Device_Data{1, 4}.ROI([2 4]));
mov=double(readBinMov([fpath{i} '/frames1.bin'],sz(2),sz(1)));
%mov=double(readBinMov_times([fpath{i} '/frames1.bin'],sz(2),sz(1),[1:30000]));
CamTrigger=find(Device_Data{1, 2}.Counter_Inputs.data(2:end)-Device_Data{1, 2}.Counter_Inputs.data(1:end-1));
mov_test=mov(:,:,150:250);
try mov_test = single(mov_test)./single(max(mov_test.data(:)));
catch disp('change to vm')

mov_test=vm(mov_test); mov_test = single(mov_test)./single(max(mov_test.data(:))); end
mov_test = movmean(mov_test,10,3);
mov_ref = squeeze(median(mov_test,3));
[mov_mc,xyField]=optical_flow_motion_correction_LBH(vm(mov(:,:,1:length(CamTrigger))),mov_ref,'normcorre');
mov_mc=vm(mov_mc);
mov_mc.transpose.savebin([fpath{i} '/mc.bin'])
%mcTrace = squeeze(mean(xyField,[1 2]));

shifts_r = squeeze(cat(3,xyField(:).shifts));
shifts_nr = cat(ndims(xyField(1).shifts)+1,xyField(:).shifts);
shifts_nr = reshape(shifts_nr,[],ndims(mov_mc)-1,size(mov_mc,3));
%shifts_x = squeeze(shifts_nr(:,2,:))';
%shifts_y = squeeze(shifts_nr(:,1,:))';
mcTrace=squeeze(shifts_nr)';

save([fpath{i} '/mcTrace.mat'],'mcTrace')
end

%moviefixsc(mov)
%hold all
%plot(dmd_current_roi{1}(:,1),dmd_current_roi{1}(:,2))
%Blue=DAQ_waves.amplitude(4,round([1:size(mov,3)]*1.27*1e-3/1e-5));
%%
i=1; DAQ_rate=0.00001;
disp([fpath{i}])
load([fpath{i} '/output_data.mat'])
Blue=Device_Data{1, 2}.buffered_tasks(1, 1).channels(1, 2).data;
a=find(Device_Data{1, 2}.Counter_Inputs.data(2:end)-Device_Data{1, 2}.Counter_Inputs.data(1:end-1));
frm_rate=(a(2)-a(1))*DAQ_rate;
sz=Device_Data{1, 4}.ROI([2 4]);
mov_mc=double(readBinMov([fpath{i} '/mc.bin'],sz(2),sz(1)));
t=[1:size(mov_mc,3)-1]*frm_rate; t(t>length(Blue)*DAQ_rate)=[];


[roi int]=clicky(mov_mc); int=int-median(int,1);
close all
figure; ax1=[];
tiledlayout(5,1)
nexttile([1 1])
imshow2(mean(mov_mc,3),[]); hold all
for r=1:length(roi); plot(roi{r}(:,1),roi{r}(:,2),'linewidth',2); hold all; end
ax1=[ax1 nexttile([3 1])]; 
plot(t,-int(1:length(t),:)+[1:length(roi)]*10); set(gca,'FontSize',18)
title(fpath_mac2window([fpath{i}]))
%ylim([min(-int(10:length(t))) max(-int(10:length(t)))])
ylim([0 (length(roi)+1)*10])
ax1=[ax1 nexttile([1 1])];
plot([1:length(Blue)]*0.00001,Blue,'linewidth',2); set(gca,'FontSize',18);
linkaxes(ax1,'x')
axis tight
set(gcf,'color','white')
saveas(gcf,[fpath{i} '/clicky_plot'])

%%



for i=1%:length(fpath)
load([fpath{i} '/output_data.mat'])
sz=Device_Data{1, 4}.ROI([2 4]);

mov_mc=double(readBinMov([fpath{i} '/mc.bin'],sz(2),sz(1)));
mov_res= mov_mc-mean(mov_mc,3);
im_G=imgaussfilt(mean(mov_mc,3),3);
[centers radii]=Cell_segment_circle_10x(im_G,0.7);
Result{i}.centers=cell_detection_manual(mean(mov_mc,3),centers,[]);
Result{i}.c_ftprnt=mask_footprint(Result{i}.centers,movmean(mov_res(:,:,1000:end),10,3),[],20);
Result{i}.traces=-(tovec(mov_mc)'*tovec(Result{i}.c_ftprnt))';
Result{i}.spike=[];
tmp=squeeze(Result{i}.traces) - movmedian(squeeze(Result{i}.traces),4,2); tmp=tmp./get_moving_threshold(tmp,1,2000);
Result{i}.spike=find_spike_bh(tmp,4,3);
tmp=squeeze(Result{i}.traces) - movmedian(squeeze(Result{i}.traces),100,2); tmp=tmp./get_moving_threshold(tmp,1,2000);
Result{i}.spike=(Result{i}.spike+find_spike_bh(tmp,3,1.5))>0;
end



%%
for i=1:length(fpath)
    show_traces_spikes(Result{i}.traces,Result{i}.spike,Result{i}.spike)
end
%%
mov_mc=double(mov_mc);
[roi sig]=clicky(mov_mc);
mov_res=mov_mc-movmean(mov_mc,300,3);
%%
%sig=-sig';
tmp=(int) - movmedian(int,70,2); tmp=tmp./get_threshold(tmp,1);
spike=find_spike_bh(tmp,3.5);
plot(int)
hold all
S=int; S(~spike)=NaN;
plot(S,'r.')

%%
ss=find(spike);
for s=1:length(ss)
    try
    sta_mov(:,:,:,s)=mov_res(:,:,ss(s)-50:ss(s)+100);
    end
end
STA=mean(sta_mov,4);
%%
for i=1:length(fpath)
    load([fpath{i} '/output_data.mat'])
    sz=Device_Data{1, 4}.ROI([2 4]);
    mov_mc=double(readBinMov([fpath{i} '/mc.bin'],sz1,sz2));
end
%% High-pass filter
mov_mc=double(mov);
im_G=imgaussfilt(mean(mov_mc,3),2);
mov_res=SeeResiduals(mov_mc,squeeze(mean(movmean(mov_mc,200,3),[1 2])));
%im_G=imgaussfilt(mean(mov_mc,3)-medfilt2(mean(mov_mc,3),[8 8]),2);
%im_G=imgaussfilt(mean(mov_mc,3)-medfilt2(mean(mov_mc,3),[16 16]),2);
[centers radii]=Cell_segment_circle_10x(im_G);
centers=cell_detection_manual(mean(mov_mc,3),centers);
%%
mcTrace = squeeze(mean(xyField,[1 2]));
mov_res = SeeResiduals(mov_res,mcTrace);
mov_res = SeeResiduals(mov_res,mcTrace.^2);
mov_res = SeeResiduals(mov_res,mcTrace(:,1).*mcTrace(:,2));
c_ftprnt=mask_footprint(centers,mov_res,[],10);
show_footprnt(c_ftprnt,mov_mc)
coord=get_coord(c_ftprnt);

traces=-(tovec(mov_mc)'*tovec(c_ftprnt>0))';
traces_hi=squeeze(traces) - movmean(squeeze(traces),400,2); 
%traces_hi=squeeze(traces) - mean(squeeze(traces),2); 
traces_hi=traces_hi./range(traces_hi,2);
%traces_hi=-traces_hi./mean(traces_hi,2);
%%
figure; scale=0.6;
tiledlayout(10,4)
ax1 = nexttile([2 2]);
colr = max(colormap(jet(size(c_ftprnt,3))),0);
imshow2(squeeze(sum(c_ftprnt.*reshape(colr,1,1,[],3),3)),[]); hold all;
text(coord(:,1)',coord(:,2)',num2str([1:size(c_ftprnt,3)]'),'color','w')

ax4 = nexttile([2 2]);
imshow2(mean(mov_mc,3),[]); colormap('gray')

linkaxes([ax1 ax4],'xy')

ax2 = nexttile([7 4]);
lines=plot(traces_hi'+[1:size(traces_hi,1)]*scale);
arrayfun(@(l,c) set(l,'Color',c{:}),lines,num2cell(jet(size(c_ftprnt,3)),2))
axis tight
set(gca,'ytick',[1:size(traces_hi,1)]*scale,'yticklabel',[1:size(traces_hi,1)])
ax3 = nexttile([1 4]);
plot(Blue)

linkaxes([ax2 ax3],'x')
saveas(gcf,[fpath 'voltage_trace_plot'])
save([fpath 'result.mat'],'c_ftprnt','traces')
