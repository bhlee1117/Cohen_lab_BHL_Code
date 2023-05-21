% Analysis on AAV expression sample and plot, in house YQ201
% 2022/10/13, Byung Hun Lee

clear
%filename='/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/20220808_YQ201_IH_bistability_first_explore/231910_YQ201_IH_M2_expression_p3';
[fname,fpath] = uigetfile('*.*','MultiSelect','on'); if fname == 0, return;end
%mov=vm(fpath);
fpath={fpath};
load(fullfile(fpath{1},'settings.mat'));
load(fullfile(fpath{1},'mcTrace.mat'));
Sz = importdata([fpath{1} '/experimental_parameters.txt']);
sz1=Sz.data(1); sz2=Sz.data(2);
mov_mc=double(readBinMov([fpath{1} '/mc.bin'],sz2,sz1));
moviefixsc(mov_mc)
hold all
plot(dmd_current_roi{1}(:,1),dmd_current_roi{1}(:,2))
%Blue=DAQ_waves.amplitude(4,round([1:size(mov_mc,3)-1]*1.27*1e-3/1e-5)+110);
Blue=DAQ_waves.amplitude(4,round([1:size(mov_mc,3)]*1.02*1e-3/1e-5));
%%

mov_test=mov(:,:,150:250);
try mov_test = single(mov_test)./single(max(mov_test.data(:)));
catch disp('change to vm')

    mov_test=vm(mov_test); mov_test = single(mov_test)./single(max(mov_test.data(:))); end
mov_test = movmean(mov_test,10,3);
mov_ref = squeeze(median(mov_test,3));
[mov_mc,xyField]=optical_flow_motion_correction_LBH((mov),mov_ref,'optic_flow');
%method: 'optic_flow' or 'normcorre'

%% High-pass filter

im_G=imgaussfilt(mean(mov_mc,3),2);
%im_G=imgaussfilt(mean(mov_mc,3)-medfilt2(mean(mov_mc,3),[8 8]),2);
%im_G=imgaussfilt(mean(mov_mc,3)-medfilt2(mean(mov_mc,3),[16 16]),2);
[centers radii]=Cell_segment_circle_2x(im_G);
centers=cell_detection_manual(mean(mov_mc,3),centers);
%%
%mcTrace = squeeze(mean(xyField,[1 2]));
mov_res=SeeResiduals(mov_mc,squeeze(mean(movmean(mov_mc,200,3),[1 2])));
mov_res = SeeResiduals(mov_res,mcTrace);
mov_res = SeeResiduals(mov_res,mcTrace.^2);
mov_res = SeeResiduals(mov_res,mcTrace(:,1).*mcTrace(:,2));
c_ftprnt=mask_footprint(centers,mov_res,[],6);
show_footprnt(c_ftprnt,mov_mc)
coord=get_coord(c_ftprnt);

traces=-(tovec(mov_res)'*tovec(c_ftprnt))';
traces_hi=squeeze(traces) - movmean(squeeze(traces),200,2);
%traces_hi=squeeze(traces) - mean(squeeze(traces),2);
noise=get_threshold(traces,1);
%traces_hi=traces_hi./range(traces_hi,2);
traces_hi=traces_hi./noise;
%%
figure; scale=10; t=[1:size(traces_hi,2)]*1*1e-3;

tiledlayout(10,4)

ax1 = nexttile([2 2]);
colr = max(colormap(jet(size(c_ftprnt,3))),0);
imshow2(squeeze(sum(c_ftprnt.*reshape(colr,1,1,[],3),3)),[]); hold all;
text(coord(:,1)',coord(:,2)',num2str([1:size(c_ftprnt,3)]'),'color','w')

ax4 = nexttile([2 2]);
imshow2(mean(mov_mc,3),[]); colormap('gray')

linkaxes([ax1 ax4],'xy')

ax2 = nexttile([7 4]);
order=[47 setdiff([1:size(traces_hi,1)],47)];
lines=plot(t,traces_hi(order,:)'+[1:size(traces_hi,1)]*scale);
arrayfun(@(l,c) set(l,'Color',c{:}),lines,num2cell(jet(size(c_ftprnt,3)),2))
axis tight
set(gca,'ytick',[1:size(traces_hi,1)]*scale,'yticklabel',order)
ax3 = nexttile([1 4]);
show_multiROI_waveform(DAQ_waves,3);


linkaxes([ax2 ax3],'x')
saveas(gcf,[fpath{1} 'voltage_trace_plot'])
save([fpath{1} 'result.mat'],'c_ftprnt','traces')

%%

writeMov('DMD2_stim',-movmean(mov_res(:,:,4000:end),30,3),[],200,1,[-5 70],[],Blue(4000:end))
%writeMov('DMD3_stim',-movmean(mov_res(:,:,8701:17400),30,3),[],200,1.25,[-5 70],[],Blue(8701:17400))
%writeMov('DMD23_stim',-movmean(mov_res(:,:,17401:end),30,3),[],200,1.25,[-5 70],[],Blue(17401:end))


