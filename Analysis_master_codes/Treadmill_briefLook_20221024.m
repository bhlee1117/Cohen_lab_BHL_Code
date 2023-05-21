% Analysis on AAV expression sample and plot, in house YQ201
% 2022/10/13, Byung Hun Lee

%clear
%filename='/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/20220808_YQ201_IH_bistability_first_explore/231910_YQ201_IH_M2_expression_p3';
[fname,fpath] = uigetfile('*.*','MultiSelect','on'); if fname == 0, return;end
mov=vm(fpath,[112000:122000]);
%mov=double(readBinMov_times(filename,304,512,[1:10000]));
moviefixsc(mov)
load(fullfile(fpath,'settings.mat'))
Blue=DAQ_waves.amplitude(4,round([1:size(mov,3)]*1.27*1e-3/1e-5));
%%

mov_test=mov(:,:,150:250);
mov_test = single(mov_test)./single(max(mov_test.data(:)));
mov_test = movmean(mov_test,10,3);
mov_ref = squeeze(median(mov_test,3));
[mov_mc,xyField]=optical_flow_motion_correction_LBH(mov,mov_ref);

%% High-pass filter
mov_res=SeeResiduals(mov_mc,squeeze(mean(movmean(mov_mc,200,3),[1 2])));
im_G=imgaussfilt(mean(mov_mc,3)-medfilt2(mean(mov_mc,3),[8 8]),2);
[centers radii]=Cell_segment_circle_2x(im_G);
centers=cell_detection_manual(mean(mov_mc,3),centers);
%%
mcTrace = squeeze(mean(xyField,[1 2]));
mov_res = SeeResiduals(mov_mc,mcTrace);
mov_res = SeeResiduals(mov_res,mcTrace.^2);
mov_res = SeeResiduals(mov_res,mcTrace(:,1).*mcTrace(:,2));
c_ftprnt=mask_footprint(centers,mov_res,[],5);
show_footprnt(c_ftprnt,mov_mc)
coord=get_coord(c_ftprnt);

traces=-(tovec(mov_mc)'*tovec(c_ftprnt>0))';
traces_hi=squeeze(traces) - movmean(squeeze(traces),200,2); 
traces_hi=traces_hi./range(traces_hi,2);
%%
figure; scale=0.7;
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
