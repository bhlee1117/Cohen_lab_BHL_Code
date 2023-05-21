% Analysis on AAV expression sample and plot, in house YQ201
% 2022/10/13, Byung Hun Lee

clear
%filename='/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/20220808_YQ201_IH_bistability_first_explore/231910_YQ201_IH_M2_expression_p3';
[fpath] = uigetfile_n_dir; %if fname == 0, return;end
% mov=vm(fpath); load(fullfile(fpath,'settings.mat'));
% %mov=double(readBinMov_times(filename,304,512,[1:10000]));
% moviefixsc(mov)
% hold all
% plot(dmd_current_roi{1}(:,1),dmd_current_roi{1}(:,2))
% Blue=DAQ_waves.amplitude(4,round([1:size(mov,3)]*1.27*1e-3/1e-5));
%%
mov=vm(fpath{1});
mov_test=mov(:,:,150:250);
mov_test = single(mov_test)./single(max(mov_test.data(:)));
mov_test = movmean(mov_test,10,3);
mov_ref = squeeze(median(mov_test,3)); mov=double(mov);

im_G=imgaussfilt(mean(mov,3)-medfilt2(mean(mov,3),[8 8]),2);
%im_G=imgaussfilt(mean(mov_mc,3)-medfilt2(mean(mov_mc,3),[16 16]),2);
[centers radii]=Cell_segment_circle_2x(im_G);
centers=cell_detection_manual(mean(mov,3),centers);
%% High-pass filter

for i=1:length(fpath)
mov=vm(fpath{i});
[mov_mc,xyField]=optical_flow_motion_correction_LBH(mov,mov_ref,'optic_flow');

mov_res=SeeResiduals(mov_mc,squeeze(mean(movmean(mov_mc,200,3),[1 2])));
mcTrace = squeeze(mean(xyField,[1 2]));
mov_res = SeeResiduals(mov_res,mcTrace);
mov_res = SeeResiduals(mov_res,mcTrace.^2);
mov_res = SeeResiduals(mov_res,mcTrace(:,1).*mcTrace(:,2));

c_ftprnt=mask_footprint(centers,mov_res,[],6);
show_footprnt(c_ftprnt,mov_mc)
coord=get_coord(c_ftprnt);
traces{i}=-(tovec(mov_mc)'*tovec(c_ftprnt>0))';
traces_hi{i}=squeeze(traces{i}) - movmean(squeeze(traces{i}),400,2); 
%traces_hi=squeeze(traces) - mean(squeeze(traces),2); 
traces_hi{i}=traces_hi{i}./range(traces_hi{i},2);
end

%%
figure; N=31;
tiledlayout(9,1);hold all
axs=[]; inten=[10:-1:2]; g=1;
for i=[1 4 2 3 5 6 7 8 9]
ax=nexttile;
axs=[axs ax];
plot(traces{g}(N,500:end))
title(['Orange: ' num2str(inten(g)) 'V'])
g=g+1;
end

linkaxes(axs,'x')

