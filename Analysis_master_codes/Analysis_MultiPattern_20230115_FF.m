% Analysis on AAV expression sample and plot, in house YQ201
% 2022/10/13, Byung Hun Lee

%clear
addpath('Y:\Lab\Labmembers\Byung Hun Lee\Code\Codes(~2021)\uigetfile_n_dir')
[fpath] = uigetfile_n_dir();

%%
for i=1:length(fpath)
load(fullfile(fpath{i},'settings.mat'));
load(fullfile(fpath{i},'mcTrace.mat'));
Sz = importdata([fpath{i} '/experimental_parameters.txt']);
sz1=Sz.data(1); sz2=Sz.data(2);
mov_mc=double(readBinMov_times([fpath{i} '/mc.bin'],sz2,sz1,[1:500]));
im_G=imgaussfilt(mean(mov_mc,3),2);
[centers radii]=Cell_segment_circle_2x(im_G,0.95);
Result{i}.centers=cell_detection_manual(mean(mov_mc,3),centers);
end
%%

for i=1:length(fpath)
load(fullfile(fpath{i},'settings.mat'));
load(fullfile(fpath{i},'mcTrace.mat'));
Sz = importdata([fpath{i} '/experimental_parameters.txt']);
sz1=Sz.data(1); sz2=Sz.data(2);
mov_mc=double(readBinMov([fpath{i} '/mc.bin'],sz2,sz1));
a=find(DAQ_waves.amplitude(1,[1:500])); frm_rate=1/((a(2)-a(1)-2)*1e-5);

if ~isempty(dmd_mask_sequence_rois)
[~, Result{i}.Blue]=show_multiROI_waveform(DAQ_waves,size(dmd_mask_sequence_rois,2),frm_rate,false);
else
Result{i}.Blue=DAQ_waves.amplitude(4,round([1:size(mov_mc,3)]*(1/frm_rate+0.00002)/1e-5));
end
%Result{i}.Blue=DAQ_waves.amplitude(4,round([1:size(mov_mc,3)]*(1/frm_rate)/1e-5));
%mov_res=SeeResiduals(mov_mc,squeeze(mean(movmean(mov_mc,200,3),[1 2])));

mov_res= mov_mc-mean(mov_mc,3);
mov_res = SeeResiduals(mov_res,mcTrace);
mov_res = SeeResiduals(mov_res,mcTrace.^2);
mov_res = SeeResiduals(mov_res,mcTrace(:,1).*mcTrace(:,2));

Result{i}.c_ftprnt=mask_footprint(Result{i}.centers,movmean(mov_res(:,:,1000:end),10,3),[],6);
%Result{i}.c_ftprnt=Result{3}.c_ftprnt;
Result{i}.ref_im=mean(mov_mc,3);
%show_footprnt(Result{i}.c_ftprnt,mov_mc)
Result{i}.coord=get_coord(Result{i}.c_ftprnt);

if ~isempty(dmd_mask_sequence_rois)
dmd_mask_sequence_rois{1}=dmd_mask_sequence_rois{1};
dmd_mask_sequence_rois{2}=dmd_mask_sequence_rois{1};

for j=1:size(dmd_mask_sequence_rois,2)
roi=dmd_pat2roi(dmd_mask_sequence_rois{j},size(Result{i}.ref_im));
Result{i}.clist_p{j}=find(roi(sub2ind(size(Result{i}.ref_im),round(Result{i}.coord(:,2)),round(Result{i}.coord(:,1)))));
end
end
Result{i}.frm_rate=frm_rate;
Result{i}.DMDPatt=dmd_mask_sequence_rois;
Result{i}.traces=-(tovec(mov_mc)'*tovec(Result{i}.c_ftprnt))';
Result{i}.traces_bin=-(tovec(mov_mc)'*tovec(Result{i}.c_ftprnt>0))';
Result{i}.traces_hi=squeeze(Result{i}.traces) - movmedian(squeeze(Result{i}.traces),150,2);
Result{i}.traces_bin_hi=squeeze(Result{i}.traces_bin) - movmedian(squeeze(Result{i}.traces_bin),200,2);
%traces_hi=squeeze(traces) - mean(squeeze(traces),2);
Result{i}.noise=get_threshold(Result{i}.traces,1);
%traces_hi=traces_hi./range(traces_hi,2);
Result{i}.traces_hi=Result{i}.traces_hi./Result{i}.noise;
Result{i}.traces_bin_hi=Result{i}.traces_bin_hi./get_threshold(Result{i}.traces_bin_hi,1);
tmp=squeeze(Result{i}.traces) - movmedian(squeeze(Result{i}.traces),70,2); tmp=tmp./get_threshold(tmp,1);
Result{i}.spike=find_spike_bh(tmp,3);
end

%%
figure; i=1; %cmap=distinguishable_colors(size(Result{i}.DMDPatt,2));
cmap=[0.9 0.7 0.13;0 0.45 0.74;0.85 0.33 0.1];
scale=10; t=[1:size(Result{i}.traces_hi,2)]/Result{i}.frm_rate;
 tr=Result{i}.traces-median(Result{i}.traces,2); fprnt=Result{i}.c_ftprnt;
 %tr=tr./range(tr,2);
 tr=tr./get_threshold(Result{i}.traces,1);
%tr=Result{i}.traces_hi; fprnt=Result{i}.c_ftprnt;
%order=[1:size(tr,1)];
[noi]=find(sum(Result{i}.spike,2)>5)';
%[~, noi]=sort(sum(Result{i}.spike,2),'descend'); noi=noi';

%order=[47 setdiff([1:size(tr,1)],47)];

%tiledlayout(10,4)
ax1 = subplot(8,2,1);%nexttile([2 2]);
% 
%ax1 = nexttile([2 2]);
colr = max(colormap(jet(size(fprnt,3))),0);
imshow2(squeeze(sum(fprnt.*reshape(colr,1,1,[],3),3)),[]); hold all;
text(Result{i}.coord(:,1)',Result{i}.coord(:,2)',num2str([1:size(fprnt,3)]'),'color','w')

%ax4 = nexttile([2 2]);
ax4 = subplot(8,2,2);
imshow2(Result{i}.ref_im,[]); hold all
line_color=zeros(size(noi,2),3); line_color2=[];
order=[];
for j=2:size(Result{i}.DMDPatt,2)
plot(Result{i}.DMDPatt{j}(:,1),Result{i}.DMDPatt{j}(:,2),'color',cmap(j,:))
order=[order; noi(find(ismember(noi,Result{i}.clist_p{j})))'];
line_color(find(ismember(noi,Result{i}.clist_p{j})),:)=repmat(cmap(j,:),sum(ismember(Result{i}.clist_p{j},noi)),1);
line_color2=[line_color2; repmat(cmap(j,:),sum(ismember(Result{i}.clist_p{j},noi)),1)];
end
order=[order; setdiff([noi],order)'];
line_color2(end+1:size(noi,2),:)=0;

linkaxes([ax1 ax4],'xy')

%ax2 = nexttile([7 4]);
ax2 = subplot(8,2,3:14);
lines=plot(t,tr(order,:)'+[1:size(noi,2)]*scale);
%arrayfun(@(l,c) set(l,'Color',c{:}),lines,num2cell(line_color,2))
try arrayfun(@(l,c) set(l,'Color',c{:}),lines,num2cell(line_color2,2))
catch arrayfun(@(l,c) set(l,'Color',c{:}),lines,num2cell(jet(length(noi)),2))
end
axis tight
hold all
S=tr; S(~Result{i}.spike)=NaN;
plot(t,S(order,:)'+[1:size(noi,2)]*scale,'r.')
set(gca,'ytick',[1:size(noi,2)]*scale,'yticklabel',order)
%title(fpath_mac2window(fpath{i}),'Interpreter', 'none')
title((fpath{i}),'Interpreter', 'none')
hold all
[a b]=unique(bwlabel(Result{i}.Blue(1,:)>0));
line([t(b(2:end))' t(b(2:end))']',repmat([0 size(noi,2)*scale],size(b,1)-1,1)','color','r')


%ax3 = nexttile([1 4]);
ax3 = subplot(8,2,15:16);
plot(t,Result{i}.Blue)
axis tight
%show_multiROI_waveform(DAQ_waves,3,Result{i}.frm_rate,true);

linkaxes([ax2 ax3],'x')
saveas(gcf,[fpath{i} 'voltage_trace_plot'])
save([fpath{1} 'result.mat'],'Result','fpath')

%%

%^writeMov('DMD2_stim',-movmean(mov_res(:,:,4000:end),30,3),[],200,1,[-5 70],[],Blue(4000:end))
%writeMov('DMD3_stim',-movmean(mov_res(:,:,8701:17400),30,3),[],200,1.25,[-5 70],[],Blue(8701:17400))
%writeMov('DMD23_stim',-movmean(mov_res(:,:,17401:end),30,3),[],200,1.25,[-5 70],[],Blue(17401:end))


