% Analysis Clustering
%
clear
[fpath_perm] = uigetfile_n_dir();

%% load footprint
fp='/Volumes/ByungHun4TB/Cohen_lab_Data/20230206_BHLm017_Treadmill_Assembly/173944_BHLm017_Ane_Const_Bp05';
load(fullfile(fp,'Result.mat'))
%%


for i=1:length(fpath_perm)
load(fullfile(fpath_perm{i},'settings.mat'));
Sz = importdata([fpath_perm{i} '/experimental_parameters.txt']);
sz1=Sz.data(1); sz2=Sz.data(2);
mov_mc=double(readBinMov([fpath_perm{i} '/mc.bin'],sz2,sz1));
a=find(DAQ_waves.amplitude(1,[1:500])); frm_rate=1/((a(2)-a(1)-2)*1e-5);

if ~isempty(dmd_mask_sequence_rois)
[~, Result_perm{i}.Blue]=show_multiROI_waveform(DAQ_waves,size(dmd_mask_sequence_rois,2),frm_rate,false);
else
Result_perm{i}.Blue=DAQ_waves.amplitude(4,round([1:size(mov_mc,3)]*(1/frm_rate+0.00002)/1e-5));
end

Result_perm{i}.ref_im=mean(mov_mc,3);
corr_im=xcorr2(Result_perm{1}.ref_im,Result{1}.ref_im); [~, arg]=max(corr_im(:));
[s1 s2]=ind2sub(size(corr_im,1),arg);
s1=sz2-s1; s2=sz1-s2;
Result_perm{i}.c_ftprnt=zeros(sz2,sz1,size(Result{1}.c_ftprnt,3));
Result_perm{i}.c_ftprnt=Result{1}.c_ftprnt;
%Result_perm{i}.c_ftprnt(1:end-abs(s1),1:end-abs(s2),:)=...
%Result{1}.c_ftprnt((sz2-s1)/2-,ceil(sz1/2-(sz1+s2)/2):floor(sz1/2+(sz1+s2)/2),:);
Result_perm{i}.coord=get_coord(Result_perm{i}.c_ftprnt);

if ~isempty(dmd_mask_sequence_rois)
%dmd_mask_sequence_rois{1}=dmd_mask_sequence_rois{1};
%dmd_mask_sequence_rois{2}=dmd_mask_sequence_rois{1};

for j=1:size(dmd_mask_sequence_rois,2) %find cells in the each dmd mask
    
    roi=dmd_pat2roi(dmd_mask_sequence_rois{j},size(Result_perm{i}.ref_im));

Result_perm{i}.clist_p{j}=find(roi(sub2ind(size(Result_perm{i}.ref_im),round(Result_perm{i}.coord(:,2)),round(Result_perm{i}.coord(:,1)))));
end
end
Result_perm{i}.frm_rate=frm_rate;
Result_perm{i}.DMDPatt=dmd_mask_sequence_rois;
Result_perm{i}.traces=-(tovec(mov_mc)'*tovec(Result_perm{i}.c_ftprnt))';
Result_perm{i}.traces_bin=-(tovec(mov_mc)'*tovec(Result_perm{i}.c_ftprnt>0))';
Result_perm{i}.traces_hi=squeeze(Result_perm{i}.traces) - movmedian(squeeze(Result_perm{i}.traces),150,2);
Result_perm{i}.traces_bin_hi=squeeze(Result_perm{i}.traces_bin) - movmedian(squeeze(Result_perm{i}.traces_bin),200,2);
%traces_hi=squeeze(traces) - mean(squeeze(traces),2);
Result_perm{i}.noise=get_threshold(Result_perm{i}.traces,1);
%traces_hi=traces_hi./range(traces_hi,2);
Result_perm{i}.traces_hi=Result_perm{i}.traces_hi./Result_perm{i}.noise;
Result_perm{i}.traces_bin_hi=Result_perm{i}.traces_bin_hi./get_threshold(Result_perm{i}.traces_bin_hi,1);

tmp=squeeze(Result_perm{i}.traces) - movmedian(squeeze(Result_perm{i}.traces),20,2); tmp=tmp./get_threshold(tmp,1);
Result_perm{i}.spike=find_spike_bh(tmp,4,1.5);
tmp=squeeze(Result_perm{i}.traces) - movmedian(squeeze(Result_perm{i}.traces),100,2); tmp=tmp./get_threshold(tmp,1);
Result_perm{i}.spike=(Result_perm{i}.spike+find_spike_bh(tmp,4,1.5))>0;
end
%%
figure; i=8; %cmap=distinguishable_colors(size(Result{i}.DMDPatt,2));
cmap=[0.9 0.7 0.13;0 0.45 0.74;0.85 0.33 0.1];
scale=10; t=[1:size(Result_perm{i}.traces_hi,2)]/Result_perm{i}.frm_rate;
 tr=Result_perm{i}.traces-median(Result_perm{i}.traces,2); fprnt=Result_perm{i}.c_ftprnt;
 %tr=tr./range(tr,2);
 tr=tr./get_threshold(Result_perm{i}.traces,1);
%tr=Result{i}.traces_hi; fprnt=Result{i}.c_ftprnt;
%order=[1:size(tr,1)];
%[noi]=find(sum(Result_perm{i}.spike,2)>0)';
[noi]=[1:size(Result_perm{i}.spike,1)];
%[noi]=find(sum(Result_perm{i}.spike,2)>5)';
%[~, noi]=sort(sum(Result{i}.spike,2),'descend'); noi=noi';

%order=[47 setdiff([1:size(tr,1)],47)];
c=1;
vs=compete_ICA(c,:);
tiledlayout(10,4)
% 
ax1 = nexttile([2 2]);
colr = max(colormap(jet(size(fprnt,3))),0);
imshow2(squeeze(sum(fprnt.*reshape(colr,1,1,[],3),3)),[]); hold all;
text(Result_perm{i}.coord(:,1)',Result_perm{i}.coord(:,2)',num2str([1:size(fprnt,3)]'),'color','w')

ax4 = nexttile([2 2]);
imshow2(Result_perm{i}.ref_im,[]); hold all
line_color=zeros(size(noi,2),3); line_color2=[];
order=[];
% for j=2:size(Result_perm{i}.DMDPatt,2)
% plot(Result_perm{i}.DMDPatt{j}(:,1),Result_perm{i}.DMDPatt{j}(:,2),'color',cmap(j,:))
% order=[order; noi(find(ismember(noi,Result_perm{i}.clist_p{j})))'];
% line_color(find(ismember(noi,Result_perm{i}.clist_p{j})),:)=repmat(cmap(j,:),sum(ismember(Result_perm{i}.clist_p{j},noi)),1);
% line_color2=[line_color2; repmat(cmap(j,:),sum(ismember(Result_perm{i}.clist_p{j},noi)),1)];
% end
plot(Result_perm{i}.DMDPatt{2}(:,1),Result_perm{i}.DMDPatt{2}(:,2),'color',cmap(2,:))
plot(Result_perm{i}.DMDPatt{3}(:,1),Result_perm{i}.DMDPatt{3}(:,2),'color',cmap(3,:))
order=[order; find(CAs(:,vs(1)))];
order=[order; find(CAs(:,vs(2)))];
line_color2=[line_color2; repmat(cmap(2,:),sum(ismember(find(CAs(:,vs(1))),noi)),1)];
line_color2=[line_color2; repmat(cmap(3,:),sum(ismember(find(CAs(:,vs(2))),noi)),1)];    
order=[order; setdiff([noi],order)'];
line_color2(end+1:size(noi,2),:)=0;

linkaxes([ax1 ax4],'xy')

ax2 = nexttile([7 4]);
lines=plot(t,tr(order,:)'+[1:size(noi,2)]*scale);
%arrayfun(@(l,c) set(l,'Color',c{:}),lines,num2cell(line_color,2))
try arrayfun(@(l,c) set(l,'Color',c{:}),lines,num2cell(line_color2,2))
catch arrayfun(@(l,c) set(l,'Color',c{:}),lines,num2cell(jet(length(noi)),2))
end
axis tight
hold all
S=tr; S(~Result_perm{i}.spike)=NaN;
plot(t,S(order,:)'+[1:size(noi,2)]*scale,'r.')
set(gca,'ytick',[1:size(noi,2)]*scale,'yticklabel',order)
title(fpath_mac2window(fpath_perm{i}),'Interpreter', 'none')
hold all
[a b]=unique(bwlabel(Result_perm{i}.Blue(1,:)>0));
line([t(b(2:end))' t(b(2:end))']',repmat([0 size(noi,2)*scale],size(b,1)-1,1)','color','r')


ax3 = nexttile([1 4]);
lines2=plot(t,Result_perm{i}.Blue);
arrayfun(@(l,c) set(l,'Color',c{:}),lines2,num2cell(cmap([2 3],:),2))
axis tight
%show_multiROI_waveform(DAQ_waves,3,Result{i}.frm_rate,true);

linkaxes([ax2 ax3],'x')
saveas(gcf,[fpath_perm{i} '/voltage_trace_plot'])