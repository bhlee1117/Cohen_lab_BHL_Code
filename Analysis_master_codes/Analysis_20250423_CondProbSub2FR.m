run_line('/Users/bhlee1117/Documents/GitHub/Cohen_lab_BHL_Code/Analysis_master_codes/Analysis_20250331_FigureMakingSeeSaw.m',[1 2])

%% Spike probability
for f=foi

nROI=size(NormalizedTrace_dirt{f},1);
nTime=size(Subthreshold{f},2);
roisD_order_ind=cellfun(@find,roisD_order{f},'UniformOutput',false);
labelvec{f}=NaN(1,size(NormalizedTrace_dirt{f},1));
for dClass=1:5
    labelvec{f}(cell2mat(roisD_order_ind(dClass,:)'))=dClass;
end

perispike_time=[-3:5];
perispike_frame=unique([tovec(find(double(allSpikeMat{f}(1,:)==1))'+perispike_time); find(CStrace{f})']);
perispike_frame(perispike_frame<=0 | perispike_frame>nTime)=[];
nonvalid_frame=find(sum(isnan(Subthreshold{f}),1)>0);

Blue_on_frame=find(imdilate(BlueStim{f}>0, [ones(1, 1), 1, ones(1, 200)]));
Badframe=unique([Blue_on_frame nonvalid_frame]);
Goodframe=setdiff([1:nTime],Badframe);

SpikeHeight=get_STA(NormalizedTrace_dirt{f},allSpikeMat{f}(1,:),0,0);
SpikeHeight=mean(SpikeHeight(labelvec{f}==2));

normSub=Subthreshold{f}(:,Goodframe)/SpikeHeight;
[~, normSub_spike]=get_STA(normSub,allSpikeMat{f}(1,:).*(BlueStim{f}==0),1,3);
[~, normSub_SS]=get_STA(normSub,allSpikeClassVecMat{f}(1,:).*(BlueStim{f}==0),1,3);
[~, normSub_CS]=get_STA(normSub,allSpikeClassVecMat{f}(2,:).*(BlueStim{f}==0),1,3);
normSub_spike=max(normSub_spike,[],3);

subthreshold_bin_edge=[-2:0.05:2];
Sub_count=[]; Sub_Spike_count=[]; Sub_SS_count=[]; Sub_CS_count=[];
for r=1:nROI
Sub_count(r,:)=histcounts(normSub(r,:),subthreshold_bin_edge);
Sub_Spike_count(r,:)=histcounts(normSub_spike(r,:),subthreshold_bin_edge);
Sub_SS_count(r,:)=histcounts(normSub_SS(r,:),subthreshold_bin_edge);
Sub_CS_count(r,:)=histcounts(normSub_CS(r,:),subthreshold_bin_edge);
end
Sub_Spike_count=movmean(Sub_Spike_count,5,2);
Sub_SS_count=movmean(Sub_SS_count,5,2);
Sub_CS_count=movmean(Sub_CS_count,5,2);

figure(f); clf; tiledlayout(3,3); ax1=[];
Sub_binEdge=subthreshold_bin_edge(1:end-1)+diff(subthreshold_bin_edge(1:2));
nexttile([1 1]);
l=plot(Sub_binEdge,Sub_count');
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(turbo(nROI),2));
xlabel('Subthreshold (Spike height)')
ylabel('Counts')
nexttile([1 1]);
l=plot(Sub_binEdge,Sub_Spike_count');
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(turbo(nROI),2));
xlabel('Subthreshold (Spike height)')
ylabel('Counts')
title('Subthreshold at spike')
nexttile([1 1]);
l=plot(Sub_binEdge,Sub_CS_count');
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(turbo(nROI),2));
xlabel('Subthreshold (Spike height)')
ylabel('Counts')
title('Subthreshold at SS')
nexttile([1 1]);
l=plot(Sub_binEdge,Sub_SS_count');
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(turbo(nROI),2));
xlabel('Subthreshold (Spike height)')
ylabel('Counts')
title('Subthreshold at 1st CS')

ax1=[ax1 nexttile([1 1])];
l=plot(Sub_binEdge,(Sub_Spike_count./Sub_count)');
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(turbo(nROI),2));
xlabel('Subthreshold (Spike height)')
ylabel('P(AP|subthreshold)')
xlim([-0.4 0.6])

ax1=[ax1 nexttile([1 1])];
SS_prob=Sub_SS_count./Sub_count;
l=plot(Sub_binEdge,(SS_prob)');
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(turbo(nROI),2));
xlabel('Subthreshold (Spike height)')
ylabel('P(SS|subthreshold)')
title('SS')
xlim([-0.4 0.6])

ax1=[ax1 nexttile([1 1])];
CS_prob=Sub_CS_count./Sub_count;
l=plot(Sub_binEdge,(CS_prob)');
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(turbo(nROI),2));
xlabel('Subthreshold (Spike height)')
ylabel('P(1st CS| subthreshold)')
title('CS')
xlim([-0.4 0.6])

show_lab=[1 2 3 4 5];
cmap_prob=hsv(5); ax2=[];
ax2=[ax2 nexttile([1 1])];
for l=show_lab
plot(Sub_binEdge,mean(SS_prob(labelvec{f}==l,:),1,'omitnan'),'color',cmap_prob(l,:),'LineWidth',2); hold all
xlabel('Subthreshold (Spike height)')
ylabel('P(AP| subthreshold)')
xlim([-0.4 0.6])
title('SS')
end
legend(label_str);

ax2=[ax2 nexttile([1 1])];
for l=show_lab
plot(Sub_binEdge,mean(CS_prob(labelvec{f}==l,:),1,'omitnan'),'color',cmap_prob(l,:),'LineWidth',2); hold all
xlabel('Subthreshold (Spike height)')
ylabel('P(AP| subthreshold)')
xlim([-0.4 0.6])
title('CS')
end
legend(label_str);
linkaxes(ax1,'x')

end
