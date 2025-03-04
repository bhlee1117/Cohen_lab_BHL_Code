clear
clc;
cd '/Volumes/BHL18TB_D2/Arranged_Data/Prism_OptopatchResult';
[~, ~, raw] = xlsread(['/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/' ...
    'Prism_OptopatchData_Arrangement.xlsx'], 'Sheet1', 'B5:K154');

save_to='/Volumes/BHL18TB_D2/Arranged_Data/Prism_OptopatchResult';
fpath=raw(:,1);
Mouse=cell2mat(raw(:,2));
NeuronInd=cell2mat(raw(:,5));
CamType=raw(:,3);
StructureData=raw(:,10);

%% Stimulation power 0.15V 0.3V 0.5V
g=1; ref_spike_ROI=5;
clear NormTrace subth_trace spike_erode_trace Blue skewBlue
for i=[148 147 146]
    load([fpath{i} '/Result.mat'])
tr_norm= Result.traces_bvMask-movprc(Result.traces_bvMask,100,30,2);
tr_norm= tr_norm./get_threshold(tr_norm,1);
spike{g}= find_spike_bh(tr_norm,5,3);
F0=tovec(Result.ref_im)'*tovec(Result.ftprnt.*double(max(Result.bvMask,[],3)==0));
NormTrace{g}=Result.traces_bvMask./F0';
NormTrace{g}=NormTrace{g}-prctile(NormTrace{g},30,2);
subth_trace{g}=get_subthreshold(NormTrace{g},spike{g}(ref_spike_ROI,:),5,15);
spike_time=tovec(find(spike{g}(ref_spike_ROI,:))'+[-5:50]);
spike_time(spike_time<1 | spike_time>size(Result.traces,2))=[];
spike_time=unique(spike_time);
spike_erode_trace{g}=zeros(1,size(Result.traces,2));
spike_erode_trace{g}(spike_time)=1;
Blue{g}=Result.Blue;
g=g+1;
end
[~, BlueOffTime]=get_blueoffTrace(zeros(1,length(Result.Blue)),Result.Blue,100);
BlueTime= {BlueOffTime,Result.Blue>0}; % Blue off: b=1; Blue on: b=2; Excitation: inh=1; inhibition: inh=2;

figure(20); clf;
tiledlayout(4,2); ax1=[]; ax2=[]; 
ax1=[ax1 nexttile([1 2])];
for j=1:3
plot(NormTrace{j}(ref_spike_ROI,:)); hold all
end
xlabel('Time (ms)')
ylabel('Voltage')
legend({'0.15V','0.3V','0.5V'})
ax1=[ax1 nexttile([1 2])];
for j=1:3
plot(subth_trace{j}(ref_spike_ROI,:)); hold all
end
xlabel('Time (ms)')
ylabel('Subthreshold')
legend({'0.15V','0.3V','0.5V'})
ax1=[ax1 nexttile([1 2])];
for j=1:3
plot(Blue{j}); hold all
end
xlabel('Time (ms)')
ylabel('Blue AOTF voltage')
legend({'0.15V','0.3V','0.5V'})
linkaxes(ax1,'x')

ax2=[ax2 nexttile([1 1])];
binEdges=[-0.05:0.001:0.05];
for j=1:3
    histogram(subth_trace{j}(ref_spike_ROI,BlueTime{1} & spike_erode_trace{j}==0),binEdges,'Normalization','probability'); hold all
end
xlabel('Subthreshold')
ylabel('PDF')
legend({'0.15V','0.3V','0.5V'})
title('Blue Off')
ax2=[ax2 nexttile([1 1])];
for j=1:3
    skewBlue{j}(1,:)=skewness(subth_trace{j}(:,BlueTime{1} & spike_erode_trace{j}==0),1,2);
    skewBlue{j}(2,:)=skewness(subth_trace{j}(:,BlueTime{2} & spike_erode_trace{j}==0),1,2);
    histogram(subth_trace{j}(ref_spike_ROI,BlueTime{2} & spike_erode_trace{j}==0),binEdges,'Normalization','probability'); hold all
end
xlabel('Subthreshold')
ylabel('PDF')
legend({'0.15V','0.3V','0.5V'})
title('Blue On')
linkaxes(ax2,'x')
%% Low stim
for i=146
    load([fpath{i} '/Result.mat'])
    load(fullfile(fpath{i},"output_data.mat"))
    sz=double(Device_Data{1, 3}.ROI([2 4]));
    mov_mc=double(readBinMov([fpath{i} '/mc_ShutterReg01.bin'],sz(2),sz(1)));
end
[blueDMDimg bluePatt]=get_blueDMDPatt(Device_Data);

bound=5;
mean_F=squeeze(mean(mov_mc(bound:end-bound,bound:end-bound,:),[1 2]));
[~, blueOff]=get_blueoffTrace(mean_F,[Result.Blue],70);
[y_fit]=expfitDM_2(find(blueOff)',mean_F(find(blueOff)),[1:size(mov_mc,3)]',1000);
bkg(1,:)=y_fit;
mov_res=mov_mc-mean(mov_mc,3);
mov_res = SeeResiduals(mov_res,Result.mc);
mov_res = SeeResiduals(mov_res,Result.mc.^2);
mov_res = SeeResiduals(mov_res,Result.mc(:,1).*Result.mc(:,end));
mov_res= SeeResiduals(mov_res,bkg,1);

mov_blueoff=mov_res(:,:,BlueTime{1} & spike_erode_trace{3}==0);
mov_blueon=mov_res(:,:,BlueTime{2} & spike_erode_trace{3}==0);


%% Calculate inter-distance
nROI=size(Result.ftprnt,3);
SkelDend = Skeletonize_dendrite(Result.ref_im,4,0.01,25);
interDendDist=[];
for i=1:nROI
    i
    for j=1:nROI
        [interDendDist(i,j), ~]=geodesic_distance(SkelDend,get_coord(Result.ftprnt(:,:,i)),get_coord(Result.ftprnt(:,:,j)));
    end
end

%% Low Stim plot
coord_1d=dim_reduce(get_coord(Result.ftprnt));
[~, dist_order]=sort(coord_1d,'descend');
som_roi=find(dist_order==1);
geodist=interDendDist(1,:)'.*sign(coord_1d-coord_1d(1));
show_footprnt_contour(Result.ftprnt(:,:,dist_order),Result.ref_im)
ref_spike_ROI=5;
tr_norm= Result.traces_bvMask-movprc(Result.traces_bvMask,100,30,2);
tr_norm= tr_norm./get_threshold(tr_norm,1);
Result.spike= find_spike_bh(tr_norm,5,3);
F0=tovec(Result.ref_im)'*tovec(Result.ftprnt.*double(max(Result.bvMask,[],3)==0));
NormTrace=Result.traces_bvMask./F0';
subth_trace=get_subthreshold(NormTrace,Result.spike(ref_spike_ROI,:),5,10);

figure(6); clf;
tiledlayout(5,1)
nexttile([1 1])
imshow2(Result.ref_im,[]); hold all
plot(bluePatt{1}([1:end 1],2),bluePatt{1}([1:end 1],1),'color',[0 0.6 1],'LineWidth',2)
ax1=[];
title_str={'No Stim.','Low Stim.','No Stim.','Low Stim.'};
for n=1:4
ax1=[ax1 nexttile([1 1])];
t=[1:2000]+2000*(n-1);
imagesc(NormTrace(dist_order,t),[-2 4]*0.01); hold all
sp=find(Result.spike(ref_spike_ROI,t));
plot(sp,ones(length(sp),1)+nROI-1,'color','k','marker','^','MarkerFaceColor',[1 0 0],'LineStyle','none')
colormap(ax1(n),turbo);
title(title_str{n})
end

%%
thres_off=[0.015 0.002]; thres_on=[0.015 0.002];
length_thres_off=10; length_thres_on=20;
AfterSpikeSub=squeeze(mean(reshape(subth_trace(:,find(Result.spike(1,:))'+[4:7]),nROI,[],4),[2 3]));
F0=imgaussfilt(Result.ref_im,4);

Subth_off=subth_trace(:,BlueTime{1} & spike_erode_trace{3}==0); Subth_off=Subth_off-median(Subth_off,2);
Subth_on=AfterSpikeSub-subth_trace(:,BlueTime{2} & spike_erode_trace{3}==0);
Subth_on=Subth_on-median(Subth_on,2);

troff=[]; tron=[]; EItrace=zeros(2,size(Result.traces,2));
for n=1:nROI
[trans_off]=detect_transient(Subth_off(n,:),thres_off,ones(1,size(Subth_off,2)));
troff(n,:)=ismember(trans_off.interval,find(trans_off.length>length_thres_off));
[trans_on]=detect_transient(Subth_on(n,:),thres_on,ones(1,size(Subth_on,2)));
tron(n,:)=ismember(trans_on.interval,find(trans_on.length>length_thres_on));
end

EItrace(1,BlueTime{1} & spike_erode_trace{3}==0)=max(troff,[],1);
EItrace(2,BlueTime{2} & spike_erode_trace{3}==0)=max(tron,[],1);
for e=1:2
    EIvec{e}=[];
    bwEI=bwlabel(EItrace(e,:));
    for b=1:max(bwEI)
        EIvec{e}(:,b)=mean(subth_trace(:,find(bwEI==b)),2);
    end
end


figure(55); clf; tiledlayout(3,2);
ax1=nexttile([1 2]);
imagesc(NormTrace(dist_order,:),[-0.02 0.05])

ax1=[ax1 nexttile([1 2])];
plot(EItrace')
colormap('Turbo')

linkaxes(ax1,'x')

nexttile([1 1])
ExcImg=std(-mov_res(bound:end-bound,bound:end-bound,EItrace(1,:)>0),0,3);
%imshow2(ExcImg./F0(bound:end-bound,bound:end-bound),[])
imshow2(imgaussfilt(ExcImg,2),[])

nexttile([1 1])
InhImg=std(mov_res(bound:end-bound,bound:end-bound,EItrace(2,:)>0),0,3);
%imshow2(InhImg./F0(bound:end-bound,bound:end-bound),[])
imshow2(imgaussfilt(InhImg,2),[])