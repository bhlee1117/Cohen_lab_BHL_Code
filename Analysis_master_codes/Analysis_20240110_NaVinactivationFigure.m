clear
clc;
cd '/Volumes/BHL18TB_D2/Arranged_Data/NaV_InactivationResult';

Blue_CAL=importdata('/Volumes/BHL18TB_D2/Arranged_Data/NaV_InactivationResult/BlueCalibration_BU.mat');
Ad_Result=importdata('/Volumes/BHL18TB_D2/20230901_PP68_Adaptation_IsoKetCpp/Result_V2CheRiffST_20240110.mat');
PC_Result=importdata('/Volumes/BHL18TB_D2/YQ601_PlaceCellResults/PF_Result_20240111.mat');
Slice_Result={'/Volumes/cohen_lab/Lab/Labmembers/Pojeong Park/Data/230731_PP070_Nav_inactivation/Cell1_PTX_50uM/173841PP070_P16_Cell1_soma_3sec_ramp+EFS_3V' ...
              '/Volumes/cohen_lab/Lab/Labmembers/Pojeong Park/Data/230731_PP070_Nav_inactivation/Cell1_PTX_50uM/174015PP070_P16_Cell1_soma_3sec_cons+EFS_3V'};

%% Optogenetic Stimulation (Pulse)
cBlue=[0 0.4 1];
cTrace=[0.1 0.1 0.1];
cCS=[0.8 0 0.2];

containString={'Awake','ket','CPP'};
for i=1:length(Ad_Result)
    sp=split(Ad_Result{i}.fpath,'_');
    cstr=find(contains(sp,'ell'));
    drug_str=find(contains(sp,containString));
    if isempty(drug_str)
        drug_ind(i)=0;
    else
        for j=1:length(containString)
        if ~isempty(find(contains(sp{drug_str},containString{j})))
        drug_ind(i)=j;
        end
        end
    end
    cellnumber(i)=str2num(cell2mat(regexp(sp{cstr},'\d+','match')));
end

Cell_list=unique(cellnumber);

for N=1:length(Cell_list)
Cell_ind=Cell_list(N);    
Scratch_list=find(cellnumber==Cell_ind);
[drug_sort, drug_sort_ind]=sort(drug_ind(Scratch_list),'ascend');
for d=1:4 % 4 states, Iso, Awake, Ket, CPP
Ad_Result_NList{N,d}=[];    
Scratch_state_list=find(drug_ind(Scratch_list)==(d-1));
for j=1:length(Scratch_state_list) %find cells match to state
n=Scratch_list(Scratch_state_list(j));

tr=Ad_Result{n}.trace;
tr=tr-prctile(tr,20);
Ad_Result_NList{N,d}.trace(j,:)=tr./get_threshold(tr,1);

Ad_Result_NList{N,d}.Blue(j,:)=Ad_Result{n}.Blue; 
Ad_Result_NList{N,d}.ref_im(j)=prctile(Ad_Result{n}.ref_im(:),90);
Ad_Result_NList{N,d}.spike(j,:)=Ad_Result{n}.spike;
Ad_Result_NList{N,d}.Waveform=sum(Ad_Result_NList{N,d}.Blue>0,2)>10000;
Ad_Result_NList{N,d}.Subth(j,:)=Ad_Result{n}.Subth;
Ad_Result_NList{N,d}.CS_trace(j,:)=Ad_Result{n}.CS_trace;

end
end
end

noi=[3 7 9]; Waveform=0; t_pulse=[3400:8400]; d=1;
show_ind=[2 1 3];
figure;
tiledlayout(4,1)

for g=1:length(noi)
N=noi(g);
tmp=find(Ad_Result_NList{N,d}.Waveform==0);
show_tr=Ad_Result_NList{N,d}.trace(show_ind(g),:);
F_zero=Ad_Result_NList{N,d}.ref_im(show_ind(g));
show_tr_CS=show_tr; show_tr_CS(Ad_Result_NList{N,d}.CS_trace(show_ind(g),:)==0)=NaN;
nexttile([1 1])
plot(show_tr(t_pulse)/F_zero,'color',cTrace); hold all; axis off
plot(show_tr_CS(t_pulse)/F_zero,'color',cCS)
line([t_pulse(end)-t_pulse(1)+50 t_pulse(end)-t_pulse(1)+50],[0.01 0.06],'color','k','linewidth',2)
end
nexttile([1 1])
plot(Ad_Result_NList{N,d}.Blue(show_ind(g),t_pulse),'color',cBlue); axis off

figure; t_Ramp=[9100:12100];
tiledlayout(4,1)
for g=1:length(noi)
N=noi(g);
tmp=find(Ad_Result_NList{N,d}.Waveform==0);
F_zero=Ad_Result_NList{N,d}.ref_im(show_ind(g));
show_tr=Ad_Result_NList{N,d}.trace(show_ind(g),:);
show_tr_CS=show_tr; show_tr_CS(Ad_Result_NList{N,d}.CS_trace(show_ind(g),:)==0)=NaN;
nexttile([1 1])
plot(show_tr(t_Ramp)/F_zero,'color',cTrace); hold all; axis off
plot(show_tr_CS(t_Ramp)/F_zero,'color',cCS)

end
nexttile([1 1])
plot(Ad_Result_NList{N,d}.Blue(show_ind(g),t_Ramp),'color',cBlue); axis off

%% Statistics
Pulse_bw=[19:23]; Ramp_bw=[24];
CSRate=[]; SpikeNumber=[]; CSNumber=[]; CSNumber_t=[]; SpikeNumber_t=[];
CSNumber_t_ramp=[]; SpikeNumber_t_ramp=[]; Blue_intensity=[];
pulsewaveform_list=cellfun(@(x) find(x.Waveform==0), Ad_Result_NList(:,1), 'UniformOutput', false);
cmap=distinguishable_colors(2);
BinTime=200; BinTime_R=150;
plotInd={[1 2],[1],[1 2 3],[1],[1 2],[1],[1 3],[5 6],[1 2],[1 2],[1]};
for n=1:size(Ad_Result_NList,1)
    BlueBW = bwlabel(Ad_Result_NList{1}.Blue(pulsewaveform_list{n}(1),:));
    [~, ind]=max(sum(Ad_Result_NList{n,1}.spike(pulsewaveform_list{n},:),2));
    ind=pulsewaveform_list{n}(ind);
    ind=pulsewaveform_list{n}(plotInd{n});
    Blue_AOTF=[];

    for s=1:length(Pulse_bw)
        StimTime=find(BlueBW==Pulse_bw(s));        
    tBin=[1:floor(length(StimTime)/BinTime)]*BinTime;
    tBin=[[1 tBin]' [tBin(1:end) length(StimTime)]'];
        Sp_trace = Ad_Result_NList{n,1}.spike(ind,StimTime);
        CS_trace = Ad_Result_NList{n,1}.spike(ind,StimTime) .* (Ad_Result_NList{n,1}.CS_trace(ind,StimTime)>0);
        SpikeNumber{n}(:,s) = sum(Sp_trace,2);
        CSNumber{n}(:,s)    = sum(CS_trace,2);
        Blue_AOTF(:,s) = max(Ad_Result_NList{n,1}.Blue(ind,StimTime),[],2);
        Blue_intensity{n}(:,s) = Blue_CAL(find_index_bh(Blue_CAL(:,1),round(Blue_AOTF(:,s),1)),2);

        for t=1:size(tBin,1)
            CSNumber_t{n}(:,t,s)    = sum(CS_trace(:,tBin(t,1):tBin(t,2)),2);
            SpikeNumber_t{n}(:,t,s)    = sum(Sp_trace(:,tBin(t,1):tBin(t,2)),2);
        end
    end
    CSRate{n}=CSNumber{n}./SpikeNumber{n};

    for s=1:length(Ramp_bw)
        StimTime=find(BlueBW==Ramp_bw(s));    
         tBin=[1:floor(length(StimTime)/BinTime_R)]*BinTime_R;
    tBin=[[1 tBin]' [tBin(1:end) length(StimTime)]'];
        Sp_trace_ramp = Ad_Result_NList{n,1}.spike(ind,StimTime);
        CS_trace_ramp = Ad_Result_NList{n,1}.spike(ind,StimTime) .* (Ad_Result_NList{n,1}.CS_trace(ind,StimTime)>0);

        for t=1:size(tBin,1)
            CSNumber_t_ramp{n}(:,t,s)    = sum(CS_trace_ramp(:,tBin(t,1):tBin(t,2)),2);
            SpikeNumber_t_ramp{n}(:,t,s)    = sum(Sp_trace_ramp(:,tBin(t,1):tBin(t,2)),2);
        end
    end

end

Blue_intensity=mean(cell2mat(Blue_intensity'),1);
CSRateMat=cell2mat(CSRate');
M = mean(CSRateMat,1,'omitnan'); S = std(CSRateMat,0,1,'omitnan');

figure(2); clf;
for d=1:size(CSRateMat)
    plot(Blue_intensity,CSRateMat(d,:),':','color','k','LineWidth',0.5); hold all
end
errorbar(Blue_intensity,M,S,'color',cmap(1,:),'LineWidth',2)
xlim([0 8])
xlabel('Blue intensity (mW/mm^2)')
ylabel('Fraction of complex spike')

%CSNumber_tMat=cell2mat(CSNumber_t');
CSNumber_tMat=cell2mat(CSNumber_t')./cell2mat(SpikeNumber_t');
t=BinTime*[0.5:size(CSNumber_tMat,2)-0.5]/1000;
figure(3); clf;
cmap_t=turbo(5);
for s=1:5
    N_t=sum(~isnan(CSNumber_tMat(:,:,s)));
    M_t=mean(CSNumber_tMat(:,:,s),1,'omitnan');
    S_t=std(CSNumber_tMat(:,:,s),0,1,'omitnan');
    errorbar(t,M_t,S_t./sqrt(N_t),'color',cmap_t(s,:),'LineWidth',2); hold all
end
xlim([0 0.6])
xlabel('Time after stimulus onset (s)')
ylabel('Fraction of complex spike')
legend({'0.3 mW/mm^2','1.4 mW/mm^2','3.2 mW/mm^2','5.4 mW/mm^2','7.7 mW/mm^2'})

figure(4); clf;

%CSNumber_tMat_ramp=cell2mat(CSNumber_t_ramp')./cell2mat(SpikeNumber_t_ramp');
t=BinTime_R*[0.5:size(cell2mat(CSNumber_t_ramp'),2)-0.5]/1000;
cmap_t_ramp=[0.1 0.1 0.1; 0.8 0 0.2];
CS_sum=sum(cell2mat(CSNumber_t_ramp'),1);
Sp_sum=sum(cell2mat(SpikeNumber_t_ramp'),1);
bar(t, Sp_sum, 'FaceColor', cmap_t_ramp(1,:)); hold all
bar(t, CS_sum, 'FaceColor', cmap_t_ramp(2,:));

xlabel('Time after stimulus onset (s)')
ylabel('Number of spikes')
legend({'Total spikes (n=19 neurons)','Complex spikes (n=19 neurons)'})
xlim([0 3.1])

%% Place cell data
% Show all the place cell with complex spiking map

cmap=[0.5 0.05 0.15;0.1 0.1 0.1]/5;
f1=figure(1); clf;
coarse_bin=3;
for i=1:length(PC_Result)
    [sz_fprntY, sz_fprntX, nNeurons]=size(PC_Result{i}.c_ftprnt);
    StimN=find(mean(PC_Result{i}.spike(:,find(PC_Result{i}.Blue>0)),2)>mean(PC_Result{i}.spike,2)*5)';
    UnstimN=setdiff([1:nNeurons],StimN);
    for n=[StimN UnstimN]
        nexttile([1 1])
        pos_bin=size(PC_Result{i}.Lap_CS,2);
        
        Lap_tmp_CS=repmat(PC_Result{i}.Lap_CS(:,:,n),1,3);
        Lap_tmp_CS=movmean(Lap_tmp_CS,coarse_bin,2);
        Lap_tmp_CS=(Lap_tmp_CS(:,pos_bin+1:2*pos_bin));
        

        Lap_tmp_SS=repmat(PC_Result{i}.Lap_FR(:,:,n),1,3);
        Lap_tmp_SS=movmean(Lap_tmp_SS,coarse_bin,2);
        Lap_tmp_SS=(Lap_tmp_SS(:,pos_bin+1:2*pos_bin));
        
        composite_LR=(squeeze(sum(cat(3,Lap_tmp_CS,Lap_tmp_SS).*reshape(cmap,1,1,[],3),3)));
        %imagesc(reshape(rescale(composite_LR(:)),[],pos_bin,3))
        imagesc(composite_LR)
 
        if ismember(n,UnstimN)
            title(['Unstimulated ' num2str(i) '-' num2str(n)])
        else
            title([num2str(i) '-' num2str(n)])
        end
    end
end


%% Filter good PCs

PC_list=[1 1;3 1;3 3;4 2;4 3;5 1;6 1;9 1;12 3;18 1;20 1;23 1;27 1;30 1;32 1;32 2;35 1;35 2];
show_PC=7; pos_bin=size(PC_Result{PC_list(show_PC,1)}.Lap_FR,2);
Lap_FR_low=movmean(repmat(PC_Result{PC_list(show_PC,1)}.Lap_FR(:,:,PC_list(show_PC,2)),1,3),5,2);
Lap_FR_low=Lap_FR_low(:,pos_bin+1:2*pos_bin);
meanPFmap=mean(Lap_FR_low,1,'omitnan');
EndLap=max(PC_Result{PC_list(show_PC,1)}.VR(8,:));
StimLap=unique(double(PC_Result{PC_list(show_PC,1)}.Blue>0).*PC_Result{PC_list(show_PC,1)}.VR(8,:));
StimLap=StimLap(2:end);
ShowLap=setdiff([1:EndLap],[StimLap 26 15]);

[~, PF]=detect_transient(meanPFmap,[10 1],[]);
[~, PF_peak]=max(meanPFmap);

show_tr=PC_Result{PC_list(show_PC,1)}.normTraces(PC_list(show_PC,2),:)-movprc(PC_Result{PC_list(show_PC,1)}.normTraces(PC_list(show_PC,2),:),300,30);
show_traces_align_Position_wCS(show_tr,PC_Result{PC_list(show_PC,1)}.CS_trace(PC_list(show_PC,2),:), ...
    median(find(PF>0))*115/pos_bin, [-2000:2000],PC_Result{PC_list(show_PC,1)}.VR,ShowLap)

Lap_tmp_CS=repmat(PC_Result{PC_list(show_PC,1)}.Lap_CS(:,:,PC_list(show_PC,2)),1,3);
Lap_tmp_CS=movmean(Lap_tmp_CS,5,2);
Lap_tmp_CS=(Lap_tmp_CS(:,pos_bin+1:2*pos_bin));

bin_rep=repmat([1:pos_bin],1,2);
PFcenter_bin=bin_rep(find(bin_rep==round(median(find(PF>0))),1,'last')-pos_bin/2+1:find(bin_rep==round(median(find(PF>0))),1,'last')+pos_bin/2);

Lap_tmp_SS=repmat(PC_Result{PC_list(show_PC,1)}.Lap_FR(:,:,PC_list(show_PC,2)),1,3);
Lap_tmp_SS=movmean(Lap_tmp_SS,5,2);
Lap_tmp_SS=(Lap_tmp_SS(:,pos_bin+1:2*pos_bin));

composite_LR=(squeeze(sum(cat(3,Lap_tmp_CS,Lap_tmp_SS).*reshape(cmap,1,1,[],3),3)));
%imagesc(reshape(rescale(composite_LR(:)),[],pos_bin,3))
figure;
imagesc(composite_LR(ShowLap,PFcenter_bin,:))

%%

[~, ~, PF_list]=xlsread(['/Volumes/BHL18TB_D2/Arranged_Data/NaV_InactivationResult/' ...
    'PlaceFieldData_Arrangement_20210112.xlsx'], 'Sheet1', 'A5:E24');
PF_list= cell2mat(PF_list);
pos_bin=150;
PF_Bin=30; Bound=5;
FR_PF_Bin=[]; CS_PF_Bin=[];
for i=1:size(PF_list,1)
    PC_fileind=PC_list(PF_list(i,1),1); PC_Nind=PC_list(PF_list(i,1),2);

    SP_map=PC_Result{PC_fileind}.Lap_FR(:,:,PC_Nind);
    CS_map=PC_Result{PC_fileind}.Lap_CS(:,:,PC_Nind);

    SP_map=repmat(SP_map,1,2); CS_map=repmat(CS_map,1,2);
    if PF_list(i,2)>PF_list(i,3)
    FR_PF= mean(SP_map(PF_list(i,4):PF_list(i,5),PF_list(i,2)-Bound:PF_list(i,3)+Bound+pos_bin),1,'omitnan');
    CS_PF= mean(CS_map(PF_list(i,4):PF_list(i,5),PF_list(i,2)-Bound:PF_list(i,3)+Bound+pos_bin),1,'omitnan');
    else
    FR_PF= mean(SP_map(PF_list(i,4):PF_list(i,5),PF_list(i,2)-Bound:PF_list(i,3)+Bound),1,'omitnan');
    CS_PF= mean(CS_map(PF_list(i,4):PF_list(i,5),PF_list(i,2)-Bound:PF_list(i,3)+Bound),1,'omitnan');
    end
OldBin=[0:size(CS_PF,2)-1];
NewBin=[0:(size(CS_PF,2)-1)/(PF_Bin-1):size(CS_PF,2)-1];
FR_PF_Bin(i,:)=interp1([0:size(FR_PF,2)-1],FR_PF,NewBin,'linear');
CS_PF_Bin(i,:)=interp1([0:size(CS_PF,2)-1],CS_PF,NewBin,'linear');
end

figure(2); clf;
CSrate_PF_Bin=movmean(CS_PF_Bin./FR_PF_Bin,3,2);
plot(rescale(NewBin),CSrate_PF_Bin','color',[0.7 0.7 0.7]); hold all
errorbar(rescale(NewBin),mean(CSrate_PF_Bin,1),std(CSrate_PF_Bin,0,1)./sqrt(size(PF_list,1)),'LineWidth',2)

%% Slice data
EndFrame=3850; ROI=[];
for f=1:length(Slice_Result)
load(fullfile(Slice_Result{f},"output_data.mat"))
sz=double(Device_Data{1, 4}.ROI([2 4]));
CamCounter=Device_Data{1, 2}.Counter_Inputs(1, 1).data;
CamTrigger=find(CamCounter(2:end)-CamCounter(1:end-1));
mov=double(readBinMov([Slice_Result{f} '/frames1.bin'],sz(2),sz(1)));

    Blue{f}=Device_Data{1, 2}.buffered_tasks(1, 1).channels(2).data(CamTrigger);
    EFS{f}=Device_Data{1, 2}.buffered_tasks(1, 1).channels(3).data(CamTrigger);

    mov=mov(:,:,1:EndFrame); Blue{f}=Blue{f}(1:EndFrame); EFS{f}=EFS{f}(1:EndFrame);
    bkg = zeros(1, size(mov,3));
    bkg(1,:)=movmedian(get_blueoffTrace(squeeze(mean(mov,[1 2])),[Blue{f}],30),200);
    
    mov_res= mov-reshape(bkg,1,1,[]);
    if isempty(ROI)
    [ROI, Som_tr{f}]=clicky(mov_res,mean(mov,3));
    [F_tr{f}]=apply_clicky(ROI,mov,'no');
    else
    [Som_tr{f}]=apply_clicky(ROI,mov_res,'no');
    [F_tr{f}]=apply_clicky(ROI,mov,'no');
    end
end

%% Plot
figure; clf; tiledlayout(8,1); 
cmap_slice=[turbo(8)];
for f=1:length(Slice_Result)
nexttile([3 1])
dFF=-(Som_tr{f})./median(F_tr{f},1); dFF=dFF-median(dFF,1);
%dFF=Som_tr{f}-median(Som_tr{f},1); dFF=-(dFF'./get_threshold(dFF',1))';
ls=plot(dFF(20:EndFrame,[1 3 4])-[1:3]/20,'color',[0.1 0.1 0.1]); axis tight off; hold all
arrayfun(@(l,c) set(l,'Color',c{:}),ls,num2cell(cmap_slice([2 5 7],:),2))
%plot(rescale2(-Som_tr{f}(20:EndFrame,[1 3 4]),1)-[1:3],'color',[0.1 0.1 0.1]); axis off; hold all
nexttile([1 1])
plot(Blue{f}(20:EndFrame),'color',[0 0.5 1],'LineWidth',1); axis tight off; hold all
plot(rescale(EFS{f}(20:EndFrame)),'color',[1 0.5 0],'LineWidth',1); axis off
end