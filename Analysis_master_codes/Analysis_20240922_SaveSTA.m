clear
clc;
[~, ~, raw] = xlsread(['/Volumes/cohen_lab/Lab/Labmembers' ...
    '/Byung Hun Lee/Data/PrismPCdata_Arrangement.xlsx'], 'Sheet1', 'C5:Q31');

ref_ROI=cellfun(@(x) (str2num(num2str(x))),raw(:,10),'UniformOutput',false);
basal_ROI=cellfun(@(x) (str2num(num2str(x))),raw(:,11),'UniformOutput',false);
apical_ROI=cellfun(@(x) (str2num(num2str(x))),raw(:,12),'UniformOutput',false);
fpath=raw(:,1)';
StructureData=raw(:,10);
BadROI=cellfun(@(x) (str2num(num2str(x))),raw(:,15),'UniformOutput',false);
EndFrame=cell2mat(raw(:,13));
ifmotionReject=cell2mat(raw(:,14));
save_figto='/Volumes/BHL_WD18TB/PP72_PlaceCellResults';
place_bin=150; time_segment=15000; overlap=200;
alignedMovFN = {'STA_Mat_SS','STA_Mat_CS','STA_Mat_dSP','STA_Mat_BS'};
bound=6;
title_str={'Basal','Apical','Peri-Soma'};
%%
f=26; load(fullfile(fpath{f},'PC_Result.mat'),'Result')
cd(fpath{f})
rois={basal_ROI{f},apical_ROI{f},ref_ROI{f}};
nROI=size(Result.normTraces,1);
nTau_bAP=[-20:20];
nTau={[-50:10],[-50:10],[-50:10],[-50:10]}; %SS, CS, dSP, Brst
% Isolated Somatic spike
som_spike=find(Result.spike(1,:)); 
ss_time=find(Result.SpClass(1,:)); brst=bwlabel((ss_time(2:end)-ss_time(1:end-1))<15);
SpClass=Result.SpClass; BS_trace=zeros(1,size(Result.traces,2));
for b=1:max(brst)
    bwn=find(brst==b);
    SpClass(1,ss_time([bwn bwn(end)+1]))=0;
    SpClass(4,ss_time([bwn(1)]))=1;
    BS_trace(1,[ss_time(bwn): ss_time(bwn(end)+1)])=b;
end
tr_ref=Result.normTraces(ref_ROI{f},:);
tr_sub=mean(tr_ref,1)-movprc(mean(tr_ref,1),200,20,2);
tr_sub=get_subthreshold(tr_sub,Result.spike(1,:),5,10);
[trans tr_trace]=detect_transient2(tr_sub,[5 1.5],Result.spike(1,:),15);

dist_order=Result.dist_order;

bAP_ref=[];
for s=som_spike
    isnearby=sum(ismember(s+nTau_bAP,som_spike))>1;
    isnearbyCS=sum(ismember(s+nTau_bAP,find(Result.CStrace)))>1;
    ispartCS=tr_trace(s)>0;
    if ~isnearby & ~isnan(s) & ~isnearbyCS & ~ispartCS
        bAP_ref=[bAP_ref s];
    end
end
sp_na=sum((bAP_ref'+nTau_bAP)<0 | (bAP_ref'+nTau_bAP)>size(Result.traces,2),2)==0;
bAP_ref=bAP_ref(sp_na);

% Isolated Somatic spike
bAP_s=[];
for s=som_spike
    isnearby=sum(ismember(s+nTau{1},som_spike))>1;
    isnearbyCS=sum(ismember(s+nTau{1},find(Result.CStrace)))>1;
    ispartCS=tr_trace(s)>0;
    if ~isnearby & ~isnan(s) & ~isnearbyCS & ~ispartCS
        bAP_s=[bAP_s s];
    end
end
sp_na=sum((bAP_s'+nTau{1})<0 | (bAP_s'+nTau{1})>size(Result.traces,2),2)==0;
bAP_s=bAP_s(sp_na);

dSP_s=[];
for s=find(Result.SpClass(3,:))
    isnearby=sum(ismember(s+nTau{3},som_spike))>1;
    isnearbyCS=sum(ismember(s+nTau{3},find(Result.CStrace)))>1;
    ispartCS=tr_trace(s)>0;
    if ~isnearby & ~isnan(s) & ~isnearbyCS & ~ispartCS
        dSP_s=[dSP_s s];
    end
end
sp_na=sum((dSP_s'+nTau{3})<0 | (dSP_s'+nTau{3})>size(Result.traces,2),2)==0;
dSP_s=dSP_s(sp_na);

% Isolated Complex spike
CS_s=[];
C_list=find(Result.SpClass(2,:));
CS_label=max(Result.spike,[],1).*bwlabel(Result.CStrace);
CS_list=[];
for g=1:length(C_list)
    s=C_list(g);
    s_tmp=Result.spike(1,:);
    s_tmp(find(CS_label==g))=0;
    isnearbyCS=max(bwlabel(Result.CStrace(s+nTau{2})))>1;
    isnearbyS=sum(ismember(s+nTau{2},find(s_tmp)))>0;
    if ~isnearbyCS & ~isnearbyS
        CS_s=[CS_s s];
        CS_list=[CS_list g];
    end
end
sp_na=sum((CS_s'+nTau{2})<0 | (CS_s'+nTau{2})>size(Result.traces,2),2)==0;
CS_s=CS_s(sp_na);

BS_s=[];
B_list=find(SpClass(4,:));
for g=1:length(B_list)
    s=B_list(g);
    
    isnearby=sum(ismember(s+nTau{4},find(BS_trace~=BS_trace(s) & Result.spike(1,:)==1)))>0;
    isnearbyCS=sum(ismember(s+nTau{4},find(Result.CStrace)))>0;
    if ~isnearbyCS & ~isnearby
        BS_s=[BS_s s];
    end
end
sp_na=sum((BS_s'+nTau{4})<0 | (BS_s'+nTau{4})>size(Result.traces,2),2)==0;
BS_s=BS_s(sp_na);


STA_SS=squeeze(mean(reshape(Result.normTraces(:,bAP_s'+nTau_bAP),nROI,[],length(nTau_bAP)),2));
STA_SS= STA_SS - mean(STA_SS(:,1:10),2);
F_ref=mean(STA_SS(:,-nTau_bAP(1)+[5:10]),2);
%F_ref=(tovec(imgaussfilt(Result.ref_im,1))'*tovec(Result.ftprnt)/Result.SpikeHeight_fit(1))';

SilentPeriod=ones(1,size(Result.traces,2));
sp_time=find(max(Result.spike,[],1))';
sp_na=sum((find(max(Result.spike,[],1))'+[-10:150])<0 | (find(max(Result.spike,[],1))'+[-10:150])>size(Result.traces,2),2)==0;
SilentPeriod(sp_time(sp_na)+[-10:150])=NaN;
t_fit=find(~isnan(SilentPeriod) & Result.Blue==0);

NormalizedTrace=(Result.normTraces)./F_ref;
% lwpass=NaN(size(NormalizedTrace));
% lwpass(:,t_fit)=NormalizedTrace(:,t_fit);
% lwpass=movmedian(lwpass,30000,2,'omitnan');
% NormalizedTrace=NormalizedTrace-lwpass;
NormalizedTrace_dirt=NormalizedTrace;
NormalizedTrace_dirt(Result.dirtTrace>0)=NaN;
%NormalizedTrace_dirt(:,17000:18000)=NaN;
NormalizedTrace_dirt(:,Result.motionReject)=NaN;
NormalizedTrace_ch=cellfun(@(x) x./F_ref,Result.norm_trace_check,'UniformOutput',false);
NormalizedTrace_ch{1}(:,Result.motionReject)=NaN; NormalizedTrace_ch{2}(:,Result.motionReject)=NaN;
NormalizedTrace_ch{1}(Result.dirtTrace>0)=NaN; NormalizedTrace_ch{2}(Result.dirtTrace>0)=NaN;

STA_CSmat=reshape(NormalizedTrace_dirt(:,CS_s'+nTau{2}),nROI,[],length(nTau{2}));
STA_SSmat=reshape(NormalizedTrace_dirt(:,bAP_s'+nTau{1}),nROI,[],length(nTau{1}));
STA_dSPmat=reshape(NormalizedTrace_dirt(:,dSP_s'+nTau{3}),nROI,[],length(nTau{3}));
STA_BSmat=reshape(NormalizedTrace_dirt(:,BS_s'+nTau{4}),nROI,[],length(nTau{4}));

STA_CSmat_ch=cellfun(@(x) reshape(x(:,CS_s'+nTau{2}),nROI,[],length(nTau{2})),NormalizedTrace_ch,'UniformOutput',false);
STA_SSmat_ch=cellfun(@(x) reshape(x(:,bAP_s'+nTau{1}),nROI,[],length(nTau{1})),NormalizedTrace_ch,'UniformOutput',false);
STA_dSPmat_ch=cellfun(@(x) reshape(x(:,dSP_s'+nTau{3}),nROI,[],length(nTau{3})),NormalizedTrace_ch,'UniformOutput',false);
STA_BSmat_ch=cellfun(@(x) reshape(x(:,BS_s'+nTau{4}),nROI,[],length(nTau{4})),NormalizedTrace_ch,'UniformOutput',false);
%%
bound=6;
%f_tmp='/Volumes/BHL18TB_D1/20240218/134705BHLm117_FOV2_VR2';
spikeframes_class={bAP_s,CS_s,dSP_s,BS_s};

    StackedSpike=[];
    alignmovlist=dir(fullfile(fpath{f},'STA_Mat*.tiff'));
    alignmov_ind=ones(1,4);
    for l=1:length(alignmovlist)
        delete(fullfile(fpath{f},alignmovlist(l).name));
    end

    for c=1:4 % Get frames of spikes of each classes
        fileInd{c}=ceil(spikeframes_class{c}/time_segment);
        frameInd{c}=mod(spikeframes_class{c},time_segment);
    end
    %line up movie segments at simple spike

    load([fpath{f} '/output_data.mat'])
    sz=double(Device_Data{1, 3}.ROI([2 4]));
    frm_end=max(Device_Data{1, 2}.Counter_Inputs(1, 1).data);
    f_seg=[[1:time_segment:frm_end] frm_end+1]; f_seg(2:end)=f_seg(2:end)-1;
    g=zeros(1,4);
    ss=[0 0 0];
    for j=unique(cell2mat(fileInd))
        j

        mov_mc=double(readBinMov([fpath{f} '/mc_ShutterReg' num2str(j,'%02d') '.bin'],sz(2),sz(1)));
        load([fpath{f} '/mcTrace' num2str(j,'%02d') '.mat']);

        if j==length(f_seg)-1
            mc=mcTrace.xymean;
        else
            mov_mc=mov_mc(:,:,1:time_segment);
            mc=mcTrace.xymean(1:time_segment,:);
        end

        mov_res= mov_mc-mean(mov_mc,3);
        bkg = zeros(2, size(mov_mc,3));
        bkg(1,:) = linspace(-1, 1, size(mov_mc,3));  % linear term
        bkg(2,:) = linspace(-1, 1, size(mov_mc,3)).^2;  % quadratic term
        mov_res = SeeResiduals(mov_res,mc);
        mov_res = SeeResiduals(mov_res,mc.^2);
        mov_res = SeeResiduals(mov_res,mc(:,1).*mc(:,end));
        mov_res= SeeResiduals(mov_res,bkg,1);

        for c=1:4
            STA_movie_align{c}=[];
            for k=find(fileInd{c}==j)
                if frameInd{c}(k)+nTau{c}(1)>1 && frameInd{c}(k)+nTau{c}(end)<size(mov_res,3)
                    mov_seg=mov_res(bound:end-bound,bound:end-bound,frameInd{c}(k)+nTau{c});
                    %mov_seg=mat2gray(mov_seg);
                    STA_movie_align{c} = cat(3,STA_movie_align{c}, mov_seg);
                    g(c)=g(c)+1;
                    StackedSpike{c}(1,g(c))=1;
                    StackedSpike{c}(2,g(c))=time_segment*(j-1)+frameInd{c}(k);
                else
                    STA_movie_align{c} = cat(3,STA_movie_align{c}, zeros(sz(2)-bound*2+1,sz(1)-bound*2+1,length(nTau{c})));
                    g(c)=g(c)+1;
                    StackedSpike{c}(1,g(c))=0;
                    StackedSpike{c}(2,g(c))=NaN;
                end
            end
            if ~isempty(STA_movie_align{c})

                % write_tif_stack(STA_movie_align{c},fullfile(f_tmp,[alignedMovFN{c} '_' num2str(alignmov_ind(c)) '.tiff']))
                % alignmovlist=dir(fullfile(f_tmp,[alignedMovFN{c} '*.tiff']));  
                write_tif_stack(STA_movie_align{c},fullfile(fpath{f},[alignedMovFN{c} '_' num2str(alignmov_ind(c)) '.tiff']))
                alignmovlist=dir(fullfile(fpath{f},[alignedMovFN{c} '*.tiff']));
                if alignmovlist(alignmov_ind(c)).bytes > 2.0*10^9
                    disp('Move on to 2nd tiff')
                    alignmov_ind(c) = alignmov_ind(c) + 1;
                end
            end
        end
        g
    end

    Result.StackedSpike=StackedSpike;
    save(fullfile(fpath{f},'PC_Result.mat'),'Result','fpath','-v7.3')

