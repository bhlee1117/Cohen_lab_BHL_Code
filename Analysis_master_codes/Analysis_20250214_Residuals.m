%% Set the path
clear
clc;
cd '/Users/bhlee1117/Documents/GitHub/Cohen_lab_BHL_Code/Analysis_master_codes';
[~, ~, raw] = xlsread(['/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/' ...
    'Prism_OptopatchData_Arrangement.xlsx'], 'Sheet1', 'C5:M164');

save_to='/Volumes/BHL18TB_D2/Arranged_Data/Prism_OptopatchResult';
fpath=raw(:,1);
Mouse=cell2mat(raw(:,2));
NeuronInd=cell2mat(raw(:,5));
CamType=raw(:,3);
StructureData=raw(:,10);
time_segment=25000;
%% Load the data

i=102; bound=6;
    load([fpath{i} '/OP_Result.mat'])
    cd(fpath{i});
    load(fullfile(fpath{i},"output_data.mat"))
    sz=double(Device_Data{1, 3}.ROI([2 4]));

    Result.Blue=Device_Data{1, 2}.buffered_tasks(1, 1).channels(1, 2).data;
    CamCounter=Device_Data{1, 2}.Counter_Inputs(1, 1).data;
    CamTrigger=find(CamCounter(2:end)-CamCounter(1:end-1));
    DMDtrigger=Device_Data{1, 2}.buffered_tasks(1, 2).channels(1, 3).data;
    DMDtrigger=DMDtrigger(CamTrigger);
    Result.DMDtrigger=[0 (DMDtrigger(2:end)-DMDtrigger(1:end-1))>0];
    Result.Blue=Result.Blue(CamTrigger);
    [Result.blueDMDimg Result.bluePatt]=get_blueDMDPatt(Device_Data,'stack');

    DMDbluetrace=(Result.Blue>0).*cumsum(Result.DMDtrigger)+1;
    DMDbluetrace(Result.Blue==0)=0;

    mov_mc=readBinMov_BHL(fpath{i},3);
    load(fullfile(fpath{i},'mcTrace01.mat'))
    mean_F=squeeze(mean(mov_mc(bound:end-bound,bound:end-bound,:),[1 2]));
    [~, blueOff]=get_blueoffTrace(mean_F,[Result.Blue],70);
    [y_fit]=expfitDM_2(find(blueOff)',mean_F(find(blueOff)),[1:size(mov_mc,3)]',1000);

mov_res= mov_mc-mean(mov_mc,3);
mov_res = SeeResiduals(mov_res,y_fit);
mov_res = SeeResiduals(mov_res,mcTrace.xymean(:,:));
mov_res = SeeResiduals(mov_res,mcTrace.xymean(:,:).^2);
mov_res = SeeResiduals(mov_res,mcTrace.xymean(:,1).*mcTrace.xymean(:,end));
dirtMov_dilate = tracking_dirt(mov_res,0.3);
% %%
% [u,s,v] = svds(tovec(mov_res(:,:,1:5000)-mean(mov_res(:,:,1:5000),3)),20);
% reshape_u=reshape(u,sz(2),sz(1),[]);
% bvMask=[];
% [~, bvMask]=get_ROI(max(abs(reshape_u),[],3),bvMask);
% 
% Result.bvMask=bvMask;
% Result.traces_bvMask=(tovec(mov_res.*double(max(Result.bvMask,[],3)==0))'*tovec(Result.ftprnt))';
 Result.dirtTrace=(tovec(dirtMov_dilate)'*tovec(Result.ftprnt))';
%% Get STAs

SomTr=Result.traces_bvMask(1,:);
SomTr=SomTr/get_threshold(SomTr,1);
Result.spike(1,:)=find_spike_bh(SomTr,5,3);
VoltageTr=Result.traces_bvMask;
VoltageTr(Result.dirtTrace>0)=NaN;
%Result.F0_PCA=get_F0PCA(VoltageTr,3);
VoltageTr=VoltageTr./Result.F0_PCA;
%VoltageTr=pcafilterTrace(VoltageTr,10);

mov_res_filt=mov_res;
mov_res_filt(dirtMov_dilate>0)=NaN;
mov_vec=tovec(-mov_res_filt);


STAmov=[]; STAtr=[];
for b=1:max(DMDbluetrace)
PattBlue=double(DMDbluetrace==b);
onsetTime=find((PattBlue(2:end)-PattBlue(1:end-1))==1)+1;
PattBlueOnset=ind2vec(length(DMDbluetrace),onsetTime,1);
[~, spMat]=get_STA(Result.spike(1,:),PattBlueOnset,30,50);
Nsp=squeeze(sum(spMat,3));
[~, m]=get_STA(mov_vec,PattBlueOnset,30,50);
STAmov(:,:,:,b)=toimg(squeeze(mean(m(:,Nsp==0,:),2,'omitnan')),sz(2),sz(1));
[~, trMat]=get_STA(VoltageTr,PattBlueOnset,30,50);
STAtr(:,:,b)=squeeze(mean(trMat(:,Nsp==0,:),2,'omitnan'));
end


%%
figure(11); clf;
STAtr=STAtr-median(STAtr(:,1:15,:),2,'omitnan');
imshow_patch(STAtr(Result.dist_order,:,:),[-0.005 0.02])
colormap(turbo)
figure(12); clf;
STAmov=STAmov-median(STAmov(:,:,1:15,:),3,'omitnan');
imshow_patch(squeeze(mean(imgaussfilt(STAmov(5:end-5,5:end-5,40:80,:),2),3,'omitnan')),[-0.5 2]);
colormap(turbo)