%% Set the path
clear
clc;
cd '/Users/bhlee1117/Documents/GitHub/Cohen_lab_BHL_Code/Analysis_master_codes';
[~, ~, raw] = xlsread(['/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/' ...
    'Prism_OptopatchData_Arrangement.xlsx'], 'Sheet1', 'C5:R205');

save_to='/Volumes/BHL18TB_D2/Arranged_Data/Prism_OptopatchResult';
fpath=raw(:,1);
Mouse=cell2mat(raw(:,2));
NeuronInd=cell2mat(raw(:,5));
CamType=raw(:,3);
StructureData=raw(:,10);
StimProtocol=raw(:,6);
TimeSegFrame=cell2mat(cellfun(@(x) (str2num(num2str(x))),raw(:,11),'UniformOutput',false));
PixelSize=cell2mat(cellfun(@(x) (str2num(num2str(x))),raw(:,12),'UniformOutput',false));
[~, unqInd] = unique([Mouse NeuronInd] ,'row');
set(0,'DefaultFigureWindowStyle','docked')
foi=[124 111 195:200];
% %% Set clicky ROIs
% for f=foi(3:end)
% load(fullfile(fpath{f},'output_data.mat'))
% load(fullfile(fpath{f},'mcTrace01.mat'))
% Blue=Device_Data{1, 2}.buffered_tasks(1, 1).channels(1, 2).data;
% CamCounter=Device_Data{1, 2}.Counter_Inputs(1, 1).data;
% CamTrigger=find(CamCounter(2:end)-CamCounter(1:end-1));
% Blue=Blue(CamTrigger);
% [blueDMDimg, bluePatt, blueCenter]=get_blueDMDPatt(Device_Data,'stack');
% sz=double(Device_Data{1, 3}.ROI([2 4]));
% DMDtrigger = Device_Data{1, 2}.buffered_tasks(1, 2).channels(1, 3).data;
% DMDtrigger = DMDtrigger(CamTrigger);
% DMDtrigger = [0 (diff(DMDtrigger) > 0)];
% DMDbluetrace=(Blue>0).*(cumsum(DMDtrigger)+1);
% DMDbluetrace(Blue==0)=0;
% mov_mc=double(readBinMov_times(fullfile(fpath{f},'mc_ShutterReg01.bin'),sz(2),sz(1),[1:1000]));
% 
% figure(10); clf;
% show_footprnt_contour(blueDMDimg,mean(mov_mc,3))
% 
% figure(11); clf;
% Result.clicky_ROI=clicky(mov_mc);
% save(fullfile(fpath{f},'OP_Result.mat'),"Result",'-v7.3')
% end
%%
nTau=[50 150];
for f=foi(1:end)
    f
    clear mov_res mov_mc mov_res_filt dirtMov_dilate
load(fullfile(fpath{f},'output_data.mat'))
load(fullfile(fpath{f},'mcTrace01.mat'))
load(fullfile(fpath{f},'OP_Result.mat'))
Blue=Device_Data{1, 2}.buffered_tasks(1, 1).channels(1, 2).data;
CamCounter=Device_Data{1, 2}.Counter_Inputs(1, 1).data;
CamTrigger=find(CamCounter(2:end)-CamCounter(1:end-1));
Blue=Blue(CamTrigger);
[blueDMDimg, bluePatt, blueCenter]=get_blueDMDPatt(Device_Data,'stack');
sz=double(Device_Data{1, 3}.ROI([2 4]));
DMDtrigger = Device_Data{1, 2}.buffered_tasks(1, 2).channels(1, 3).data;
DMDtrigger = DMDtrigger(CamTrigger);
DMDtrigger = [0 (diff(DMDtrigger) > 0)];
DMDbluetrace=(Blue>0).*(cumsum(DMDtrigger)+1);
DMDbluetrace(Blue==0)=0;

mov_mc=double(readBinMov_times(fullfile(fpath{f},'mc_ShutterReg01.bin'),sz(2),sz(1),[1:length(CamTrigger)]));
dirtMov_dilate = double(readBinMov_times(fullfile(fpath{f},'dirt_mov.bin'),sz(2),sz(1),[1:length(CamTrigger)]));
avgImg=mean(mov_mc,3); nTime=size(mov_mc,3);

[~, t_fit]=get_blueoffTrace([1:nTime],Blue>0,20,20);
t_fit=find(t_fit);
meanF = squeeze(mean(mov_mc,[1 2]));
bleaching_fit=expfitDM_2(t_fit',meanF(t_fit),[1:nTime]',[10000]);
mcTrace.xymean=mcTrace.xymean(1:size(mov_mc,3),:);
Result.ref_im=avgImg;

bkg = zeros(1, size(mov_mc,3));
bkg(1,:) = bleaching_fit;  % linear term
mov_res=mov_mc-mean(mov_mc,3);
mov_res = SeeResiduals(mov_res,mcTrace.xymean);
mov_res = SeeResiduals(mov_res,mcTrace.xymean.^2);
mov_res = SeeResiduals(mov_res,mcTrace.xymean(:,1).*mcTrace.xymean(:,end));
mov_res = SeeResiduals(mov_res,bkg,1);

mov_res_filt=mov_res;
mov_res_filt(dirtMov_dilate>0)=NaN;

Result.Blue=Blue; Result.DMDbluetrace=DMDbluetrace;
figure(10); clf;
show_footprnt_contour(blueDMDimg,mean(mov_mc,3))

figure(11); clf;
Result.blueDMDimg=blueDMDimg;
[Result.clicky_tr]=apply_clicky(Result.clicky_ROI,mov_res_filt);
Result.clicky_tr_raw=apply_clicky(Result.clicky_ROI,mov_mc);
Result.clicky_tr=Result.clicky_tr';
Result.clicky_tr_raw=Result.clicky_tr_raw';

DMDBlueSTAmat=[]; STAmov=[];
clicky_tr_filt=movmean(Result.clicky_tr,10,2);
% clicky_tr_filt=get_bandstop(clicky_tr_filt,1000,[50 75]);
% clicky_tr_filt=get_bandstop(clicky_tr_filt,1000,[9 11]);
% clicky_tr_filt(8500:9000)=NaN;
for b=1:max(DMDbluetrace)
[~, DMDBlueSTAmat{b}]=get_STA(clicky_tr_filt,[0 diff(DMDbluetrace==b)==1],nTau(1),nTau(2));
[STAmov{b}]=get_STA(tovec(mov_res_filt),[0 diff(DMDbluetrace==b)==1],nTau(1),nTau(2));
STAmov{b}=reshape(STAmov{b},sz(2),sz(1),[]);
DMDBlueSTAmat{b}=squeeze(DMDBlueSTAmat{b});
end
DMDBlueSTA=-cell2mat(cellfun(@(x) mean(x,1,'omitnan'),DMDBlueSTAmat,'UniformOutput',false)');
DMDBlueSTA=DMDBlueSTA-median(DMDBlueSTA(:,end-100:end),2);
Result.STAmov=STAmov;
Result.DMDBlueSTAmat=DMDBlueSTAmat;
Result.DMDBlueSTA=DMDBlueSTA;

save(fullfile(fpath{f},'OP_Result.mat'),"Result",'-v7.3')
disp(['Result save to ' fullfile(fpath{f},'OP_Result.mat')])
end

%% Concatenate data
GridBlueResult=[];
g=1; tic;
for f=foi
    load(fullfile(fpath{f},'OP_Result.mat'))
    GridBlueResult{g}=Result;
    GridBlueResult{g}.PixelSize=PixelSize(f);
    g=g+1;
end
toc;
%%
figure(19); clf;
for g=1:length(GridBlueResult)
    nexttile([1 1]);
    imagesc(GridBlueResult{g}.DMDBlueSTA)
    title(g)
end
%%
 time2show=[10:30];
 DecayDat=[]; OrangeDat=[];
for ind2show=[1:6 8]
clickyROI_center=mean(GridBlueResult{ind2show}.clicky_ROI{1}(1:end-1,:),1);
blueDMDpattcenter=get_coord(GridBlueResult{ind2show}.blueDMDimg);
distMat=distance_mat(blueDMDpattcenter,clickyROI_center).*GridBlueResult{ind2show}.PixelSize;
[~, c]=min(distMat); %distMat(1:c-1)=-distMat(1:c-1); 
ref_im=(GridBlueResult{ind2show}.ref_im-100)-medfilt2(GridBlueResult{ind2show}.ref_im-100,[80 80]);
[OrangeProfile kymoROI]=polyLineKymo_lineInput(ref_im,5,30,blueDMDpattcenter);
kymoROI_center=cell2mat(cellfun(@(x) mean(x(1:end-1,:),1),kymoROI,'UniformOutput',false)');
distMat_orange=distance_mat(kymoROI_center,clickyROI_center).*GridBlueResult{ind2show}.PixelSize;
%[~, c]=min(distMat_orange); distMat_orange(1:c-1)=-distMat_orange(1:c-1);
OrangeProfile=rescale(OrangeProfile);
dVpolyLine=mean(GridBlueResult{ind2show}.DMDBlueSTA(:,nTau(1)+time2show),2);
[maxdV, m]=max(dVpolyLine);
if c ~= m
    shift=c-m;
    if shift>0
    distMat=distMat(1+shift:end);
    dVpolyLine=dVpolyLine(1:end-shift);
    else
    distMat=distMat(1:end+shift);
    dVpolyLine=dVpolyLine(1-shift:end);    
    end
end
dVpolyLine=dVpolyLine./maxdV;
distMat=distMat(1:length(dVpolyLine)); 
distMat=distMat-min(distMat);
distMat_orange=distMat_orange-min(distMat_orange);
GridBlueResult{ind2show}.DMDBlueSTA=GridBlueResult{ind2show}.DMDBlueSTA-median(GridBlueResult{ind2show}.DMDBlueSTA(:,1:40),2);
DecayDat=[DecayDat; {[distMat dVpolyLine]}];
OrangeDat=[OrangeDat; {[distMat_orange OrangeProfile']}];
end

show_rep=5;
figure(21); clf; tiledlayout(4,2);
cmap_blue=winter(size(GridBlueResult{show_rep}.blueDMDimg,3));
blueDMDpattcenter=get_coord(GridBlueResult{show_rep}.blueDMDimg);
nexttile([1 2])
imshow2(GridBlueResult{show_rep}.ref_im,[]); hold all
plot(blueDMDpattcenter(:,1),blueDMDpattcenter(:,2),'r');
nexttile([1 2])
imshow2(GridBlueResult{show_rep}.ref_im,[]); hold all
plot(GridBlueResult{show_rep}.clicky_ROI{1}(:,1),GridBlueResult{show_rep}.clicky_ROI{1}(:,2),'r');
drawScaleBar(100/GridBlueResult{show_rep}.PixelSize,'Horizontal')
nexttile([1 2])
imshow2(GridBlueResult{show_rep}.ref_im,[]); hold all
for b=1:size(GridBlueResult{show_rep}.DMDBlueSTA,1)
    bd=bwboundaries(GridBlueResult{show_rep}.blueDMDimg(:,:,b)>0);
plot(bd{1}(:,2),bd{1}(:,1),'color',cmap_blue(b,:),'LineWidth',1.5)
end
nexttile([1 1]);
DMDBlueSTD=cell2mat(cellfun(@(x) std(x,0,1,'omitnan')./sqrt(size(x,1)),GridBlueResult{show_rep}.DMDBlueSTAmat,'UniformOutput',false)');
l=[];
for b=1:size(DMDBlueSTD,1)
l(b)=errorbar_shade([-nTau(1):nTau(2)],GridBlueResult{show_rep}.DMDBlueSTA(b,:),DMDBlueSTD(b,:),cmap_blue(b,:)); hold all;
%plot([-nTau(1):nTau(2)],GridBlueResult{show_rep}.DMDBlueSTA(b,:),'color',cmap_blue(b,:)); hold all;
end
box off;
xlabel('Time (ms)');
ylabel('\DeltaF'); 
legend(l,num2str(DecayDat{show_rep}(:,1),2))

nexttile([1 1]);
exp_decayconst=[];
for g=1:length(DecayDat)
    plot(DecayDat{g}(:,1),DecayDat{g}(:,2),'color',[0.7 0.7 0.7]); hold all
     [~, params, R2] = fitExpDecay_1param(DecayDat{g}(:,1), DecayDat{g}(:,2), [0:50]);
exp_decayconst(g,:)=[params(1) R2];
end
hold all
[M S Xc N]=binning_data(DecayDat,[0:3:40]);
errorbar(Xc,M,S./sqrt(sum(cellfun(@sum,N),2))','r','linewidth',1.5);
xlim([0 30]); ylim([0 1]);
box off;
xlabel('Stimulation distance from ROI (\mum)');
ylabel('Normalized \DeltaF');
% plot(DecayDat{show_rep}(:,1),DecayDat{show_rep}(:,2)); hold all
% plot(OrangeDat{show_rep}(:,1),OrangeDat{show_rep}(:,2),'r')
