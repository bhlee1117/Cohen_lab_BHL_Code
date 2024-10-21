clear
clc;
cd '/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/FromBackup/PP72_PlaceCellResults';
[~, ~, raw] = xlsread(['/Volumes/cohen_lab/Lab/Labmembers' ...
    '/Byung Hun Lee/Data/PrismPCdata_Arrangement.xlsx'], 'Sheet1', 'C5:P23');

% [~, ~, NeuronsToUse]=xlsread(['/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/' ...
%     'PlaceCellData_Arrangement.xlsx'], 'Sheet1', 'L8:M46');
%
% NeuronsToUse=cellfun(@(x) (str2num(num2str(x))),NeuronsToUse,'UniformOutput',false);
ref_ROI=cellfun(@(x) (str2num(num2str(x))),raw(:,10),'UniformOutput',false);
fpath=raw(:,1)';
StructureData=raw(:,9);
EndFrame=cell2mat(raw(:,13));
ifmotionReject=cell2mat(raw(:,14));
save_figto='/Volumes/BHL_WD18TB/PP72_PlaceCellResults';
place_bin=150; time_segment=15000; overlap=200;
alignedMovFN = {'STA_Mat_SS','STA_Mat_CS','STA_Mat_dSP'};
bound=6;
save_at='/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Updates/2024/20240408_Movs_Figs/';
%%
f=18; load(fullfile(fpath{f},'PC_Result.mat'),'Result')
basal_ROI=cellfun(@(x) (str2num(num2str(x))),raw(:,11),'UniformOutput',false);
apical_ROI=cellfun(@(x) (str2num(num2str(x))),raw(:,12),'UniformOutput',false);
rois={basal_ROI{f},apical_ROI{f}};

nROI=size(Result.normTraces,1);
nTau={[-130:20],[-130:100],[-30:20]}; %SS, CS, dSP
nTau_bAP=[-20:20];

som_spike=find(Result.spike(1,:));
tr_ref=Result.normTraces(ref_ROI{f},:);
tr_sub=mean(tr_ref,1)-movprc(mean(tr_ref,1),200,20,2);
tr_sub=get_subthreshold(tr_sub,Result.spike(1,:),5,10);
[trans tr_trace]=detect_transient2(tr_sub,[5 1.5],Result.spike(1,:),15);

SkelDend = Skeletonize_dendrite(Result.ref_im,4,0.02,25);
interDendDist=[];
for i=1:size(Result.normTraces,1)
    for j=1:size(Result.normTraces,1)
        [interDendDist(i,j), ~]=geodesic_distance(SkelDend,get_coord(Result.ftprnt(:,:,i)),get_coord(Result.ftprnt(:,:,j)));
    end
end
coord_1d=dim_reduce(get_coord(Result.ftprnt));
[~, dist_order]=sort(coord_1d,'ascend');
som_roi=find(dist_order==1);
geodist=interDendDist(1,:)'.*sign(coord_1d-coord_1d(1));
show_footprnt_contour(Result.ftprnt(:,:,dist_order),Result.ref_im)

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
SS_list = ismember(find(Result.SpClass(1,:)), bAP_s);

% Reference bAPs for F_ref
bAP_ref=[];
for s=som_spike
    isnearby=sum(ismember(s+nTau_bAP,som_spike))>1;
    isnearbyCS=sum(ismember(s+nTau_bAP,find(Result.CStrace)))>1;
    ispartCS=tr_trace(s)>0;
    if ~isnearby & ~isnan(s) & ~isnearbyCS & ~ispartCS
        bAP_ref=[bAP_ref s];
    end
end

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

dSP=[]; roi_dsp=[3 4 5 6 7];
for s=find(Result.SpClass(3,:))
    isnearby=sum(ismember(s+nTau{3},find(Result.SpClass(1,:))))>1;
    isnearbyCS=sum(ismember(s+nTau{3},find(Result.CStrace)))>1;
    ispartCS=tr_trace(s)>0;
    if ~isnearby & ~isnan(s) & ~isnearbyCS & ~ispartCS
        dSP=[dSP s];
    end
end
dSP_roi=dSP(ismember(dSP, find(max(Result.spike(roi_dsp,:)))));
dSP_list = ismember(find(Result.SpClass(3,:)), dSP_roi);


prc_normTr=Result.normTraces;
%prc_normTr=Result.normTraces-movprc(Result.normTraces,500,35,2);
STA_SSmat=reshape(prc_normTr(:,bAP_s'+nTau{1}),nROI,[],length(nTau{1}));
STA_SS=squeeze(mean(reshape(prc_normTr(:,bAP_ref'+nTau_bAP),nROI,[],length(nTau_bAP)),2));
STA_SS_abs=squeeze(mean(reshape(Result.traces(:,bAP_ref'+nTau_bAP),nROI,[],length(nTau_bAP)),2))./sum(tovec(Result.ftprnt))';
F_ROI=tovec(Result.ftprnt)'*tovec(Result.ref_im);
F_ref=mean(STA_SS(:,-nTau_bAP(1)+[10:14]),2);
STA_SSmat=STA_SSmat./F_ref;
prc_normTrCS=Result.normTraces;
%prc_normTrCS=Result.normTraces-movprc(Result.normTraces,500,35,2);
STA_CSmat=reshape(prc_normTrCS(:,CS_s'+nTau{2}),nROI,[],length(nTau{2}));
STA_CS_abs=squeeze(mean(reshape(Result.traces(:,CS_s'+nTau_bAP),nROI,[],length(nTau_bAP)),2))./sum(tovec(Result.ftprnt))';
STA_CSmat=STA_CSmat./F_ref;
CSpike=Result.spike.*Result.CStrace;
STA_CSpikemat=reshape(CSpike(:,CS_s'+nTau{2}),nROI,[],length(nTau{2}));
STA_dSPmat=reshape(prc_normTrCS(:,dSP_roi'+nTau{3}),nROI,[],length(nTau{3}));
STA_dSPmat=STA_dSPmat./F_ref;
%%
f=18; load(fullfile(fpath{f},'PC_Result.mat'),'Result')
StructureStack=(double(tiffreadVolume(StructureData{f})));
StructureStack(StructureStack<100)=NaN;
StructureStack=StructureStack-100;
illumination_field=imgaussfilt_NaN(max(StructureStack,[],3),150);
StructureStack_norm=max(StructureStack,[],3)./illumination_field;
StructureStack_norm(StructureStack_norm<0.7)=prctile(StructureStack_norm(:),10);

figure(3); clf;
imshow2(max(StructureStack_norm,[],3),[])
g=1; ROIrmv=[];
while g
    h = drawpolygon('Color','r');
    if size(h.Position,1)==1 %no more ROI
        g=0;
    else
        ROIrmv=[ROIrmv; {h.Position}];
        hold all
        plot(h.Position(:,1),h.Position(:,2))
    end
end
ROIrmvmask=roi2mask(ROIrmv,size(StructureStack_norm,1),size(StructureStack_norm,2));
close(figure(3));
StructureStack_norm(ROIrmvmask)=0;

[RegImg,tformReg]=imReg(Result.ref_im,StructureStack_norm);
Result.tform=tformReg;
Result.Struckmax=StructureStack_norm;
save(fullfile(fpath{f},'PC_Result.mat'),"Result",'-v7.3')

%% Representative CS kymos
figure(27); clf; cax=[-2 6];
show_cs=[4 5 6 9 17 44 45 49 33 71 51 53 79 92 100 77 76 104 108 137 148];
noi=[setdiff(1:nROI,[2 4])];
for i=show_cs
    nexttile([1 1])
    imagesc(squeeze(STA_CSmat(dist_order(noi),i,50:250)),cax)
    title(num2str(i))
end
figure(28); clf;
for i=[2:2:42]
    nexttile([1 1])
    imagesc(squeeze(STA_SSmat(dist_order(noi),i,:)),cax)
    title(num2str(i))
end

figure(29); clf; g=1;
for i=[6 7 8 9 12 14 15 28 29 22 31 34 42 41 48 35 36 63 64 99]
    nexttile([1 1])
    imagesc(nTau{3},[1:length(noi)],squeeze(STA_dSPmat(dist_order(noi),i,:)),cax)
    title(num2str(g))
    g=g+1;
end
%%
figure(3); clf; ax1=[]; ax2=[];
cmap=distinguishable_colors(6);
tiledlayout(2,3)
ax1=[ax1 nexttile([1 1])];
imagesc(nTau{1},[1:1:nROI],squeeze(mean(STA_SSmat(dist_order,:,:),2)),[-1 3])
xlabel('Peri-spike time (ms)')
ylabel('ROI ordered along basal-distal axis')
ax1=[ax1 nexttile([1 1])];
imagesc(nTau{2},[1:1:nROI],squeeze(mean(STA_CSmat(dist_order,:,:),2)),[-1 3])
xlabel('Peri-spike time (ms)')
ylabel('ROI ordered along basal-distal axis')
ax1=[ax1 nexttile([1 1])];
imagesc(nTau{3},[1:1:nROI],squeeze(mean(STA_dSPmat(dist_order,:,:),2)),[-1 3])
%set(gca,'ytick',[1:2:nROI],'YTickLabel',geodist(dist_order([1:2:nROI])))
xlabel('Peri-spike time (ms)')
ylabel('ROI ordered along basal-distal axis')
colormap(turbo)

ax2=[ax2 nexttile([1 1])];
h1=errorbar_shade(nTau{1},squeeze(mean(mean(STA_SSmat(rois{1},:,:),1),2))',squeeze(std(mean(STA_SSmat(rois{1},:,:),1),0,2))',cmap(1,:)); hold all
h2=errorbar_shade(nTau{1},squeeze(mean(mean(STA_SSmat(rois{2},:,:),1),2))',squeeze(std(mean(STA_SSmat(rois{2},:,:),1),0,2))',cmap(2,:));
axis tight
xlabel('Peri-spike time (ms)')
legend([h1 h2],{'Basal','Apical'})
ax2=[ax2 nexttile([1 1])];
h1=errorbar_shade(nTau{2},squeeze(mean(mean(STA_CSmat(rois{1},:,:),1),2))',squeeze(std(mean(STA_CSmat(rois{1},:,:),1),0,2))',cmap(1,:)); hold all
h2=errorbar_shade(nTau{2},squeeze(mean(mean(STA_CSmat(rois{2},:,:),1),2))',squeeze(std(mean(STA_CSmat(rois{2},:,:),1),0,2))',cmap(2,:));
legend([h1 h2],{'Basal','Apical'})
axis tight
xlabel('Peri-spike time (ms)')
ax2=[ax2 nexttile([1 1])];
h1= errorbar_shade(nTau{3},squeeze(mean(mean(STA_dSPmat(rois{1},:,:),1),2))',squeeze(std(mean(STA_dSPmat(rois{1},:,:),1),0,2))',cmap(1,:)); hold all
h2= errorbar_shade(nTau{3},squeeze(mean(mean(STA_dSPmat(rois{2},:,:),1),2))',squeeze(std(mean(STA_dSPmat(rois{2},:,:),1),0,2))',cmap(2,:));
legend([h1 h2],{'Basal','Apical'})
axis tight
xlabel('Peri-spike time (ms)')
%linkaxes(ax1,'x')
linkaxes(ax2,'y')

%%
s_list={SS_list,CS_list,dSP_list};
for spclass_ind=1:3
    %load aligned movie somatic spike
    alignmovlist=dir(fullfile(fpath{f},[alignedMovFN{spclass_ind} '*.tiff']));
    AlignMov=[];
    for l=1:length(alignmovlist)
        l
        AlignMov=cat(3,AlignMov,readtiff(fullfile(fpath{f},alignmovlist(l).name)));
    end
    sz_align=size(AlignMov)
    AlignMov=double(reshape(AlignMov,sz_align(1),sz_align(2),length(nTau{spclass_ind}),[]));
    %AlignMov=AlignMov-median(AlignMov(:,:,:,:),3);
    %AlignMov=AlignMov-mean(AlignMov(:,:,1:6,:),3);
    STAmovie{spclass_ind}=mat2gray(mean(AlignMov(:,:,:,s_list{spclass_ind}),4))*range(STA_SS_abs(1,:));
    STAmovie{spclass_ind}=-STAmovie{spclass_ind};
    STAmovie{spclass_ind}=STAmovie{spclass_ind}-prctile(STAmovie{spclass_ind},30,3);
end
%%
bound=6; n_pc=[12 15 5];
for i=1:3
    STAmovie{i}=pcafilt(STAmovie{i},n_pc(i));
end
F_refImg=mean(STAmovie{1}(:,:,-nTau{1}(1)+[6:10]),3); F_refImg(F_refImg<0.1)=prctile(F_refImg(:),30);
F_refImg=imgaussfilt(F_refImg,3);
%F_refImg=imgaussfilt(Result.ref_im(bound:end-bound,bound:end-bound),3);
Rfixed = imref2d([size(Result.ref_im,1) size(Result.ref_im,2)]);
inverseTform = invert(Result.tform);
revertedStruct = imwarp(Result.Struckmax, inverseTform,'OutputView',Rfixed);
%revertedStruct = imwarp(Result.Struckmax, inverseTform,'OutputView',Rfixed);
revertedStruct(revertedStruct==0)=prctile(revertedStruct(:),30);
revertedStruct=mat2gray(revertedStruct);
revertedStruct=revertedStruct(bound:end-bound,bound:end-bound);
revertedStruct_filt=revertedStruct-imgaussfilt(revertedStruct,15);
revertedStruct(revertedStruct==1)=0; revertedStruct(revertedStruct<0.05)=0;
SkelDend = Skeletonize_dendrite(Result.ref_im,10,0.01,25); SkelDend=SkelDend(bound:end-bound,bound:end-bound);
for spclass_ind=1:3
    STAnorm_sub=STAmovie{spclass_ind};%-mean(STAmovie{spclass_ind}(:,:,1:5),3);
    STAmovie_norm{spclass_ind}=imgaussfilt3(STAnorm_sub./F_refImg,[1 1 0.1]);%.*SkelDend(bound:end-bound,bound:end-bound);

    colorSTA=grs2rgb(tovec(STAmovie_norm{spclass_ind}),colormap('jet'),-1,5);
    colorSTA=reshape(colorSTA,size(STAmovie_norm{spclass_ind},1),size(STAmovie_norm{spclass_ind},2),size(STAmovie_norm{spclass_ind},3),[]);
    colorSTA=permute(colorSTA,[1 2 4 3]);

    colorSTA2=grs2rgb(tovec(STAmovie_norm{spclass_ind}.*revertedStruct),colormap('jet'),-0.5,2);
    colorSTA2=reshape(colorSTA2,size(STAmovie_norm{spclass_ind},1),size(STAmovie_norm{spclass_ind},2),size(STAmovie_norm{spclass_ind},3),[]);
    colorSTA2=permute(colorSTA2,[1 2 4 3]);
    %STAmovie_normStr{spclass_ind}=colorSTA.*revertedStruct*5;%.*SkelDend(bound:end-bound,bound:end-bound);
    STAmovie_normStr{spclass_ind}=colorSTA2.*revertedStruct*5;%.*SkelDend(bound:end-bound,bound:end-bound);
end

sptype={'SS','CS','dSP'};
cax=[-0.1 0.3]*30; crop_roi=[21.5 57.5 328 128]; toi={[1:51],[71:251],[1:51]};
for sp_class=1:3;
    %mov_show=STAmovie_norm{sp_class};
    mov_show=imrotate(STAmovie_normStr{sp_class},-14);
    mov_show=mov_show(56:180,22:350,:,toi{sp_class});
    writeMov4d(['STA_dFFStruct_grsrgbMask' sptype{sp_class}],mov_show,[nTau{sp_class}(toi{sp_class})],10,1,cax)
end
%%
figure(24); clf; ax1=[]; crop_roi=[21.5000   57.5000  328.0000  128.0000];
sptype={'Simple Spike','Complex Spike','dSpike'};
for i=1:3
    ax1=[ax1 nexttile([1 1])];
    stamov_sub=STAmovie_normStr{i};
    mean_subth=mean(stamov_sub(:,:,:,-nTau{i}(1)+[-3:-1]),4);
    mean_subth=imgaussfilt(mean_subth,1);
    mean_subth=imcrop(imrotate(mean_subth,-14,"bilinear"),crop_roi);
    imagesc(mean_subth,[-0.1 0.3]*10)
    axis equal tight
    title(sptype{i})
end
colormap(jet)
linkaxes(ax1,'xy')

%%
figure(24); clf; ax1=[]; cax=[25 25 25];
load('revertStructure_rmvMask.mat')
sptype={'Simple Spike','Complex Spike','dSpike'};
revertedStruct_th=revertedStruct;
revertedStruct_th(revertedStruct_th<0.07)=0;
revertedStruct_th(ROIrmvmask)=0;
revertedStruct_th=revertedStruct_th - medfilt2(revertedStruct_th,[30 30]); revertedStruct_th(revertedStruct_th<0)=0;
for i=1:2
    ax1=[ax1 nexttile([1 1])];
    stamov_sub=STAmovie{i};
    peak_amp=max(stamov_sub(:,:,-nTau{i}(1):end),[],3).*double(revertedStruct>0.06);
    peak_amp=imgaussfilt(peak_amp,1);
    peakcolor=grs2rgb(peak_amp,colormap('jet'),-3,cax(i));
    peakAmpMap=peakcolor.*revertedStruct_th*10;
    peakAmpMap=imcrop(imrotate(peakAmpMap,-14,"bilinear"),crop_roi);
    imagesc(peakAmpMap)
    axis equal tight
    title(sptype{i})
end
ax1=[ax1 nexttile([1 1])];
stamov_sub=STAmovie{3};
    peak_amp=mean(stamov_sub(:,:,-nTau{3}(1)+[0:2]),3).*double(revertedStruct>0.06);
    peak_amp=imgaussfilt(peak_amp,1);
    peakcolor=grs2rgb(peak_amp,colormap('jet'),-3,cax(3));
    peakAmpMap=peakcolor.*revertedStruct_th*10;
    peakAmpMap=imcrop(imrotate(peakAmpMap,-14,"bilinear"),crop_roi);
    imagesc(peakAmpMap)
    axis equal tight
    title(sptype{3})
colormap(jet)
linkaxes(ax1,'xy')
%%
figure(21); clf;
cmap=distinguishable_colors(6);
markers = {'o', 's', 'd', '^', 'v', '>', '<', 'p', 'h', '+', '*'}; m=1;
for f=[1 4 11 14 18]
    load(fullfile(fpath{f},'PC_Result.mat'),'Result')
    run_time=Result.VR(end,:)>0.005;
    for i=1:3
        h{i}=plot([1 2]-0.1+0.05*i,1000*[mean(Result.SpClass(i,run_time)) mean(Result.SpClass(i,~run_time))],markers{m},'color',cmap(i,:),'markersize',10); hold all
    end
    m=m+1;
end
legend([h{1} h{2} h{3}],{'SS','CS','dSP'})
xlim([0.5 2.5])
set(gca,'xtick',[1 2],'xticklabel',{'Walking', 'Resting'})
ylabel('Number of spikes per second')

%%
clear
clc;
[~, ~, raw] = xlsread(['/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/' ...
    'Prism_OptopatchData_Arrangement.xlsx'], 'Sheet1', 'B5:K96');

save_to='/Volumes/BHL18TB_D1/Arranged_Data/Prism_OptopatchResult';
fpath=raw(:,1);
Mouse=cell2mat(raw(:,2));
NeuronInd=cell2mat(raw(:,5));
CamType=raw(:,3);
StructureData=raw(:,10);

%%
nTau2=[-40:120];

for i=[75]
    cd(fpath{i})
    load('Result.mat')
    Result.normTrace=Result.traces./get_threshold(Result.traces,1);
    Result.spike=find_spike_bh(Result.normTrace-movprc(Result.normTrace,200,30,2),5,4);

    Blue=Result.Blue;
    blueOff = Blue == 0;
    blueOff2 = imerode(blueOff, [ones(1,20), zeros(1, 20)]);
    Blue_di=~blueOff2;
    bwBlue_di=bwlabel(Blue_di);

    ref_ROI=find(sum(Result.spike(1:5,:).*bwBlue_di,2)==max(sum(Result.spike(1:5,:).*bwBlue_di,2)),1);
    nROI=size(Result.normTrace,1);
    tr=Result.normTrace(ref_ROI,:); t=[1:length(tr)];
    spike=Result.spike(ref_ROI,:);

    sp_pulse=[];
    for b=1:max(bwBlue_di)
        t_tmp=find(bwBlue_di==b);
        if max(bwBlue_di)<50
            sp_pulse=[sp_pulse find(spike(t_tmp))+t_tmp(1)-1];
        else
            sp_pulse=[sp_pulse find(spike(t_tmp),1,'first')+t_tmp(1)-1];
        end
    end

    sTau=sp_pulse(1:end-1)'+nTau2;
    spikeMat=reshape(Result.normTrace(:,sTau),nROI,length(sp_pulse)-1,[]);
    STAtrace=squeeze(mean(spikeMat,2));
    STAtrace=STAtrace-prctile(STAtrace, 5, 2);
    Result.spikeMat = spikeMat;
    Result.STAtrace = STAtrace;

    load(fullfile(fpath{i},"output_data.mat"))
    load([fpath{i} '/mcTrace' num2str(1,'%02d') '.mat']);
    switch char(CamType(i))
        case 'flash'
            sz=double(Device_Data{1, 4}.ROI([2 4]));
        case 'fusion'
            sz=double(Device_Data{1, 3}.ROI([2 4]));
    end

    CamCounter=Device_Data{1, 2}.Counter_Inputs(1, 1).data;
    CamTrigger=find(CamCounter(2:end)-CamCounter(1:end-1));

    mov_mc=double(readBinMov_times([fpath{i} '/mc_ShutterReg' num2str(1,'%02d') '.bin'],sz(2),sz(1),[1:length(CamTrigger)]));
    mov_res= mov_mc-mean(mov_mc,3);
    bkg = zeros(1, size(mov_mc,3));
    bkg(1,:)=movmedian(get_blueoffTrace(squeeze(mean(mov_mc,[1 2])),Blue,30),3000,'omitnan');
    mov_res = SeeResiduals(mov_res,Result.mc);
    mov_res = SeeResiduals(mov_res,Result.mc.^2);
    mov_res = SeeResiduals(mov_res,Result.mc(:,1).*Result.mc(:,end));
    mov_res= SeeResiduals(mov_res,bkg,1);

    STA_tmp=reshape(mov_res(:,:,sTau),sz(2),sz(1),[],length(nTau2));
    Result.STAmovie=squeeze(mean(STA_tmp,3));
end
moviefixsc(Result.STAmovie)
save(['Result'],"Result",'-v7.3')
%%
rois={[16 20],[5:7]};
somspfpath=73; ddspfpath=75;
load(fullfile(fpath{somspfpath},'Result.mat'))
STAmov_som=Result.STAmovie;
STAtr_som=Result.STAtrace;
F_ref=(mean(STAtr_som(:,-nTau2(1)+[5:7]),2));
STAtr_som=STAtr_som./F_ref;
STA_somSpMat=Result.spikeMat./F_ref;
SomBlue=Result.BlueDMDimg;
load(fullfile(fpath{ddspfpath},'Result.mat'))
STAmov_dd=Result.STAmovie;
STAtr_dd=Result.STAtrace;
STAtr_dd=STAtr_dd./F_ref;
STA_ddSpMat=Result.spikeMat./F_ref;
ddBlue=Result.BlueDMDimg;
nROI=size(Result.ftprnt,3);
coord_1d=dim_reduce(get_coord(Result.ftprnt));
[~, dist_order]=sort(coord_1d,'descend');

%%
figure(36); clf; ax2=[]; ax3=[]; noi=setdiff([1:nROI],24);
cmap=distinguishable_colors(6);
tiledlayout(4,2)
ax3=[ax3 nexttile([1 1])];
blue_bd=cell2mat(bwboundaries(SomBlue));
imshow2(Result.ref_im,[]); hold all
plot(blue_bd(:,2),blue_bd(:,1),'color',[0 0.5 1])
ax3=[ax3 nexttile([1 1])];
blue_bd=cell2mat(bwboundaries(ddBlue));
imshow2(Result.ref_im,[]); hold all
plot(blue_bd(:,2),blue_bd(:,1),'color',cmap(1,:))

ax3=[ax3 nexttile([1 1])];
roi_bd=cell2mat(bwboundaries(max(Result.ftprnt(:,:,rois{1}),[],3)));
imshow2(Result.ref_im,[]); hold all
plot(roi_bd(:,2),roi_bd(:,1),'color',[0 0.5 1])
ax3=[ax3 nexttile([1 1])];
roi_bd=cell2mat(bwboundaries(max(Result.ftprnt(:,:,rois{2}),[],3)));
imshow2(Result.ref_im,[]); hold all
plot(roi_bd(:,2),roi_bd(:,1),'color',cmap(2,:))

ax1=[nexttile([1 1])];
imagesc(nTau2,[1:1:nROI],STAtr_som(dist_order(noi),:,:),[-0.5 2])
xlabel('Peri-spike time (ms)')
ylabel('ROI ordered along basal-distal axis')
ax4=[nexttile([1 1])];
imagesc(nTau2,[1:1:nROI],STAtr_dd(dist_order(noi),:,:),[-0.5 2])
xlabel('Peri-spike time (ms)')
ylabel('ROI ordered along basal-distal axis')
colormap(ax1,turbo)
colormap(ax4,turbo)

ax2=[ax2 nexttile([1 1])];
h1=errorbar_shade(nTau2,squeeze(mean(mean(STA_somSpMat(rois{1},:,:),1),2))',squeeze(std(mean(STA_somSpMat(rois{1},:,:),1),0,2))',cmap(1,:)); hold all
h2=errorbar_shade(nTau2,squeeze(mean(mean(STA_somSpMat(rois{2},:,:),1),2))',squeeze(std(mean(STA_somSpMat(rois{2},:,:),1),0,2))',cmap(2,:));
axis tight
xlabel('Peri-spike time (ms)')
legend([h1 h2],{'Basal','Distal'})
ax2=[ax2 nexttile([1 1])];
h1=errorbar_shade(nTau2,squeeze(mean(mean(STA_ddSpMat(rois{1},:,:),1),2))',squeeze(std(mean(STA_ddSpMat(rois{1},:,:),1),0,2))',cmap(1,:)); hold all
h2=errorbar_shade(nTau2,squeeze(mean(mean(STA_ddSpMat(rois{2},:,:),1),2))',squeeze(std(mean(STA_ddSpMat(rois{2},:,:),1),0,2))',cmap(2,:));
legend([h1 h2],{'Basal','Distal'})
axis tight
xlabel('Peri-spike time (ms)')
linkaxes(ax2,'y')

%%
bound=6;
i=73; load(fullfile(fpath{i},'Result.mat'))
STAmovie{1}=-Result.STAmovie;
STAmovie{1}=STAmovie{1}-mean(STAmovie{1}(:,:,[1:10]),3);
STAmovie{1}=pcafilt(STAmovie{1},7);
i=75; load(fullfile(fpath{i},'Result.mat'))
STAmovie{2}=-Result.STAmovie;
STAmovie{2}=STAmovie{2}-mean(STAmovie{2}(:,:,[1:10]),3);
STAmovie{2}=pcafilt(STAmovie{2},7);

F_refImg=mean(STAmovie{1}(:,:,-nTau2(1)+[5:10]),3);
F_refImg=imgaussfilt(F_refImg,3);
%F_refImg=Result.ref_im;
Rfixed = imref2d([size(Result.ref_im,1) size(Result.ref_im,2)]);
inverseTform = invert(Result.tform);
revertedStruct = imwarp(Result.Structure, inverseTform,'OutputView',Rfixed);
%revertedStruct = imwarp(Result.Struckmax, inverseTform,'OutputView',Rfixed);
revertedStruct(revertedStruct==0)=prctile(revertedStruct(:),30);
revertedStruct=mat2gray(revertedStruct);
%revertedStruct=revertedStruct(bound:end-bound,bound:end-bound);
revertedStruct_filt=revertedStruct-imgaussfilt(revertedStruct,15);
%SkelDend = Skeletonize_dendrite(Result.ref_im,10,0.01,25);
for spclass_ind=1:2
    STAnorm_sub=STAmovie{spclass_ind};
    STAmovie_norm{spclass_ind}=imgaussfilt3(STAnorm_sub./F_refImg,[1.5 1.5 0.1]);%.*SkelDend(bound:end-bound,bound:end-bound);
    colorSTA2=grs2rgb(tovec(STAmovie_norm{spclass_ind}),colormap('jet'),-0.5,1.8);
    colorSTA2=reshape(colorSTA2,size(STAmovie_norm{spclass_ind},1),size(STAmovie_norm{spclass_ind},2),size(STAmovie_norm{spclass_ind},3),[]);
    colorSTA2=permute(colorSTA2,[1 2 4 3]);
    STAmovie_normStr{spclass_ind}=colorSTA2.*revertedStruct*3;%.*SkelDend(bound:end-bound,bound:end-bound);
end

sptype={'SomStim','ddStim'};
cax=[-0.1 0.3]*30; crop_roi=[94.5100    6.5100  339.9800  158.9800];
for sp_class=1:2;
    %mov_show=STAmovie_norm{sp_class};
    mov_show=imrotate(STAmovie_normStr{sp_class}(6:165,94:400,:,:),90);
    %writeMov4d(['STA_dFFStructgrsrgb' sptype{sp_class}],mov_show,[nTau2],10,1,cax)
end

%% Subthreshold figure;
figure(21); clf;

nexttile([1 1]);
imshow2(mean(STAmovie_normStr{1}(:,:,:,35:40),4),[]); hold all

blue_bd=cell2mat(bwboundaries(SomBlue));
plot(blue_bd(:,2),blue_bd(:,1),'color',[0 0.5 1])
nexttile([1 1]);
imshow2(mean(STAmovie_normStr{2}(:,:,:,35:40),4),[]); hold all
blue_bd=cell2mat(bwboundaries(ddBlue));
plot(blue_bd(:,2),blue_bd(:,1),'color',cmap(1,:))



%% Generate SNAPT movie
for i=[73]
    load(fullfile(fpath{i},'Result.mat'))
    mask=max(Result.Structure_bin,[],3)>0.01;
    maskSTA=max(-Result.STAmovie,[],3)./Result.ref_im>0.05;
    StrImg=max(Result.Structure,[],3);
    STAmovie=mat2gray(-Result.STAmovie);
    STAmovie=STAmovie-prctile(STAmovie,10,3);
    STAmovie=mat2gray(STAmovie(:,:,-nTau2(1)-5:-nTau2(1)+10));
    tformReg=Result.tform;
    [Result.SNAPT Result.dtimg]=generate_SNAPTmov(mat2gray(STAmovie),mask,StrImg,tformReg);
    revSNAPT = imwarp(Result.SNAPT, inverseTform,'OutputView',Rfixed);
    [yR xR zR]=size(Result.Structure);
    bluePatt = bwboundaries(Result.BlueDMDimg);
    figure(20); clf;
    v = VideoWriter([fpath{i} '/SNAPT_movie'],'MPEG-4');
    %v = VideoWriter([fpath{i} '/SNAPT_movie'],'Uncompressed AVI');

    open(v);
    subframeT = 0.025; % ms
    initialT = -2; % ms
    finalT = 2; % ms
    times = initialT:subframeT:finalT;

    for j = 1:length(times)
        clf;
        set(gca,'units','pixels','position',[50 50 500 400])
        imshow(revSNAPT(:,50:end-50,:,j),[])
        pbaspect([size(double(revSNAPT(:,50:end-50,:,j)),2) size(double(revSNAPT(:,:,:,j)),1) 1]),colormap(gray)
        hold all
        plot(bluePatt{1}(:,2)-50,bluePatt{1}(:,1),'color',[0 0.6 1],'linewidth',2)
        axis off
        text(2,20,[num2str(times(j)+0.95) 'ms'], 'FontSize', 20, 'color', [0.99 0.99 0.99])% the value 1. is to adjust timing by eyes
        pause(0.1)
        set(gcf,'color','w')    % Sets background to white
        frame = getframe(gcf);
        writeVideo(v,frame);
        pause(0.1);
    end;
    close(v);
end
save(fullfile(fpath{i},'Result.mat'),'Result')
%%
i=73;
load(fullfile(fpath{i},'Result.mat'))

nROI=size(Result.normTrace,1);
figure(25); clf;
nexttile([1 1])
dtimg=Result.dtimg;
dtimg(isnan(dtimg))=prctile(dtimg(:),30);
dtimg=imgaussfilt(dtimg-prctile(dtimg(:),0.001),5);
dtimg=grs2rgb((double(dtimg)),colormap("jet"));
delayImg=dtimg.*revertedStruct*5;
imagesc(delayImg(15:end-15,70:end-100,:))
axis equal tight
title('bAP delay')
colorbarHandle = colorbar;
colorbarHandle.Ticks = [0 0.5 1]; % Set the tick positions
colorbarHandle.TickLabels = {'0', '1.5 ms', '3 ms'}; % Set the tick labels

nexttile([1 1])
dff=imgaussfilt(STAmovie_norm{1}(15:end-15,70:end-100,-nTau2(1)+1),5);
dffImg=grs2rgb(dff,colormap('jet'));
dffImg=dffImg.*revertedStruct(15:end-15,70:end-100)*5;
imagesc(dffImg)
title('Spike amplitude')
axis equal tight
colormap(jet)

i=73;
load(fullfile(fpath{i},'Result.mat'))
nexttile([1 1])
coord_1d=dim_reduce(get_coord(Result.ftprnt));
STAtr=Result.STAtrace-mean(Result.STAtrace(:,1:15),2);
F_ref=mean(STAtr(:,-nTau2(1)+[3:8]),2);
STAtr=STAtr./F_ref;
[~, dist_order]=sort(coord_1d,'descend');
l=plot(nTau2,STAtr(dist_order,:)');
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(turbo(nROI),2))
xlabel('Peri-spike time (ms)')

%%


%%
i=46; nTau3=[-20:20];
load(fullfile(fpath{i},'Result.mat'))
nROI=size(Result.traces,1);
tr=Result.traces-movprc(Result.traces,50,25,2);
tr=tr./get_threshold(tr,1);
Result.normTrace=Result.traces./get_threshold(Result.traces,1);
sp=find_spike_bh(tr,4,3);
coord_1d=dim_reduce(get_coord(Result.ftprnt));
[~, dist_order]=sort(coord_1d,'descend');
dist_order=dist_order([1:4 6 5 7:nROI]);
sp_list=find(sp(2,:)); rmv_ind=find(sum((sp_list'+nTau3)<1,2)>1 | sum((sp_list'+nTau3)>size(Result.traces,2),2)>1);
sp_list(rmv_ind)=[];
STA_SS=squeeze(mean(reshape(Result.normTrace(:,sp_list'+nTau3),nROI,[],length(nTau3)),2));
STA_SS=STA_SS-mean(STA_SS(:,1:5),2);
F_ref3=mean(STA_SS(:,-nTau3(1)+[5:8]),2);
dFFtrace=Result.normTrace./F_ref3;
%dFFtrace=dFFtrace-movprc(dFFtrace,500,20,2);
show_traces_spikes(dFFtrace(dist_order,:),sp(dist_order,:),Result.Blue);
figure(27); clf;
show_footprnt_contour(Result.ftprnt(:,:,dist_order),Result.ref_im)
figure(26); clf;
l=plot(STA_SS(dist_order,:)'./F_ref3(dist_order)');
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(turbo(nROI),2))
xlabel('Peri-spike time (ms)')

%%

i=46; nTau3=[-20:20];
load(fullfile(fpath{i},'Result.mat'))
nROI=size(Result.traces,1);
tr=Result.traces-movprc(Result.traces,50,25,2);
tr=tr./get_threshold(tr,1);
Result.normTrace=Result.traces./get_threshold(Result.traces,1);
sp=find_spike_bh(tr,4,3);
coord_1d=dim_reduce(get_coord(Result.ftprnt));
[~, dist_order]=sort(coord_1d,'descend');
dist_order=dist_order([1:4 6 5 7:nROI]);
sp_list=find(sp(2,:)); rmv_ind=find(sum((sp_list'+nTau3)<1,2)>1 | sum((sp_list'+nTau3)>size(Result.traces,2),2)>1);
sp_list(rmv_ind)=[];
STA_SS=squeeze(mean(reshape(Result.normTrace(:,sp_list'+nTau3),nROI,[],length(nTau3)),2));
STA_SS=STA_SS-mean(STA_SS(:,1:5),2);
F_ref3=mean(STA_SS(:,-nTau3(1)+[5:8]),2);
dFFtrace=Result.normTrace./F_ref3;
%dFFtrace=dFFtrace-movprc(dFFtrace,500,20,2);
show_traces_spikes(dFFtrace(dist_order,:),sp(dist_order,:),Result.Blue);
figure(27); clf;
show_footprnt_contour(Result.ftprnt(:,:,dist_order),Result.ref_im)
figure(26); clf;
l=plot(STA_SS(dist_order,:)'./F_ref3(dist_order)');
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(turbo(nROI),2))
xlabel('Peri-spike time (ms)')

%%

i=86; nTau3=[-20:20];
load(fullfile(fpath{i},'Result.mat'))
nROI=size(Result.traces,1);
tr=Result.traces-movprc(Result.traces,50,25,2);
tr=tr./get_threshold(tr,1);
Result.normTrace=Result.traces./get_threshold(Result.traces,1);
sp=find_spike_bh(tr,4,3);
coord_1d=dim_reduce(get_coord(Result.ftprnt));
[~, dist_order]=sort(coord_1d,'ascend');
%dist_order=dist_order([1:4 6 5 7:nROI]);
sp_list=find(sp(2,:)); rmv_ind=find(sum((sp_list'+nTau3)<1,2)>1 | sum((sp_list'+nTau3)>size(Result.traces,2),2)>1);
sp_list(rmv_ind)=[];
STA_SS=squeeze(mean(reshape(Result.normTrace(:,sp_list'+nTau3),nROI,[],length(nTau3)),2));
STA_SS=STA_SS-mean(STA_SS(:,1:5),2);
F_ref3=mean(STA_SS(:,-nTau3(1)+[3:6]),2);
dFFtrace=Result.normTrace./F_ref3;
%dFFtrace=dFFtrace-movprc(dFFtrace,500,20,2);
show_traces_spikes(dFFtrace(dist_order,:),sp(dist_order,:),Result.Blue);
figure(27); clf;
show_footprnt_contour(Result.ftprnt(:,:,dist_order),Result.ref_im)
figure(26); clf;
l=plot(STA_SS(dist_order,:)'./F_ref3(dist_order)');
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(turbo(nROI),2))
xlabel('Peri-spike time (ms)')

%%

i=73; nTau3=[-20:20];
load(fullfile(fpath{i},'Result.mat'))
nROI=size(Result.traces,1);
tr=Result.traces-movprc(Result.traces,50,25,2);
tr=tr./get_threshold(tr,1);
Result.normTrace=Result.traces./get_threshold(Result.traces,1);
sp=find_spike_bh(tr,4,3);
coord_1d=dim_reduce(get_coord(Result.ftprnt));
[~, dist_order]=sort(coord_1d,'descend');
%dist_order=dist_order([1:4 6 5 7:nROI]);
sp_list=find(sp(2,:)); rmv_ind=find(sum((sp_list'+nTau3)<1,2)>1 | sum((sp_list'+nTau3)>size(Result.traces,2),2)>1);
sp_list(rmv_ind)=[];
STA_SS=squeeze(mean(reshape(Result.normTrace(:,sp_list'+nTau3),nROI,[],length(nTau3)),2));
STA_SS=STA_SS-mean(STA_SS(:,1:10),2);
F_ref3=mean(STA_SS(:,-nTau3(1)+[3:6]),2);
dFFtrace=Result.normTrace./F_ref3;
%dFFtrace=dFFtrace-movprc(dFFtrace,500,20,2);
show_traces_spikes(dFFtrace(dist_order,:),sp(dist_order,:),Result.Blue);
figure(27); clf;
show_footprnt_contour(Result.ftprnt(:,:,dist_order),Result.ref_im)
figure(26); clf;
l=plot(nTau3,STA_SS(dist_order,:)'./F_ref3(dist_order)');
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(turbo(nROI),2))
xlabel('Peri-spike time (ms)')

%%
SomRP_path=[15 22 54 72 86 92 4];
DenRP_path=[13 23 55 74 88 91 3];
catpath=[SomRP_path DenRP_path];
for f=1:length(catpath)
    i=catpath(f);
    load(fullfile(fpath{i},'Result'))
    coord_1d=dim_reduce(get_coord(Result.ftprnt));
    coord_1d=abs(coord_1d-coord_1d(1));
    Result.normTraces=Result.traces-prctile(Result.traces,25,2);
    Result.normTraces=Result.normTraces./get_threshold(Result.normTraces,1);
    ref_ind=find(coord_1d<80);
    tr_ref=mean(Result.normTraces(ref_ind,:),1);
    tr=Result.normTraces;
    nROI=size(Result.normTraces,1);
    sp=find_spike_bh(tr-movmedian(tr,50,2),4,3);
    sp_soma=find_spike_bh(tr_ref-movprc(tr_ref,200,30,2),4,2);

    tr_sub=mean(tr_ref,1)-movprc(mean(tr_ref,1),200,45,2);
    tr_sub=get_subthreshold(tr_sub,sp_soma,5,10);
    

    [trans tr_trace]=detect_transient2(tr_sub,[2 1],sp_soma,20);
    transcand=cell2mat(cellfun(@(x) length(x)>2,trans.ISI,'UniformOutput',false));
    meanISI_frnt=cellfun(@(x) mean(x(1:2)),trans.ISI(transcand));
    meanISI_first3=zeros(1,length(trans.length));
    meanISI_first3(transcand)=meanISI_frnt;

    CS_ind=find(trans.spike_number>1 & meanISI_first3<20);
    CS_trace=ismember(tr_trace,CS_ind);
    CS_spike=sp_soma.*bwlabel(CS_trace);
    [~, CS_spike_time]=unique(CS_spike);

    sp_total=max([sp_soma; sp(2:end,:)],[],1);
    bAP_ind=zeros(1,size(tr,2));
    bAP_ind(unique(find(sp_soma)'+[-1:2]))=1;

    SpikeClassMat=zeros(3,size(tr,2));
    SpikeClassMat(1,:)=sp_soma.*(1-CS_trace); %bAPs
    SpikeClassMat(2,CS_spike_time(2:end))=1; %Complex spikes
    SpikeClassMat(3,:)=sp_total.*(1-bAP_ind); %dSpikes

    Result.spike=[sp_soma; sp(2:end,:)];
    Result.SpClass=SpikeClassMat;
    Result.CStrace=CS_trace;
    
    show_traces_spikes(Result.normTraces,Result.spike,[Result.SpClass(1:2,:); Result.Blue]);
    save(fullfile(fpath{i},'Result'),'Result')
end

%%

% SomRP_path=[15 22 54 72 86 92 4];
% DenRP_path=[13 23 55 74 88 91 3];

SomRP_path=[15 22 54 72 86 92 4];
DenRP_path=[13 23 55 74 88 91 3];
path_cat=[SomRP_path; DenRP_path];
rheo_bin=[0:1.3:10];
pulseFR=[]; RampFR=[]; pulseBlue_rb=[];
for f=1:size(path_cat,2)
    for r=1:2
    i=path_cat(r,f);
    load(fullfile(fpath{i},'Result'))
    bwblue=bwlabel(Result.Blue);
    Ramp_period=find(bwblue==max(bwblue));
    Rheobase_frame=find(max(Result.spike(:,Ramp_period),[],1))+Ramp_period(1)-1;
    Rheobase_blue=Result.Blue(round(mean(Rheobase_frame(1:3))));
    %Rheobase_blue=1;

    firstpulst_bw=max(bwblue)-5;
    for b=1:5
    pulse_period=find(bwblue==firstpulst_bw+b-1);    
    pulseBlue_rb(b,f,r)=median(Result.Blue(pulse_period))/Rheobase_blue;
    pulseFR(b,:,f,r)=sum(Result.SpClass(1:2,pulse_period),2)';
    end

    RampFR(:,:,f,r)=NaN(length(rheo_bin),2);
    Ramp_period=find(bwblue==max(bwblue));    
    binnedRheo_ramp=ceil(Result.Blue(Ramp_period)/Rheobase_blue/(rheo_bin(2)-rheo_bin(1)));
    for b=binnedRheo_ramp
    bin_frm=find(binnedRheo_ramp==b);
    RampFR(b,:,f,r)=sum(Result.SpClass(1:2,bin_frm+Ramp_period(1)-1),2)';
    end
    end
end

figure(27); clf; cmap=distinguishable_colors(6); bin_width=1.3;
tiledlayout(6,2);
load(fullfile(fpath{SomRP_path(5)},'Result'));
tr_somStim=rescale(Result.normTraces(1,:));
tr_cssomStim=Result.CStrace;
load(fullfile(fpath{DenRP_path(6)},'Result'));
tr_ddStim=rescale(Result.normTraces(1,:));
tr_csddStim=Result.CStrace;
ax1=nexttile([3 2]);
plot(tr_somStim,'color',cmap(1,:)); hold all
show_cs=tr_somStim.*tr_cssomStim; show_cs(tr_cssomStim==0)=NaN;
plot(show_cs,'color',[0.7 0 0.2]); hold all
plot(tr_ddStim+1,'color',cmap(1,:)); hold all
show_cs=tr_ddStim.*tr_csddStim; show_cs(tr_csddStim==0)=NaN;
axis off
text(-500,0.7,'Soma Stimulation','FontSize',13)
text(-500,1.7,'Distal dendrite Stimulation','FontSize',13)
text(9500,2,'Complex spikes','FontSize',13,'color',[0.7 0 0.2])
plot(show_cs+1,'color',[0.7 0 0.2]); hold all
ax2=nexttile([1 2]);
plot(Result.Blue)
linkaxes([ax1 ax2],'x')

nexttile([2 1])
CS_fracPulseMatSom=squeeze(pulseFR(:,2,:,1)./sum(pulseFR(:,:,:,1),2));
CS_fracPulseMatDD=squeeze(pulseFR(:,2,:,2)./sum(pulseFR(:,:,:,2),2));
markers = {'o', 's', 'd', '^', 'v', '>', '<', 'p', 'h', '+', '*'};
for f=1:size(pulseFR,3)
plot(squeeze(pulseBlue_rb(:,f,1))',CS_fracPulseMatSom(:,f),'marker',markers{f},'color',cmap(1,:),'linestyle','none'); hold all
plot(squeeze(pulseBlue_rb(:,f,2))',CS_fracPulseMatDD(:,f),'marker',markers{f},'color',cmap(2,:),'linestyle','none'); hold all
end
    clear Ms Ss
for b=1:max(max(ceil(pulseBlue_rb(:,:,1)/bin_width)))
    Ms(b)=mean(CS_fracPulseMatSom(b==ceil(pulseBlue_rb(:,:,1)/bin_width)),'omitnan');
    Ss(b)=std(CS_fracPulseMatSom(b==ceil(pulseBlue_rb(:,:,1)/bin_width)),'omitnan');
end
clear Md Sd
for b=1:max(max(ceil(pulseBlue_rb(:,:,2)/bin_width)))
    Md(b)=mean(CS_fracPulseMatDD(b==ceil(pulseBlue_rb(:,:,2)/bin_width)),'omitnan');
    Sd(b)=std(CS_fracPulseMatDD(b==ceil(pulseBlue_rb(:,:,2)/bin_width)),'omitnan');
end
errorbar([1:max(max(ceil(pulseBlue_rb(:,:,1)/bin_width)))],Ms,Ss,'color',cmap(1,:),'LineWidth',1.5); hold all
errorbar([1:max(max(ceil(pulseBlue_rb(:,:,2)/bin_width)))],Md,Sd,'color',cmap(2,:),'LineWidth',1.5);

xlabel('Optical Rheobase')
ylabel('Fraction of complex spike')
title('Pulse Stimulation')

nexttile([2 1])
CS_fracRampMatSom=squeeze(RampFR(:,2,:,1)./sum(RampFR(:,:,:,1),2));
CS_fracRampMatDD=squeeze(RampFR(:,2,:,2)./sum(RampFR(:,:,:,2),2));
% for f=1:size(RampFR,3)
% plot(rheo_bin,CS_fracRampMatSom,'color',cmap(1,:)); hold all
% plot(rheo_bin,CS_fracRampMatDD,'color',cmap(2,:))
% end
errorbar(rheo_bin,mean(CS_fracRampMatSom,2,'omitnan'),std(CS_fracRampMatSom,0,2,'omitnan'),'color',cmap(1,:),'LineWidth',1.5); hold all
errorbar(rheo_bin,mean(CS_fracRampMatDD,2,'omitnan'),std(CS_fracRampMatDD,0,2,'omitnan'),'color',cmap(2,:),'LineWidth',1.5)
xlabel('Optical Rheobase')
ylabel('Fraction of complex spike')
title('Ramp Stimulation')




