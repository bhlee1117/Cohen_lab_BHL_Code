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

%% Stim
for i=145
    load([fpath{i} '/Result.mat'])
    load(fullfile(fpath{i},"output_data.mat"))
    sz=double(Device_Data{1, 3}.ROI([2 4]));
    mov_mc=double(readBinMov([fpath{i} '/mc_ShutterReg01.bin'],sz(2),sz(1)));
    CamCounter=Device_Data{1, 2}.Counter_Inputs(1, 1).data;
    CamTrigger=find(CamCounter(2:end)-CamCounter(1:end-1));
    DMDtrigger=Device_Data{1, 2}.buffered_tasks(1, 2).channels(1, 3).data;
    DMDtrigger=DMDtrigger(CamTrigger);
end
[blueDMDimg bluePatt]=get_blueDMDPatt(Device_Data,'stack');

%%
bound=5;
mean_F=squeeze(mean(mov_mc(bound:end-bound,bound:end-bound,:),[1 2]));
[~, ~, blueOff]=get_blueoffTrace(mean_F,[Result.Blue],50,0);
[y_fit]=expfitDM_2(find(blueOff)',mean_F(find(blueOff)),[1:size(mov_mc,3)]',[1000 100]);
bkg(1,:)=y_fit;
figure(6); clf;
plot(bkg)
hold all
plot(mean_F)
mov_res=mov_mc-mean(mov_mc,3);
mov_res = SeeResiduals(mov_res,Result.mc);
mov_res = SeeResiduals(mov_res,Result.mc.^2);
mov_res = SeeResiduals(mov_res,Result.mc(:,1).*Result.mc(:,end));
mov_res= SeeResiduals(mov_res,bkg,1);
%%
coord_1d=dim_reduce(get_coord(Result.ftprnt));
[~, dist_order]=sort(coord_1d,'descend');
F0=tovec(Result.ref_im)'*tovec(Result.ftprnt.*double(max(Result.bvMask,[],3)==0));
NormTrace=Result.traces_bvMask./F0';

tr_norm= Result.traces_bvMask-movprc(Result.traces_bvMask,100,30,2);
tr_norm= tr_norm./get_threshold(tr_norm,1);
spike = find_spike_bh(tr_norm,4,1.5);
spike_time=tovec(find(spike(1,:))'+[-5:20]);
spike_time(spike_time<1 | spike_time>size(Result.traces,2))=[];
spike_time=unique(spike_time);
spike_erode_trace=zeros(1,size(Result.traces,2));
spike_erode_trace(spike_time)=1;

nTau=[-30:90]; nROI=size(NormTrace,1);
bwDMDtrigger=cumsum(double([0 DMDtrigger(2:end)-DMDtrigger(1:end-1)]>0))+1;
bwBlue=Result.Blue.*bwDMDtrigger;
BlueOnset=[0 Result.Blue(2:end)-Result.Blue(1:end-1)]>0;
STA_Patt_Mov=[]; STA_Patt_Mat=[];

for b=1:max(bwBlue)
    Patt_blueonset{b}=find(bwBlue==b & BlueOnset);
    isSpike=sum(ismember(Patt_blueonset{b}'+nTau,find(spike_erode_trace)),2)>0;
    rshMov=reshape(mov_res(:,:,Patt_blueonset{b}'+nTau),sz(2),sz(1),[],length(nTau));
    rshMat=reshape(NormTrace(:,Patt_blueonset{b}'+nTau),nROI,[],length(nTau));
    STA_Patt_Mov(:,:,:,b)=squeeze(mean(rshMov(:,:,~isSpike,:),3));
    STA_Patt_Mat(:,:,b)=squeeze(mean(rshMat(:,~isSpike,:),2));
end
STA_Patt_Mov=STA_Patt_Mov-mean(STA_Patt_Mov(:,:,1:20,:),3);
%%
figure(10); clf;
overlap_BlueDMD=tovec(Result.Structure_bin)'*tovec(blueDMDimg);
nPatt=size(blueDMDimg,3); noi=[20 21]; dTrdt=[];
dTr=squeeze(mean(STA_Patt_Mat(noi,:,:),1));
dTr_soma=squeeze(mean(STA_Patt_Mat(1,:,:),1));
for p=1:nPatt
dTrdt(:,p)=get_slope(dTr(:,p),5);    
end
dTr=movmean(dTr,7,1);
dTr=dTr-mean(dTr(1:25,:),1);
nexttile([1 1])
l=plot(nTau,dTr); 
%l=plot(nTau,dTr-prctile(dTr,40,1)); 
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(turbo(nPatt),2));
xlabel('Time from stimulation onset (ms)')
ylabel('\DeltaF/F')
title(['Stim. triggered average trace ' newline 'of stimulated branch'])

nexttile([1 1])
l=plot(nTau,movmean(dTrdt,10,1)); 
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(turbo(nPatt),2));
xlabel('Time from stimulation onset (ms)')
ylabel('dV/dt')
title(['dV/dt of Stim. TA'])

nexttile([1 1])
plot(overlap_BlueDMD,max(dTr(-nTau(1)+1:-nTau(1)+40,:),[],1))
ylabel('Amplitude of Stim-TA trace')
yyaxis right
plot(overlap_BlueDMD,sum(dTr(-nTau(1)+1:-nTau(1)+40,:),1))
ylabel('AUC of Stim-TA trace')
xlabel('Area in blue stimulation')

nexttile([1 1])
plot(overlap_BlueDMD,max(dTrdt(-nTau(1)+1:-nTau(1)+40,:),[],1))
ylabel('max dV/dt')
xlabel('Area in blue stimulation')
%%
figure(11); clf;
tiledlayout(4,2); ax1=[];
for p=1:nPatt
STA_filt=imgaussfilt3(squeeze(-STA_Patt_Mov(bound:end-bound,bound:end-bound,:,p)).*Result.Structure_bin(bound:end-bound,bound:end-bound),[1.5 1.5 1]);
ax1=[ax1 nexttile([1 1])];
show_im=max(STA_filt(:,:,-nTau(1)+5:-nTau(1)+25),[],3);
show_im=imgaussfilt(show_im,3);
%imshow2(show_im,[2 15])
 show_im=grs2rgb(show_im,turbo,2,15);
 imshow2(show_im.*mat2gray(Result.Structure(bound:end-bound,bound:end-bound))*1.5,[]); hold all
 plot(bluePatt{p}(:,2)-bound,bluePatt{p}(:,1)-bound,'color',[0 .5 1],'LineWidth',1.5)
title(['Stimulation area = ' num2str(overlap_BlueDMD(p)*1.17^2,3),'\mum^2'])
end
linkaxes(ax1,'xy')
%%
[~, max_frm]=max(dTr,[],1); pixelsize=6.5/10*180/100;
max_img=[]; kymoImg=[];
for p=1:8
STA_PattFilt=imgaussfilt3(-STA_Patt_Mov(:,:,:,p),[2 2 2]);    
max_img(:,:,p)=STA_PattFilt(:,:,max_frm(p));
end

[max_imgTr kymoRoi]=polyLineKymo3(max_img,7,12,Result.Structure);
kymoROI_cent=cellfun(@mean,kymoRoi,'UniformOutput',false);
KymoROI_inBlue=double(cell2mat(cellfun(@(x) squeeze(blueDMDimg(round(x(2)),round(x(1)),:)),kymoROI_cent,'UniformOutput',false)));
KymoROI_inBlue(KymoROI_inBlue==0)=NaN;
figure(23); clf; tiledlayout(2,4); kymo_noi=setdiff([1:length(kymoRoi)],[12 13]);
for p=1:8
nexttile([1 1])
STA_PattFilt=imgaussfilt3(-STA_Patt_Mov(:,:,:,p),[2 2 2]);    
kymoImg(:,:,p)=apply_clicky(kymoRoi,STA_PattFilt,'no');
imagesc(nTau,[1:length(kymo_noi)],kymoImg(:,kymo_noi,p)',[-5 10]); hold all
plot(KymoROI_inBlue(p,kymo_noi),[1:length(kymo_noi)],'ro')
xlabel('Stimulation onset time (ms)')
ylabel('Kymo ROI #')
end

figure(22); clf;
tiledlayout(2,2);
nexttile([1 2])
imshow2(Result.Structure,[]); hold all; cmap_kymo=jet(length(kymoRoi));
for p=1:length(kymoRoi)
    plot(kymoRoi{p}(:,1),kymoRoi{p}(:,2),'color',cmap_kymo(p,:))
end

kymoROI_cent_offset=cell2mat(kymoROI_cent');
kymoROI_cent_offset=(kymoROI_cent_offset-kymoROI_cent_offset(1,:))*pixelsize;
kymoDist=sqrt(sum(kymoROI_cent_offset.^2,2));

nexttile([1 1]); poi=3;
fit_ind=find((KymoROI_inBlue(poi,kymo_noi)==1));
fit_ind=[max(fit_ind):length(kymo_noi)];
x=kymoDist(kymo_noi(fit_ind));
[y_fit lengthConst]=expfitDM_2(x,kymoImg(max_frm(poi),kymo_noi(fit_ind),poi)',[min(x):max(x)]',[50]);
plot(kymoDist(kymo_noi),kymoImg(max_frm(poi),kymo_noi,poi)); hold all
plot([min(x):max(x)],y_fit,'r')
title(['Length Const: ' num2str(lengthConst) '\mum'])
xlabel('Distance (\mum)')

nexttile([1 1]); poi=4;
fit_ind=find((KymoROI_inBlue(poi,kymo_noi)==1));
fit_ind=[max(fit_ind):length(kymo_noi)];
x=kymoDist(kymo_noi(fit_ind));
[y_fit lengthConst]=expfitDM_2(x,kymoImg(max_frm(poi),kymo_noi(fit_ind),poi)',[min(x):max(x)]',[50]);
plot(kymoDist(kymo_noi),kymoImg(max_frm(poi),kymo_noi,poi)); hold all
plot([min(x):max(x)],y_fit,'r')
title(['Length Const: ' num2str(lengthConst) '\mum'])
xlabel('Distance (\mum)')
%% pca analysis
bound=7;
mov_res_mask=mov_res(:,:,:).*double(max(Result.bvMask,[],3)==0);
subMov=tovec(imresize(imgaussfilt3(mov_res_mask(bound:end-bound,bound:end-bound,:),[1 1 0.1]),1/4));
subMov=subMov-mean(subMov,2);
covMat=subMov'*subMov;

[V, D] = eig(covMat);
D = diag(D);
D = D(end:-1:1);
V = V(:,end:-1:1);
vSign = sign(max(V) - max(-V));  % make the largest value always positive
V = V.*vSign;

nPCs=28;
eigImg=toimg(tovec(mov_res(:,:,:))*V(:,1:nPCs),size(mov_res,1),size(mov_res,2));
figure(4); clf;
for n=1:nPCs
    nexttile([1 1])
    imshow2(eigImg(:,:,n),[])
    title(['PC #', num2str(n), ' Fraction : ' num2str(D(n)/sum(D),2)])
end

[V_ics, mixmat, sepmat]=sorted_ica(V(:,1:nPCs),10);
icsImg=toimg(tovec(mov_res(:,:,:))*V_ics,size(mov_res,1),size(mov_res,2));
figure(5); clf;
for n=1:size(V_ics,2)
    nexttile([1 1])
    %show_footprnt_contour(Result.bvMask,icsImg(:,:,n))
    imshow2(icsImg(:,:,n),[])
    title(['ICS #', num2str(n)])
end
colormap('gray')
%%