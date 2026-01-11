clear
clc;
cd '/Volumes/BHL18TB_D2/Arranged_Data/Prism_OptopatchResult';
[~, ~, raw] = xlsread(['/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/' ...
    'Prism_OptopatchData_Arrangement.xlsx'], 'Sheet1', 'C5:S200');

save_to='/Volumes/BHL18TB_D2/Arranged_Data/Prism_OptopatchResult';
fpath=raw(:,1);
Mouse=cell2mat(raw(:,2));
NeuronInd=cell2mat(raw(:,5));
CamType=raw(:,3);
StructureData=raw(:,10);
PixelSize=cell2mat(raw(:,12));

%% Low stim
f=100;
load([fpath{f} '/OP_Result.mat'])
load(fullfile(fpath{f},"output_data.mat"))
sz=size(Result.ref_im');
mov_mc=[];
mov_mc=cat(3,mov_mc,double(readBinMov([fpath{f} '/mc_ShutterReg' num2str(1,'%02d') '.bin'],sz(2),sz(1))));
load([fpath{f} '/mcTrace' num2str(1,'%02d') '.mat']);
mc=mcTrace.xymean;
t_fit=[1:size(mov_mc,3)];
[bleaching_fit] = expfitDM_2(t_fit',squeeze(mean(mov_mc,[1 2])),t_fit',[100000 10000]);
bv_trace=tovec(mov_mc)'*tovec(Result.bvMask);
bv_trace=SeeResiduals(permute(bv_trace,[2 3 1]),Result.normTraces(1:5,:)');
bv_trace=squeeze(bv_trace)';
bv_trace=movmean(bv_trace,5,1);

mov_res= mov_mc-median(mov_mc,3);
mov_res = SeeResiduals(mov_res,mc);
mov_res = SeeResiduals(mov_res,mc.^2);
mov_res = SeeResiduals(mov_res,mc(:,1).*mc(:,end));
mov_res = SeeResiduals(mov_res,bleaching_fit,1);
mov_res = SeeResiduals(mov_res,bv_trace,1);
mov_res = imresize(mov_res,0.7);
sz2=size(mov_res);

mov_res=-mov_res;
mov_res=maskEdge(mov_res,5,0);
mov_res_sub=get_subthreshold(tovec(mov_res),Result.spike(1,:),11,25);
mov_res_sub=reshape(mov_res_sub,sz2(1),sz2(2),[]);
mov_res_sub(isnan(mov_res_sub))=0;
mov_res_sub=maskBinary(mov_res_sub,imresize(max(Result.bvMask,[],3),0.7),NaN);
mov_res_sub=fillNaNmovie(mov_res_sub);
mov_res_sub_filt=imgaussfilt_NaN(mov_res_sub,1.5);

DirtMov=double(readBinMov([fpath{f} '/DirtMov.bin'],sz(2),sz(1)));
DirtMov=imresize(DirtMov,0.7);
DirtFrame=squeeze(sum(DirtMov,[1 2]));
%mov_res=tovec(mov_res);
%F0img=get_F0img(toimg(mov_res,sz(2),sz(1)));'
%%
[nROI nTime]=size(Result.normTraces);
normTr=Result.normTraces./Result.F0_PCA;
subTr=get_subthreshold(normTr,Result.spike(1,:),7,25);
periSpikeTime=unique(find(Result.spike(1,:)>0)'+[-3:30]);
periSpikeTime(periSpikeTime<1 | periSpikeTime>nTime)=[];
SilentTime=setdiff([1:nTime],periSpikeTime);
subTr_silent=subTr; subTr_silent(:,periSpikeTime)=NaN;
%subTr_silent=subTr_silent-movmedian(subTr_silent,300,2,'omitnan');

Dsign=ones(1,size(Result.interDendDist,2));
Dsign(Result.dist_order(1:find(Result.dist_order==1)-1))=-1;
dendaxis=Result.interDendDist.*Dsign;
dbin=[-300 -50 50 270 400];
labelvec=zeros(length(dbin)-1,nROI);
for b=1:length(dbin)-1
    labelvec(b,dendaxis>dbin(b) & dendaxis<dbin(b+1))=1;
end
normTr_label=(pcafilterTrace(normTr,1:10)'*labelvec')'./sum(labelvec,2);
%normTr_label2=(pcafilterTrace(normTr,1:10)'*labelvec')'./sum(labelvec,2);
subTr_label=(subTr'*labelvec')'./sum(labelvec,2);
subTr_label_silent=(subTr_silent'*labelvec')'./sum(labelvec,2);

peaktroughMat=[];
for n=1:size(subTr_label,1)
    % [~, peaks] = findpeaks(subTr_label(n,:),'MinPeakHeight', 1,'MinPeakDistance',50, ...
    %     'MinPeakProminence', 0.5);  % Prominence can also be tuned
    [~, peaks] = findpeaks(movmean(subTr_label(n,:)-movprc(subTr_label(n,:),1000,30,2),50),'MinPeakHeight', 0.7,'MinPeakDistance',50, ...
         'MinPeakProminence', 0.5);  % Prominence can also be tuned
    peaks(ismember(peaks,periSpikeTime))=[];
    subTr_tmp=subTr_label_silent(n,:);
    subTr_tmp(Result.Blue>0)=NaN;
    subTr_tmp=movmedian(subTr_tmp,1000,'omitnan');
    localAmps=subTr_label(n,peaks)-subTr_tmp(peaks);

    [~, troughs] = findpeaks(movmean(median(tovec(subTr_label_silent(:,Result.Blue>0)),'omitnan')-subTr_label(n,:),50),'MinPeakHeight', 0.7,'MinPeakDistance',50, ...
        'MinPeakProminence', 0.35);  % Prominence can also be tuned
    troughs(ismember(troughs,periSpikeTime))=[];
    subTr_tmp2=subTr_label_silent(n,:);
    subTr_tmp2(Result.Blue==0)=NaN;
    subTr_tmp2=movmedian(subTr_tmp2,1000,'omitnan');
    localAmps2=subTr_label(n,troughs)-subTr_tmp2(troughs);

    peaktroughMat=[peaktroughMat; [repmat([n Result.interDendDist(1,n)' 1],length(peaks),1) peaks' subTr(peaks)' localAmps' Result.Blue(peaks)']];
    peaktroughMat=[peaktroughMat; [repmat([n Result.interDendDist(1,n)' 0],length(troughs),1) troughs' subTr(troughs)' localAmps2' Result.Blue(troughs)']];
end
FieldName={'ROI','Distance','IsPeak','frame','Amplitude','LocalAmplitude','IsBlue'};
peaktroughMat=array2table(peaktroughMat,'VariableNames',FieldName);

figure(2); clf;
imagesc(normTr(Result.dist_order,:),[-2 5]);
hold all
%plot(find(Result.spike(1,:))-1,4,'ro')
plot(peaktroughMat.frame(peaktroughMat.IsPeak>0 & ismember(peaktroughMat.frame,find(Result.Blue==0)))-2,5,'ro','LineWidth',2) %Peak
plot(peaktroughMat.frame(peaktroughMat.IsPeak==0 & ismember(peaktroughMat.frame,find(Result.Blue>0)))-3,4.5,'ko','LineWidth',2) %Trough
plot(Result.Blue-4,'b')
colormap(turbo);
%% Filter movie
% nTau_EXCINHpatt=[50];
% strImg=imresize(Result.Structure,0.7);
% binImg=imresize(highpassDendrite(Result.ref_im,Result.SWC)>5,0.7);

% BlueOn_silent=setdiff(find(Result.Blue>0),periSpikeTime);
% BlueOndF=median(mov_res_sub_filt(:,:,BlueOn_silent),3,'omitnan');
% inhibitory_frame=peaktroughMat.frame(peaktroughMat.IsPeak==0 & peaktroughMat.IsBlue>0);
% excitatory_frame=peaktroughMat.frame(peaktroughMat.IsPeak>0 & peaktroughMat.IsBlue==0);
% 
% [~, inhibitory_dFMat, inhibitory_frame]=get_STA(tovec(BlueOndF-mov_res_sub_filt),inhibitory_frame,nTau_EXCINHpatt,nTau_EXCINHpatt);
% inhibitory_dFMat=reshape(permute(inhibitory_dFMat,[1 3 2]),sz2(1),sz2(2),[]);
% [~, excitatory_dFMat, excitatory_frame]=get_STA(tovec(mov_res_sub_filt),excitatory_frame,nTau_EXCINHpatt,nTau_EXCINHpatt);
% excitatory_dFMat=reshape(permute(excitatory_dFMat,[1 3 2]),sz2(1),sz2(2),[]);
% 
% [~,~,u_all] = get_eigvector(tovec(mov_res_sub_filt(:,:,DirtFrame==0))',40);
% [~,~,u_exc] = get_eigvector(tovec(imgaussfilt(excitatory_dFMat,1))',40);
% [~,~,u_inh] = get_eigvector(tovec(imgaussfilt(inhibitory_dFMat,1))',40);
% u_all=reshape(u_all,sz2(1),sz2(2),[]); 
% u_inh=reshape(u_inh,sz2(1),sz2(2),[]); u_exc=reshape(u_exc,sz2(1),sz2(2),[]);

% template2use=cat(3,u_all(:,:,1:2),u_inh(:,:,1:2),u_exc(:,:,1:2));
% mov_res_sub_filt2=pcafilt_template(mov_res_sub_filt,template2use);

%F0img=std(mov_res_sub_filt2,0,3);
F0img=get_F0img_PCA(mov_res_sub_filt);
% F0img2=get_F0img(toimg(mov_res,sz2(1),sz2(2)));
%F0img=sqrt(mean(mov_res_sub_filt(:,:,1:end-1).*mov_res_sub_filt(:,:,2:end),3));

%%
nTau_EXCINHpatt=[50];
strImg=imresize(Result.Structure,0.7);
binImg=imresize(highpassDendrite(Result.ref_im,Result.SWC)>5,0.7);
cmap=gen_colormap([0 0.2 1; 0 0 0; 1 0 0],256);

BlueOn_silent=setdiff(find(Result.Blue>0),periSpikeTime);
BlueOndF=prctile(mov_res_sub_filt2(:,:,BlueOn_silent),45,3);
inhibitory_frame=peaktroughMat.frame(peaktroughMat.IsPeak==0 & peaktroughMat.IsBlue>0);
excitatory_frame=peaktroughMat.frame(peaktroughMat.IsPeak>0 & peaktroughMat.IsBlue==0);

[~, inhibitory_dFMat, inhibitory_frame]=get_STA(tovec(BlueOndF-mov_res_sub_filt),inhibitory_frame,nTau_EXCINHpatt,nTau_EXCINHpatt);
%inhibitory_dFMat=prctile(inhibitory_dFMat,80,3)-inhibitory_dFMat;
[~, inhibitory_dirtMat]=get_STA(DirtFrame',inhibitory_frame,nTau_EXCINHpatt,nTau_EXCINHpatt);
inhibitory_dirtMat=squeeze(inhibitory_dirtMat);
inhInd2use=find(sum(inhibitory_dirtMat(:,nTau_EXCINHpatt+1+[-3:3]),2)<100000);
inhibitory_frame=inhibitory_frame(inhInd2use);
[inhibitory_frame isort]=sort(inhibitory_frame);
%inhibitory_dF=reshape(permute(inhibitory_dF,[1 3 2]),sz(2),sz(1),[]);

[~, excitatory_dFMat, excitatory_frame]=get_STA(tovec(mov_res_sub_filt),excitatory_frame,nTau_EXCINHpatt,nTau_EXCINHpatt);
%excitatory_dFMat=excitatory_dFMat-prctile(excitatory_dFMat,30,3);
[~, excitatory_dirtMat]=get_STA(DirtFrame',excitatory_frame,nTau_EXCINHpatt,nTau_EXCINHpatt);
excitatory_dirtMat=squeeze(excitatory_dirtMat);
excInd2use=find(sum(excitatory_dirtMat(:,nTau_EXCINHpatt+1+[-3:3]),2)<100000);
excitatory_frame=excitatory_frame(excInd2use);
[excitatory_frame esort]=sort(excitatory_frame);
%excitatory_dF=reshape(permute(excitatory_dF,[1 3 2]),sz(2),sz(1),[]);

inhibitory_dF=reshape(mean(inhibitory_dFMat(:,inhInd2use,nTau_EXCINHpatt+1+[-3:3]),3),sz2(1),sz2(2),[]);
excitatory_dF=reshape(mean(excitatory_dFMat(:,excInd2use,nTau_EXCINHpatt+1+[-3:3]),3),sz2(1),sz2(2),[]);

inhibitory_dF=inhibitory_dF(:,:,isort);
excitatory_dF=excitatory_dF(:,:,esort);
% [~, DInh, eigImgInh]=get_eigvector(tovec(imgaussfilt(inhibitory_dF,1))',10);
% [~, DExc, eigImgExc]=get_eigvector(tovec(imgaussfilt(excitatory_dF,1))',10);
% [V D E]=get_eigvector(tovec(mov_res_sub)',10);

% Show footprints
F0img=maskEdge(F0img,5,NaN);
F0img_filt=medfilt2nan(F0img,[3 3]); 

Ftprnts_resz=imresize(Result.ftprnt,0.7);
inhibitory_dF_norm=tovec(maskEdge(inhibitory_dF,5,NaN)); 
inhibitory_dF_norm(tovec(imresize(Result.Structure_bin,0.7)==0),:)=NaN;
inhibitory_dF_norm=toimg(inhibitory_dF_norm,sz2(1),sz2(2));
inhibitory_dF_norm=(inhibitory_dF_norm)./F0img_filt;
inhibitory_dF_norm_colored=[]; excitatory_dF_norm_colored=[];
dFF_range_inh=[-1 1]*1.5; dFF_range_exc=[-1 1]*2.5;
for p=1:size(inhibitory_dF_norm,3)
    p_tmp=inhibitory_dF_norm(:,:,p);
    p_tmp=medfilt2nan(p_tmp,[15 15]);
    inhibitory_dF_norm(:,:,p)=-maskEdge(imgaussfilt_NaN(p_tmp,2),5,NaN);
    inhibitory_dF_norm_colored(:,:,:,p)= grs2rgb(double(inhibitory_dF_norm(:,:,p).*binImg), cmap ,dFF_range_inh(1),dFF_range_inh(2)).*(strImg>0.01);
    inhibitory_dF_norm_colored(:,:,:,p) = grs2rgb(double(strImg), colormap("gray")) + inhibitory_dF_norm_colored(:,:,:,p);
end
inhibitory_dF_norm_colored=maskEdge(inhibitory_dF_norm_colored,5,NaN);

excitatory_dF_norm=tovec(maskEdge(excitatory_dF,5,NaN)); 
excitatory_dF_norm(tovec(imresize(Result.Structure_bin,0.7)==0),:)=NaN;
excitatory_dF_norm=toimg(excitatory_dF_norm,sz2(1),sz2(2));
excitatory_dF_norm=(excitatory_dF_norm)./F0img_filt;
for p=1:size(excitatory_dF_norm,3)
    p_tmp=excitatory_dF_norm(:,:,p);
    p_tmp=medfilt2nan(p_tmp,[5 5]);
    excitatory_dF_norm(:,:,p)=maskEdge(imgaussfilt_NaN(p_tmp,2),5,NaN);
    excitatory_dF_norm_colored(:,:,:,p)= grs2rgb(double(excitatory_dF_norm(:,:,p).*binImg), cmap ,dFF_range_exc(1),dFF_range_exc(2)).*(strImg>0.01);
    excitatory_dF_norm_colored(:,:,:,p) = grs2rgb(double(strImg), colormap("gray")) + excitatory_dF_norm_colored(:,:,:,p);
end
excitatory_dF_norm_colored=maskEdge(excitatory_dF_norm_colored,5,NaN);

figure(22); clf;
for p=1:size(excitatory_dF_norm_colored,4)
nexttile([1 1]);
imshow2(excitatory_dF_norm_colored(:,:,:,p),[]); title([num2str(excitatory_frame(p)),', p= ' num2str(p)]);
end
figure(23); clf;
for p=1:size(inhibitory_dF_norm_colored,4)
nexttile([1 1]);
imshow2(inhibitory_dF_norm_colored(:,:,:,p),[]); title([num2str(inhibitory_frame(p)),', p= ' num2str(p)]);
end
%%
%show Inhibitory footprints

showpatt=[5 7]; xrange2show=[35:600]; somacoord=get_coord(Ftprnts_resz(:,:,1));
dax=([1:size(inhibitory_dF_norm,2)]-somacoord(1))*PixelSize(f);
offset=0; cmap_tr=gen_colormap(Plasma,4); cmap_sub=cmap_tr./3*2;
g=1; INH_str=[];
BlueOnSteadyLevel=prctile(subTr_label_silent(:,Result.Blue>0),40,2);
BlueOffSteadyLevel=prctile(subTr_label_silent(:,Result.Blue==0),50,2);
figure(55); clf; tiledlayout(1,length(showpatt)+1,'Padding','compact'); 
t2show=[-50:50]; ax2=[]; catfprnt=[];
szfprnt=size(inhibitory_dF_norm_colored);
for p=showpatt;%1:size(inhibitory_dF_norm_colored,4)
    catfprnt=[catfprnt ones(szfprnt(1),50,3) fliplr(inhibitory_dF_norm_colored(:,xrange2show,:,p))];
    INH_str{g}=['INH. #' num2str(g)];

    ax2=[ax2 nexttile(g,[1 1])];
    l_sub=plot(t2show,subTr_label(:,inhibitory_frame(p)+t2show)'+[1:4]*offset,'linewidth',1.5); hold all
    arrayfun(@(l,c) set(l,'Color',c{:}),l_sub,num2cell(cmap_tr,2));
    axis off;

    nexttile((length(showpatt)+1),[1 1]);
    inhfprnt=inhibitory_dF_norm(:,:,p); inhfprnt(strImg==0)=NaN;
    plot(-dax(xrange2show),mean(inhfprnt(:,xrange2show),1,'omitnan'),'linewidth',1.5); box off; axis tight; hold all;
    xlabel('Distance (\mum)'); ylabel('Mean subthreshold (Z score)');
    g=g+1;
end
nexttile(length(showpatt),[1 1]);
drawScaleBar(10,'horizontal','color',[0 0 0],'Linewidth',3,'Position',[47 1.4]);
drawScaleBar(1,'vertical','color',[0 0 0],'Linewidth',3,'Position',[47 0.4]);
linkaxes(ax2); ylim([-0.5 1.5]);
nexttile((length(showpatt)+1),[1 1]); legend(INH_str,'Box','off','Location','best');
nexttile(1,[1 1]); legend(({'Basal','Soma','Apical','Distal'}),'Box','off','Location','best','NumColumns',2);
set_fontsize(13);
figure(65); clf;
nexttile(1,[1 length(showpatt)]);
imshow2(catfprnt,[]); 
cb=colorbar; cb.Label.String='Z score'; colormap(gen_colormap([0 0.5 1;1 1 1;1 0 0],256)); cb.Ticks=[0 1]; cb.TickLabels=dFF_range_inh;
set_fontsize(15);

showpatt2=[13 18];
offset=0; g=1; EXC_str=[];
figure(56); clf; tiledlayout(1,length(showpatt2)+1,'Padding','compact'); 
t2show=[-50:50]; ax2=[]; catfprnt=[];
szfprnt=size(excitatory_dF_norm_colored);
for p=showpatt2;%1:size(inhibitory_dF_norm_colored,4)
    catfprnt=[catfprnt ones(szfprnt(1),50,3) fliplr(excitatory_dF_norm_colored(:,xrange2show,:,p))];
    EXC_str{g}=['EXC. #' num2str(g)];

    ax2=[ax2 nexttile(g,[1 1])];
    l_sub=plot(t2show,subTr_label(:,excitatory_frame(p)+t2show)'+[1:4]*offset,'linewidth',1.5); hold all
    arrayfun(@(l,c) set(l,'Color',c{:}),l_sub,num2cell(cmap_tr,2));
    axis off;

    nexttile((length(showpatt2)+1),[1 1]);
    inhfprnt=excitatory_dF_norm(:,:,p); inhfprnt(strImg==0)=NaN;
    plot(-dax(xrange2show),mean(inhfprnt(:,xrange2show),1,'omitnan'),'linewidth',1.5); box off; axis tight; hold all;
    xlabel('Distance (\mum)'); ylabel('Mean subthreshold (Z score)');
    g=g+1;
end
nexttile(length(showpatt2),[1 1]);
drawScaleBar(10,'horizontal','color',[0 0 0],'Linewidth',3,'Position',[47 1.4]);
drawScaleBar(1,'vertical','color',[0 0 0],'Linewidth',3,'Position',[47 0.4]);
linkaxes(ax2); ylim([-0.25 1.8]);
nexttile((length(showpatt2)+1),[1 1]); legend(EXC_str,'Box','off','Location','best');
%nexttile(1,[1 1]); legend(({'Basal','Soma','Apical','Distal'}),'Box','off','Location','best','NumColumns',2);
set_fontsize(13);
figure(66); clf;
nexttile(1,[1 length(showpatt2)]);
imshow2(catfprnt,[]); 
cb=colorbar; cb.Label.String='Z score'; colormap(gen_colormap([0 0.5 1;1 1 1;1 0 0],256)); cb.Ticks=[0 1]; cb.TickLabels=dFF_range_exc;
set_fontsize(15);

figure(57); clf; tiledlayout(1,4);
offset=6; cmap_tr=gen_colormap(Plasma,4); cmap_sub=cmap_tr./2;
Ftprnt_lab=toimg(tovec(Result.ftprnt>0)*labelvec',sz(2),sz(1));
t=[1:nTime];
SWC=Result.SWC;
SWC(:,4)=5; SWC(1,4)=30;
% nexttile([1 1]);
% showScaleScatter([1:4], SWC, Ftprnt_lab ,gen_colormap(Plasma,256),[1 4]);
% drawScaleBar(100/PixelSize(f),'vertical')
nexttile([1 4]);
l=plot(t,fliplr(normTr_label')+[1:4]*offset); hold all
l_sub=plot(t,fliplr(subTr_label')+[1:4]*offset);
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(flipud(cmap_tr),2));
arrayfun(@(l,c) set(l,'Color',c{:}),l_sub,num2cell(flipud(cmap_sub),2));
hold all
%plot(find(Result.spike(1,:)),offset*4.8,'k.')
plot([excitatory_frame(showpatt2); excitatory_frame(showpatt2)],[5 offset*5],'r--','linewidth',1.5)
plot([inhibitory_frame(showpatt); inhibitory_frame(showpatt)],[5 offset*5],'b--','linewidth',1.5)
plot(Result.Blue*5,'color',[0 0.5 1],'linewidth',1.5)
axis off; xlim([100 9999]);
drawScaleBar(1000,'Horizontal')


%% Low stim (Cell2)
f=126;
load([fpath{f} '/OP_Result.mat'])
load(fullfile(fpath{f},"output_data.mat"))
sz=size(Result.ref_im');
mov_mc=[];
mov_mc=cat(3,mov_mc,double(readBinMov([fpath{f} '/mc_ShutterReg' num2str(1,'%02d') '.bin'],sz(2),sz(1))));
load([fpath{f} '/mcTrace' num2str(1,'%02d') '.mat']);
mc=mcTrace.xymean;
t_fit=[1:size(mov_mc,3)];
[bleaching_fit] = expfitDM_2(t_fit',squeeze(mean(mov_mc,[1 2])),t_fit',[100000 10000]);
bv_trace=tovec(mov_mc)'*tovec(Result.bvMask);
bv_trace=SeeResiduals(permute(bv_trace,[2 3 1]),Result.normTraces(1:5,:)');
bv_trace=squeeze(bv_trace)';
bv_trace=movmean(bv_trace,5,1);

mov_res= mov_mc-median(mov_mc,3);
mov_res = SeeResiduals(mov_res,mc);
mov_res = SeeResiduals(mov_res,mc.^2);
mov_res = SeeResiduals(mov_res,mc(:,1).*mc(:,end));
mov_res = SeeResiduals(mov_res,bleaching_fit,1);
mov_res = SeeResiduals(mov_res,bv_trace,1);
mov_res = imresize(mov_res,0.7);
sz2=size(mov_res);

mov_res=-mov_res;
mov_res=maskEdge(mov_res,5,0);
mov_res_sub=get_subthreshold(tovec(mov_res),Result.spike(1,:),11,25);
mov_res_sub=reshape(mov_res_sub,sz2(1),sz2(2),[]);
mov_res_sub(isnan(mov_res_sub))=0;
mov_res_sub=maskBinary(mov_res_sub,imresize(max(Result.bvMask,[],3),0.7),NaN);
mov_res_sub=fillNaNmovie(mov_res_sub);
mov_res_sub_filt=imgaussfilt_NaN(mov_res_sub,1.5);

DirtMov=double(readBinMov([fpath{f} '/DirtMov.bin'],sz(2),sz(1)));
DirtMov=imresize(DirtMov,0.7);
DirtFrame=squeeze(sum(DirtMov,[1 2]));
%mov_res=tovec(mov_res);
%F0img=get_F0img(toimg(mov_res,sz(2),sz(1)));'
%%
[nROI nTime]=size(Result.normTraces);
normTr=Result.normTraces./Result.F0_PCA;
subTr=get_subthreshold(normTr,Result.spike(1,:),7,25);
periSpikeTime=unique(find(Result.spike(1,:)>0)'+[-3:30]);
periSpikeTime(periSpikeTime<1 | periSpikeTime>nTime)=[];
SilentTime=setdiff([1:nTime],periSpikeTime);
subTr_silent=subTr; subTr_silent(:,periSpikeTime)=NaN;
%subTr_silent=subTr_silent-movmedian(subTr_silent,300,2,'omitnan');

Dsign=ones(1,size(Result.interDendDist,2));
Dsign(Result.dist_order(1:find(Result.dist_order==1)-1))=-1;
dendaxis=Result.interDendDist.*Dsign;
dbin=[-300 -50 50 200 400];
labelvec=zeros(length(dbin)-1,nROI);
for b=1:length(dbin)-1
    labelvec(b,dendaxis>dbin(b) & dendaxis<dbin(b+1))=1;
end
normTr_label=(pcafilterTrace(normTr,1:10)'*labelvec')'./sum(labelvec,2);
%normTr_label2=(pcafilterTrace(normTr,1:10)'*labelvec')'./sum(labelvec,2);
subTr_label=(subTr'*labelvec')'./sum(labelvec,2);
subTr_label_silent=(subTr_silent'*labelvec')'./sum(labelvec,2);

peaktroughMat=[];
for n=1:size(subTr_label,1)
    % [~, peaks] = findpeaks(subTr_label(n,:),'MinPeakHeight', 1,'MinPeakDistance',50, ...
    %     'MinPeakProminence', 0.5);  % Prominence can also be tuned
    [~, peaks] = findpeaks(movmean(subTr_label(n,:)-movprc(subTr_label(n,:),1000,30,2),50),'MinPeakHeight', 0.7,'MinPeakDistance',50, ...
         'MinPeakProminence', 0.5);  % Prominence can also be tuned
    peaks(ismember(peaks,periSpikeTime))=[];
    subTr_tmp=subTr_label_silent(n,:);
    subTr_tmp(Result.Blue>0)=NaN;
    subTr_tmp=movmedian(subTr_tmp,1000,'omitnan');
    localAmps=subTr_label(n,peaks)-subTr_tmp(peaks);

    [~, troughs] = findpeaks(movmean(median(tovec(subTr_label_silent(:,Result.Blue>0)),'omitnan')-subTr_label(n,:),50),'MinPeakHeight', 0.7,'MinPeakDistance',50, ...
        'MinPeakProminence', 0.35);  % Prominence can also be tuned
    troughs(ismember(troughs,periSpikeTime))=[];
    subTr_tmp2=subTr_label_silent(n,:);
    subTr_tmp2(Result.Blue==0)=NaN;
    subTr_tmp2=movmedian(subTr_tmp2,1000,'omitnan');
    localAmps2=subTr_label(n,troughs)-subTr_tmp2(troughs);

    peaktroughMat=[peaktroughMat; [repmat([n Result.interDendDist(1,n)' 1],length(peaks),1) peaks' subTr(peaks)' localAmps' Result.Blue(peaks)']];
    peaktroughMat=[peaktroughMat; [repmat([n Result.interDendDist(1,n)' 0],length(troughs),1) troughs' subTr(troughs)' localAmps2' Result.Blue(troughs)']];
end
FieldName={'ROI','Distance','IsPeak','frame','Amplitude','LocalAmplitude','IsBlue'};
peaktroughMat=array2table(peaktroughMat,'VariableNames',FieldName);

figure(2); clf;
imagesc(normTr(Result.dist_order,:),[-2 5]);
hold all
%plot(find(Result.spike(1,:))-1,4,'ro')
plot(peaktroughMat.frame(peaktroughMat.IsPeak>0 & ismember(peaktroughMat.frame,find(Result.Blue==0)))-2,5,'ro','LineWidth',2) %Peak
plot(peaktroughMat.frame(peaktroughMat.IsPeak==0 & ismember(peaktroughMat.frame,find(Result.Blue>0)))-3,4.5,'ko','LineWidth',2) %Trough
plot(Result.Blue-4,'b')
colormap(turbo);
%% Filter movie
% nTau_EXCINHpatt=[50];
% strImg=imresize(Result.Structure,0.7);
% binImg=imresize(highpassDendrite(Result.ref_im,Result.SWC)>5,0.7);

% BlueOn_silent=setdiff(find(Result.Blue>0),periSpikeTime);
% BlueOndF=median(mov_res_sub_filt(:,:,BlueOn_silent),3,'omitnan');
% inhibitory_frame=peaktroughMat.frame(peaktroughMat.IsPeak==0 & peaktroughMat.IsBlue>0);
% excitatory_frame=peaktroughMat.frame(peaktroughMat.IsPeak>0 & peaktroughMat.IsBlue==0);
% 
% [~, inhibitory_dFMat, inhibitory_frame]=get_STA(tovec(BlueOndF-mov_res_sub_filt),inhibitory_frame,nTau_EXCINHpatt,nTau_EXCINHpatt);
% inhibitory_dFMat=reshape(permute(inhibitory_dFMat,[1 3 2]),sz2(1),sz2(2),[]);
% [~, excitatory_dFMat, excitatory_frame]=get_STA(tovec(mov_res_sub_filt),excitatory_frame,nTau_EXCINHpatt,nTau_EXCINHpatt);
% excitatory_dFMat=reshape(permute(excitatory_dFMat,[1 3 2]),sz2(1),sz2(2),[]);
% 
% [~,~,u_all] = get_eigvector(tovec(mov_res_sub_filt(:,:,DirtFrame==0))',40);
% [~,~,u_exc] = get_eigvector(tovec(imgaussfilt(excitatory_dFMat,1))',40);
% [~,~,u_inh] = get_eigvector(tovec(imgaussfilt(inhibitory_dFMat,1))',40);
% u_all=reshape(u_all,sz2(1),sz2(2),[]); 
% u_inh=reshape(u_inh,sz2(1),sz2(2),[]); u_exc=reshape(u_exc,sz2(1),sz2(2),[]);

% template2use=cat(3,u_all(:,:,1:2),u_inh(:,:,1:2),u_exc(:,:,1:2));
% mov_res_sub_filt2=pcafilt_template(mov_res_sub_filt,template2use);

%F0img=std(mov_res_sub_filt2,0,3);
F0img=get_F0img_PCA(mov_res_sub_filt);
% F0img2=get_F0img(toimg(mov_res,sz2(1),sz2(2)));
%F0img=sqrt(mean(mov_res_sub_filt(:,:,1:end-1).*mov_res_sub_filt(:,:,2:end),3));

%%
nTau_EXCINHpatt=[50];
strImg=imresize(Result.Structure,0.7);
binImg=imresize(highpassDendrite(Result.ref_im,Result.SWC)>5,0.7);
cmap=gen_colormap([0 0.2 1; 0 0 0; 1 0 0],256);

BlueOn_silent=setdiff(find(Result.Blue>0),periSpikeTime);
BlueOndF=prctile(mov_res_sub_filt(:,:,BlueOn_silent),50,3);
inhibitory_frame=peaktroughMat.frame(peaktroughMat.IsPeak==0 & peaktroughMat.IsBlue>0);
excitatory_frame=peaktroughMat.frame(peaktroughMat.IsPeak>0 & peaktroughMat.IsBlue==0);

[~, inhibitory_dFMat, inhibitory_frame]=get_STA(tovec(BlueOndF-mov_res_sub_filt),inhibitory_frame,nTau_EXCINHpatt,nTau_EXCINHpatt);
%inhibitory_dFMat=prctile(inhibitory_dFMat,80,3)-inhibitory_dFMat;
[~, inhibitory_dirtMat]=get_STA(DirtFrame',inhibitory_frame,nTau_EXCINHpatt,nTau_EXCINHpatt);
inhibitory_dirtMat=squeeze(inhibitory_dirtMat);
inhInd2use=find(sum(inhibitory_dirtMat(:,nTau_EXCINHpatt+1+[-3:3]),2)<100000);
inhibitory_frame=inhibitory_frame(inhInd2use);
[inhibitory_frame isort]=sort(inhibitory_frame);
%inhibitory_dF=reshape(permute(inhibitory_dF,[1 3 2]),sz(2),sz(1),[]);

[~, excitatory_dFMat, excitatory_frame]=get_STA(tovec(mov_res_sub_filt),excitatory_frame,nTau_EXCINHpatt,nTau_EXCINHpatt);
%excitatory_dFMat=excitatory_dFMat-prctile(excitatory_dFMat,30,3);
[~, excitatory_dirtMat]=get_STA(DirtFrame',excitatory_frame,nTau_EXCINHpatt,nTau_EXCINHpatt);
excitatory_dirtMat=squeeze(excitatory_dirtMat);
excInd2use=find(sum(excitatory_dirtMat(:,nTau_EXCINHpatt+1+[-3:3]),2)<100000);
excitatory_frame=excitatory_frame(excInd2use);
[excitatory_frame esort]=sort(excitatory_frame);
%excitatory_dF=reshape(permute(excitatory_dF,[1 3 2]),sz(2),sz(1),[]);

inhibitory_dF=reshape(mean(inhibitory_dFMat(:,inhInd2use,nTau_EXCINHpatt+1+[-3:3]),3),sz2(1),sz2(2),[]);
excitatory_dF=reshape(mean(excitatory_dFMat(:,excInd2use,nTau_EXCINHpatt+1+[-3:3]),3),sz2(1),sz2(2),[]);

inhibitory_dF=inhibitory_dF(:,:,isort);
excitatory_dF=excitatory_dF(:,:,esort);
% [~, DInh, eigImgInh]=get_eigvector(tovec(imgaussfilt(inhibitory_dF,1))',10);
% [~, DExc, eigImgExc]=get_eigvector(tovec(imgaussfilt(excitatory_dF,1))',10);
% [V D E]=get_eigvector(tovec(mov_res_sub)',10);

% Show footprints
F0img=maskEdge(F0img,5,NaN);
F0img_filt=medfilt2nan(F0img,[3 3]); 

Ftprnts_resz=imresize(Result.ftprnt,0.7);
inhibitory_dF_norm=tovec(maskEdge(inhibitory_dF,5,NaN)); 
inhibitory_dF_norm(tovec(imresize(Result.Structure_bin,0.7)==0),:)=NaN;
inhibitory_dF_norm=toimg(inhibitory_dF_norm,sz2(1),sz2(2));
inhibitory_dF_norm=(inhibitory_dF_norm)./F0img_filt;
inhibitory_dF_norm_colored=[]; excitatory_dF_norm_colored=[];
dFF_range_inh=[-1 1]*2.5; dFF_range_exc=[-1 1]*3.5;
for p=1:size(inhibitory_dF_norm,3)
    p_tmp=inhibitory_dF_norm(:,:,p);
    p_tmp=medfilt2nan(p_tmp,[15 15]);
    inhibitory_dF_norm(:,:,p)=-maskEdge(imgaussfilt_NaN(p_tmp,2),5,NaN);
    inhibitory_dF_norm_colored(:,:,:,p)= grs2rgb(double(inhibitory_dF_norm(:,:,p).*binImg), cmap ,dFF_range_inh(1),dFF_range_inh(2)).*(strImg>0.01);
    inhibitory_dF_norm_colored(:,:,:,p) = grs2rgb(double(strImg), colormap("gray")) + inhibitory_dF_norm_colored(:,:,:,p);
end
inhibitory_dF_norm_colored=maskEdge(inhibitory_dF_norm_colored,5,NaN);

excitatory_dF_norm=tovec(maskEdge(excitatory_dF,5,NaN)); 
excitatory_dF_norm(tovec(imresize(Result.Structure_bin,0.7)==0),:)=NaN;
excitatory_dF_norm=toimg(excitatory_dF_norm,sz2(1),sz2(2));
excitatory_dF_norm=(excitatory_dF_norm)./F0img_filt;
for p=1:size(excitatory_dF_norm,3)
    p_tmp=excitatory_dF_norm(:,:,p);
    p_tmp=medfilt2nan(p_tmp,[5 5]);
    excitatory_dF_norm(:,:,p)=maskEdge(imgaussfilt_NaN(p_tmp,2),5,NaN);
    excitatory_dF_norm_colored(:,:,:,p)= grs2rgb(double(excitatory_dF_norm(:,:,p).*binImg), cmap ,dFF_range_exc(1),dFF_range_exc(2)).*(strImg>0.01);
    excitatory_dF_norm_colored(:,:,:,p) = grs2rgb(double(strImg), colormap("gray")) + excitatory_dF_norm_colored(:,:,:,p);
end
excitatory_dF_norm_colored=maskEdge(excitatory_dF_norm_colored,5,NaN);

figure(22); clf;
for p=1:size(excitatory_dF_norm_colored,4)
nexttile([1 1]);
imshow2(excitatory_dF_norm_colored(:,:,:,p),[]); title([num2str(excitatory_frame(p)),', p= ' num2str(p)]);
end
sgtitle('Excitatory footprints');
figure(23); clf;
for p=1:size(inhibitory_dF_norm_colored,4)
nexttile([1 1]);
imshow2(inhibitory_dF_norm_colored(:,:,:,p),[]); title([num2str(inhibitory_frame(p)),', p= ' num2str(p)]);
end
sgtitle('Inhibitory footprints');
%% show footprints (Figure S)

showpatt=[1 11]; xrange2show=[30:380]; somacoord=get_coord(Ftprnts_resz(:,:,1));
dax=([1:size(inhibitory_dF_norm,2)]-somacoord(1))*PixelSize(f);
offset=0; cmap_tr=gen_colormap(Plasma,4); cmap_sub=cmap_tr./3*2;
g=1; INH_str=[];
BlueOnSteadyLevel=prctile(subTr_label_silent(:,Result.Blue>0),40,2);
BlueOffSteadyLevel=prctile(subTr_label_silent(:,Result.Blue==0),50,2);
figure(55); clf; tiledlayout(1,length(showpatt)+1,'Padding','compact'); 
t2show=[-50:50]; ax2=[]; catfprnt=[];
szfprnt=size(inhibitory_dF_norm_colored);
for p=showpatt;%1:size(inhibitory_dF_norm_colored,4)
    catfprnt=[catfprnt ones(szfprnt(1),50,3) fliplr(inhibitory_dF_norm_colored(:,xrange2show,:,p))];
    INH_str{g}=['INH. #' num2str(g)];

    ax2=[ax2 nexttile(g,[1 1])];
    l_sub=plot(t2show,subTr_label(:,inhibitory_frame(p)+t2show)'+[1:4]*offset,'linewidth',1.5); hold all
    arrayfun(@(l,c) set(l,'Color',c{:}),l_sub,num2cell(cmap_tr,2));
    axis off;

    nexttile((length(showpatt)+1),[1 1]);
    inhfprnt=inhibitory_dF_norm(:,:,p); inhfprnt(strImg==0)=NaN;
    plot(-dax(xrange2show),mean(inhfprnt(:,xrange2show),1,'omitnan'),'linewidth',1.5); box off; axis tight; hold all;
    xlabel('Distance (\mum)'); ylabel('Mean subthreshold (Z score)');
    g=g+1;
end
nexttile(length(showpatt),[1 1]);
drawScaleBar(10,'horizontal','color',[0 0 0],'Linewidth',3,'Position',[47 1.4]);
drawScaleBar(1,'vertical','color',[0 0 0],'Linewidth',3,'Position',[47 0.4]);
linkaxes(ax2); ylim([-0.6 2.5]);
nexttile((length(showpatt)+1),[1 1]); legend(INH_str,'Box','off','Location','best');
set_fontsize(13);
figure(65); clf;
nexttile(1,[1 length(showpatt)]);
imshow2(catfprnt,[]); 
cb=colorbar; cb.Label.String='Z score'; colormap(gen_colormap([0 0.5 1;1 1 1;1 0 0],256)); cb.Ticks=[0 1]; cb.TickLabels=dFF_range_inh;
set_fontsize(15);

showpatt2=[16 28];
offset=0; g=1; EXC_str=[];
figure(56); clf; tiledlayout(1,length(showpatt2)+1,'Padding','compact'); 
t2show=[-50:50]; ax2=[]; catfprnt=[];
szfprnt=size(excitatory_dF_norm_colored);
for p=showpatt2;%1:size(inhibitory_dF_norm_colored,4)
    catfprnt=[catfprnt ones(szfprnt(1),50,3) fliplr(excitatory_dF_norm_colored(:,xrange2show,:,p))];
    EXC_str{g}=['EXC. #' num2str(g)];

    ax2=[ax2 nexttile(g,[1 1])];
    l_sub=plot(t2show,subTr_label(:,excitatory_frame(p)+t2show)'+[1:4]*offset,'linewidth',1.5); hold all
    arrayfun(@(l,c) set(l,'Color',c{:}),l_sub,num2cell(cmap_tr,2));
    axis off;

    nexttile((length(showpatt2)+1),[1 1]);
    inhfprnt=excitatory_dF_norm(:,:,p); inhfprnt(strImg==0)=NaN;
    plot(-dax(xrange2show),mean(inhfprnt(:,xrange2show),1,'omitnan'),'linewidth',1.5); box off; axis tight; hold all;
    xlabel('Distance (\mum)'); ylabel('Mean subthreshold (Z score)');
    g=g+1;
end
nexttile(length(showpatt2),[1 1]);
drawScaleBar(10,'horizontal','color',[0 0 0],'Linewidth',3,'Position',[47 1.4]);
drawScaleBar(1,'vertical','color',[0 0 0],'Linewidth',3,'Position',[47 0.4]);
linkaxes(ax2); ylim([-0.2 2.3]);
nexttile((length(showpatt2)+1),[1 1]); legend(EXC_str,'Box','off','Location','best');
nexttile(1,[1 1]); legend(({'Basal','Soma','Apical','Distal'}),'Box','off','Location','best','NumColumns',2);
set_fontsize(13);
figure(66); clf;
nexttile(1,[1 length(showpatt2)]);
imshow2(catfprnt,[]); 
cb=colorbar; cb.Label.String='Z score'; colormap(gen_colormap([0 0.5 1;1 1 1;1 0 0],256)); cb.Ticks=[0 1]; cb.TickLabels=dFF_range_exc;
set_fontsize(15);

figure(57); clf; tiledlayout(1,4);
offset=6; cmap_tr=gen_colormap(Plasma,4); cmap_sub=cmap_tr./2;
Ftprnt_lab=toimg(tovec(Result.ftprnt>0)*labelvec',sz(2),sz(1));
t=[1:nTime];
SWC=Result.SWC;
SWC(:,4)=5; SWC(1,4)=30;
% nexttile([1 1]);
% showScaleScatter([1:4], SWC, Ftprnt_lab ,gen_colormap(Plasma,256),[1 4]);
% drawScaleBar(100/PixelSize(f),'vertical')
nexttile([1 4]);
l=plot(t,fliplr(normTr_label')+[1:4]*offset); hold all
l_sub=plot(t,fliplr(subTr_label')+[1:4]*offset);
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(flipud(cmap_tr),2));
arrayfun(@(l,c) set(l,'Color',c{:}),l_sub,num2cell(flipud(cmap_sub),2));
hold all
%plot(find(Result.spike(1,:)),offset*4.8,'k.')
plot([excitatory_frame(showpatt2); excitatory_frame(showpatt2)],[5 offset*5],'r--','linewidth',1.5)
plot([inhibitory_frame(showpatt); inhibitory_frame(showpatt)],[5 offset*5],'b--','linewidth',1.5)
plot(Result.Blue*5,'color',[0 0.5 1],'linewidth',1.5)
axis off; xlim([100 9999]);
drawScaleBar(1000,'Horizontal')


%% Low stim (Cell3)
f=129;
load([fpath{f} '/OP_Result.mat'])
load(fullfile(fpath{f},"output_data.mat"))
sz=size(Result.ref_im');
mov_mc=[];
mov_mc=cat(3,mov_mc,double(readBinMov([fpath{f} '/mc_ShutterReg' num2str(1,'%02d') '.bin'],sz(2),sz(1))));
load([fpath{f} '/mcTrace' num2str(1,'%02d') '.mat']);
mc=mcTrace.xymean;
t_fit=[1:size(mov_mc,3)];
[bleaching_fit] = expfitDM_2(t_fit',squeeze(mean(mov_mc,[1 2])),t_fit',[100000 10000]);
bv_trace=tovec(mov_mc)'*tovec(Result.bvMask);
bv_trace=SeeResiduals(permute(bv_trace,[2 3 1]),Result.normTraces(1:5,:)');
bv_trace=squeeze(bv_trace)';
bv_trace=movmean(bv_trace,5,1);

mov_res= mov_mc-median(mov_mc,3);
mov_res = SeeResiduals(mov_res,mc);
mov_res = SeeResiduals(mov_res,mc.^2);
mov_res = SeeResiduals(mov_res,mc(:,1).*mc(:,end));
mov_res = SeeResiduals(mov_res,bleaching_fit,1);
mov_res = SeeResiduals(mov_res,bv_trace,1);
mov_res = imresize(mov_res,0.7);
sz2=size(mov_res);

mov_res=-mov_res;
mov_res=maskEdge(mov_res,5,0);
mov_res_sub=get_subthreshold(tovec(mov_res),Result.spike(1,:),11,25);
mov_res_sub=reshape(mov_res_sub,sz2(1),sz2(2),[]);
mov_res_sub(isnan(mov_res_sub))=0;
mov_res_sub=maskBinary(mov_res_sub,imresize(max(Result.bvMask,[],3),0.7),NaN);
mov_res_sub=fillNaNmovie(mov_res_sub);
mov_res_sub_filt=imgaussfilt_NaN(mov_res_sub,1.5);

DirtMov=double(readBinMov([fpath{f} '/DirtMov.bin'],sz(2),sz(1)));
DirtMov=imresize(DirtMov,0.7);
DirtFrame=squeeze(sum(DirtMov,[1 2]));
%DirtFrame=zeros(1,size(mov_res,3));
%mov_res=tovec(mov_res);
%F0img=get_F0img(toimg(mov_res,sz(2),sz(1)));'
%%
[nROI nTime]=size(Result.normTraces);
normTr=Result.normTraces./Result.F0_PCA;
subTr=get_subthreshold(normTr,Result.spike(1,:),7,25);
periSpikeTime=unique(find(Result.spike(1,:)>0)'+[-3:30]);
periSpikeTime(periSpikeTime<1 | periSpikeTime>nTime)=[];
SilentTime=setdiff([1:nTime],periSpikeTime);
subTr_silent=subTr; subTr_silent(:,periSpikeTime)=NaN;
%subTr_silent=subTr_silent-movmedian(subTr_silent,300,2,'omitnan');

Dsign=ones(1,size(Result.interDendDist,2));
Dsign(Result.dist_order(1:find(Result.dist_order==1)-1))=-1;
dendaxis=Result.interDendDist.*Dsign;
dbin=[-300 -50 50 200 400];
labelvec=zeros(length(dbin)-1,nROI);
for b=1:length(dbin)-1
    labelvec(b,dendaxis>dbin(b) & dendaxis<dbin(b+1))=1;
end
normTr_label=(pcafilterTrace(normTr,1:10)'*labelvec')'./sum(labelvec,2);
%normTr_label2=(pcafilterTrace(normTr,1:10)'*labelvec')'./sum(labelvec,2);
subTr_label=(subTr'*labelvec')'./sum(labelvec,2);
subTr_label_silent=(subTr_silent'*labelvec')'./sum(labelvec,2);

peaktroughMat=[];
for n=1:size(subTr_label,1)
    % [~, peaks] = findpeaks(subTr_label(n,:),'MinPeakHeight', 1,'MinPeakDistance',50, ...
    %     'MinPeakProminence', 0.5);  % Prominence can also be tuned
    [~, peaks] = findpeaks(movmean(subTr_label(n,:)-movprc(subTr_label(n,:),1000,30,2),50),'MinPeakHeight', 0.7,'MinPeakDistance',50, ...
         'MinPeakProminence', 0.5);  % Prominence can also be tuned
    peaks(ismember(peaks,periSpikeTime))=[];
    subTr_tmp=subTr_label_silent(n,:);
    subTr_tmp(Result.Blue>0)=NaN;
    subTr_tmp=movmedian(subTr_tmp,1000,'omitnan');
    localAmps=subTr_label(n,peaks)-subTr_tmp(peaks);

    [~, troughs] = findpeaks(movmean(median(tovec(subTr_label_silent(:,Result.Blue>0)),'omitnan')-subTr_label(n,:),50),'MinPeakHeight', 0.7,'MinPeakDistance',50, ...
        'MinPeakProminence', 0.35);  % Prominence can also be tuned
    troughs(ismember(troughs,periSpikeTime))=[];
    subTr_tmp2=subTr_label_silent(n,:);
    subTr_tmp2(Result.Blue==0)=NaN;
    subTr_tmp2=movmedian(subTr_tmp2,1000,'omitnan');
    localAmps2=subTr_label(n,troughs)-subTr_tmp2(troughs);

    peaktroughMat=[peaktroughMat; [repmat([n Result.interDendDist(1,n)' 1],length(peaks),1) peaks' subTr(peaks)' localAmps' Result.Blue(peaks)']];
    peaktroughMat=[peaktroughMat; [repmat([n Result.interDendDist(1,n)' 0],length(troughs),1) troughs' subTr(troughs)' localAmps2' Result.Blue(troughs)']];
end
FieldName={'ROI','Distance','IsPeak','frame','Amplitude','LocalAmplitude','IsBlue'};
peaktroughMat=array2table(peaktroughMat,'VariableNames',FieldName);

figure(2); clf;
imagesc(normTr(Result.dist_order,:),[-2 5]);
hold all
%plot(find(Result.spike(1,:))-1,4,'ro')
plot(peaktroughMat.frame(peaktroughMat.IsPeak>0 & ismember(peaktroughMat.frame,find(Result.Blue==0)))-2,5,'ro','LineWidth',2) %Peak
plot(peaktroughMat.frame(peaktroughMat.IsPeak==0 & ismember(peaktroughMat.frame,find(Result.Blue>0)))-3,4.5,'ko','LineWidth',2) %Trough
plot(Result.Blue-4,'b')
colormap(turbo);
%% Filter movie
% nTau_EXCINHpatt=[50];
% strImg=imresize(Result.Structure,0.7);
% binImg=imresize(highpassDendrite(Result.ref_im,Result.SWC)>5,0.7);

% BlueOn_silent=setdiff(find(Result.Blue>0),periSpikeTime);
% BlueOndF=median(mov_res_sub_filt(:,:,BlueOn_silent),3,'omitnan');
% inhibitory_frame=peaktroughMat.frame(peaktroughMat.IsPeak==0 & peaktroughMat.IsBlue>0);
% excitatory_frame=peaktroughMat.frame(peaktroughMat.IsPeak>0 & peaktroughMat.IsBlue==0);
% 
% [~, inhibitory_dFMat, inhibitory_frame]=get_STA(tovec(BlueOndF-mov_res_sub_filt),inhibitory_frame,nTau_EXCINHpatt,nTau_EXCINHpatt);
% inhibitory_dFMat=reshape(permute(inhibitory_dFMat,[1 3 2]),sz2(1),sz2(2),[]);
% [~, excitatory_dFMat, excitatory_frame]=get_STA(tovec(mov_res_sub_filt),excitatory_frame,nTau_EXCINHpatt,nTau_EXCINHpatt);
% excitatory_dFMat=reshape(permute(excitatory_dFMat,[1 3 2]),sz2(1),sz2(2),[]);
% 
% [~,~,u_all] = get_eigvector(tovec(mov_res_sub_filt(:,:,DirtFrame==0))',40);
% [~,~,u_exc] = get_eigvector(tovec(imgaussfilt(excitatory_dFMat,1))',40);
% [~,~,u_inh] = get_eigvector(tovec(imgaussfilt(inhibitory_dFMat,1))',40);
% u_all=reshape(u_all,sz2(1),sz2(2),[]); 
% u_inh=reshape(u_inh,sz2(1),sz2(2),[]); u_exc=reshape(u_exc,sz2(1),sz2(2),[]);

% template2use=cat(3,u_all(:,:,1:2),u_inh(:,:,1:2),u_exc(:,:,1:2));
% mov_res_sub_filt2=pcafilt_template(mov_res_sub_filt,template2use);

%F0img=std(mov_res_sub_filt2,0,3);
F0img=get_F0img_PCA(mov_res_sub_filt);
% F0img2=get_F0img(toimg(mov_res,sz2(1),sz2(2)));
%F0img=sqrt(mean(mov_res_sub_filt(:,:,1:end-1).*mov_res_sub_filt(:,:,2:end),3));

%%
nTau_EXCINHpatt=[50];
strImg=imresize(Result.Structure,0.7);
binImg=imresize(highpassDendrite(Result.ref_im,Result.SWC)>5,0.7);
cmap=gen_colormap([0 0.2 1; 0 0 0; 1 0 0],256);

BlueOn_silent=setdiff(find(Result.Blue>0),periSpikeTime);
BlueOndF=prctile(mov_res_sub_filt(:,:,BlueOn_silent),50,3);
inhibitory_frame=peaktroughMat.frame(peaktroughMat.IsPeak==0 & peaktroughMat.IsBlue>0);
excitatory_frame=peaktroughMat.frame(peaktroughMat.IsPeak>0 & peaktroughMat.IsBlue==0);

[~, inhibitory_dFMat, inhibitory_frame]=get_STA(tovec(BlueOndF-mov_res_sub_filt),inhibitory_frame,nTau_EXCINHpatt,nTau_EXCINHpatt);
%inhibitory_dFMat=prctile(inhibitory_dFMat,80,3)-inhibitory_dFMat;
[~, inhibitory_dirtMat]=get_STA(DirtFrame',inhibitory_frame,nTau_EXCINHpatt,nTau_EXCINHpatt);
inhibitory_dirtMat=squeeze(inhibitory_dirtMat);
inhInd2use=find(sum(inhibitory_dirtMat(:,nTau_EXCINHpatt+1+[-3:3]),2)<100000);
inhibitory_frame=inhibitory_frame(inhInd2use);
[inhibitory_frame isort]=sort(inhibitory_frame);
%inhibitory_dF=reshape(permute(inhibitory_dF,[1 3 2]),sz(2),sz(1),[]);

[~, excitatory_dFMat, excitatory_frame]=get_STA(tovec(mov_res_sub_filt),excitatory_frame,nTau_EXCINHpatt,nTau_EXCINHpatt);
%excitatory_dFMat=excitatory_dFMat-prctile(excitatory_dFMat,30,3);
[~, excitatory_dirtMat]=get_STA(DirtFrame',excitatory_frame,nTau_EXCINHpatt,nTau_EXCINHpatt);
excitatory_dirtMat=squeeze(excitatory_dirtMat);
excInd2use=find(sum(excitatory_dirtMat(:,nTau_EXCINHpatt+1+[-3:3]),2)<100000);
excitatory_frame=excitatory_frame(excInd2use);
[excitatory_frame esort]=sort(excitatory_frame);
%excitatory_dF=reshape(permute(excitatory_dF,[1 3 2]),sz(2),sz(1),[]);

inhibitory_dF=reshape(mean(inhibitory_dFMat(:,inhInd2use,nTau_EXCINHpatt+1+[-3:3]),3),sz2(1),sz2(2),[]);
excitatory_dF=reshape(mean(excitatory_dFMat(:,excInd2use,nTau_EXCINHpatt+1+[-3:3]),3),sz2(1),sz2(2),[]);

inhibitory_dF=inhibitory_dF(:,:,isort);
excitatory_dF=excitatory_dF(:,:,esort);
% [~, DInh, eigImgInh]=get_eigvector(tovec(imgaussfilt(inhibitory_dF,1))',10);
% [~, DExc, eigImgExc]=get_eigvector(tovec(imgaussfilt(excitatory_dF,1))',10);
% [V D E]=get_eigvector(tovec(mov_res_sub)',10);

% Show footprints
F0img=maskEdge(F0img,5,NaN);
F0img_filt=medfilt2nan(F0img,[3 3]); 

Ftprnts_resz=imresize(Result.ftprnt,0.7);
inhibitory_dF_norm=tovec(maskEdge(inhibitory_dF,5,NaN)); 
inhibitory_dF_norm(tovec(imresize(Result.Structure_bin,0.7)==0),:)=NaN;
inhibitory_dF_norm=toimg(inhibitory_dF_norm,sz2(1),sz2(2));
inhibitory_dF_norm=(inhibitory_dF_norm)./F0img_filt;
inhibitory_dF_norm_colored=[]; excitatory_dF_norm_colored=[];
dFF_range_inh=[-1 1]*2.5; dFF_range_exc=[-1 1]*3.5;
for p=1:size(inhibitory_dF_norm,3)
    p_tmp=inhibitory_dF_norm(:,:,p);
    p_tmp=medfilt2nan(p_tmp,[15 15]);
    inhibitory_dF_norm(:,:,p)=-maskEdge(imgaussfilt_NaN(p_tmp,2),5,NaN);
    inhibitory_dF_norm_colored(:,:,:,p)= grs2rgb(double(inhibitory_dF_norm(:,:,p).*binImg), cmap ,dFF_range_inh(1),dFF_range_inh(2)).*(strImg>0.01);
    inhibitory_dF_norm_colored(:,:,:,p) = grs2rgb(double(strImg), colormap("gray")) + inhibitory_dF_norm_colored(:,:,:,p);
end
inhibitory_dF_norm_colored=maskEdge(inhibitory_dF_norm_colored,5,NaN);

excitatory_dF_norm=tovec(maskEdge(excitatory_dF,5,NaN)); 
excitatory_dF_norm(tovec(imresize(Result.Structure_bin,0.7)==0),:)=NaN;
excitatory_dF_norm=toimg(excitatory_dF_norm,sz2(1),sz2(2));
excitatory_dF_norm=(excitatory_dF_norm)./F0img_filt;
for p=1:size(excitatory_dF_norm,3)
    p_tmp=excitatory_dF_norm(:,:,p);
    p_tmp=medfilt2nan(p_tmp,[5 5]);
    excitatory_dF_norm(:,:,p)=maskEdge(imgaussfilt_NaN(p_tmp,2),5,NaN);
    excitatory_dF_norm_colored(:,:,:,p)= grs2rgb(double(excitatory_dF_norm(:,:,p).*binImg), cmap ,dFF_range_exc(1),dFF_range_exc(2)).*(strImg>0.01);
    excitatory_dF_norm_colored(:,:,:,p) = grs2rgb(double(strImg), colormap("gray")) + excitatory_dF_norm_colored(:,:,:,p);
end
excitatory_dF_norm_colored=maskEdge(excitatory_dF_norm_colored,5,NaN);

figure(22); clf;
for p=1:size(excitatory_dF_norm_colored,4)
nexttile([1 1]);
imshow2(excitatory_dF_norm_colored(:,:,:,p),[]); title([num2str(excitatory_frame(p)),', p= ' num2str(p)]);
end
sgtitle('Excitatory footprints');
figure(23); clf;
for p=1:size(inhibitory_dF_norm_colored,4)
nexttile([1 1]);
imshow2(inhibitory_dF_norm_colored(:,:,:,p),[]); title([num2str(inhibitory_frame(p)),', p= ' num2str(p)]);
end
sgtitle('Inhibitory footprints');
%% show footprints (Figure S)

showpatt=[2 7]; xrange2show=[10:365]; somacoord=get_coord(Ftprnts_resz(:,:,1));
dax=([1:size(inhibitory_dF_norm,2)]-somacoord(1))*PixelSize(f);
offset=0; cmap_tr=gen_colormap(Plasma,4); cmap_sub=cmap_tr./3*2;
g=1; INH_str=[];
BlueOnSteadyLevel=prctile(subTr_label_silent(:,Result.Blue>0),40,2);
BlueOffSteadyLevel=prctile(subTr_label_silent(:,Result.Blue==0),50,2);
figure(55); clf; tiledlayout(1,length(showpatt)+1,'Padding','compact'); 
t2show=[-100:100]; ax2=[]; catfprnt=[];
szfprnt=size(inhibitory_dF_norm_colored);
for p=showpatt;%1:size(inhibitory_dF_norm_colored,4)
    catfprnt=[catfprnt ones(szfprnt(1),50,3) fliplr(inhibitory_dF_norm_colored(:,xrange2show,:,p))];
    INH_str{g}=['INH. #' num2str(g)];

    ax2=[ax2 nexttile(g,[1 1])];
    l_sub=plot(t2show,subTr_label(:,inhibitory_frame(p)+t2show)'+[1:4]*offset,'linewidth',1.5); hold all
    arrayfun(@(l,c) set(l,'Color',c{:}),l_sub,num2cell(cmap_tr,2));
    axis off;

    nexttile((length(showpatt)+1),[1 1]);
    inhfprnt=inhibitory_dF_norm(:,:,p); inhfprnt(strImg==0)=NaN;
    plot(-dax(xrange2show),mean(inhfprnt(:,xrange2show),1,'omitnan'),'linewidth',1.5); box off; axis tight; hold all;
    xlabel('Distance (\mum)'); ylabel('Mean subthreshold (Z score)');
    g=g+1;
end
nexttile(length(showpatt),[1 1]);
drawScaleBar(10,'horizontal','color',[0 0 0],'Linewidth',3,'Position',[47 1.4]);
drawScaleBar(1,'vertical','color',[0 0 0],'Linewidth',3,'Position',[47 0.4]);
linkaxes(ax2); ylim([-0.6 2.5]);
nexttile((length(showpatt)+1),[1 1]); legend(INH_str,'Box','off','Location','best');
set_fontsize(13);
figure(65); clf;
nexttile(1,[1 length(showpatt)]);
imshow2(catfprnt,[]); 
cb=colorbar; cb.Label.String='Z score'; colormap(gen_colormap([0 0.5 1;1 1 1;1 0 0],256)); cb.Ticks=[0 1]; cb.TickLabels=dFF_range_inh;
set_fontsize(15);

showpatt2=[7 10];
offset=0; g=1; EXC_str=[];
figure(56); clf; tiledlayout(1,length(showpatt2)+1,'Padding','compact'); 
t2show=[-50:50]; ax2=[]; catfprnt=[];
szfprnt=size(excitatory_dF_norm_colored);
for p=showpatt2;%1:size(inhibitory_dF_norm_colored,4)
    catfprnt=[catfprnt ones(szfprnt(1),50,3) fliplr(excitatory_dF_norm_colored(:,xrange2show,:,p))];
    EXC_str{g}=['EXC. #' num2str(g)];

    ax2=[ax2 nexttile(g,[1 1])];
    l_sub=plot(t2show,subTr_label(:,excitatory_frame(p)+t2show)'+[1:4]*offset,'linewidth',1.5); hold all
    arrayfun(@(l,c) set(l,'Color',c{:}),l_sub,num2cell(cmap_tr,2));
    axis off;

    nexttile((length(showpatt2)+1),[1 1]);
    inhfprnt=excitatory_dF_norm(:,:,p); inhfprnt(strImg==0)=NaN;
    plot(-dax(xrange2show),mean(inhfprnt(:,xrange2show),1,'omitnan'),'linewidth',1.5); box off; axis tight; hold all;
    xlabel('Distance (\mum)'); ylabel('Mean subthreshold (Z score)');
    g=g+1;
end
nexttile(length(showpatt2),[1 1]);
drawScaleBar(10,'horizontal','color',[0 0 0],'Linewidth',3,'Position',[47 1.4]);
drawScaleBar(1,'vertical','color',[0 0 0],'Linewidth',3,'Position',[47 0.4]);
linkaxes(ax2); ylim([-0.8 2.6]);
nexttile((length(showpatt2)+1),[1 1]); legend(EXC_str,'Box','off','Location','best');
nexttile(1,[1 1]); legend(({'Basal','Soma','Apical','Distal'}),'Box','off','Location','best','NumColumns',2);
set_fontsize(13);
figure(66); clf;
nexttile(1,[1 length(showpatt2)]);
imshow2(catfprnt,[]); 
cb=colorbar; cb.Label.String='Z score'; colormap(gen_colormap([0 0.5 1;1 1 1;1 0 0],256)); cb.Ticks=[0 1]; cb.TickLabels=dFF_range_exc;
set_fontsize(15);

figure(57); clf; tiledlayout(1,4);
offset=6; cmap_tr=gen_colormap(Plasma,4); cmap_sub=cmap_tr./2;
Ftprnt_lab=toimg(tovec(Result.ftprnt>0)*labelvec',sz(2),sz(1));
t=[1:nTime];
SWC=Result.SWC;
SWC(:,4)=5; SWC(1,4)=30;
% nexttile([1 1]);
% showScaleScatter([1:4], SWC, Ftprnt_lab ,gen_colormap(Plasma,256),[1 4]);
% drawScaleBar(100/PixelSize(f),'vertical')
nexttile([1 4]);
l=plot(t,fliplr(normTr_label')+[1:4]*offset); hold all
l_sub=plot(t,fliplr(subTr_label')+[1:4]*offset);
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(flipud(cmap_tr),2));
arrayfun(@(l,c) set(l,'Color',c{:}),l_sub,num2cell(flipud(cmap_sub),2));
hold all
%plot(find(Result.spike(1,:)),offset*4.8,'k.')
plot([excitatory_frame(showpatt2); excitatory_frame(showpatt2)],[5 offset*5],'r--','linewidth',1.5)
plot([inhibitory_frame(showpatt); inhibitory_frame(showpatt)],[5 offset*5],'b--','linewidth',1.5)
plot(Result.Blue*5,'color',[0 0.5 1],'linewidth',1.5)
axis off; xlim([100 9999]);
drawScaleBar(1000,'Horizontal')