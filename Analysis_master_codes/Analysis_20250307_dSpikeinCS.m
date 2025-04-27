
clear
clc;
cd '/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/Statistics_Optopatch_Prism';
[~, ~, raw] = xlsread(['/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/' ...
    'Prism_OptopatchData_Arrangement.xlsx'], 'Sheet1', 'C5:Q175');

save_dir='/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Figures/invivoPrism/FigureOptopatch';
fpath=raw(:,1);
Mouse=cell2mat(raw(:,2));
NeuronInd=cell2mat(raw(:,5));
CamType=raw(:,3);
StructureData=raw(:,10);
StimROI=raw(:,6);
StimWfn=raw(:,7);
isGoodCell=cell2mat(raw(:,11));
PixelSize=cell2mat(cellfun(@(x) (str2num(num2str(x))),raw(:,12),'UniformOutput',false));
refROI=cellfun(@(x) (str2num(num2str(x))),raw(:,14),'UniformOutput',false);
maintrunkROI=cellfun(@(x) (str2num(num2str(x))),raw(:,15),'UniformOutput',false);
place_bin=150; time_segment=15000; overlap=200;
alignedMovFN = {'STA_Mat_SS','STA_Mat_CS','STA_Mat_dSP'};
bound=6;
title_str={'Basal','Apical','Peri-Soma'};
[~, unqInd] = unique([Mouse NeuronInd] ,'row');
set(0,'DefaultFigureWindowStyle','docked')
% %% get F0_img
% 
% bound=7; i=54;
% SameCellInd=find(Mouse==Mouse(i) & NeuronInd==NeuronInd(i));
% mov_res_small=[];
% for j=SameCellInd'
% j
%     load([fpath{j} '/OP_Result.mat'])
% 
%     mov_mc=readBinMov_BHL(fpath{j},3);
%     mean_F=squeeze(mean(mov_mc(bound:end-bound,bound:end-bound,:),[1 2]));
%     [~, blueOff]=get_blueoffTrace(mean_F,[Result.Blue],70);
%     [y_fit]=expfitDM_2(find(blueOff)',mean_F(find(blueOff)),[1:size(mov_mc,3)]',1000);
% 
%     mov_res= mov_mc-mean(mov_mc,3);
%     mov_res = SeeResiduals(mov_res,y_fit);
%     mov_res = SeeResiduals(mov_res,Result.mc(:,:));
%     mov_res = SeeResiduals(mov_res,Result.mc(:,:).^2);
%     mov_res = SeeResiduals(mov_res,Result.mc(:,1).*Result.mc(:,end));
% 
%     if isfield(Result,'dirtTrace')
%     validFrm=find(sum(isnan(Result.dirtTrace),1)==0);
%     else
%     validFrm=[1:size(mov_res,3)];
%     end
%     mov_res=movmean(mov_res,15,3);
%     mov_res=mov_res(:,:,validFrm);
%     mov_res_small=cat(3,mov_res_small,imresize(mov_res,1/4));
% end
% 
% [V D]=get_eigvector(tovec(mov_res_small),10);
% F0img=sqrt(sum((V.^2).*D(1:10)',2));
% F0img=toimg(F0img,size(mov_res_small,1),size(mov_res_small,2));
% F0img=imresize(F0img,4);   
% 
% for j=SameCellInd'
%     load([fpath{j} '/OP_Result.mat'])
%     Result.F0PCAimg=F0img;
%     save([fpath{j} '/OP_Result.mat'],'Result','-v7.3')
% end


%%
figure(21); clf; tiledlayout(11,1); ax1=[]; ax3=[];
dendriteaxis_bin=[-100:25:300]; cax=[-0.009 0.025]; t_show=[2.25 12.25]; cmap=[194 102 56; 118 85 157;]/256;

blue_cat=[]; normTr_cat=[]; spike_cat=[];
for i=[73 75];
ax1=[ax1 nexttile([2 1])];
load(fullfile(fpath{i},'OP_Result.mat'))
taxis=[1:size(Result.normTraces,2)]/1000;

%noi=maintrunkROI{i};
noi=[1:size(Result.ftprnt,3)];
noi_dist=ismember(Result.dist_order,noi);

Dsign=ones(1,size(Result.interDendDist,2));
Dsign(Result.dist_order(1:find(Result.dist_order==1)-1))=-1;
dendaxis=Result.interDendDist(1,Result.dist_order(noi_dist)).*Dsign(Result.dist_order(noi_dist));

normResidue=SeeResiduals(permute(Result.normTraces,[3 1 2]),Result.mc);
normResidue=SeeResiduals(permute(Result.normTraces,[3 1 2]),Result.mc.^2);
normResidue=SeeResiduals(permute(Result.normTraces,[3 1 2]),Result.mc(:,1).*Result.mc(:,end));
normResidue=squeeze(normResidue);
%F0PCA=get_F0PCA(normResidue);

normTr=normResidue./Result.F0_PCA;
normTr=pcafilterTrace(normTr,[1:7]);
normTr=normTr(Result.dist_order(noi_dist),:);
normTr_cat=[normTr_cat normTr];
blue_cat=[blue_cat Result.Blue];
spike_cat=[spike_cat Result.spike];

[~, sortind]=sort(dendaxis,'ascend');
[Xq, Yq] = meshgrid([taxis], dendriteaxis_bin);
normTr_interp = interp2([taxis], dendaxis(sortind), normTr(sortind,:), Xq, Yq, 'linear');
pcolor(taxis,dendriteaxis_bin,normTr_interp);
caxis(cax);
colormap(turbo);
shading interp
cb=colorbar;
cb.Label.String='\DeltaF/F';
set(gca,'XTick',t_show,'XTickLabel',t_show-t_show(1))
ylim([0 220])

ax3=[ax3 nexttile([3 1])];
plot(taxis,mean(normTr(1,:),1,'omitnan'),'color',cmap(1,:)); hold all
plot(taxis,mean(normTr(end,:),1,'omitnan'),'color',cmap(2,:));
ylim([-0.015 0.03])
axis off
end

ax2=nexttile([1 1]);
plot(taxis, Result.Blue)
linkaxes([ax1 ax2 ax3],'x')
xlim(t_show)
axis off

%%
figure(22); clf; tiledlayout(2,1)
ax1=nexttile([1 1]);
imagesc(normTr_cat,[-0.005 0.02]);
colormap(turbo);
ax2=nexttile([1 1]);
imagesc(pcafilterTrace(normTr_cat,[1:2]),[-0.005 0.02]);
colormap(turbo);
linkaxes([ax1 ax2])

bwBlue=bwlabel(blue_cat);
spvec=spike_cat;

[~, normTrMat spind]=get_STA(normTr_cat(1,:),spike_cat(1,:),1,1);
[~, shift]=max(normTrMat,[],3);
spTime=spind+shift-2;
spvec_shifted=ind2vec(size(normTrMat,2),spTime,1);

BlueonSpike=[];
for b=1:max(bwBlue)-1
    sptime_tmp=find(spvec_shifted(1,find(bwBlue==b,1)+[0:50])>0,1);
    BlueonSpike=[BlueonSpike sptime_tmp+find(bwBlue==b,1)-1];
end
BlueonSpike_vec=ind2vec(size(normTr_cat,2),BlueonSpike,1);
[~, BlueonsetSpikeMat]=get_STA(normTr_cat,BlueonSpike_vec,20,50);

normTr_cat_hi=normTr_cat-movmedian(normTr_cat,250,2);
[V D eigTrace]=get_eigvector(normTr_cat_hi,10);
[~, eigTraceMat]=get_STA(eigTrace(:,1:2)',BlueonSpike_vec,3,50);
figure(23); clf;
tiledlayout(2,1);
ax1=nexttile([1 1]);
imagesc(normTr_cat,[cax])
colormap(turbo);
ax2=nexttile([1 1]);
plot(eigTrace(:,1:2))
linkaxes([ax1 ax2],'x')

spikewdSpike=find(max(squeeze(eigTraceMat(2,:,:))'<-0.01 & squeeze(eigTraceMat(1,:,:))'>0.065,[],1)>0);
dSpike=find_spike_bh(normTr_cat_hi(31,:),0.015,0.001);
blue_cat_dilate=ind2vec(size(normTr_cat,2),tovec(BlueonSpike'+[0:100]),1);
dSpike=dSpike.*blue_cat_dilate;
[~, dSpikeMat spind]=get_STA(normTr_cat,dSpike,20,20);

%%
mov_res_cat=[];
for i=[73 75]
    load(fullfile(fpath{i},'OP_Result.mat'))
    mov=readBinMov_BHL(fpath{i});
    mov_res= mov-mean(mov,3);
    bkg = zeros(2, size(mov,3));
    bkg(1,:) = linspace(-1, 1, size(mov,3));  % linear term
    bkg(2,:) = linspace(-1, 1, size(mov,3)).^2;  % quadratic term
    mov_res = SeeResiduals(mov_res,Result.mc);
    mov_res = SeeResiduals(mov_res,Result.mc.^2);
    mov_res = SeeResiduals(mov_res,Result.mc(:,1).*Result.mc(:,end));
    mov_res= SeeResiduals(mov_res,bkg,1);
    mov_res_cat=cat(3,mov_res_cat,imresize(mov_res,1/2));
end

dSpikeTAmovie=get_STA(tovec(mov_res_cat-movmedian(mov_res_cat,250,3)),dSpike,100,30);
dSpikeTAmovie=-toimg(dSpikeTAmovie,size(mov_res_cat,1),size(mov_res_cat,2));
dSpikeTAmovie=dSpikeTAmovie-mean(dSpikeTAmovie(:,:,1:10),3);

filevec=zeros(1,size(blue_cat,2)); filevec(1:19999)=1;
STAmovie=get_STA(tovec(mov_res_cat-movmedian(mov_res_cat,250,3)),BlueonSpike_vec.*filevec,100,30);
STAmovie=-toimg(STAmovie,size(mov_res_cat,1),size(mov_res_cat,2));
STAmovie=STAmovie-mean(STAmovie(:,:,1:10),3);

dSpikeTAmovie_norm=imresize(dSpikeTAmovie,2)./Result.F0PCAimg;
STAmovie_norm=imresize(STAmovie,2)./Result.F0PCAimg;

revertedStruct= Result.Structure;
revertedStruct(revertedStruct==0)=prctile(revertedStruct(:),30);
revertedStruct=mat2gray(revertedStruct);

crop_roi=[94 6  300  159];
cax_sub=[-0.005 0.03];
dSpikeTAmovie_norm_filt=imgaussfilt(dSpikeTAmovie_norm,2).*(Result.Structure_bin>0.01);
%dSpikeTAmovie_norm_filt=pcafilt(dSpikeTAmovie_norm_filt,10);
colorSTA=grs2rgb(tovec(dSpikeTAmovie_norm_filt),colormap('turbo'),cax_sub(1),cax_sub(2));
colorSTA=reshape(colorSTA,size(dSpikeTAmovie_norm_filt,1),size(dSpikeTAmovie_norm_filt,2),size(dSpikeTAmovie_norm_filt,3),[]);
colorSTA=permute(colorSTA,[1 2 4 3]);

dSpikemov_colored=colorSTA.*revertedStruct*3;%.*SkelDend(bound:end-bound,bound:end-bound);
writeMov4d(fullfile(save_dir,['dSpikeTA_in_CS']),dSpikemov_colored(crop_roi(2):crop_roi(2)+crop_roi(4),crop_roi(1):crop_roi(1)+crop_roi(3),:,50:end),[-51:30],10,1,cax_sub)

STAmovie_norm_filt=imgaussfilt(STAmovie_norm,2).*(Result.Structure_bin>0.01);
colorSTA2=grs2rgb(tovec(STAmovie_norm_filt),colormap('turbo'),cax_sub(1),cax_sub(2));
colorSTA2=reshape(colorSTA2,size(STAmovie_norm_filt,1),size(STAmovie_norm_filt,2),size(STAmovie_norm_filt,3),[]);
colorSTA2=permute(colorSTA2,[1 2 4 3]);
STAmov_colored=colorSTA2.*revertedStruct*3;%.*SkelDend(bound:end-bound,bound:end-bound);

writeMov4d(fullfile(save_dir,['STAmov']),STAmov_colored(crop_roi(2):crop_roi(2)+crop_roi(4),crop_roi(1):crop_roi(1)+crop_roi(3),:,50:end),[-51:30],10,1,cax_sub)

figure(24); clf;
imshow2([STAmov_colored(crop_roi(2):crop_roi(2)+crop_roi(4),crop_roi(1):crop_roi(1)+crop_roi(3),:,101); ...
         dSpikemov_colored(crop_roi(2):crop_roi(2)+crop_roi(4),crop_roi(1):crop_roi(1)+crop_roi(3),:,101)])
cb=colorbar;
colormap(turbo);
cb.Label.String='\DeltaF/F';
cb.Ticks=[0 1];
cb.TickLabels=cax_sub;

figure(25); clf;
[~, SpMat]=get_STA(spvec_shifted,ind2vec(size(normTr_cat,2),BlueonSpike,1),20,100);
[~, dSPMat]=get_STA(dSpike,ind2vec(size(normTr_cat,2),BlueonSpike,1),20,100);
SpMatsom=squeeze(SpMat(1,1:find(BlueonSpike>19999,1)-1,:));
SpMatdend=squeeze(SpMat(1,find(BlueonSpike>19999,1):end,:));
dSpMatsom=squeeze(dSPMat(1,1:find(BlueonSpike>19999,1)-1,:));
dSpMatdend=squeeze(dSPMat(1,find(BlueonSpike>19999,1):end,:));
[~, SptimeSom]=find(SpMatsom>0);
[~, SptimeDen]=find(SpMatdend>0);
[~, dSptimeSom]=find(dSpMatsom>0);
[~, dSptimeDen]=find(dSpMatdend>0);

ax2=nexttile([1 1]); t_edege=[-10:2:60];
histogram(SptimeSom-20,t_edege,'FaceColor',[0.6 0.6 0.6]); hold all
histogram(dSptimeSom-20,t_edege,'FaceColor',[1 0 0]); hold all
xlabel('Time from first spike (ms)')
ylabel('# of spike')
title('Soma stimulation')
legend({'Somatic spike','dSpike'})
ax1=nexttile([1 1]);
histogram(SptimeDen-20,t_edege,'FaceColor',[0.6 0.6 0.6]); hold all
histogram(dSptimeDen-20,t_edege,'FaceColor',[1 0 0]); hold all
xlabel('Time from first spike (ms)')
ylabel('# of spike')
title('Distal dendrite stimulation')
legend({'Somatic spike','dSpike'})
linkaxes([ax1 ax2])
%% Subthreshold figure;
bound=6;
STAmovie=[]; g=1; BlueonSpike=[]; BlueBoundary=[];
for i=[73 75]
    load(fullfile(fpath{i},'OP_Result.mat'))
    mov=readBinMov_BHL(fpath{i});
    mov_res= mov-mean(mov,3);
    bkg = zeros(2, size(mov,3));
    bkg(1,:) = linspace(-1, 1, size(mov,3));  % linear term
    bkg(2,:) = linspace(-1, 1, size(mov,3)).^2;  % quadratic term
    mov_res = SeeResiduals(mov_res,Result.mc);
    mov_res = SeeResiduals(mov_res,Result.mc.^2);
    mov_res = SeeResiduals(mov_res,Result.mc(:,1).*Result.mc(:,end));
    mov_res= SeeResiduals(mov_res,bkg,1);

    bwBlue=bwlabel(Result.Blue);
    spvec=Result.spike;
    [~, normTrMat spind]=get_STA(Result.normTraces(1,:),Result.spike(1,:),1,1);
    [~, shift]=max(normTrMat,[],3);
    spTime=spind+shift-2;
    spvec_shifted=ind2vec(size(Result.normTraces,2),spTime,1);

    BlueonSpike=[];
    for b=1:max(bwBlue)
        sptime_tmp=find(spvec_shifted(1,find(bwBlue==b,1)+[0:50])>0,1);
        BlueonSpike=[BlueonSpike sptime_tmp+find(bwBlue==b,1)-1];
    end

    statmp = get_STA(tovec(mov_res),ind2vec(size(mov_res,3),BlueonSpike,1),50,50);
    STAmovie(:,:,:,g)=toimg(statmp,size(mov_res,1),size(mov_res,2));

    BlueBoundary{g}=(Result.blueDMDimg);
    g=g+1;
end

%mov_filt=imresize(pcafilt(imresize(mov_res,1/5),3),5);
%lag1mean=abs(sqrt(mean(mov_filt(:,:,2:end).*mov_filt(:,:,1:end-1),3)));
mov_res_shrink=imresize(movmedian(mov_res,10,3),1/5);
[V D]=get_eigvector(tovec(mov_res_shrink),10);
F0img=sqrt(sum((V.^2).*D(1:10)',2));
F0img=toimg(F0img,size(mov_res_shrink,1),size(mov_res_shrink,2));
F0img=imresize(F0img,5);

STAmov_norm=-STAmovie./F0img;
STAmov_norm=STAmov_norm-mean(STAmov_norm(:,:,[1:15],:),3);

Rfixed = imref2d([size(Result.ref_im,1) size(Result.ref_im,2)]);
inverseTform = invert(Result.tform);
%revertedStruct = imwarp(Result.Structure, inverseTform,'OutputView',Rfixed);
revertedStruct= Result.Structure;
revertedStruct(revertedStruct==0)=prctile(revertedStruct(:),30);
revertedStruct=mat2gray(revertedStruct);
%%
STAmovie_normStr=[];
crop_roi=[94 6  300  159];
STAmov_norm_crop=STAmov_norm(crop_roi(2):crop_roi(2)+crop_roi(4),crop_roi(1):crop_roi(1)+crop_roi(3),:,:);
revertedStruct_crop=revertedStruct(crop_roi(2):crop_roi(2)+crop_roi(4),crop_roi(1):crop_roi(1)+crop_roi(3));
cax_sub=[0.004 0.01];
cax_sp=[0.004 0.022];
STAmov_norm_crop_filt=[];
for spclass_ind=1:2
    %STAmovie_norm{spclass_ind}=imgaussfilt3(STAnorm_sub./F_refImg,[1.5 1.5 0.1]);%.*SkelDend(bound:end-bound,bound:end-bound);
    STAmov_norm_crop_filt{spclass_ind}=imgaussfilt(pcafilt(STAmov_norm_crop(:,:,:,spclass_ind),15),4);
    colorSTA=grs2rgb(tovec(STAmov_norm_crop_filt{spclass_ind}),colormap('turbo'),cax_sub(1),cax_sub(2));
    colorSTA=reshape(colorSTA,size(STAmov_norm_crop,1),size(STAmov_norm_crop,2),size(STAmov_norm_crop,3),[]);
    colorSTA=permute(colorSTA,[1 2 4 3]);

    colorSTA2=grs2rgb(tovec(STAmov_norm_crop_filt{spclass_ind}),colormap('turbo'),cax_sp(1),cax_sp(2));
    colorSTA2=reshape(colorSTA2,size(STAmov_norm_crop,1),size(STAmov_norm_crop,2),size(STAmov_norm_crop,3),[]);
    colorSTA2=permute(colorSTA2,[1 2 4 3]);

    STAmovie_normStr{spclass_ind}=colorSTA.*revertedStruct_crop*3;%.*SkelDend(bound:end-bound,bound:end-bound);
    STAmovie_normStr2{spclass_ind}=colorSTA2.*revertedStruct_crop*3;%.*SkelDend(bound:end-bound,bound:end-bound);
end

% sptype={'SomStim','ddStim'};
% cax=[-0.005,0.02];
% writeMov4d(fullfile(save_dir,['STA_dFFStructgrsrgb_ShortPulse']),[imrotate(STAmovie_normStr{1},90) imrotate(STAmovie_normStr{2},90)],[-50:50],10,1,cax)

figure(21); clf; tiledlayout(1,6)
t_show=[41:50]; t_show_spike=[51:53]; ax1=[];
for stimROI=1:2
subimage=grs2rgb(mean(STAmov_norm_crop_filt{stimROI}(:,:,t_show),3),colormap(turbo),cax_sub(1),cax_sub(2));
subimage=subimage.*revertedStruct_crop*3;
spimage=grs2rgb(max(STAmov_norm_crop_filt{stimROI}(:,:,t_show_spike),[],3),colormap(turbo),cax_sp(1),cax_sp(2));
spimage=spimage.*revertedStruct_crop*3;

ax1=[ax1 nexttile([1 1])];
imshow2(imrotate([Result.Structure(crop_roi(2):crop_roi(2)+crop_roi(4),crop_roi(1):crop_roi(1)+crop_roi(3))],90),[]); hold all
Blbd=bwboundaries(imrotate(BlueBoundary{stimROI}(crop_roi(2):crop_roi(2)+crop_roi(4),crop_roi(1):crop_roi(1)+crop_roi(3)),90));
plot(Blbd{1}(:,2),Blbd{1}(:,1),'color',[0 0.6 1])
nexttile([1 1]);
imshow2(imrotate([subimage],90),[]); hold all
cb=colorbar;
cb.Ticks=[0 0.5 1];
cb.TickLabels=[cax_sub(1) mean(cax_sub) cax_sub(2)];
cb.Label.String = '\DeltaF/F';

nexttile([1 1]);
imshow2(imrotate([spimage],90),[]); hold all
cb=colorbar;
cb.Ticks=[0 0.5 1];
cb.TickLabels=[cax_sp(1) mean(cax_sp) cax_sp(2)];
cb.Label.String = '\DeltaF/F';
%plot(BlueBoundary{1}(:,2),BlueBoundary{1}(:,1),'color',[0 0.5 1])
end
colormap(turbo);
colormap(ax1(1),gray);
colormap(ax1(2),gray);  

%% Concatenate data
StimROI_Ind={'Soma','Distal Dend','WF'};
StimWfn_Ind={'Ramp Stim','Short Pulse'};
CatResult=[];
%foi_somDend=[1 22 14 18 10 26 25 32 46];
foi_somDend=[1 22 14 18 10 26 25 32 46 47 48 43 34 37 35];
g2=1;
for i=unqInd(foi_somDend)' %Rois

    SameCellInd=find(Mouse==Mouse(i) & NeuronInd==NeuronInd(i));
    isSoma = ~cellfun(@isempty,strfind(StimROI(SameCellInd), StimROI_Ind{1}));
    isDend = ~cellfun(@isempty,strfind(StimROI(SameCellInd), StimROI_Ind{2}));
    isWF = ~cellfun(@isempty,strfind(StimROI(SameCellInd), StimROI_Ind{3}));
    isRamp = ~cellfun(@isempty,strfind(StimWfn(SameCellInd), StimWfn_Ind{1}));
    isSP   = ~cellfun(@isempty,strfind(StimWfn(SameCellInd), StimWfn_Ind{2}));

    ROIwvf_ind=[isSoma isDend isWF isRamp isSP];
    validind=find(sum(ROIwvf_ind,2)>=2 & isGoodCell(SameCellInd));
    patterns = [1 0 0 1 0; 1 0 0 0 1; 0 1 0 1 0; 0 1 0 0 1; 0 0 1 1 0; 0 0 1 0 1];
    values = [1; 2; 3; 4; 5; 6]; % 1=soma, ramp; 2=soma, sp; 3=dend, ramp; 4=dend, sp; 5=WF, ramp; 6=WF, sp;

    g=ones(1,6);
    for j=1:length(validind)
        f2read=SameCellInd(validind(j))
        load(fullfile(fpath{f2read},'OP_Result.mat'),'Result');
        wfn = values(find(ismember(patterns, ROIwvf_ind(validind(j),:), 'rows')));
        CatResult{wfn,g(wfn),g2}=Result;
        CatResult{wfn,g(wfn),g2}.fpath=fpath{f2read};
        CatResult{wfn,g(wfn),g2}.pixelsize=PixelSize(f2read);
        CatResult{wfn,g(wfn),g2}.maintrunkROI=maintrunkROI{f2read};
        %CatResult{wfn,g(wfn),g2}.F0_PCA=get_F0PCA(get_subthreshold(CatResult{wfn,g(wfn),g2}.normTraces,CatResult{wfn,g(wfn),g2}.spike(1,:),7,15));
        %CatResult{wfn,g(wfn),g2}.F0_PCA=get_F0PCA(CatResult{wfn,g(wfn),g2}.normTraces,3);
        g(wfn)=g(wfn)+1;
    end
    g2=g2+1;
end

%% CA2 neuron example
nn=4; cax=[-0.005 0.02];
normTr=[]; dSpike=[]; SpMat=[]; dSpMat=[]; bAPMat=[]; bAP_spTime=[]; spike_order=[];
for region_index=1:2 %1:soma 2:dendrite
nTime=size(CatResult{region_index*2,1,nn}.normTraces,2);
normTr{region_index}=CatResult{region_index*2,1,nn}.normTraces./CatResult{region_index*2,1,nn}.F0_PCA;
normTr{region_index}=normTr{region_index}(CatResult{region_index*2,1,nn}.dist_order,:);
somaROI=find(abs(CatResult{region_index*2,1,nn}.interDendDist(1,CatResult{region_index*2,1,nn}.dist_order))<30);
apicalROI=[size(normTr{region_index},1)-2:size(normTr{region_index},1)];
normTr_hi=normTr{region_index}-movmedian(normTr{region_index},250,2);
dSpike=find_spike_bh(mean(normTr_hi(apicalROI,:),1,'omitnan'),0.015,0.002,5);

bwBlue=bwlabel(CatResult{region_index*2,1,nn}.Blue);
[~, normTrMat spind]=get_STA(mean(normTr_hi(somaROI,:),1,'omitnan'),CatResult{region_index*2,1,nn}.spike(1,:),1,1);
[~, shift]=max(normTrMat,[],3);
spTime=spind+shift-2;
spvec_shifted=ind2vec(nTime,spTime,1);

BlueonSpike=[];
for b=1:max(bwBlue)-1
    sptime_tmp=find(spvec_shifted(1,find(bwBlue==b,1)+[0:50])>0,1);
    BlueonSpike=[BlueonSpike sptime_tmp+find(bwBlue==b,1)-1];
end

blue_cat_dilate=ind2vec(nTime,tovec(BlueonSpike'+[0:100]),1);
dSpike=dSpike.*blue_cat_dilate;
spvec_shifted=spvec_shifted.*blue_cat_dilate;
spike_order{region_index}=get_BurstOrder(spvec_shifted,50);

[~, SpMat{region_index}]=get_STA(spvec_shifted,ind2vec(nTime,BlueonSpike,1),20,100);
[~, dSpMat{region_index}]=get_STA(dSpike,ind2vec(nTime,BlueonSpike,1),20,100);
[~, bAPMat{region_index} bAP_spTime{region_index}]=get_STA(normTr{region_index},spvec_shifted,2,4);
end

%bAPMat_cat=cell2mat(cellfun(@(x) reshape(permute(x,[1 3 2]),size(x,1),[]),bAPMat(1),'UniformOutput',false));

[~, Sptime]=cellfun(@(x) find(squeeze(x)),SpMat,'UniformOutput',false);
[~, dSptime]=cellfun(@(x) find(squeeze(x)),dSpMat,'UniformOutput',false);

figure(25); clf;
ax2=nexttile([1 1]); t_edege=[-10:2:60];
histogram(Sptime{1}-20,t_edege,'FaceColor',[0.6 0.6 0.6]); hold all
histogram(dSptime{1}-20,t_edege,'FaceColor',[1 0 0]); hold all
xlabel('Time from first spike (ms)')
ylabel('# of spike')
title('Somatic dendrite stimulation')
legend({'Somatic spike','dSpike'})
ax1=nexttile([1 1]);
histogram(Sptime{2}-20,t_edege,'FaceColor',[0.6 0.6 0.6]); hold all
histogram(dSptime{2}-20,t_edege,'FaceColor',[1 0 0]); hold all
xlabel('Time from first spike (ms)')
ylabel('# of spike')
title('Distal dendrite stimulation')
legend({'Somatic spike','dSpike'})
linkaxes([ax1 ax2])

figure(27); clf; tiledlayout(2,3)
bAPMat_cat=max(bAPMat{1},[],3);
branch_ind={[15 16 18], [12 14 17]};
[~ ,s]=sort(mean(bAPMat_cat(branch_ind{1},:),1,'omitnan'));
spord=spike_order{1}(bAP_spTime{1});
cmap=hsv(max(spord));
cmap=cmap([3 1 4 2],:);

ax1=nexttile([1 3]);
imagesc(bAPMat_cat([1 2 3 6 8 11 13 15 16 18 12 14 17],s))
colormap(turbo);
cb=colorbar;
cb.Label.String=['bAP amplitude' newline '(\DeltaF/F)'];
xlabel('Spike ID, sorted by upper dendrite bAP amplitude')

ax1=nexttile([1 1]);
scatter(mean(bAPMat_cat(branch_ind{1},:),1,'omitnan'),mean(bAPMat_cat(branch_ind{2},:),1,'omitnan'),15,cmap(spord,:),'filled')
xlabel('bAP amplitude of upper dendrite')
ylabel('bAP amplitude of lower dendrite')
axis equal; grid on;
title('Soma stimulation')

ax1=nexttile([1 1]);
h=histogram(mean(bAPMat_cat(branch_ind{1},:),1,'omitnan'),25); hold all
histogram(mean(bAPMat_cat(branch_ind{2},:),1,'omitnan'),h.BinEdges);
legend({'Upper','Lower'})
xlabel('bAP amplitude')
ylabel('Counts')
grid on; box off;

bAPMat_cat2=max(bAPMat{2},[],3);
[~ ,s]=sort(mean(bAPMat_cat2(branch_ind{2},:),1,'omitnan'));
spord=spike_order{2}(bAP_spTime{2});
ax1=nexttile([1 1]);
scatter(mean(bAPMat_cat2(branch_ind{1},:),1,'omitnan'),mean(bAPMat_cat2(branch_ind{2},:),1,'omitnan'),15,cmap(spord,:),'filled')
xlabel('bAP amplitude of upper dendrite')
ylabel('bAP amplitude of lower dendrite')
axis equal; grid on;
title('Dendrite stimulation')

%% CA2 neuron, in movie domain
STAmovie=[]; BlueBoundary=[]; g=1;
for i=[57]
    load(fullfile(fpath{i},'OP_Result.mat'))
    mov=readBinMov_BHL(fpath{i});
    mov_res= mov-mean(mov,3);
    bkg = zeros(2, size(mov,3));
    bkg(1,:) = linspace(-1, 1, size(mov,3));  % linear term
    bkg(2,:) = linspace(-1, 1, size(mov,3)).^2;  % quadratic term
    mov_res = SeeResiduals(mov_res,Result.mc);
    mov_res = SeeResiduals(mov_res,Result.mc.^2);
    mov_res = SeeResiduals(mov_res,Result.mc(:,1).*Result.mc(:,end));
    mov_res= SeeResiduals(mov_res,bkg,1);

    bwBlue=bwlabel(Result.Blue);
    spvec=Result.spike;
    [~, normTrMat spind]=get_STA(Result.normTraces(1,:),Result.spike(1,:),1,1);
    [~, shift]=max(normTrMat,[],3);
    spTime=spind+shift-2;
    spvec_shifted=ind2vec(size(Result.normTraces,2),spTime,1);

    BlueonSpike=[];
    for b=1:max(bwBlue)
        sptime_tmp=find(spvec_shifted(1,find(bwBlue==b,1)+[0:50])>0,1);
        BlueonSpike=[BlueonSpike sptime_tmp+find(bwBlue==b,1)-1];
    end

    statmp = get_STA(tovec(mov_res),ind2vec(size(mov_res,3),BlueonSpike,1),50,50);
    STAmovie(:,:,:,g)=toimg(statmp,size(mov_res,1),size(mov_res,2));
g=g+1;
end

mov_res=mov_res.*double(max(Result.bvMask,[],3)==0);
subMov=tovec(imresize(mov_res,1/2));
[V D eigTrace]=get_eigvector(subMov,30);
mov_res_residual=SeeResiduals(mov_res,V(:,1:2));

figure(28); clf; tiledlayout(2,2);

for r=[7 14]
tr=tovec(mov_res)'*tovec(Result.ftprnt(:,:,r));
tr_res=tovec(mov_res_residual)'*tovec(Result.ftprnt(:,:,r));
corrimg_res=toimg(corr(tovec(mov_res)',tr),size(mov_res,1),size(mov_res,2));
corrimg_residual=toimg(corr(tovec(mov_res_residual)',tr_res),size(mov_res,1),size(mov_res,2));
bd=cell2mat(bwboundaries(bwlabel(Result.ftprnt(:,:,r)>100)));

nexttile([1 1])
imshow2(corrimg_res,[-0.05 0.3]); hold all
%plot(bd(:,2),bd(:,1),'r--')
nexttile([1 1])
imshow2(corrimg_residual,[-0.05 0.3]); hold all
%plot(bd(:,2),bd(:,1),'r--')
end
colormap('turbo')

figure(29); clf; cax_sub=[-0.005 0.021]; cax_diff=[-0.007 0.005];
bAP_vec=ind2vec(nTime,bAP_spTime{1},1);
[~ ,s]=sort(mean(bAPMat_cat(branch_ind{1},:),1,'omitnan'));
[~ ,s2]=sort(mean(bAPMat_cat(17,:),1,'omitnan'));
[~, bAP_STAmov]=get_STA(tovec(mov_res),bAP_vec,10,10);

bAP_STAmov_branch1=toimg(squeeze(mean(bAP_STAmov(:,s(end-5:end),:),2,'omitnan')),size(mov_res,1),size(mov_res,2));
bAP_STAmov_branch1=-bAP_STAmov_branch1./Result.F0PCAimg.*double(max(Result.ftprnt>100,[],3));
bAP_STAmov_branch1=imgaussfilt(bAP_STAmov_branch1,3);

bAP_STAmov_branch2=toimg(squeeze(mean(bAP_STAmov(:,s2(end-3:end),:),2,'omitnan')),size(mov_res,1),size(mov_res,2));
bAP_STAmov_branch2=-bAP_STAmov_branch2./Result.F0PCAimg.*double(max(Result.ftprnt>100,[],3));
bAP_STAmov_branch2=imgaussfilt(bAP_STAmov_branch2,3);

colorSTA=grs2rgb(tovec(bAP_STAmov_branch1),colormap('turbo'),cax_sub(1),cax_sub(2));
colorSTA=reshape(colorSTA,size(bAP_STAmov_branch1,1),size(bAP_STAmov_branch1,2),size(bAP_STAmov_branch1,3),[]);
colorSTA=permute(colorSTA,[1 2 4 3]);

colorSTA2=grs2rgb(tovec(bAP_STAmov_branch2),colormap('turbo'),cax_sub(1),cax_sub(2));
colorSTA2=reshape(colorSTA2,size(bAP_STAmov_branch2,1),size(bAP_STAmov_branch2,2),size(bAP_STAmov_branch2,3),[]);
colorSTA2=permute(colorSTA2,[1 2 4 3]);

bAP_STAmov_branch1_colored=colorSTA.*mat2gray(Result.ref_im)*3;%.*SkelDend(bound:end-bound,bound:end-bound);
writeMov4d(fullfile(save_dir,['CA2_branch1_STA']),bAP_STAmov_branch1_colored,[-11:10],10,1,cax_sub)

bAP_STAmov_branch2_colored=colorSTA2.*mat2gray(Result.ref_im)*3;%.*SkelDend(bound:end-bound,bound:end-bound);
writeMov4d(fullfile(save_dir,['CA2_branch2_STA']),bAP_STAmov_branch2_colored,[-11:10],10,1,cax_sub)

diffMov=bAP_STAmov_branch2-bAP_STAmov_branch1;
colorSTA3=grs2rgb(tovec(diffMov),colormap('turbo'),cax_diff(1),cax_diff(2));
colorSTA3=reshape(colorSTA3,size(diffMov,1),size(diffMov,2),size(diffMov,3),[]);
colorSTA3=permute(colorSTA3,[1 2 4 3]);

bAP_STAmov_branch3_colored=colorSTA3.*mat2gray(Result.ref_im)*3;%.*SkelDend(bound:end-bound,bound:end-bound);
writeMov4d(fullfile(save_dir,['CA2_branch2-1_STA']),bAP_STAmov_branch3_colored,[-11:20],10,1,cax_sub)

figure(30); clf; t_ref=14;
nexttile([1 1])
imshow2(bAP_STAmov_branch1_colored(:,:,:,t_ref),[])
nexttile([1 1])
imshow2(bAP_STAmov_branch2_colored(:,:,:,t_ref),[])
nexttile([1 1])
imshow2(bAP_STAmov_branch3_colored(:,:,:,t_ref),[])