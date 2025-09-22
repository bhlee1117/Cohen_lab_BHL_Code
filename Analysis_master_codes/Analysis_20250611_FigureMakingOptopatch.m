
clear
clc;
cd '/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/Statistics_Optopatch_Prism';
[~, ~, raw] = xlsread(['/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/' ...
    'Prism_OptopatchData_Arrangement.xlsx'], 'Sheet1', 'C5:Q192');

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
TimeSegFrame=cell2mat(cellfun(@(x) (str2num(num2str(x))),raw(:,11),'UniformOutput',false));
place_bin=150; time_segment=15000; overlap=200;
alignedMovFN = {'STA_Mat_SS','STA_Mat_CS','STA_Mat_dSP'};
bound=6;
title_str={'Basal','Apical','Peri-Soma'};
[~, unqInd] = unique([Mouse NeuronInd] ,'row');
set(0,'DefaultFigureWindowStyle','docked')
%
% %% Concatenate data
% StimROI_Ind={'Soma','Distal Dend','WF'};
% StimWfn_Ind={'Ramp Stim','Short Pulse'};
% OpResult=[];
% %foi_somDend=[1 22 14 18 10 26 25 32 46];
% %foi_somDend=[1 22 14 18 10 26 25 32 46 47 48 43 34 37 35];
% foi_somDend=[1 22 14 18 10 26 25 32 46 47 48 43 34 37 35 2:4 6:8 12 13 17 20 29 31 39 41 42 44 49 50 53 58];
% g2=1;
% for i=unqInd(foi_somDend)' %Neuron Index
%
%     SameCellInd=find(Mouse==Mouse(i) & NeuronInd==NeuronInd(i));
%     isSoma = ~cellfun(@isempty,strfind(StimROI(SameCellInd), StimROI_Ind{1}));
%     isDend = ~cellfun(@isempty,strfind(StimROI(SameCellInd), StimROI_Ind{2}));
%     isWF = ~cellfun(@isempty,strfind(StimROI(SameCellInd), StimROI_Ind{3}));
%     isRamp = ~cellfun(@isempty,strfind(StimWfn(SameCellInd), StimWfn_Ind{1}));
%     isSP   = ~cellfun(@isempty,strfind(StimWfn(SameCellInd), StimWfn_Ind{2}));
%
%     ROIwvf_ind=[isSoma isDend isWF isRamp isSP];
%     validind=find(sum(ROIwvf_ind,2)>=2 & isGoodCell(SameCellInd));
%     patterns = [1 0 0 1 0; 1 0 0 0 1; 0 1 0 1 0; 0 1 0 0 1; 0 0 1 1 0; 0 0 1 0 1];
%     values = [1; 2; 3; 4; 5; 6]; % 1=soma, ramp; 2=soma, sp; 3=dend, ramp; 4=dend, sp; 5=WF, ramp; 6=WF, sp;
%
%     g=ones(1,6);
%     for j=1:length(validind)
%         f2read=SameCellInd(validind(j))
%         load(fullfile(fpath{f2read},'OP_Result.mat'),'Result');
%         wfn = values(find(ismember(patterns, ROIwvf_ind(validind(j),:), 'rows')));
%         if ~isempty(wfn)
%             OpResult{wfn,g(wfn),g2}=Result;
%             OpResult{wfn,g(wfn),g2}.fpath=fpath{f2read};
%             OpResult{wfn,g(wfn),g2}.fileInd=f2read;
%             OpResult{wfn,g(wfn),g2}.pixelsize=PixelSize(f2read);
%             OpResult{wfn,g(wfn),g2}.maintrunkROI=maintrunkROI{f2read};
%             OpResult{wfn,g(wfn),g2}.NeuronIndex=find(unqInd(foi_somDend)==i);
%             OpResult{wfn,g(wfn),g2}.MouseNumber=Mouse(f2read);
%             Dsign=ones(1,size(OpResult{wfn,g(wfn),g2}.interDendDist,2));
%             Dsign(OpResult{wfn,g(wfn),g2}.dist_order(1:find(OpResult{wfn,g(wfn),g2}.dist_order==1)-1))=-1;
%             OpResult{wfn,g(wfn),g2}.dendaxis=OpResult{wfn,g(wfn),g2}.interDendDist(1,:).*Dsign * OpResult{wfn,g(wfn),g2}.pixelsize;
%
%             if wfn==1 || wfn==3 || wfn==5
%                 bwBlue=bwlabel(OpResult{wfn,g(wfn),g2}.Blue>0);
%                 bluePeriod=regionprops(bwBlue,'Area');
%                 bluePeriod = cat(1, bluePeriod.Area);
%                 blueRampN=find(bluePeriod>2500);
%                 RampSpike=find(OpResult{wfn,g(wfn),g2}.spike(1,:).*(bwBlue==blueRampN));
%                 if ~isempty(RampSpike)
%                     Rheobase=mean(OpResult{wfn,g(wfn),g2}.Blue(RampSpike(1:3)));
%                     OpResult{wfn,g(wfn),g2}.Rheobase=Rheobase;
%                     OpResult{wfn,g(wfn),g2}.RheobaseBlue=OpResult{wfn,g(wfn),g2}.Blue/Rheobase;
%                 else
%                     PulseSpike=find(OpResult{wfn,g(wfn),g2}.spike(1,:).*(bwBlue>0));
%                     Rheobase=mean(OpResult{wfn,g(wfn),g2}.Blue(PulseSpike(1:3)));
%                     OpResult{wfn,g(wfn),g2}.Rheobase=Rheobase;
%                     OpResult{wfn,g(wfn),g2}.RheobaseBlue=OpResult{wfn,g(wfn),g2}.Blue/Rheobase;
%                 end
%             end
%             %OpResult{wfn,g(wfn),g2}.F0_PCA=get_F0PCA(get_subthreshold(OpResult{wfn,g(wfn),g2}.normTraces,OpResult{wfn,g(wfn),g2}.spike(1,:),7,15));
%             %OpResult{wfn,g(wfn),g2}.F0_PCA=get_F0PCA(OpResult{wfn,g(wfn),g2}.normTraces,3);
%             g(wfn)=g(wfn)+1;
%         end
%     end
%     g2=g2+1;
% end
%
% %% Show Kymos of each stimulation & Save/Load Data
% KymoCell=[]; OpResultLabel=[];
% [Protocol Repeat Neuron] = ind2sub(size(OpResult), find(cellfun(@(x) ~isempty(x),OpResult)));
% g=1;
% for img=1:length(Protocol)
%     RstTemp=OpResult{Protocol(img),Repeat(img),Neuron(img)};
%     bwBlue=bwlabel(RstTemp.Blue>0);
%     bluePeriod=regionprops(bwBlue,'Area');
%     bluePeriod=cat(1,bluePeriod.Area);
%     blueValid=find(bluePeriod<100);
%     bwBlue=bwlabel(ismember(bwBlue,blueValid));
%     for b=1:max(bwBlue)
%         KymoFrame=find(bwBlue==b);
%         KymoCell{g} = RstTemp.normTraces./RstTemp.F0_PCA;
%         KymoCell{g} = KymoCell{g}(RstTemp.dist_order,KymoFrame);
%         OpResultLabel(g,:)=[Protocol(img),Repeat(img),Neuron(img),b];
%         g=g+1;
%     end
% end
% Kymo2use = InteractivelabelImages(KymoCell);
% fieldName={'Stim_Region_Waveform','Repeat','Neuron','BlueStim_Order'};
% OpResultLabel=array2table(OpResultLabel,'VariableNames',fieldName);

%% Load the data
% save(fullfile('/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/Statistics_Optopatch_Prism/Kymo2Use.mat'),'Kymo2use','OpResultLabel','OpResult','-v7.3')
%save(fullfile('/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/Statistics_Optopatch_Prism/OpResult.mat'),'OpResult','-v7.3')
load(fullfile('/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/Statistics_Optopatch_Prism/Kymo2Use.mat'))
load(fullfile('/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/Statistics_Optopatch_Prism/OpResult.mat'))
OpResultLabel=table2array(OpResultLabel);
%% Show representative neuron
figure(20); clf; tiledlayout(2,1,'padding','compact');
nexttile([1 1])
i=72;
load(fullfile(fpath{i},'OP_Result.mat'))
imshow2(Result.Structure',[0 4]); hold all
plot(Result.bluePatt{1, 1}(:,1),Result.bluePatt{1, 1}(:,2),'color',[0 0.6 1])
view([180 90])

nexttile([1 1])
i=74;
load(fullfile(fpath{i},'OP_Result.mat'))
imshow2(Result.Structure',[0 4]); hold all
plot(Result.bluePatt{1, 1}(:,1),Result.bluePatt{1, 1}(:,2),'color',[0 0.6 1])
set(gca,'ydir','reverse')
view([180 90])


%% Show STA trace
i=[73]; nTau=[-15:15];
load(fullfile(fpath{i},'OP_Result.mat'))
taxis=[1:size(Result.normTraces,2)]/1000;

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
normTr=pcafilterTrace(normTr,[1:15]);
normTr=normTr(Result.dist_order(noi_dist),:);

STAtr=get_STA(normTr,Result.spike(1,:),-nTau(1),nTau(end));
STAtr=STAtr-median(STAtr(:,1:5),2);
figure(31); clf;
l=plot(nTau,STAtr');
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(turbo(33),2));

% %% Stimulation subthreshold figure;
% bound=6;
% STAmovie=[]; g=1; BlueonSpike=[]; BlueBoundary=[];
% for i=[73 75]
%     load(fullfile(fpath{i},'OP_Result.mat'))
%     mov=readBinMov_BHL(fpath{i});
%     mov_res= mov-mean(mov,3);
%     bkg = zeros(2, size(mov,3));
%     bkg(1,:) = linspace(-1, 1, size(mov,3));  % linear term
%     bkg(2,:) = linspace(-1, 1, size(mov,3)).^2;  % quadratic term
%     mov_res = SeeResiduals(mov_res,Result.mc);
%     mov_res = SeeResiduals(mov_res,Result.mc.^2);
%     mov_res = SeeResiduals(mov_res,Result.mc(:,1).*Result.mc(:,end));
%     mov_res= SeeResiduals(mov_res,bkg,1);
%
%     bwBlue=bwlabel(Result.Blue);
%     spvec=Result.spike;
%     [~, normTrMat spind]=get_STA(Result.normTraces(1,:),Result.spike(1,:),1,1);
%     [~, shift]=max(normTrMat,[],3);
%     spTime=spind+shift-2;
%     spvec_shifted=ind2vec(size(Result.normTraces,2),spTime,1);
%
%     BlueonSpike=[];
%     for b=1:max(bwBlue)
%         sptime_tmp=find(spvec_shifted(1,find(bwBlue==b,1)+[0:50])>0,1);
%         BlueonSpike=[BlueonSpike sptime_tmp+find(bwBlue==b,1)-1];
%     end
%
%     statmp = get_STA(tovec(mov_res),ind2vec(size(mov_res,3),BlueonSpike,1),50,50);
%     STAmovie(:,:,:,g)=toimg(statmp,size(mov_res,1),size(mov_res,2));
%
%     BlueBoundary{g}=(Result.blueDMDimg);
%     g=g+1;
% end
%
% mov_filt=imresize(pcafilt(imresize(mov_res,1/5),3),5);
% lag1mean=abs(sqrt(mean(mov_filt(:,:,2:end).*mov_filt(:,:,1:end-1),3)));
% mov_res_shrink=imresize(movmedian(mov_res,10,3),1/5);
% [V D]=get_eigvector(tovec(mov_res_shrink),10);
% F0img=sqrt(sum((V.^2).*D(1:10)',2));
% F0img=toimg(F0img,size(mov_res_shrink,1),size(mov_res_shrink,2));
% F0img=imresize(F0img,5);
% %F0img=get_F0img(mov_res);
%
% STAmov_norm=-STAmovie./F0img;
% STAmov_norm=STAmov_norm-mean(STAmov_norm(:,:,[1:15],:),3);
%
% Rfixed = imref2d([size(Result.ref_im,1) size(Result.ref_im,2)]);
% inverseTform = invert(Result.tform);
% %revertedStruct = imwarp(Result.Structure, inverseTform,'OutputView',Rfixed);
% revertedStruct= Result.Structure;
% revertedStruct(revertedStruct==0)=prctile(revertedStruct(:),30);
% revertedStruct=mat2gray(revertedStruct);
%
% STAmovie_normStr=[];
% crop_roi=[94 6  300  159];
% STAmov_norm_crop=STAmov_norm(crop_roi(2):crop_roi(2)+crop_roi(4),crop_roi(1):crop_roi(1)+crop_roi(3),:,:);
% revertedStruct_crop=revertedStruct(crop_roi(2):crop_roi(2)+crop_roi(4),crop_roi(1):crop_roi(1)+crop_roi(3));
% cax_sub=[0.004 0.01];
% cax_sp=[0.004 0.025];
% STAmov_norm_crop_filt=[];
% for spclass_ind=1:2
%     %STAmovie_norm{spclass_ind}=imgaussfilt3(STAnorm_sub./F_refImg,[1.5 1.5 0.1]);%.*SkelDend(bound:end-bound,bound:end-bound);
%     STAmov_norm_crop_filt{spclass_ind}=imgaussfilt(pcafilt(STAmov_norm_crop(:,:,:,spclass_ind),15),4);
%     colorSTA=grs2rgb(tovec(STAmov_norm_crop_filt{spclass_ind}),colormap('turbo'),cax_sub(1),cax_sub(2));
%     colorSTA=reshape(colorSTA,size(STAmov_norm_crop,1),size(STAmov_norm_crop,2),size(STAmov_norm_crop,3),[]);
%     colorSTA=permute(colorSTA,[1 2 4 3]);
%
%     colorSTA2=grs2rgb(tovec(STAmov_norm_crop_filt{spclass_ind}),colormap('turbo'),cax_sp(1),cax_sp(2));
%     colorSTA2=reshape(colorSTA2,size(STAmov_norm_crop,1),size(STAmov_norm_crop,2),size(STAmov_norm_crop,3),[]);
%     colorSTA2=permute(colorSTA2,[1 2 4 3]);
%
%     STAmovie_normStr{spclass_ind}=colorSTA.*revertedStruct_crop*3;%.*SkelDend(bound:end-bound,bound:end-bound);
%     STAmovie_normStr2{spclass_ind}=colorSTA2.*revertedStruct_crop*3;%.*SkelDend(bound:end-bound,bound:end-bound);
% end
%
% % sptype={'SomStim','ddStim'};
% % cax=[-0.005,0.02];
% % writeMov4d(fullfile(save_dir,['STA_dFFStructgrsrgb_ShortPulse']),[imrotate(STAmovie_normStr{1},90) imrotate(STAmovie_normStr{2},90)],[-50:50],10,1,cax)
%
% figure(21); clf; tiledlayout(2,3)
% t_show=[37:41]; t_show_spike=[51:53]; ax1=[];
% for stimROI=1:2
%     subimage=grs2rgb(mean(STAmov_norm_crop_filt{stimROI}(:,:,t_show),3),colormap(turbo),cax_sub(1),cax_sub(2));
%     subimage=subimage.*revertedStruct_crop*3;
%     spimage=grs2rgb(max(STAmov_norm_crop_filt{stimROI}(:,:,t_show_spike),[],3),colormap(turbo),cax_sp(1),cax_sp(2));
%     spimage=spimage.*revertedStruct_crop*3;
%
%     ax1=[ax1 nexttile([1 1])];
%     imshow2(imrotate([Result.Structure(crop_roi(2):crop_roi(2)+crop_roi(4),crop_roi(1):crop_roi(1)+crop_roi(3))],90),[]); hold all
%     Blbd=bwboundaries(imrotate(BlueBoundary{stimROI}(crop_roi(2):crop_roi(2)+crop_roi(4),crop_roi(1):crop_roi(1)+crop_roi(3)),90));
%     plot(Blbd{1}(:,2),Blbd{1}(:,1),'color',[0 0.6 1])
%     nexttile([1 1]);
%     imshow2(imrotate([subimage],90),[]); hold all
%     cb=colorbar;
%     cb.Ticks=[0 0.5 1];
%     cb.TickLabels=[cax_sub(1) mean(cax_sub) cax_sub(2)];
%     cb.Label.String = '\DeltaF/F';
%
%     nexttile([1 1]);
%     imshow2(imrotate([spimage],90),[]); hold all
%     cb=colorbar;
%     cb.Ticks=[0 0.5 1];
%     cb.TickLabels=[cax_sp(1) mean(cax_sp) cax_sp(2)];
%     cb.Label.String = '\DeltaF/F';
%     %plot(BlueBoundary{1}(:,2),BlueBoundary{1}(:,1),'color',[0 0.5 1])
% end
% colormap(turbo);
% colormap(ax1(1),gray);
% colormap(ax1(2),gray);
% %%
% for n=1:size(OpResult,3)
%     [nonemptyr nonemptyc]=find(cellfun(@(x) ~isempty(x), OpResult(:,:,n)));
%     catTr=[];
%     if ~isempty(nonemptyr)
%     for e=1:length(nonemptyr)
%         normsubTrtmp=get_subthreshold(OpResult{nonemptyr(e),nonemptyc(e),n}.normTraces,OpResult{nonemptyr(e),nonemptyc(e),n}.spike(1,:),7,17);
%         catTr=[catTr normsubTrtmp];
%     end
%     [F0PCA nPC]=get_F0PCA(catTr);
%     if nPC==1
%     [F0PCA nPC]=get_F0PCA(catTr,5);
%     end
%     for e=1:length(nonemptyr)
%        OpResult{nonemptyr(e),nonemptyc(e),n}.F0_PCA=F0PCA;
%     end
%     end
% end
%% Short pulses Soma vs Dend
nTau=[60 100];
SP_STA=[]; distFromSoma=[]; SP_STA_std=[];
dendriteaxis_bin=[-280:50:600];
dendriteaxis_bin2=[-300 -50 50 140 460];
perisomadist=[-70 70]; cmap=[0 0 0; 1 0 0];  cmap_light=[0.7 0.7 0.7; 1 0.7 0.7];
prespikesubtime=[-10:-3]; dax=[];
LabelMat_SP=[]; AlignedbAPall_SP=[]; dax_all=[]; g=1;
for n=1:size(OpResult,3) %neuron
    for wvf=[2 4] %stim region
        for rep=1:size(OpResult,2) %repeat
            if ~isempty(OpResult{wvf,rep,n})
                normTr=OpResult{wvf,rep,n}.normTraces./OpResult{wvf,rep,n}.F0_PCA;
                %noi=OpResult{wvf,rep,n}.maintrunkROI; % only main trunks
                noi=[1:size(OpResult{wvf,rep,n}.normTraces,1)]; % all ROIs
                dax=OpResult{wvf,rep,n}.dendaxis(noi);
                MouseInd=Mouse(OpResult{wvf,rep,n}.fileInd); fileInd=OpResult{wvf,rep,n}.fileInd;
                spikeheight=get_STA(mean(normTr(dax>dendriteaxis_bin2(2) & dax <= dendriteaxis_bin2 (3),:),1,'omitnan'),OpResult{wvf,rep,n}.spike(1,:),0,0);
                %normTr=normTr./spikeheight;

                nROI=length(noi);
                nTime=size(OpResult{wvf,rep,n}.normTraces,2);
                noi_dist=find(ismember(OpResult{wvf,rep,n}.dist_order,noi));
                dOrder=OpResult{wvf,rep,n}.dist_order(noi_dist);

                [~, blueOff2, blueOff3] = get_blueoffTrace(zeros(1,nTime), OpResult{wvf,rep,n}.Blue>0, 25, 15);
                blueDialated=(1-blueOff2);  bwBlue=bwlabel(blueDialated>0);
                [bb blueonset]=unique(bwBlue);  blueonset(bb==0)=[];
                som_spike=max(OpResult{wvf,rep,n}.spike([1],:),[],1).*(blueDialated>0);
                som_spike_blueIndexed=max(som_spike,[],1).*(bwBlue);

                SubLarge=get_subthreshold(normTr, max([som_spike; bwBlue>0],[],1),50,500); %subthreshold
                SubTr=get_subthreshold(normTr, max([som_spike],[],1),7,17); %subthreshold
                normTrSub=normTr-SubLarge;

                [STAbAP, AlignedbAP_SP, Sptime]=get_STA(normTrSub,som_spike>0,nTau(1),nTau(2));
                AlignedbAP_SP=permute(AlignedbAP_SP,[1 3 2]);
                AlignedbAP_SP_preSubSubtracted=AlignedbAP_SP-mean(AlignedbAP_SP(:,nTau(1)+prespikesubtime,:),2,'omitnan');
                AlignedbAPall_SP=[AlignedbAPall_SP; squeeze(num2cell(AlignedbAP_SP(dOrder,:,:),[1 2]))];
                dax_all{g,1}=dax; dax_all{g,2}=dax(OpResult{wvf,rep,n}.dist_order);
                sptimefromblueonset=[]; sporderinblue=[]; BluePulseN=[];
                for s=1:length(Sptime)
                    sptimefromblueonset(s)=Sptime(s)-blueonset(som_spike_blueIndexed(Sptime(s)));
                    som_spike_blueIndexed_sub=find(som_spike_blueIndexed==bwBlue(Sptime(s)));
                    sporderinblue(s)=find(som_spike_blueIndexed_sub==Sptime(s));
                    BluePulseN(s)=bwBlue(Sptime(s));
                end

                [~, tmax]=max(STAbAP(:,nTau(1)+[-2:4]),[],2);
                tmax=tmax+nTau(1)-3;
                [AUCbAP, AUCrawbAP_SP]=get_AUC(AlignedbAP_SP,repmat(tmax,1,1,size(AlignedbAP_SP,3)),2,3);
                [~, AUCrawbAP_SP_short, ~, AmprawbAP_SP_short]=get_AUC(AlignedbAP_SP,repmat(tmax,1,1,size(AlignedbAP_SP,3)),1,1);
                [~, AUCrawbAP_SP_subtracted]=get_AUC(AlignedbAP_SP_preSubSubtracted,repmat(tmax,1,1,size(AlignedbAP_SP,3)),2,3);

                AUCbAP_SP_cell = num2cell(AUCrawbAP_SP, 1); % Change this line to use AUC vs AUC raw for Tx rate
                AUCbAP_SP_cell = cellfun(@(x) [dax' x(:,1)],AUCbAP_SP_cell,'UniformOutput',false);
                AUCshortbAP_SP_cell = num2cell(AUCrawbAP_SP_short, 1);
                AUCshortbAP_SP_cell = cellfun(@(x) [dax' x(:,1)],AUCshortbAP_SP_cell,'UniformOutput',false);
                AmpbAP_SP_cell = num2cell(AmprawbAP_SP_short, 1);
                AmpbAP_SP_cell = cellfun(@(x) [dax' x(:,1)],AmpbAP_SP_cell,'UniformOutput',false);
                KinkbAP_SP_cell = num2cell(AUCrawbAP_SP_subtracted, 1);
                KinkbAP_SP_cell = cellfun(@(x) [dax' x(:,1)],KinkbAP_SP_cell,'UniformOutput',false);

                [~, AUCbAPbin_SP_cell, dcenter, dbin_index] = zscore_binning(AUCbAP_SP_cell, dendriteaxis_bin2);
                [~, AUCshortbAPbin_SP_cell] = zscore_binning(AUCshortbAP_SP_cell, dendriteaxis_bin2);
                [~, AmpbAPbin_SP_cell] = zscore_binning(AmpbAP_SP_cell, dendriteaxis_bin2);
                [~, KinkbAPbin_SP_cell] = zscore_binning(KinkbAP_SP_cell, dendriteaxis_bin2);
                emptycell=cellfun(@isempty,AUCbAPbin_SP_cell);

                preSub_SP_Soma=squeeze(mean(AlignedbAP_SP(find(dbin_index{2,1}),nTau(1)+prespikesubtime,:),[1 2],'omitnan'));
                preSub_SP_Apical=squeeze(mean(AlignedbAP_SP(find(dbin_index{4,1}),nTau(1)+prespikesubtime,:),[1 2],'omitnan'));
                trnsmit_SP_bAPAUC=NaN(length(Sptime),1);
                trnsmit_SP_bAPAUC(find(~emptycell))=cellfun(@(x,y) x(4,2)/y(2,2),AUCbAPbin_SP_cell(~emptycell),AUCshortbAPbin_SP_cell(~emptycell)); %Divided by short AUC

                bAPAUC_SP_Apical=cellfun(@(x) x(4,2),AUCbAPbin_SP_cell(~emptycell));
                bAPAUC_SP_Soma=cellfun(@(x) x(2,2),AUCshortbAPbin_SP_cell(~emptycell));
                bAPKink_SP_Soma=cellfun(@(x) x(2,2),KinkbAPbin_SP_cell(~emptycell));
                bAPKink_SP_Apical=cellfun(@(x) x(4,2),KinkbAPbin_SP_cell(~emptycell));
                bAPAmp_SP_Soma=cellfun(@(x) x(2,2),AmpbAPbin_SP_cell(~emptycell));
                bAPAmp_SP_Apical=cellfun(@(x) x(4,2),AmpbAPbin_SP_cell(~emptycell));
                ISI=[NaN diff(Sptime)]';

                catlab=[repmat(n,length(Sptime),1) Sptime' sptimefromblueonset' ...
                    BluePulseN' repmat(MouseInd,length(Sptime),1) repmat(fileInd,length(Sptime),1) repmat(wvf,length(Sptime),1) repmat(rep,length(Sptime),1)...
                    sporderinblue' trnsmit_SP_bAPAUC bAPAUC_SP_Apical bAPAUC_SP_Soma ISI preSub_SP_Apical preSub_SP_Soma bAPKink_SP_Soma bAPKink_SP_Apical, bAPAmp_SP_Soma, bAPAmp_SP_Apical];
                LabelMat_SP=[LabelMat_SP; catlab];  %Session number, Neuron #, Spike Frame, Stimulation Pattern, Rheobase, Spike Order, Mouse#
                g=g+1;
            end
        end
    end
end
fieldName={'NeuronID','SpikeFrame','SpikeFrameBlueonset','StimOrder','MouseID','FileID','StimRegion','RepeatID',...
    'SpikeOrder','TransmissionRatio','AUC_apical','AUC_Soma','ISI','PreSub_apical','PreSub_soma','Kink_soma','Kink_apical','Amp_Soma','Amp_apical'};
LabelMat_SP=array2table(LabelMat_SP,'VariableNames',fieldName);

%% Show spike order and TX ratio (figure 4)

timebin=[0 5 20 40 60 100];
maxTXratio=6; stimstr={'Soma','Dendrite'}; g=1; cmap_sporder=gray(10); cmap_sporder=cmap_sporder(1:2:10,:);
figure(14); clf; tiledlayout(2,1);
for wvf=[2 4]
    nexttile([1 1]);
    % [M S t_center N SubSet]=binning_data({[LabelMat_SP.SpikeFrameBlueonset(show_ind) LabelMat_SP.TransmissionRatio(show_ind)]},timebin);
    % %p=get_pValue(SubSet,0);
    % scatter_heatmap(LabelMat_SP.SpikeFrameBlueonset(show_ind),LabelMat_SP.TransmissionRatio(show_ind),100,200);
    sp2show_sporder=[];
    for sporder=1:4
        % show_ind=find(LabelMat_SP.StimRegion==wvf & abs(LabelMat_SP.TransmissionRatio)<maxTXratio & LabelMat_SP.SpikeOrder==sporder); %Pulse & Soma
        show_ind=find(LabelMat_SP.NeuronID==5 & LabelMat_SP.StimRegion==wvf & abs(LabelMat_SP.TransmissionRatio)<maxTXratio & LabelMat_SP.SpikeOrder==sporder); %Pulse & Soma
        sp2show_sporder{sporder}=LabelMat_SP.TransmissionRatio(show_ind);
    end
    p=Violin_wPoints(sp2show_sporder,cmap_sporder(1:4,:));
    % hold all;
    % errorbar(t_center, M, S ./ sqrt(cellfun(@sum,N)'), 'Color', [1 0 0], 'LineWidth', 1.5);
    %xlim([0 100]); ylim([0 4]);
    %drawPValueLines(p,0,'TextYOffset',0.1,'XCoord',t_center);
    xlabel('Time after blue onset (ms)'); ylabel('Transmission ratio'); box off;
    title([stimstr{g} ' Stim. (Pulse)']); caxis([0 0.0112]);
    colormap([gen_colormap([0 0 0; parula(10)])])
    g=g+1;
end

figure(17); clf; g=1; cmap_stimregion=[0 0 0; 0.5 0.5 0.5];
tiledlayout(1,3,'Padding','tight');
sp2show_sporder_meanN=NaN(max(LabelMat_SP.NeuronID),3,4);
for wvf=[2 4]
    for sporder=1:3
        for n=unique(LabelMat_SP.NeuronID)'
            show_ind=find(LabelMat_SP.NeuronID==n & LabelMat_SP.StimRegion==wvf & abs(LabelMat_SP.TransmissionRatio)<maxTXratio & LabelMat_SP.SpikeOrder==sporder); %Pulse & Soma
            if isempty(show_ind)
                sp2show_sporder_meanN(n,sporder,wvf)=NaN;
            else
                sp2show_sporder_meanN(n,sporder,wvf)=mean(LabelMat_SP.TransmissionRatio(show_ind),'omitnan');
            end
        end
        length(unique(LabelMat_SP.MouseID(find(LabelMat_SP.StimRegion==wvf & abs(LabelMat_SP.TransmissionRatio)<maxTXratio & LabelMat_SP.SpikeOrder==sporder))))
    end
end
for sporder=1:3
    nexttile([1 1]);
    p=Boxplot_wPoints2(permute(sp2show_sporder_meanN(:,sporder,[2 4]),[1 3 2]),cmap_stimregion);
    p123_soma=get_pValue(sp2show_sporder_meanN(:,[1 2 3],[2]),1);
    drawPValueLines(p,0,'TextYOffset',0.1,'XCoord',[1:2]); box off;
    set(gca,'XTick',[]);
    xlabel(counting_string(sporder)); ylabel('Transmission ratio'); box off;
end

%% Kink vs Amp vs AUC
figure(15); clf;
show_ind_somaStim=find(LabelMat_SP.StimRegion==2 & abs(LabelMat_SP.TransmissionRatio)<maxTXratio & LabelMat_SP.SpikeOrder==1 & LabelMat_SP.NeuronID==1);
show_ind_dendStim=find(LabelMat_SP.StimRegion==4 & abs(LabelMat_SP.TransmissionRatio)<maxTXratio & LabelMat_SP.SpikeOrder==1 & LabelMat_SP.NeuronID==1);
nexttile([1 1])
scatter(LabelMat_SP.PreSub_soma(show_ind_somaStim),LabelMat_SP.Kink_soma(show_ind_somaStim),'filled','MarkerFaceAlpha',0.5); hold all
scatter(LabelMat_SP.PreSub_soma(show_ind_dendStim),LabelMat_SP.Kink_soma(show_ind_dendStim),'filled','MarkerFaceAlpha',0.5); hold all
xlabel('Prespike sub. (soma)'); ylabel('Kink amplitude (soma)'); legend({'Soma stim.','Apical stim.'});
nexttile([1 1])
scatter(LabelMat_SP.PreSub_apical(show_ind_somaStim),LabelMat_SP.Kink_apical(show_ind_somaStim),'filled','MarkerFaceAlpha',0.5); hold all
scatter(LabelMat_SP.PreSub_apical(show_ind_dendStim),LabelMat_SP.Kink_apical(show_ind_dendStim),'filled','MarkerFaceAlpha',0.5); hold all
xlabel('Prespike sub. (apical)'); ylabel('Kink amplitude (apical)'); legend({'Soma stim.','Apical stim.'});

nexttile([1 1])
scatter(LabelMat_SP.PreSub_apical(show_ind_somaStim),LabelMat_SP.AUC_apical(show_ind_somaStim),'filled','MarkerFaceAlpha',0.5); hold all
scatter(LabelMat_SP.PreSub_apical(show_ind_dendStim),LabelMat_SP.AUC_apical(show_ind_dendStim),'filled','MarkerFaceAlpha',0.5); hold all
xlabel('Prespike sub. (apical)'); ylabel('AUC (apical)'); legend({'Soma stim.','Apical stim.'});

nexttile([1 1])
scatter(LabelMat_SP.PreSub_apical(show_ind_somaStim),LabelMat_SP.Amp_apical(show_ind_somaStim),'filled','MarkerFaceAlpha',0.5); hold all
scatter(LabelMat_SP.PreSub_apical(show_ind_dendStim),LabelMat_SP.Amp_apical(show_ind_dendStim),'filled','MarkerFaceAlpha',0.5); hold all
xlabel('Prespike sub. (apical)'); ylabel('Amplitude (apical)'); legend({'Soma stim.','Apical stim.'});

nexttile([1 1])
scatter(LabelMat_SP.Kink_soma(show_ind_somaStim),LabelMat_SP.Kink_apical(show_ind_somaStim),'filled','MarkerFaceAlpha',0.5); hold all
scatter(LabelMat_SP.Kink_soma(show_ind_dendStim),LabelMat_SP.Kink_apical(show_ind_dendStim),'filled','MarkerFaceAlpha',0.5); hold all
xlabel('Kink amplitude (soma)'); ylabel('Kink amplitude (apical)'); legend({'Soma stim.','Apical stim.'});

nexttile([1 1])
scatter(LabelMat_SP.Amp_Soma(show_ind_somaStim),LabelMat_SP.AUC_apical(show_ind_somaStim),'filled','MarkerFaceAlpha',0.5); hold all
scatter(LabelMat_SP.Amp_Soma(show_ind_dendStim),LabelMat_SP.AUC_apical(show_ind_dendStim),'filled','MarkerFaceAlpha',0.5); hold all
xlabel('Amplitude (soma)'); ylabel('AUC (apical)'); legend({'Soma stim.','Apical stim.'});

nexttile([1 1])
scatter(LabelMat_SP.PreSub_soma(show_ind_somaStim),LabelMat_SP.PreSub_apical(show_ind_somaStim),'filled','MarkerFaceAlpha',0.5); hold all
scatter(LabelMat_SP.PreSub_soma(show_ind_dendStim),LabelMat_SP.PreSub_apical(show_ind_dendStim),'filled','MarkerFaceAlpha',0.5); hold all
xlabel('Prespike sub. (soma)'); ylabel('Prespike sub. (apical)'); legend({'Soma stim.','Apical stim.'});

nexttile([1 1])
scatter(LabelMat_SP.AUC_apical(show_ind_somaStim),LabelMat_SP.Amp_apical(show_ind_somaStim),'filled','MarkerFaceAlpha',0.5); hold all
scatter(LabelMat_SP.AUC_apical(show_ind_dendStim),LabelMat_SP.Amp_apical(show_ind_dendStim),'filled','MarkerFaceAlpha',0.5); hold all
xlabel('AUC (apical)'); ylabel('Amplitude (apical)'); legend({'Soma stim.','Apical stim.'});

nexttile([1 1])
scatter(LabelMat_SP.Kink_apical(show_ind_somaStim),LabelMat_SP.AUC_apical(show_ind_somaStim),'filled','MarkerFaceAlpha',0.5); hold all
scatter(LabelMat_SP.Kink_apical(show_ind_dendStim),LabelMat_SP.AUC_apical(show_ind_dendStim),'filled','MarkerFaceAlpha',0.5); hold all
xlabel('Kink amplitude (apical)'); ylabel('AUC (apical)'); legend({'Soma stim.','Apical stim.'});
%% Show Representative STA (figure 4)
colormap(turbo); cmap_ExTr=gen_colormap(Plasma,10);
figure(14); clf;
figure(15); clf; tiledlayout(4,1,'padding','compact');
figure(16); clf; tiledlayout(2,5,'padding','compact');
show_neuron=5; ax2=[]; cax=[-0.5 3.5]; t2show=[-22 -10 1 2 3];
SWCpoints=OpResult{1,1,show_neuron}.SWC; SWCpoints(1,3)=70;

dax=OpResult{1,1,show_neuron}.dendaxis;
ROIvec=OpResult{1,1,5}.maintrunkROI; ftprint2show=OpResult{1,1,5}.ftprnt(:,:,ROIvec); ftprintall=OpResult{1,1,5}.ftprnt;
nROI=size(OpResult{1,1,5}.ftprnt,3); ftprint2show_trace=get_coord(ftprint2show);
somROI=1; dendROI=7;

figure(14);
nexttile([1 1])
boundary_s=bwboundaries(ftprint2show(:,:,somROI)); boundary_ap=bwboundaries(ftprint2show(:,:,dendROI));
boundary_blue=OpResult{2,1,show_neuron}.bluePatt;
boundary_blue2=OpResult{4,1,show_neuron}.bluePatt;
showScaleScatter([1:nROI], SWCpoints, ftprintall,repmat([0 0 0],256,1)); hold all
plot(boundary_s{1}(:,1),boundary_s{1}(:,2),'color',cmap_ExTr(1,:),'linewidth',2);
plot(boundary_ap{1}(:,1),boundary_ap{1}(:,2),'color',cmap_ExTr(end,:),'linewidth',2);
plot(boundary_blue{1}(:,1),boundary_blue{1}(:,2),'color',[0 0.5 1],'linewidth',2);
plot(boundary_blue2{1}(:,1),boundary_blue2{1}(:,2),'color',[0 0.5 1],'linewidth',2);
plot(ftprint2show_trace(:,2),ftprint2show_trace(:,1),'r','linewidth',2);
drawScaleBar(100/OpResult{1,1,show_neuron}.pixelsize,'vertical')
view([180 -90])

for wvf=[2 4]

    figure(15);
    ax2=[ax2 nexttile([1 1])];
    sp2mean=LabelMat_SP.SpikeOrder==1 & LabelMat_SP.StimRegion==4 & LabelMat_SP.NeuronID==5;
    STA2showall=mean(cat(3,AlignedbAPall_SP{LabelMat_SP.SpikeOrder==1 & LabelMat_SP.StimRegion==wvf & LabelMat_SP.NeuronID==5}),3,'omitnan');
    STA2show=STA2showall(ismember(OpResult{1,1,show_neuron}.dist_order,ROIvec)',:);
    STA2show_interp=interp1(dax(ROIvec),STA2show,linspace(dax(ROIvec(1)),dax(ROIvec(end)),10));
    imagesc([-nTau(1):nTau(2)],ROIvec,STA2show_interp,cax);
    set_kymoYtick(dax(ROIvec));

    % ax2=[ax2 nexttile([1 1])];
    % plot([-nTau(1):nTau(2)],STA2show(1,:),'color',cmap_ExTr(1,:)); hold all;
    % plot([-nTau(1):nTau(2)],STA2show(end,:),'color',cmap_ExTr(end,:));
    % ylim(cax); xlim([-50 50]);
    % axis off;
    % drawScaleBar(1,'vertical');
    % drawScaleBar(50,'horizontal');
    xlabel('Peri-spike time (ms)')

    figure(16);
    for t=1:length(t2show)
        nexttile([1 1]);
        V2show=STA2showall(:,t2show(t)+nTau(1));
        %V2showcolor=vec2cmap(V2show,turbo,cax);
        showScaleScatter(V2show, SWCpoints, ftprintall(:,:,OpResult{1,1,show_neuron}.dist_order),'turbo',cax);
        axis tight;
        view([180 -90])
    end
end
figure(15);
linkaxes(ax2,'x'); xlim([-50 50]);

%% Show Short pulse stimulation voltage trace (Raw trace)
figure(21); clf; tiledlayout(12,6,'padding','tight'); ax1=[]; ax3=[];
dendriteaxis_bin=[-100:25:300]; cax=[-0.3 2]*0.01; t_show=[3 13]; cmap=[118 85 157; 194 102 56]/256;
SomaROI=[1]; ApicalROI=[6 7]; g=1;
ShowT_zoom=[315]+[-30:60];
for i=[73 75];
    load(fullfile(fpath{i},'OP_Result.mat'))
    taxis=[1:size(Result.normTraces,2)]/1000;

    noi=maintrunkROI{i};
    noi_dist=ismember(Result.dist_order,noi);

    Dsign=ones(1,size(Result.interDendDist,2));
    Dsign(Result.dist_order(1:find(Result.dist_order==1)-1))=-1;
    dendaxis=Result.interDendDist(1,Result.dist_order(noi_dist)).*Dsign(Result.dist_order(noi_dist));

    normResidue=permute(Result.normTraces,[3 1 2]);
    normResidue=SeeResiduals_tiling(normResidue,Result.mc,1,5000);
    normResidue=SeeResiduals_tiling(normResidue,Result.mc.^2,1,5000);
    normResidue=SeeResiduals_tiling(normResidue,Result.mc(:,1).*Result.mc(:,end),1,5000);
    normResidue=squeeze(normResidue);

    %F0PCA=get_F0PCA(normResidue);

    %normTr=normResidue./F0PCA;
    normTr=normResidue./Result.F0_PCA;
    %normTr=normTr-movprc(normTr,200,30,2);
    %subTr=get_blueoffTrace(normTr(1,:),Result.Blue,50,70);
    %normTr=normTr-subTr;
    %normTr=normTr-movmedian(normTr,100,2);
    normTr=pcafilterTrace(normTr,[1:15]);
    normTr=normTr(Result.dist_order(noi_dist),:);
    %normTr=get_bandstop(normTr,1000,[49 80]);

    ax3=[ax3 nexttile((g-1)*30+1,[5 4])];
    plot(taxis,mean(normTr(SomaROI,:),1,'omitnan'),'color',cmap(1,:),'linewidth',1.5); hold all
    plot(taxis,mean(normTr(ApicalROI,:),1,'omitnan'),'color',cmap(2,:),'linewidth',1.5);
    %ylim([-0.015 0.03])
    axis off

    ax1=[ax1 nexttile((g-1)*30+5,[2 2])];
    [~, sortind]=sort(dendaxis,'ascend');
    [Xq, Yq] = meshgrid([taxis], dendriteaxis_bin);
    normTr_interp = interp2([taxis], dendaxis(sortind), normTr(sortind,:), Xq, Yq, 'linear');
    normTr_interpShow=normTr_interp(:,ShowT_zoom);
    normTr_interpShow=normTr_interpShow-mean(normTr_interpShow(:,1:10),2);
    pcolor(taxis(ShowT_zoom),dendriteaxis_bin,normTr_interpShow);
    caxis(cax); colormap(turbo(256)); shading flat;
    cb=colorbar;
    cb.Label.String='\DeltaF/F';
    set(gca,'XTick',t_show,'XTickLabel',t_show-t_show(1),'YDir','reverse')
    ylim([0 250])

    ax1=[ax1 nexttile((g-1)*30+17,[3 2])];
    normTr_show=[mean(normTr(SomaROI,ShowT_zoom),1,'omitnan'); mean(normTr(ApicalROI,ShowT_zoom),1,'omitnan')];
    normTr_show=normTr_show-mean(normTr_show(:,1:10),2);
    l=plot(taxis(ShowT_zoom),normTr_show','linewidth',1.5);
    arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(cmap,2)); axis tight off;
    ylim(cax)
    g=g+1;
end

ax2=nexttile([2 4]);
plot(taxis, Result.Blue);
linkaxes([ax2 ax3],'x');  xlim(t_show); axis off;
ax1=[ax1 nexttile(65,[2 2])];
plot(taxis(ShowT_zoom), Result.Blue(ShowT_zoom)); axis tight off;
linkaxes([ax1],'x')
%% Pulse & Ramp
nTau=[-20 20]; nTauPeriSp=[-3:7]; nTauPeriSp2=[-3:3];
g=ones(1,4); SP_STA=[]; distFromSoma=[]; LabelMat=[]; AlignedbAPall=[];
dendriteaxis_bin=[-200 -100 -50 50 100:50:500];
dendriteaxis_bin2=[-300 -50 50 140 460];
timebin=[0 30 35 70 100 200:100:500];
Rheobin=[1 1.3 1.6 2 3 4 5];
perisomadist=[-30 30]; cmap=[0 0 0; 1 0 0];  cmap_light=[0.7 0.7 0.7; 1 0.7 0.7];
Kymos=[]; AUCbAPall=[];
%OpResultLabel=table2array(OpResultLabel);

SpikeAmp=[]; SpikeAUC=[]; dax=[]; g=1;
for wvf=[1:4]; % wvf=1 : soma stimulation, wvf=3 : Dendrite stimulation, wvf =5 : WF stimulation
    BlueStimN=OpResultLabel(find(Kymo2use'>0 & OpResultLabel.Stim_Region_Waveform==wvf),[3 4 2]); %Neuron#, bwBlue, Repeat
    for n=1:size(BlueStimN,1) % session
        cax=[];
        Neuron=BlueStimN.Neuron(n);
        rep=BlueStimN.Repeat(n);
        BluePulseN=BlueStimN.BlueStim_Order(n);
        noi=[1:size(OpResult{wvf,rep,Neuron}.normTraces,1)];
        dax{n,wvf}=OpResult{wvf,rep,Neuron}.dendaxis(noi);
        perisomaInd=find(dax{n,wvf}'>perisomadist(1) & dax{n,wvf}'<perisomadist(2));

        if ~isempty(OpResult{wvf,rep,Neuron})
            MouseInd=OpResult{wvf,rep,Neuron}.MouseNumber;
            fileInd=OpResult{wvf,rep,Neuron}.fileInd;
            normTr=OpResult{wvf,rep,Neuron}.normTraces./OpResult{wvf,rep,Neuron}.F0_PCA; %Read traces
            %noi=OpResult{wvf,rep,Neuron}.maintrunkROI;
            normTr=normTr(noi,:);
            nROI=length(noi); nTime=size(OpResult{wvf,rep,Neuron}.normTraces,2); %set ROI and time
            dOrder=OpResult{wvf,rep,Neuron}.dist_order(noi);
            Classvec = get_Class2index(OpResult{wvf,rep,Neuron}.SpClass([1 2],:));
            som_spike = Classvec .* OpResult{wvf,rep,Neuron}.spike(1,:);
            bwBlue=bwlabel(OpResult{wvf,rep,Neuron}.Blue);
            bwBlue((bwBlue-(max(bwBlue)-5))<0)=0; bwBlue=bwlabel(bwBlue>0);

            [STAbAP]=get_STA(normTr,max(som_spike>0,[],1),-nTau(1),nTau(2));
            [~, tmax]=max(STAbAP(:,-nTau(1)+[-2:4]),[],2);
            tmax=tmax-nTau(1)-3;

            orderVector = get_BurstOrder(max(som_spike,[],1), 20);
            SpikeinBlue=max(som_spike,[],1).*(bwBlue==BluePulseN);
            SpikeinBlue_ind=find(SpikeinBlue);
            BlueOnset=find(bwBlue==BluePulseN,1);
            SubLarge=get_subthreshold(normTr, max([som_spike; bwBlue>0],[],1),200,1000); %subthreshold
            %SubLarge=get_subthreshold(normTr, max([som_spike],[],1), 7,50);
            normTrSub=normTr-SubLarge;
            Kymos{n,wvf}=normTr(OpResult{wvf,rep,Neuron}.dist_order,unique(find(bwBlue==BluePulseN)'+[-50:50]));

            [~, AlignedbAP, Sptime]=get_STA(normTrSub,SpikeinBlue>0,-nTau(1),nTau(2));
            RheoBaseAmp=OpResult{wvf,rep,Neuron}.RheobaseBlue(Sptime);
            IsCS=OpResult{wvf,rep,Neuron}.CStrace(Sptime);
            BurstOrder=orderVector(Sptime);
            Sptime_inBlue=Sptime-BlueOnset+1;
            AlignedbAP=permute(AlignedbAP,[1 3 2]);
            mat2cat=squeeze(num2cell(AlignedbAP(OpResult{wvf,rep,Neuron}.dist_order,:,:), [1 2]))';
            AlignedbAPall=[AlignedbAPall mat2cat];
            [AUCbAP, AUCrawbAP kinkAmp kink_raw]=get_AUC(AlignedbAP,repmat(tmax,1,1,size(AlignedbAP,3)),2,3);
            [AUCbAP_short, AUCrawbAP_short]=get_AUC(AlignedbAP,repmat(tmax,1,1,size(AlignedbAP,3)),1,1);

            AUCbAP_cell = num2cell(AUCrawbAP, 1); % Change this line to use AUC vs AUC raw for Tx rate
            AUCbAP_cell = cellfun(@(x) [dax{n,wvf}' x(:,1)],AUCbAP_cell,'UniformOutput',false);
            AUCbAPall=[AUCbAPall AUCbAP_cell];
            AmpbAP_cell = num2cell(AUCrawbAP_short, 1);
            AmpbAP_cell = cellfun(@(x) [dax{n,wvf}' x(:,1)],AmpbAP_cell,'UniformOutput',false);

            [~, AUCbAPbin_cell, dcenter] = zscore_binning(AUCbAP_cell, dendriteaxis_bin2);
            [~, AmpbAPbin_cell, dcenter] = zscore_binning(AmpbAP_cell, dendriteaxis_bin2);
            emptycell=cellfun(@isempty,AUCbAPbin_cell);

            centroid_bAPAUC=sum(dax{n,wvf}'.*AUCrawbAP,1,'omitnan')./sum(AUCrawbAP,1,'omitnan');
            trnsmit_bAPAUC=NaN(length(Sptime),1);
            trnsmit_bAPAUC(find(~emptycell))=cellfun(@(x,y) x(4,2)/y(2,2),AUCbAPbin_cell(~emptycell),AmpbAPbin_cell(~emptycell)); %Divided by short AUC
            trnsmit_bAPAUC2=NaN(length(Sptime),1);
            trnsmit_bAPAUC2(find(~emptycell))=cellfun(@(x) x(4,2)/x(2,2),AUCbAPbin_cell(~emptycell)); %Divided by same length AUC
            bAPAUC_Apical=cellfun(@(x) x(4,2),AUCbAPbin_cell(~emptycell));
            bAPAUC_Soma=cellfun(@(x) x(2,2),AmpbAPbin_cell(~emptycell));
            ISI=[NaN diff(Sptime)]';

            [SpType, ~]=find(som_spike(:,Sptime));

            %AUCbAP=zscore(AUCbAP,[],2);
            %AUCrawbAP=AUCrawbAP./mean(AUCrawbAP(perisomaInd,:),1,'omitnan');
            ztemp = num2cell(AUCrawbAP, 1);
            ztemp = cellfun(@(x) [dax{n,wvf}' x(:,1)],ztemp,'UniformOutput',false);
            SpikeAUC=[SpikeAUC ztemp];
            dax{n,wvf}=dax{n,wvf}(OpResult{wvf,rep,Neuron}.dist_order);

            catlab=[repmat(n,length(Sptime),1) repmat(Neuron,length(Sptime),1) Sptime_inBlue' ...
                repmat(BluePulseN,length(Sptime),1) RheoBaseAmp' [1:length(Sptime)]'  repmat(MouseInd,length(Sptime),1) repmat(fileInd,length(Sptime),1) repmat(wvf,length(Sptime),1) ...
                IsCS' BurstOrder' repmat(rep,length(Sptime),1) SpType trnsmit_bAPAUC bAPAUC_Apical bAPAUC_Soma ISI];
            LabelMat=[LabelMat; catlab];  %Session number, Neuron #, Spike Frame, Stimulation Pattern, Rheobase, Spike Order, Mouse#
        end
    end
end
fieldName={'SessionID','NeuronID','SpikeFrame','StimOrder','Rheobase','SpikeOrder','MouseID','FileID','StimRegion',...
    'InCS','BurstOrder','RepeatID','SpikeType','TransmissionRatio','AUC_apical','Amp_Soma','ISI'};
LabelMat=array2table(LabelMat,'VariableNames',fieldName);
%% Show TX ratio over time (figure 5)
timebin=[0 30 35 85 110 200:100:500];
Rheobin=[1 1.3 1.6 2 3 4 5];
maxTXratio=7; stimstr={'Soma','Dendrite'}; g=1;
figure(23); clf; tiledlayout(1,2);
for wvf=[1]
    nexttile([1 1]);
    show_ind=find(LabelMat.StimOrder~=6 & LabelMat.StimRegion==wvf & abs(LabelMat.TransmissionRatio)<maxTXratio); %Pulse & Soma
    [M S t_center N SubSet]=binning_data({[LabelMat.SpikeFrame(show_ind) LabelMat.TransmissionRatio(show_ind)]},timebin);
    p=get_pValue(SubSet,0);
    scatter_density(LabelMat.SpikeFrame(show_ind),LabelMat.TransmissionRatio(show_ind),40);
    %scatter_heatmap(LabelMat.SpikeFrame(show_ind),LabelMat.TransmissionRatio(show_ind),100,200);
    hold all;
    errorbar(t_center, M, S ./ sqrt(cellfun(@sum,N)'), 'Color', [1 0 0], 'LineWidth', 1.5);
    xlim([0 500]); ylim([0 4]);
    %set(gca,'xscale','log')
    %drawPValueLines(p,0,'TextYOffset',0.1,'XCoord',t_center);
    xlabel('Time after blue onset (ms)'); ylabel('Transmission ratio'); box off;
    title([stimstr{g} ' Stim. (Pulse)']); caxis([0 0.0012]);
    colormap([gen_colormap([0 0 0; parula(10)])])

    nexttile([1 1]);
    show_ind=find(LabelMat.StimOrder==6 & LabelMat.StimRegion==wvf & abs(LabelMat.TransmissionRatio)<maxTXratio); %Pulse & Soma
    [M S Rheo_center N SubSet]=binning_data({[LabelMat.Rheobase(show_ind) LabelMat.TransmissionRatio(show_ind)]},Rheobin);
    pR=get_pValue(SubSet,0);
    scatter_density(LabelMat.Rheobase(show_ind),LabelMat.TransmissionRatio(show_ind),40);
    %scatter_heatmap(LabelMat.Rheobase(show_ind),LabelMat.TransmissionRatio(show_ind),500,200);
    hold all;
    errorbar(Rheo_center, M, S ./ sqrt(cellfun(@sum,N)'), 'Color', [1 0 0], 'LineWidth', 1.5);
    xlim([0.5 5]); ylim([0 4]); box off;
    %drawPValueLines(pR,0,'TextYOffset',0.1,'XCoord',Rheo_center);
    xlabel('Rheobase'); ylabel('Transmission ratio')
    title([stimstr{g} ' Stim. (Ramp)']); caxis([0 0.195]);
    colormap([gen_colormap([parula(10)])])
    g=g+1;
end

% Show TX ratio, spike order
TXratio_sporder=[]; g2=1;
for stimregion2show=[1 3];
    for sporder=1:13
        g=1;
        [unqLab, refind]=unique([LabelMat.SessionID LabelMat.StimRegion],'row');
        for n=refind(unqLab(:,2)==stimregion2show,1)' %soma targeted
            sessioninterested=unique(LabelMat.SessionID(n));
            stimregioninterested=unique(LabelMat.StimRegion(n));
            ind2norm=find(LabelMat.SessionID==sessioninterested & LabelMat.StimRegion==stimregioninterested);
            show_ind=find(LabelMat.StimOrder~=6 & LabelMat.StimRegion==stimregioninterested & abs(LabelMat.TransmissionRatio)<maxTXratio & LabelMat.SpikeOrder==sporder & LabelMat.SessionID==sessioninterested); %Pulse & Soma
            % TXratio_sporder(sporder,g,1)=mean(LabelMat.TransmissionRatio(show_ind)./median(LabelMat.TransmissionRatio(ind2norm),'omitnan'),'omitnan');
            % TXratio_sporder(sporder,g,2)=std(LabelMat.TransmissionRatio(show_ind)./median(LabelMat.TransmissionRatio(ind2norm),'omitnan'),'omitnan');
            TXratio_sporder(sporder,g,1)=mean(LabelMat.TransmissionRatio(show_ind),'omitnan');
            TXratio_sporder(sporder,g,2)=std(LabelMat.TransmissionRatio(show_ind),'omitnan');
            g=g+1;
        end
    end
    g2=g2+1;
end
figure(24); clf; cmap_order=gray(10); show_sporder=[1:7];
nexttile([1 1]);
TXratio_sporder(TXratio_sporder>5 | TXratio_sporder<0)=NaN;
TXratio_sporder_mean=mean(TXratio_sporder(show_sporder,:,1)',1,'omitnan');
TXratio_sporder_sem=std(TXratio_sporder(show_sporder,:,1)',0,1,'omitnan')./sqrt(sum(~isnan(TXratio_sporder(show_sporder,:,1)),2)');
%plot([1:4],TXratio_sporder(show_sporder,:,1)','color',[0.7 0.7 0.7]); hold all;
%errorbar([1:4],TXratio_sporder_mean,TXratio_sporder_sem,'r');
Violin_wPoints(TXratio_sporder(show_sporder,:,1)',cmap_order(show_sporder,:));
p_values = get_pValue(TXratio_sporder(show_sporder,:,1)', 1);
drawPValueLines(p_values,0.1);
xlim([show_sporder(1)-0.5 show_sporder(end)+0.5]); ylim([0 6]);
set(gca,'xtick',[1:5],'XTickLabel',counting_string(1:5))
ylabel('TX ratio');  xlabel('Evoked spike order');

%% Representative Pulse stim (Figure 5a)
% show example bAP amplitude/AUC image
ShowSession=[60 1]; ExampleRep=1;
ExampleSession=ShowSession(1,1); wvf=ShowSession(1,2);
ExampleNeuron=unique(LabelMat.NeuronID(find(LabelMat.SessionID==ExampleSession & LabelMat.StimRegion==wvf)));
ROIvec=setdiff(OpResult{1,1,ExampleNeuron}.maintrunkROI, []); 
ftprint2show=OpResult{1,1,ExampleNeuron}.ftprnt(:,:,ROIvec); ftprintall=OpResult{1,1,ExampleNeuron}.ftprnt;
nROI=size(OpResult{1,1,ExampleNeuron}.ftprnt,3); ftprint2show_trace=get_coord(ftprint2show);
somROI=1; dendROI=[9]; caxshow=[-2.5 6.5];
SWCpoints=OpResult{1,1,ExampleNeuron}.SWC; SWCpoints(1,3)=70; cmap_ExTr=gen_colormap(Plasma,10);

figure(25); clf;
nexttile([1 1])
boundary_s=bwboundaries(max(ftprint2show(:,:,somROI),[],3)); boundary_ap=bwboundaries(max(ftprint2show(:,:,dendROI),[],3));
boundary_blue=bwboundaries(OpResult{1,1,ExampleNeuron}.BlueDMDimg);
showScaleScatter([1:nROI], SWCpoints, ftprintall,repmat([0 0 0],256,1)); hold all
plot(boundary_s{1}(:,1),boundary_s{1}(:,2),'color',cmap_ExTr(1,:),'linewidth',2);
plot(boundary_ap{1}(:,1),boundary_ap{1}(:,2),'color',cmap_ExTr(end,:),'linewidth',2);
plot(boundary_blue{1}(:,1),boundary_blue{1}(:,2),'color',[0 0.5 1],'linewidth',2);
plot(ftprint2show_trace(:,2),ftprint2show_trace(:,1),'r','linewidth',2);
drawScaleBar(100/OpResult{1,1,ExampleNeuron}.pixelsize,'vertical')
view([210 90])

nROI=size(OpResult{1,1,ExampleNeuron}.normTraces,1);
ftprint=OpResult{1,1,ExampleNeuron}.ftprnt(:,:,OpResult{1,1,ExampleNeuron}.dist_order)>0;
dax=OpResult{1,1,ExampleNeuron}.dendaxis;


normTr=OpResult{wvf,1,ExampleNeuron}.normTraces./OpResult{wvf,1,ExampleNeuron}.F0_PCA;
normTr=normTr(ROIvec,:);
normTr=pcafilterTrace(normTr,[1 2 4 6 7 8 9]);
sublarge=get_subthreshold(normTr,OpResult{wvf,ExampleRep,ExampleNeuron}.Blue>0,300,1000);
Tr2show=[mean(normTr(somROI,:),1); mean(normTr(dendROI,:),1)];
normTr_interp=interp1(dax(ROIvec),normTr,linspace(dax(ROIvec(1)),dax(ROIvec(end)),10));
dax2show=linspace(dax(ROIvec(1)),dax(ROIvec(end)),10);

kymosSub2show=get_subthreshold(normTr,OpResult{wvf,ExampleRep,ExampleNeuron}.spike(1,:),7,19);

figure(27); clf; tiledlayout(5,1); % show example bAP amplitude/AUC image
ax1=nexttile([2 1]);
time2show=[2980:3550];
kymos2show=normTr_interp(:,time2show);
kymos2show=kymos2show-median(kymos2show(:,1:35),2);
t_ax=[1:size(kymos2show,2)]-20;
imagesc(t_ax, [1:size(kymos2show,1)], kymos2show,caxshow); colormap(ax1,'turbo'); cb=colorbar; cb.Label.String = 'Z score';
set_kymoYtick(dax2show)

ax2=nexttile([3 1]);
l=plot(t_ax,Tr2show([1 2],time2show)','linewidth',1.5);
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(cmap_ExTr([1 10],:),2));
ylabel('Z score'); xlabel('Time (ms)');
legend({'Soma','Apical dendrite'}); axis tight off;
linkaxes([ax1 ax2],'x'); box off; ylim(caxshow);

ss=1; t_start=[3000 3290]; t_showlength=90;
t_ax=[1:t_showlength+1];
figure(30); clf;
tiledlayout(5,length(t_start),'padding','compact')
for t=1:length(t_start)
    ax1=nexttile(t,[2 1]);
    kymos2show=normTr_interp(:,t_start(t):t_start(t)+t_showlength);
    %kymos2show=pcafilterTrace(kymos2show,[1:5]);
    %kymos2show=kymos2show(ShowROI,:);
    imagesc(t_ax, [1:size(kymos2show,1)], kymos2show,caxshow); colormap(ax1,'turbo'); cb=colorbar; cb.Label.String = 'Z score';
    Trcat=Tr2show(:,t_start(t):t_start(t)+t_showlength);
    set_kymoYtick(dax2show)
    ax2=nexttile(t+length(t_start)*2,[3 1]);
    l=plot(t_ax,Trcat([1 2],:)','linewidth',1.5);
    arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(cmap_ExTr([1 10],:),2));
    ylabel('\DeltaF/F'); xlabel('Time (ms)');
    axis tight off;
    ylim(caxshow);
end

Splist=find(LabelMat.SessionID==ExampleSession & LabelMat.StimRegion==1 & LabelMat.NeuronID==ExampleNeuron);
showSP={[1],[2 3 4],[13 14 15]}; excludeROI=[15 35]; cax_AUC=[3 15];
AUC2show=cell2mat(cellfun(@(x) x(:,2),AUCbAPall(Splist),'UniformOutput',false));
ftprintall2show=ftprintall;
ftprintall2show(:,:,excludeROI)=[];
AUC2show=AUC2show;%./AUC2show(1,:);
AUC2show(excludeROI,:)=[];
%tiledlayout(1,9,'Padding','tight'); 
Fcat=[];
for s=1:length(showSP)    
    figure(31+s); clf;
    %nexttile([1 1])
    showvec=mean(AUC2show(:,showSP{s}),2,'omitnan');
    showScaleScatter(showvec, SWCpoints, ftprintall2show,'turbo',cax_AUC);
    axis tight;
    view([210 90])
end
% 
% figure(32); clf;
% set(gcf,'Color','w')
% imshow2(Fcat,[]);
% colormap(turbo); 
% cb=colorbar; 
% cb.Ticks=[0 1]; caxis([0 1]);
% cb.TickLabels=num2cell(cax_AUC)';
% cb.Label.String = 'AUC (A.U.)';


%%
SomRP_path=[15 54 72 86 92 4 76]; %22
DenRP_path=[13 55 74 88 91 3 79]; %23
validresultind=find(cell2mat(cellfun(@(x) ~isempty(x),OpResult,'UniformOutput',false)));
finterest=cell2mat(cellfun(@(x) x.fileInd,OpResult(validresultind),'UniformOutput',false));
finterest_SomRP=validresultind(ismember(finterest,SomRP_path));
finterest_DenRP=validresultind(ismember(finterest,DenRP_path));
[finterest_SomRPw, finterest_SomRPr, finterest_SomRPn]=ind2sub(size(OpResult),finterest_SomRP);
[finterest_DenRPw, finterest_DenRPr, finterest_DenRPn]=ind2sub(size(OpResult),finterest_DenRP);

%path_cat=[SomRP_path; DenRP_path];
rheo_bin=[0:1.3:10];
CSfraction=NaN(size(finterest_SomRP,1),2);
for n=1:size(finterest_SomRP,1) %neuron
    resultdat=OpResult{finterest_SomRPw(n),finterest_SomRPr(n),finterest_SomRPn(n)};
    spikeInfoInd=(LabelMat.NeuronID==finterest_SomRPn(n) & LabelMat.RepeatID==finterest_SomRPr(n) & LabelMat.StimRegion==finterest_SomRPw(n));
    CSfraction(n,1)=sum(LabelMat.InCS(spikeInfoInd))/sum(spikeInfoInd);

    spikeInfoInd=(LabelMat.NeuronID==finterest_DenRPn(n) & LabelMat.RepeatID==finterest_DenRPr(n) & LabelMat.StimRegion==finterest_DenRPw(n));
    CSfraction(n,2)=sum(LabelMat.InCS(spikeInfoInd))/sum(spikeInfoInd);
end

showneuron=4;
figure(33); clf; scaleoffset=1.5; tiledlayout(1,5);
nexttile([1,4]);
result2show=OpResult{finterest_SomRPw(showneuron),finterest_SomRPr(showneuron),finterest_SomRPn(showneuron)};
ROIvec=result2show.maintrunkROI;
normTR=result2show.normTraces./get_F0PCA(result2show.normTraces);
[S]=get_STA(normTR,result2show.spike(1,:),0,0);
normTR=normTR./S(1);
somROI=1; dendROI=10;
normTR2show_S=mean(normTR(somROI,:),1,'omitnan'); normTR2show_D=mean(normTR(dendROI,:),1,'omitnan');
normTR2show_S=get_bandstop(normTR2show_S,1000,[70 140]);
normTR2show_S_CS=normTR2show_S; normTR2show_S_CS(result2show.CStrace==0)=NaN;
normTR2show_D_CS=normTR2show_D; normTR2show_D_CS(result2show.CStrace==0)=NaN;
plot(normTR2show_S,'color',[0.4 0.4 0.4]); hold all
plot(normTR2show_S_CS,'color',[1 0 0]); hold all

result2show=OpResult{finterest_DenRPw(showneuron),finterest_DenRPr(showneuron),finterest_DenRPn(showneuron)};
ROIvec=result2show.maintrunkROI;
normTR=result2show.normTraces./result2show.F0_PCA;
[S]=get_STA(normTR,result2show.spike(1,:),0,0);
normTR=normTR./S(1);
somROI=1; dendROI=10;
normTR2show_S=mean(normTR(somROI,:),1,'omitnan'); normTR2show_D=mean(normTR(dendROI,:),1,'omitnan');
normTR2show_S=get_bandstop(normTR2show_S,1000,[20 50]);
normTR2show_S_CS=normTR2show_S; normTR2show_S_CS(result2show.CStrace==0)=NaN;
normTR2show_D_CS=normTR2show_D; normTR2show_D_CS(result2show.CStrace==0)=NaN;
plot(normTR2show_S+scaleoffset,'color',[0.4 0.4 0.4]); hold all
plot(normTR2show_S_CS+scaleoffset,'color',[1 0 0]); hold all

nexttile([1 1])
Boxplot_wPoints2(CSfraction,[0 0 0])
%%


    % Ramp_period=find(bwblue==max(bwblue));
    % Rheobase_frame=find(max(Result.spike(:,Ramp_period),[],1))+Ramp_period(1)-1;
    % Rheobase_blue=Result.Blue(round(mean(Rheobase_frame(1:3))));
    % %Rheobase_blue=1;
    % 
    % firstpulst_bw=max(bwblue)-5;
    % for b=1:5
    % pulse_period=find(bwblue==firstpulst_bw+b-1);    
    % pulseBlue_rb(b,f,r)=median(Result.Blue(pulse_period))/Rheobase_blue;
    % pulseFR(b,:,f,r)=sum(Result.SpClass(1:2,pulse_period),2)';
    % end
    % 
    % RampFR(:,:,f,r)=NaN(length(rheo_bin),2);
    % Ramp_period=find(bwblue==max(bwblue));    
    % binnedRheo_ramp=ceil(Result.Blue(Ramp_period)/Rheobase_blue/(rheo_bin(2)-rheo_bin(1)));
    % for b=binnedRheo_ramp
    % bin_frm=find(binnedRheo_ramp==b);
    % RampFR(b,:,f,r)=sum(Result.SpClass(1:2,bin_frm+Ramp_period(1)-1),2)';
    % end
    % end
% 
% %% Update of dendStim/SomaStim (2025.03.20)
% 
% 
% SomRP_path=[15 22 54 72 86 92 4];
% DenRP_path=[13 23 55 74 88 91 3];
% path_cat=[SomRP_path; DenRP_path];
% rheo_bin=[0:1.3:10];
% pulseFR=[]; RampFR=[]; pulseBlue_rb=[];
% SpikeMatInfo=[];
% for f=1:size(path_cat,2)
%     for r=1:2 %soma and dendrite
%     SpikeMatInfo{f,r}=[];   
%     i=path_cat(r,f);
%     load(fullfile(fpath{i},'OP_Result'))
%     bwblue=bwlabel(Result.Blue);
%     Ramp_period=(bwblue==max(bwblue));
%     Rheobase_frame=find(max(Result.spike(:,Ramp_period),[],1))+find(Ramp_period,1)-1;
%     Rheobase_blue=Result.Blue(round(mean(Rheobase_frame(1:3))));
%     %Rheobase_blue=1;
% 
%     SP_trace=[Result.spike(1,:).*(1-Result.CStrace);Result.spike(1,:).*Result.CStrace];
%     SP_trace=SP_trace.*(1-Ramp_period).*double(Result.Blue>0);
%     BlueOnset=find((double(Result.Blue(2:end)>0)-double(Result.Blue(1:end-1)>0))>0);
% 
%     for stype=1:2 %SS and CS
%     Sptime=find(SP_trace(stype,:));
%     SptimeMod=mod(Sptime-BlueOnset(1),BlueOnset(2)-BlueOnset(1));
%     SpikeMatInfo{f,r}=[SpikeMatInfo{f,r}; [Result.Blue(Sptime)'/Rheobase_blue SptimeMod' ones(length(Sptime),1)*stype]];
%     end
%     end
% end
% 
% 
% 
% 
% figure(15); clf; cmap=[0.4 0.4 0.4; 1 0 0]; tiledlayout(1,2);
% RheoEdge=[0 1 15];
% TimeEdge=[0:100:500];
% % for r=1:2
% %     nexttile([1 1])
% % ShowDat=cell2mat(SpikeMatInfo(:,r));
% % for stype=1:2
% %     scatter(ShowDat(ShowDat(:,3)==stype,2),ShowDat(ShowDat(:,3)==stype,1),10,cmap(stype,:),'o','filled'); hold all
% % end
% % end
% 
% BinRheoRatio=cellfun(@(x) histcounts(x(x(:,3)==2,1),RheoEdge)./histcounts(x(:,1),RheoEdge),SpikeMatInfo,'UniformOutput',false);
% BinRheoRatio_weak=cellfun(@(x) x(1),BinRheoRatio,'UniformOutput',false);
% BinRheoRatio_strong=cellfun(@(x) x(2),BinRheoRatio,'UniformOutput',false);
% nexttile([1 1]);
% Boxplot_wPoints(cell2mat(BinRheoRatio_weak(:,1)),cell2mat(BinRheoRatio_weak(:,2)),[1 0 0],'pairwise')
% set(gca,'XTick',[1 2],'XTickLabel',{'Soma Stim.','Dend. Stim.'})
% ylabel('Fraction of complex spike')
% nexttile([1 1]);
% p=Boxplot_wPoints(cell2mat(BinRheoRatio_strong(:,1)),cell2mat(BinRheoRatio_strong(:,2)),[1 0 0],'pairwise')
% set(gca,'XTick',[1 2],'XTickLabel',{'Soma Stim.','Dend. Stim.'})
% ylabel('Fraction of complex spike')
% 
% 
% 
%     firstpulst_bw=max(bwblue)-5;
%     for b=1:5
%     pulse_period=find(bwblue==firstpulst_bw+b-1);    
%     pulseBlue_rb(b,f,r)=median(Result.Blue(pulse_period))/Rheobase_blue;
%     pulseFR(b,:,f,r)=sum(Result.SpClass(1:2,pulse_period),2)';
%     end
% 
%     RampFR(:,:,f,r)=NaN(length(rheo_bin),2);
%     Ramp_period=find(bwblue==max(bwblue));    
%     binnedRheo_ramp=ceil(Result.Blue(Ramp_period)/Rheobase_blue/(rheo_bin(2)-rheo_bin(1)));
%     for b=binnedRheo_ramp
%     bin_frm=find(binnedRheo_ramp==b);
%     RampFR(b,:,f,r)=sum(Result.SpClass(1:2,bin_frm+Ramp_period(1)-1),2)';
%     end
%     end
% end
% 
% 
% figure(27); clf; cmap=distinguishable_colors(6); bin_width=1.3;
% tiledlayout(6,2);
% load(fullfile(fpath{SomRP_path(5)},'OP_Result'));
% tr_somStim=rescale(Result.normTraces(1,:));
% tr_cssomStim=Result.CStrace;
% load(fullfile(fpath{DenRP_path(6)},'OP_Result'));
% tr_ddStim=rescale(Result.normTraces(1,:));
% tr_csddStim=Result.CStrace;
% ax1=nexttile([3 2]);
% plot(tr_somStim,'color','k'); hold all
% show_cs=tr_somStim.*tr_cssomStim; show_cs(tr_cssomStim==0)=NaN;
% plot(show_cs,'color',[0.7 0 0.2]); hold all
% plot(tr_ddStim+1,'color','k'); hold all
% show_cs=tr_ddStim.*tr_csddStim; show_cs(tr_csddStim==0)=NaN;
% axis off
% text(-500,0.7,'Soma Stimulation','FontSize',13)
% text(-500,1.7,'Distal dendrite Stimulation','FontSize',13)
% text(9500,2,'Complex spikes','FontSize',13,'color',[0.7 0 0.2])
% plot(show_cs+1,'color',[0.7 0 0.2]); hold all
% ax2=nexttile([1 2]);
% plot(Result.Blue)
% linkaxes([ax1 ax2],'x')
% 
% nexttile([2 1])
% CS_fracPulseMatSom=squeeze(pulseFR(:,2,:,1)./sum(pulseFR(:,:,:,1),2));
% CS_fracPulseMatDD=squeeze(pulseFR(:,2,:,2)./sum(pulseFR(:,:,:,2),2));
% markers = {'o', 's', 'd', '^', 'v', '>', '<', 'p', 'h', '+', '*'};
% for f=1:size(pulseFR,3)
% plot(squeeze(pulseBlue_rb(:,f,1))',CS_fracPulseMatSom(:,f),'marker',markers{f},'color',cmap(1,:),'linestyle','none'); hold all
% plot(squeeze(pulseBlue_rb(:,f,2))',CS_fracPulseMatDD(:,f),'marker',markers{f},'color',cmap(2,:),'linestyle','none'); hold all
% end
%     clear Ms Ss
% for b=1:max(max(ceil(pulseBlue_rb(:,:,1)/bin_width)))
%     Ms(b)=mean(CS_fracPulseMatSom(b==ceil(pulseBlue_rb(:,:,1)/bin_width)),'omitnan');
%     Ss(b)=std(CS_fracPulseMatSom(b==ceil(pulseBlue_rb(:,:,1)/bin_width)),'omitnan');
% end
% clear Md Sd
% for b=1:max(max(ceil(pulseBlue_rb(:,:,2)/bin_width)))
%     Md(b)=mean(CS_fracPulseMatDD(b==ceil(pulseBlue_rb(:,:,2)/bin_width)),'omitnan');
%     Sd(b)=std(CS_fracPulseMatDD(b==ceil(pulseBlue_rb(:,:,2)/bin_width)),'omitnan');
% end
% errorbar([1:max(max(ceil(pulseBlue_rb(:,:,1)/bin_width)))],Ms,Ss,'color',cmap(1,:),'LineWidth',1.5); hold all
% errorbar([1:max(max(ceil(pulseBlue_rb(:,:,2)/bin_width)))],Md,Sd,'color',cmap(2,:),'LineWidth',1.5);
% 
% xlabel('Optical Rheobase')
% ylabel('Fraction of complex spike')
% title('Pulse Stimulation')
% 
% nexttile([2 1])
% CS_fracRampMatSom=squeeze(RampFR(:,2,:,1)./sum(RampFR(:,:,:,1),2));
% CS_fracRampMatDD=squeeze(RampFR(:,2,:,2)./sum(RampFR(:,:,:,2),2));
% % for f=1:size(RampFR,3)
% % plot(rheo_bin,CS_fracRampMatSom,'color',cmap(1,:)); hold all
% % plot(rheo_bin,CS_fracRampMatDD,'color',cmap(2,:))
% % end
% errorbar(rheo_bin,mean(CS_fracRampMatSom,2,'omitnan'),std(CS_fracRampMatSom,0,2,'omitnan'),'color',cmap(1,:),'LineWidth',1.5); hold all
% errorbar(rheo_bin,mean(CS_fracRampMatDD,2,'omitnan'),std(CS_fracRampMatDD,0,2,'omitnan'),'color',cmap(2,:),'LineWidth',1.5)
% xlabel('Optical Rheobase')
% ylabel('Fraction of complex spike')
% title('Ramp Stimulation')
% 
















%% Representative Pulse ramp stim, soma & apical dend stim (Figure 4)
% show example bAP amplitude/AUC image

ExampleNeuron=6; 
ROIvec=setdiff(OpResult{1,1,ExampleNeuron}.maintrunkROI, []); 
ftprint2show=OpResult{1,1,ExampleNeuron}.ftprnt(:,:,ROIvec); ftprintall=OpResult{1,1,ExampleNeuron}.ftprnt;
nROI=size(OpResult{1,1,ExampleNeuron}.ftprnt,3); ftprint2show_trace=get_coord(ftprint2show);
caxshow=[-2.5 6.5];
%SWCpoints=OpResult{1,1,ExampleNeuron}.SWC; SWCpoints(1,3)=70; 
cmap_ExTr=gen_colormap(Plasma,10);
nROI=size(OpResult{1,1,ExampleNeuron}.normTraces,1);
ftprint=OpResult{1,1,ExampleNeuron}.ftprnt(:,:,OpResult{1,1,ExampleNeuron}.dist_order)>0;
dax=OpResult{1,1,ExampleNeuron}.dendaxis;
% 
% figure(25); clf;
% nexttile([1 1])
% boundary_s=bwboundaries(max(ftprint2show(:,:,somROI),[],3)); boundary_ap=bwboundaries(max(ftprint2show(:,:,dendROI),[],3));
% boundary_blue=bwboundaries(OpResult{1,1,ExampleNeuron}.BlueDMDimg);
% showScaleScatter([1:nROI], SWCpoints, ftprintall,repmat([0 0 0],256,1)); hold all
% plot(boundary_s{1}(:,1),boundary_s{1}(:,2),'color',cmap_ExTr(1,:),'linewidth',2);
% plot(boundary_ap{1}(:,1),boundary_ap{1}(:,2),'color',cmap_ExTr(end,:),'linewidth',2);
% plot(boundary_blue{1}(:,1),boundary_blue{1}(:,2),'color',[0 0.5 1],'linewidth',2);
% plot(ftprint2show_trace(:,2),ftprint2show_trace(:,1),'r','linewidth',2);
% drawScaleBar(100/OpResult{1,1,ExampleNeuron}.pixelsize,'vertical')
% view([210 90])

figure(29); clf;
tiledlayout(5,1); % show example bAP amplitude/AUC image
for wvf=[1 3]
normTr=OpResult{wvf,1,ExampleNeuron}.normTraces./OpResult{wvf,1,ExampleNeuron}.F0_PCA;
normTr=normTr(ROIvec,:);
%normTr=pcafilterTrace(normTr,[1 2 4 6 7 8 9]);
%sublarge=get_subthreshold(normTr,OpResult{wvf,ExampleRep,ExampleNeuron}.Blue>0,300,1000);
Tr2show=[mean(normTr(1,:),1); mean(normTr(end,:),1)];
Tr2showCS=Tr2show;
Tr2showCS(:,OpResult{wvf,1,ExampleNeuron}.CStrace==0)=NaN;
normTr_interp=interp1(dax(ROIvec),normTr,linspace(dax(ROIvec(1)),dax(ROIvec(end)),10));

ax1=nexttile([2 1]);
plot(Tr2show'+[0 5],'color',[0.4 0.4 0.4]); hold all
plot(Tr2showCS'+[0 5],'color',[1 0 0]); hold all
end


%%

for wvf=[2 4]

    figure(15);
    ax2=[ax2 nexttile([1 1])];
    sp2mean=LabelMat_SP.SpikeOrder==1 & LabelMat_SP.StimRegion==4 & LabelMat_SP.NeuronID==5;
    STA2showall=mean(cat(3,AlignedbAPall_SP{LabelMat_SP.SpikeOrder==1 & LabelMat_SP.StimRegion==wvf & LabelMat_SP.NeuronID==5}),3,'omitnan');
    STA2show=STA2showall(ismember(OpResult{1,1,show_neuron}.dist_order,ROIvec)',:);
    
    imagesc([-nTau(1):nTau(2)],ROIvec,STA2show_interp,cax);
    set_kymoYtick(dax(ROIvec));

    ax2=[ax2 nexttile([1 1])];
    plot([-nTau(1):nTau(2)],STA2show(1,:),'color',cmap_ExTr(1,:)); hold all;
    plot([-nTau(1):nTau(2)],STA2show(end,:),'color',cmap_ExTr(end,:));
    ylim(cax); xlim([-50 50]);
    axis off;
    drawScaleBar(1,'vertical');
    drawScaleBar(50,'horizontal');
    xlabel('Peri-spike time (ms)')

    figure(16);
    for t=1:length(t2show)
        nexttile([1 1]);
        V2show=STA2showall(:,t2show(t)+nTau(1));
        %V2showcolor=vec2cmap(V2show,turbo,cax);
        showScaleScatter(V2show, SWCpoints, ftprintall(:,:,OpResult{1,1,show_neuron}.dist_order),'turbo',cax);
        axis tight;
        view([180 -90])
    end
end
figure(15);
linkaxes(ax2,'x'); xlim([-50 50]);


figure(25); clf; cmap_ExTr=gen_colormap(Plasma,10);
for ss=1:size(ShowSession,1)
    nexttile([1 1])
    ExampleSession=ShowSession(ss,1); wvf=ShowSession(ss,2);
    ExampleNeuron=unique(LabelMat.NeuronID(find(LabelMat.SessionID==ExampleSession & LabelMat.StimRegion==wvf)));
    ExampleRep=unique(LabelMat.RepeatID(find(LabelMat.SessionID==ExampleSession & LabelMat.StimRegion==wvf)));
    BlueBin=OpResult{wvf,1,ExampleNeuron}.BlueDMDimg;
    SWCpoints=OpResult{1,1,ExampleNeuron}.SWC;
    SWCpoints(1,3)=50;
    ROIvec=[1:nROI]; %ROIvec(ShowROI([35 37 39]))=2; ROIvec(ShowROI([3]))=1;
    ftprint_s=max(ftprint(:,:,ShowROI(SomaROI)),[],3); ftprint_ap=max(ftprint(:,:,ShowROI(ApicalROI)),[],3);
    boundary_s=bwboundaries(ftprint_s); boundary_ap=bwboundaries(ftprint_ap); boundary_blue=bwboundaries(BlueBin);
    showScaleScatter(ROIvec, SWCpoints, ftprint,repmat([0 0 0],256,1)); hold all
    plot(boundary_s{1}(:,1),boundary_s{1}(:,2),'color',cmap_ExTr(1,:),'linewidth',2);
    plot(boundary_ap{1}(:,1),boundary_ap{1}(:,2),'color',cmap_ExTr(end,:),'linewidth',2);
    plot(boundary_blue{1}(:,1),boundary_blue{1}(:,2),'color',[0 0.5 1],'linewidth',2);
    drawScaleBar(100/OpResult{1,1,ExampleNeuron}.pixelsize,'vertical')
    view([180 90])
end


%%

for wvf=[2 4]

    figure(15);
    ax2=[ax2 nexttile([1 1])];
    sp2mean=LabelMat_SP.SpikeOrder==1 & LabelMat_SP.StimRegion==4 & LabelMat_SP.NeuronID==5;
    STA2showall=mean(cat(3,AlignedbAPall_SP{LabelMat_SP.SpikeOrder==1 & LabelMat_SP.StimRegion==wvf & LabelMat_SP.NeuronID==5}),3,'omitnan');
    STA2show=STA2showall(ismember(OpResult{1,1,show_neuron}.dist_order,ROIvec)',:);
    
    imagesc([-nTau(1):nTau(2)],ROIvec,STA2show_interp,cax);
    set_kymoYtick(dax(ROIvec));

    ax2=[ax2 nexttile([1 1])];
    plot([-nTau(1):nTau(2)],STA2show(1,:),'color',cmap_ExTr(1,:)); hold all;
    plot([-nTau(1):nTau(2)],STA2show(end,:),'color',cmap_ExTr(end,:));
    ylim(cax); xlim([-50 50]);
    axis off;
    drawScaleBar(1,'vertical');
    drawScaleBar(50,'horizontal');
    xlabel('Peri-spike time (ms)')

    figure(16);
    for t=1:length(t2show)
        nexttile([1 1]);
        V2show=STA2showall(:,t2show(t)+nTau(1));
        %V2showcolor=vec2cmap(V2show,turbo,cax);
        showScaleScatter(V2show, SWCpoints, ftprintall(:,:,OpResult{1,1,show_neuron}.dist_order),'turbo',cax);
        axis tight;
        view([180 -90])
    end
end
figure(15);
linkaxes(ax2,'x'); xlim([-50 50]);


figure(25); clf; cmap_ExTr=gen_colormap(Plasma,10);
for ss=1:size(ShowSession,1)
    nexttile([1 1])
    ExampleSession=ShowSession(ss,1); wvf=ShowSession(ss,2);
    ExampleNeuron=unique(LabelMat.NeuronID(find(LabelMat.SessionID==ExampleSession & LabelMat.StimRegion==wvf)));
    ExampleRep=unique(LabelMat.RepeatID(find(LabelMat.SessionID==ExampleSession & LabelMat.StimRegion==wvf)));
    BlueBin=OpResult{wvf,1,ExampleNeuron}.BlueDMDimg;
    SWCpoints=OpResult{1,1,ExampleNeuron}.SWC;
    SWCpoints(1,3)=50;
    ROIvec=[1:nROI]; %ROIvec(ShowROI([35 37 39]))=2; ROIvec(ShowROI([3]))=1;
    ftprint_s=max(ftprint(:,:,ShowROI(SomaROI)),[],3); ftprint_ap=max(ftprint(:,:,ShowROI(ApicalROI)),[],3);
    boundary_s=bwboundaries(ftprint_s); boundary_ap=bwboundaries(ftprint_ap); boundary_blue=bwboundaries(BlueBin);
    showScaleScatter(ROIvec, SWCpoints, ftprint,repmat([0 0 0],256,1)); hold all
    plot(boundary_s{1}(:,1),boundary_s{1}(:,2),'color',cmap_ExTr(1,:),'linewidth',2);
    plot(boundary_ap{1}(:,1),boundary_ap{1}(:,2),'color',cmap_ExTr(end,:),'linewidth',2);
    plot(boundary_blue{1}(:,1),boundary_blue{1}(:,2),'color',[0 0.5 1],'linewidth',2);
    drawScaleBar(100/OpResult{1,1,ExampleNeuron}.pixelsize,'vertical')
    view([180 90])
end

%%

%%
figure(41); clf;
ShowInd=LabelMat.StimRegion==1;
TransmissionRate_order=[];
for ord=1:5
    TransmissionRate_order{ord}=transmissionRate(ShowInd & LabelMat.SpikeOrder==ord);
end
p_values = Violin_wPoints(TransmissionRate_order, autumn(5));
%%
figure(222);
Neurons=unique(LabelMat.NeuronID);
for n=Neurons'
    clf;
    tiledlayout('vertical');
    nexttile([1 1])
    imshow2(OpResult{1,1,n}.ref_im,[])
    title([num2str(n)]);
    for wvf=[3]
        maxr=sum(1-cellfun(@isempty,OpResult(wvf,:,n)));
        for r=1:maxr
            ax2=[ax2 nexttile([1 1])];
            ntr=OpResult{wvf,r,n}.normTraces./OpResult{wvf,r,n}.F0_PCA;
            %ntr=ntr(OpResult{wvf,r,n}.dist_order,:);
            %imagesc(ntr,[-0.01 0.02]);
            plot(ntr(1,:))
            title(['stim region: ' num2str(wvf)]);
        end
    end
    colormap(turbo);
    %    linkaxes(ax2,'x');
    ss=input('yes?\n');
end
%% Ramp (Soma Stim)
nTau=[-20 20]; nTauPeriSp=[-3:7]; nTauPeriSp2=[-3:3]; conductionspeed=170; %um/ms
g=ones(1,4); SP_STA=[]; distFromSoma=[]; LabelMat=[];
dendriteaxis_bin=[-200 -100 -50 50 100:50:500];
dendriteaxis_bin2=[-200 -50 50 150 600];
timebin=[0 20 40 60 100 200:100:500];
perisomadist=[-50 50]; cmap=[0 0 0; 1 0 0];  cmap_light=[0.7 0.7 0.7; 1 0.7 0.7];
%BlueStimN=[2 5 1;3 3 1;3 4 1;4 2 1;4 3 1;4 4 1;5 3 1;5 5 1;6 4 1;6 5 1;6 3 2;6 4 2;6 5 2;7 1 1;7 2 1;7 3 1;];
BlueStimN=[2 6 1;3 6 1;4 6 1;5 6 1;6 6 1;6 6 2;7 6 1];

SpikeAmp=[]; SpikeAUC=[]; g=1;
for n=1:size(BlueStimN,1) % neuron
    cax=[];
    Neuron=BlueStimN(n,1);
    rep=BlueStimN(n,3);
    BluePulseN=BlueStimN(n,2);
    wvf=1;
    if ~isempty(OpResult{wvf,rep,Neuron})

        normTr=OpResult{wvf,rep,Neuron}.normTraces./OpResult{wvf,rep,Neuron}.F0_PCA;
        noi=[1:size(OpResult{wvf,rep,Neuron}.normTraces,1)];
        %noi=OpResult{wvf,rep,Neuron}.maintrunkROI;
        normTr=normTr(noi,:);

        nROI=length(noi);
        nTime=size(OpResult{wvf,rep,Neuron}.normTraces,2);
        dOrder=OpResult{wvf,rep,Neuron}.dist_order(noi);
        %[~, dOrder]=sort(OpResult{wvf,rep,Neuron}.dist_order(noi),'ascend');
        Classvec = get_Class2index(OpResult{wvf,rep,Neuron}.SpClass([1 2],:));
        som_spike = Classvec .* OpResult{wvf,rep,Neuron}.spike(1,:);

        bwBlue=bwlabel(OpResult{wvf,rep,Neuron}.Blue);
        bwBlue((bwBlue-(max(bwBlue)-5))<0)=0;
        bwBlue=bwlabel(bwBlue>0);

        [STAbAP]=get_STA(normTr,max(som_spike>0,[],1),-nTau(1),nTau(2));
        [~, tmax]=max(STAbAP(:,-nTau(1)+[-2:4]),[],2);
        tmax=tmax-nTau(1)-3;

        SpikeinBlue=max(som_spike,[],1).*(bwBlue==BluePulseN);
        SpikeinBlue_ind=find(SpikeinBlue);
        BlueOnset=find(bwBlue==BluePulseN,1);
        SubLarge=get_subthreshold(normTr, max([som_spike; bwBlue>0],[],1), 400,3000);
        normTrSub=normTr-SubLarge;

        Rheobase=mean(OpResult{wvf,rep,Neuron}.Blue(SpikeinBlue_ind(1:2)));


        Dsign=ones(1,size(OpResult{wvf,rep,Neuron}.interDendDist,2));
        Dsign(OpResult{wvf,rep,Neuron}.dist_order(1:find(OpResult{wvf,rep,Neuron}.dist_order==1)-1))=-1;
        dendaxis=OpResult{wvf,rep,Neuron}.interDendDist(1,:).*Dsign;
        dendaxis=dendaxis(noi);
        distFromSoma=dendaxis*OpResult{wvf,rep,Neuron}.pixelsize;
        perisomaInd=find(distFromSoma'>perisomadist(1) & distFromSoma'<perisomadist(2));

        SpikeMat=permute(reshape(normTr(:,SpikeinBlue_ind'+nTauPeriSp),nROI,[],length(nTauPeriSp)),[1 3 2]); %1: ROI, 2:time, 3:event
        %SpikeAmp{g}=[[NaN SpikeinBlue_ind-BlueOnset]; [distFromSoma' squeeze(max(SpikeMat,[],2))]];
        RheoBase_trace=OpResult{wvf,rep,Neuron}.Blue(SpikeinBlue_ind)/Rheobase;
        SpikeAmp{g}=[[NaN RheoBase_trace]; [distFromSoma' squeeze(max(SpikeMat,[],2))]]; %RheoBase
        SpikeAmp{g}(2:end,2:end)=SpikeAmp{g}(2:end,2:end)./mean(SpikeAmp{g}(perisomaInd+1,2:end),[1]);
        SpikeAmp{g}(2:end,:)=SpikeAmp{g}(dOrder+1,:);
        [~, reorder_dist]=sort(SpikeAmp{g}(2:end,1),'ascend');
        SpikeAmp{g}(2:end,:)=SpikeAmp{g}(reorder_dist,:);

        SPdelay=round(abs(distFromSoma)/conductionspeed);
        sub=[];
        for r=1:size(normTr,1)
            sp_vec=ind2vec(size(normTr,2),find(som_spike)+SPdelay(r),1,0);
            [~, sub(r,:)]=get_subthreshold(normTr(r,:),sp_vec,2*(2+SPdelay(r))+1,15);
        end
        normTr_sub=normTr-sub;

        SpikeMat_sub=permute(reshape(normTr_sub(:,SpikeinBlue_ind'+nTauPeriSp),nROI,[],length(nTauPeriSp)),[1 3 2]); %1: ROI, 2:time, 3:event
        AUC=squeeze(sum(SpikeMat_sub,2,'omitnan'));
        %
        % for r=1:size(SpikeMat,1)
        %     for s=1:size(SpikeMat,3)
        %     AUC(r,s)=get_AUC(squeeze(SpikeMat(r,:,s)),-nTauPeriSp(1)+1+SPdelay(r),3,3);
        %     end
        % end

        %SpikeAUC{g}=[[NaN SpikeinBlue_ind-BlueOnset]; [distFromSoma' AUC]];
        SpikeAUC{g}=[[NaN RheoBase_trace]; [distFromSoma' AUC]]; %RheoBase
        SpikeAUC{g}(2:end,2:end)=SpikeAUC{g}(2:end,2:end)./mean(SpikeAUC{g}(perisomaInd+1,2:end),[1]);
        SpikeAUC{g}(2:end,:)=SpikeAUC{g}(dOrder+1,:);

        g=g+1;
    end
end

figure(14); clf; %tiledlayout(1,2);
nexttile([1 1]);
[binnedZ binX binY]=show3Dbinning(SpikeAmp, 11, 8, 'image'); hold all
shading flat
allX=[]; allY=[]; allZ=[];
for n = 1:numel(SpikeAmp)
    timeAxis = SpikeAmp{n}(1, 2:end);
    xAxis = SpikeAmp{n}(2:end, 1);
    zData = SpikeAmp{n}(2:end, 2:end);
    [X, Y] = meshgrid(timeAxis, xAxis);
    allX = [allX; X(:)];
    allY = [allY; Y(:)];
    allZ = [allZ; zData(:)];
end
%scatter3(allX,allY,allZ,20,[0.7 0.7 0.7],'filled')
set(gca, 'YDir', 'reverse');
colormap("turbo")
%zlim([0.3 1.6])
xlabel('Optical Rheobase')
ylabel('Distance from soma (\mum)')
title('Normalized Spike Amplitude')


%% Representative pulse/Ramp

normTr=OpResult{1,1,7}.normTraces./OpResult{1,1,7}.F0_PCA;
normTr=normTr(OpResult{1,1,7}.dist_order,:);
normTr=pcafilterTrace(normTr,[1 2 3 4 9 10]);
cax=[-0.005 0.02];

Dsign=ones(1,size(OpResult{1,1,7}.interDendDist,2));
Dsign(OpResult{1,1,7}.dist_order(1:find(OpResult{1,1,7}.dist_order==1)-1))=-1;
dendaxis=OpResult{1,1,7}.interDendDist(1,:).*Dsign;
dendaxis=dendaxis(OpResult{1,1,7}.dist_order);

SomaROI=[6 10 16];
distalROI=[29 31 33];
dendriteaxis_bin=[-100:10:220];
taxis=[-99:600];
[~, sortind]=sort(dendaxis,'ascend');
[Xq, Yq] = meshgrid([taxis], dendriteaxis_bin);
normTr_interp = interp2([taxis], dendaxis(sortind), normTr(sortind,2901:3600), Xq, Yq, 'linear');
Sptime=find(OpResult{1,1,7}.spike(1,2901:3600)>0)-100;

figure(14); clf; tiledlayout(4,1);
ax1=nexttile([1 1]);
% pcolor(Xq,Yq,normTr_interp); hold all
plot(taxis,mean(normTr(SomaROI,2901:3600),1)); hold all
plot(taxis,mean(normTr(distalROI,2901:3600),1)); hold all
%plot(Sptime,normTr(SomaROI,Sptime+3000),'kv','markersize',3); axis tight off;
%set(gca,'YDir','reverse')
%caxis([cax])
%shading flat
%colormap(turbo);
ax2=nexttile([1 1]);
plot(taxis,OpResult{1,1,7}.Blue(2901:3600),'color',[0 0.6 1]); axis off;
linkaxes([ax1 ax2],'x')

N=7; wvf=3; rep=2;
SomaROI=[6 10 16];
distalROI=[29 31 33];
normTr=OpResult{wvf,rep,N}.normTraces./OpResult{wvf,rep,N}.F0_PCA;
normTr=normTr(OpResult{wvf,rep,N}.dist_order,:);
%normTr=pcafilterTrace(normTr,[1 2 3 4 9 10]);
cax=[-0.005 0.02];

Dsign=ones(1,size(OpResult{wvf,rep,N}.interDendDist,2));
Dsign(OpResult{wvf,rep,N}.dist_order(1:find(OpResult{wvf,rep,N}.dist_order==1)-1))=-1;
dendaxis=OpResult{wvf,rep,N}.interDendDist(1,:).*Dsign;
dendaxis=dendaxis(OpResult{wvf,rep,N}.dist_order);

dendriteaxis_bin=[-100:10:220];
taxis=[1:9999]; time_show=[1:9999];
[~, sortind]=sort(dendaxis,'ascend');
[Xq, Yq] = meshgrid([taxis], dendriteaxis_bin);
normTr_interp = interp2([taxis], dendaxis(sortind), normTr(sortind,time_show), Xq, Yq, 'linear');
Sptime=find(OpResult{wvf,rep,N}.spike(1,time_show)>0)-100;
%figure(15); clf; tiledlayout(2,1);
ax3=nexttile([1 1]);
plot(taxis,mean(normTr(SomaROI,time_show),1)); hold all
plot(taxis,mean(normTr(distalROI,time_show),1)); hold all
%pcolor(Xq,Yq,normTr_interp); hold all
%plot(Sptime,-110,'k','marker','v','markersize',3); axis tight off;
%set(gca,'YDir','reverse')
%caxis([cax]);
%shading flat
%colormap(turbo);
ax4=nexttile([1 1]);
plot(taxis,OpResult{wvf,rep,N}.Blue(time_show),'color',[0 0.6 1]); axis off;
linkaxes([ax4 ax3],'x')

%%
figure(15); clf;
% show dStim, CS examples
BlueStimN=[6 1;6 2;7 1]; scale=0.03;
wvf=3;
for n=1:size(BlueStimN,1) % neuron
    cax=[];
    Neuron=BlueStimN(n,1);
    rep=BlueStimN(n,2);

    normTr=OpResult{wvf,rep,Neuron}.normTraces./OpResult{wvf,rep,Neuron}.F0_PCA;
    tax=[1:size(normTr,2)];
    CStr=normTr(1,:);
    CStr(OpResult{wvf,rep,Neuron}.CStrace==0)=NaN;
    plot(tax,normTr(1,:)+scale*(n-1),'k');
    hold all
    plot(tax,CStr+scale*(n-1),'color',[0.7 0 0.1]);
end
%% Ramp (Soma Stim)
nTau=[-70:70]; nTauPeriSp=[-3:7]; nTauPeriSp2=[-3:3]; conductionspeed=170; %um/ms
g=ones(1,4); SP_STA=[]; distFromSoma=[];
dendriteaxis_bin=[-80:80:600];
perisomadist=[-50 50]; cmap=[0 0 0; 1 0 0];  cmap_light=[0.7 0.7 0.7; 1 0.7 0.7];
%BlueStimN=[2 5 1;3 3 1;3 4 1;4 2 1;4 3 1;4 4 1;5 3 1;5 5 1;6 4 1;6 5 1;6 3 2;6 4 2;6 5 2;7 1 1;7 2 1;7 3 1;];
BlueStimN=[2 6 1;3 6 1;4 6 1;5 6 1;6 6 1;6 6 2;7 6 1];

SpikeAmp=[]; SpikeAUC=[]; g=1;
for n=1:size(BlueStimN,1) % neuron
    cax=[];
    Neuron=BlueStimN(n,1);
    rep=BlueStimN(n,3);
    BluePulseN=BlueStimN(n,2);
    wvf=1;
    if ~isempty(OpResult{wvf,rep,Neuron})

        normTr=OpResult{wvf,rep,Neuron}.normTraces./OpResult{wvf,rep,Neuron}.F0_PCA;
        %noi=[1:size(OpResult{wvf,rep,Neuron}.normTraces,1)];
        noi=OpResult{wvf,rep,Neuron}.maintrunkROI;
        normTr=normTr(noi,:);

        nROI=length(noi);
        nTime=size(OpResult{wvf,rep,Neuron}.normTraces,2);
        [~, dOrder]=sort(OpResult{wvf,rep,Neuron}.dist_order(noi),'ascend');
        som_spike=max(OpResult{wvf,rep,Neuron}.SpClass([1 2],:),[],1);
        bwBlue=bwlabel(OpResult{wvf,rep,Neuron}.Blue);
        bwBlue((bwBlue-(max(bwBlue)-5))<0)=0;
        bwBlue=bwlabel(bwBlue>0);
        SpikeinBlue=som_spike.*(bwBlue==BluePulseN);
        SpikeinBlue_ind=find(SpikeinBlue);
        BlueOnset=find(bwBlue==BluePulseN,1);
        Rheobase=mean(OpResult{wvf,rep,Neuron}.Blue(SpikeinBlue_ind(1:2)));


        normTr=normTr-prctile(normTr(:,BlueOnset+[-500:-200]),30,2);

        Dsign=ones(1,size(OpResult{wvf,rep,Neuron}.interDendDist,2));
        Dsign(OpResult{wvf,rep,Neuron}.dist_order(1:find(OpResult{wvf,rep,Neuron}.dist_order==1)-1))=-1;
        dendaxis=OpResult{wvf,rep,Neuron}.interDendDist(1,:).*Dsign;
        dendaxis=dendaxis(noi);
        distFromSoma=dendaxis*OpResult{wvf,rep,Neuron}.pixelsize;
        perisomaInd=find(distFromSoma'>perisomadist(1) & distFromSoma'<perisomadist(2));

        SpikeMat=permute(reshape(normTr(:,SpikeinBlue_ind'+nTauPeriSp),nROI,[],length(nTauPeriSp)),[1 3 2]); %1: ROI, 2:time, 3:event
        %SpikeAmp{g}=[[NaN SpikeinBlue_ind-BlueOnset]; [distFromSoma' squeeze(max(SpikeMat,[],2))]];
        RheoBase_trace=OpResult{wvf,rep,Neuron}.Blue(SpikeinBlue_ind)/Rheobase;
        SpikeAmp{g}=[[NaN RheoBase_trace]; [distFromSoma' squeeze(max(SpikeMat,[],2))]]; %RheoBase
        SpikeAmp{g}(2:end,2:end)=SpikeAmp{g}(2:end,2:end)./mean(SpikeAmp{g}(perisomaInd+1,2:end),[1]);
        SpikeAmp{g}(2:end,:)=SpikeAmp{g}(dOrder+1,:);
        [~, reorder_dist]=sort(SpikeAmp{g}(2:end,1),'ascend');
        SpikeAmp{g}(2:end,:)=SpikeAmp{g}(reorder_dist,:);

        SPdelay=round(abs(distFromSoma)/conductionspeed);
        sub=[];
        for r=1:size(normTr,1)
            sp_vec=ind2vec(size(normTr,2),find(som_spike)+SPdelay(r),1,0);
            [~, sub(r,:)]=get_subthreshold(normTr(r,:),sp_vec,2*(2+SPdelay(r))+1,15);
        end
        normTr_sub=normTr-sub;

        SpikeMat_sub=permute(reshape(normTr_sub(:,SpikeinBlue_ind'+nTauPeriSp),nROI,[],length(nTauPeriSp)),[1 3 2]); %1: ROI, 2:time, 3:event
        AUC=squeeze(sum(SpikeMat_sub,2,'omitnan'));
        %
        % for r=1:size(SpikeMat,1)
        %     for s=1:size(SpikeMat,3)
        %     AUC(r,s)=get_AUC(squeeze(SpikeMat(r,:,s)),-nTauPeriSp(1)+1+SPdelay(r),3,3);
        %     end
        % end

        %SpikeAUC{g}=[[NaN SpikeinBlue_ind-BlueOnset]; [distFromSoma' AUC]];
        SpikeAUC{g}=[[NaN RheoBase_trace]; [distFromSoma' AUC]]; %RheoBase
        SpikeAUC{g}(2:end,2:end)=SpikeAUC{g}(2:end,2:end)./mean(SpikeAUC{g}(perisomaInd+1,2:end),[1]);
        SpikeAUC{g}(2:end,:)=SpikeAUC{g}(dOrder+1,:);

        g=g+1;
    end
end

figure(14); clf; %tiledlayout(1,2);
nexttile([1 1]);
[binnedZ binX binY]=show3Dbinning(SpikeAmp, 11, 8, 'image'); hold all
shading flat
allX=[]; allY=[]; allZ=[];
for n = 1:numel(SpikeAmp)
    timeAxis = SpikeAmp{n}(1, 2:end);
    xAxis = SpikeAmp{n}(2:end, 1);
    zData = SpikeAmp{n}(2:end, 2:end);
    [X, Y] = meshgrid(timeAxis, xAxis);
    allX = [allX; X(:)];
    allY = [allY; Y(:)];
    allZ = [allZ; zData(:)];
end
%scatter3(allX,allY,allZ,20,[0.7 0.7 0.7],'filled')
set(gca, 'YDir', 'reverse');
colormap("turbo")
%zlim([0.3 1.6])
xlabel('Optical Rheobase')
ylabel('Distance from soma (\mum)')
title('Normalized Spike Amplitude')
%view(0,90)
%
% nexttile([1 1]);
% [binnedZ binX binY]=show3Dbinning(SpikeAUC, 11, 8, 'image'); hold all
% shading flat
% allX=[]; allY=[]; allZ=[];
%   for n = 1:numel(SpikeAUC)
%         timeAxis = SpikeAUC{n}(1, 2:end);
%         xAxis = SpikeAUC{n}(2:end, 1);
%         zData = SpikeAUC{n}(2:end, 2:end);
%         [X, Y] = meshgrid(timeAxis, xAxis);
%         allX = [allX; X(:)];
%         allY = [allY; Y(:)];
%         allZ = [allZ; zData(:)];
%     end
% %scatter3(allX,allY,allZ,20,[0.8 0.8 0.8],'filled')
% set(gca, 'YDir', 'reverse');
% colormap("turbo")
% %zlim([0.3 1.4])
% xlabel('Optical Rheobase')
% ylabel('Distance from soma (\mum)')
% title('Normalized Spike AUC')
% %view(0,90)

% %% Pulse (WF Stim)
% nTau=[-70:70]; nTauPeriSp=[-3:7]; nTauPeriSp2=[-3:3]; conductionspeed=170; %um/ms
% g=ones(1,4); SP_STA=[]; distFromSoma=[];
% dendriteaxis_bin=[-80:80:600];
% perisomadist=[-50 50]; cmap=[0 0 0; 1 0 0];  cmap_light=[0.7 0.7 0.7; 1 0.7 0.7];
% BlueStimN=[10 3 1;10 4 1;10 5 1;11 4 1;11 5 1;12 3 1;12 4 1;12 5 1;13 4 1;13 5 1; 14 5 1;15 4 1;15 5 1]; % Neuron number, Blue Pulse # , N th session
% %BlueStimN=[2 6 1;3 6 1;4 6 1;5 6 1;6 6 1;6 6 2;7 6 1];
%
% SpikeAmp=[]; SpikeAUC=[]; g=1;
% for n=1:size(BlueStimN,1) % neuron
%     cax=[];
%     Neuron=BlueStimN(n,1);
%     rep=BlueStimN(n,3);
%     BluePulseN=BlueStimN(n,2);
%     wvf=5;
%     if ~isempty(OpResult{wvf,rep,Neuron})
%
%         normTr=OpResult{wvf,rep,Neuron}.normTraces./OpResult{wvf,rep,Neuron}.F0_PCA;
%         %noi=[1:size(OpResult{wvf,rep,Neuron}.normTraces,1)];
%         noi=OpResult{wvf,rep,Neuron}.maintrunkROI;
%         normTr=normTr(noi,:);
%
%         nROI=length(noi);
%         nTime=size(OpResult{wvf,rep,Neuron}.normTraces,2);
%         [~, dOrder]=sort(OpResult{wvf,rep,Neuron}.dist_order(noi),'ascend');
%         som_spike=max(OpResult{wvf,rep,Neuron}.SpClass([1 2],:),[],1);
%         bwBlue=bwlabel(OpResult{wvf,rep,Neuron}.Blue);
%         bwBlue((bwBlue-(max(bwBlue)-5))<0)=0;
%         bwBlue=bwlabel(bwBlue>0);
%         SpikeinBlue=som_spike.*(bwBlue==BluePulseN);
%         SpikeinBlue_ind=find(SpikeinBlue);
%         BlueOnset=find(bwBlue==BluePulseN,1);
%         normTr=normTr-prctile(normTr(:,BlueOnset+[-500:-200]),30,2);
%
%         Dsign=ones(1,size(OpResult{wvf,rep,Neuron}.interDendDist,2));
%         Dsign(OpResult{wvf,rep,Neuron}.dist_order(1:find(OpResult{wvf,rep,Neuron}.dist_order==1)-1))=-1;
%         dendaxis=OpResult{wvf,rep,Neuron}.interDendDist(1,:).*Dsign;
%         dendaxis=dendaxis(noi);
%         distFromSoma=dendaxis*OpResult{wvf,rep,Neuron}.pixelsize;
%         perisomaInd=find(distFromSoma'>perisomadist(1) & distFromSoma'<perisomadist(2));
%
%         SpikeMat=permute(reshape(normTr(:,SpikeinBlue_ind'+nTauPeriSp),nROI,[],length(nTauPeriSp)),[1 3 2]); %1: ROI, 2:time, 3:event
%         SpikeAmp{g}=[[NaN SpikeinBlue_ind-BlueOnset]; [distFromSoma' squeeze(max(SpikeMat,[],2))]];
%         SpikeAmp{g}(2:end,2:end)=SpikeAmp{g}(2:end,2:end)./mean(SpikeAmp{g}(perisomaInd+1,2:end),[1]);
%         SpikeAmp{g}(2:end,:)=SpikeAmp{g}(dOrder+1,:);
%
%         SPdelay=round(abs(distFromSoma)/conductionspeed);
%         sub=[];
%         for r=1:size(normTr,1)
%             sp_vec=ind2vec(size(normTr,2),find(som_spike)+SPdelay(r),1,0);
%         [~, sub(r,:)]=get_subthreshold(normTr(r,:),sp_vec,2*(2+SPdelay(r))+1,15);
%         end
%         normTr_sub=normTr-sub;
%
%         SpikeMat_sub=permute(reshape(normTr_sub(:,SpikeinBlue_ind'+nTauPeriSp),nROI,[],length(nTauPeriSp)),[1 3 2]); %1: ROI, 2:time, 3:event
%         AUC=squeeze(sum(SpikeMat_sub,2,'omitnan'));
%         %
%         % for r=1:size(SpikeMat,1)
%         %     for s=1:size(SpikeMat,3)
%         %     AUC(r,s)=get_AUC(squeeze(SpikeMat(r,:,s)),-nTauPeriSp(1)+1+SPdelay(r),3,3);
%         %     end
%         % end
%
%         SpikeAUC{g}=[[NaN SpikeinBlue_ind-BlueOnset]; [distFromSoma' AUC]];
%         SpikeAUC{g}(2:end,2:end)=SpikeAUC{g}(2:end,2:end)./mean(SpikeAUC{g}(perisomaInd+1,2:end),[1]);
%         SpikeAUC{g}(2:end,:)=SpikeAUC{g}(dOrder+1,:);
%
%         g=g+1;
%     end
% end
%
% figure(13); clf; tiledlayout(1,2);
% nexttile([1 1]);
% [binnedZ binX binY]=show3Dbinning(SpikeAmp, 9, 6, 'image'); hold all
% allX=[]; allY=[]; allZ=[];
%   for n = 1:numel(SpikeAmp)
%         timeAxis = SpikeAmp{n}(1, 2:end);
%         xAxis = SpikeAmp{n}(2:end, 1);
%         zData = SpikeAmp{n}(2:end, 2:end);
%         [X, Y] = meshgrid(timeAxis, xAxis);
%         allX = [allX; X(:)];
%         allY = [allY; Y(:)];
%         allZ = [allZ; zData(:)];
%     end
% scatter3(allX,allY,allZ,20,[0.7 0.7 0.7],'filled')
% set(gca, 'YDir', 'reverse');
% colormap("turbo")
% %zlim([0.3 1.6])
% xlabel('Time after blue onset (ms)')
% ylabel('Distance from soma (\mum)')
% title('Normalized Spike Amplitude')
% %view(0,90)
%
% nexttile([1 1]);
% [binnedZ binX binY]=show3Dbinning(SpikeAUC, 9, 7, 'image'); hold all
% allX=[]; allY=[]; allZ=[];
%   for n = 1:numel(SpikeAUC)
%         timeAxis = SpikeAUC{n}(1, 2:end);
%         xAxis = SpikeAUC{n}(2:end, 1);
%         zData = SpikeAUC{n}(2:end, 2:end);
%         [X, Y] = meshgrid(timeAxis, xAxis);
%         allX = [allX; X(:)];
%         allY = [allY; Y(:)];
%         allZ = [allZ; zData(:)];
%     end
% %scatter3(allX,allY,allZ,20,[0.8 0.8 0.8],'filled')
% set(gca, 'YDir', 'reverse');
% colormap("turbo")
% %zlim([0.3 1.4])
% xlabel('Time after blue onset (ms)')
% ylabel('Distance from soma (\mum)')
% title('Normalized Spike AUC')
% %view(0,90)

% %% Ramp (WF Stim)
% nTau=[-70:70]; nTauPeriSp=[-2:3]; nTauPeriSp2=[-3:3]; conductionspeed=170; %um/ms
% g=ones(1,4); SP_STA=[]; distFromSoma=[];
% dendriteaxis_bin=[-80:80:600];
% perisomadist=[-10 10]; cmap=[0 0 0; 1 0 0];  cmap_light=[0.7 0.7 0.7; 1 0.7 0.7];
% BlueStimN=[10 6 1;12 6 1;13 6 1;14 6 1];
%
% SpikeAmp=[]; SpikeAUC=[]; g=1;
% for n=1:size(BlueStimN,1) % neuron
%     cax=[];
%     Neuron=BlueStimN(n,1);
%     rep=BlueStimN(n,3);
%     BluePulseN=BlueStimN(n,2);
%     wvf=5;
%     if ~isempty(OpResult{wvf,rep,Neuron})
%
%         normTr=OpResult{wvf,rep,Neuron}.normTraces./OpResult{wvf,rep,Neuron}.F0_PCA;
%         %noi=[1:size(OpResult{wvf,rep,Neuron}.normTraces,1)];
%         noi=OpResult{wvf,rep,Neuron}.maintrunkROI;
%         normTr=normTr(noi,:);
%
%         nROI=length(noi);
%         nTime=size(OpResult{wvf,rep,Neuron}.normTraces,2);
%         [~, dOrder]=sort(OpResult{wvf,rep,Neuron}.dist_order(noi),'ascend');
%         som_spike=max(OpResult{wvf,rep,Neuron}.SpClass([1 2],:),[],1);
%         bwBlue=bwlabel(OpResult{wvf,rep,Neuron}.Blue);
%         bwBlue((bwBlue-(max(bwBlue)-5))<0)=0;
%         bwBlue=bwlabel(bwBlue>0);
%         SpikeinBlue=som_spike.*(bwBlue==BluePulseN);
%         SpikeinBlue_ind=find(SpikeinBlue);
%         BlueOnset=find(bwBlue==BluePulseN,1);
%         Rheobase=mean(OpResult{wvf,rep,Neuron}.Blue(SpikeinBlue_ind(1:2)));
%
%         normTr=normTr-prctile(normTr(:,BlueOnset+[-500:-200]),30,2);
%
%         Dsign=ones(1,size(OpResult{wvf,rep,Neuron}.interDendDist,2));
%         Dsign(OpResult{wvf,rep,Neuron}.dist_order(1:find(OpResult{wvf,rep,Neuron}.dist_order==1)-1))=-1;
%         dendaxis=OpResult{wvf,rep,Neuron}.interDendDist(1,:).*Dsign;
%         dendaxis=dendaxis(noi);
%         distFromSoma=dendaxis*OpResult{wvf,rep,Neuron}.pixelsize;
%         perisomaInd=find(distFromSoma'>perisomadist(1) & distFromSoma'<perisomadist(2));
%
%         SpikeMat=permute(reshape(normTr(:,SpikeinBlue_ind'+nTauPeriSp),nROI,[],length(nTauPeriSp)),[1 3 2]); %1: ROI, 2:time, 3:event
%         %SpikeAmp{g}=[[NaN SpikeinBlue_ind-BlueOnset]; [distFromSoma' squeeze(max(SpikeMat,[],2))]];
%         RheoBase_trace=OpResult{wvf,rep,Neuron}.Blue(SpikeinBlue_ind)/Rheobase;
%         SpikeAmp{g}=[[NaN RheoBase_trace]; [distFromSoma' squeeze(max(SpikeMat,[],2))]]; %RheoBase
%         SpikeAmp{g}(2:end,2:end)=SpikeAmp{g}(2:end,2:end)./mean(SpikeAmp{g}(perisomaInd+1,2:end),[1]);
%         %SpikeAmp{g}(2:end,:)=SpikeAmp{g}(dOrder+1,:);
%
%         SPdelay=round(abs(distFromSoma)/conductionspeed);
%         sub=[];
%         for r=1:size(normTr,1)
%             sp_vec=ind2vec(size(normTr,2),find(som_spike)+SPdelay(r),1,0);
%         [~, sub(r,:)]=get_subthreshold(normTr(r,:),sp_vec,2*(2+SPdelay(r))+1,15);
%         end
%         normTr_sub=normTr-sub;
%
%         SpikeMat_sub=permute(reshape(normTr_sub(:,SpikeinBlue_ind'+nTauPeriSp),nROI,[],length(nTauPeriSp)),[1 3 2]); %1: ROI, 2:time, 3:event
%         AUC=squeeze(sum(SpikeMat_sub,2,'omitnan'));
%         %
%         % for r=1:size(SpikeMat,1)
%         %     for s=1:size(SpikeMat,3)
%         %     AUC(r,s)=get_AUC(squeeze(SpikeMat(r,:,s)),-nTauPeriSp(1)+1+SPdelay(r),3,3);
%         %     end
%         % end
%
%         %SpikeAUC{g}=[[NaN SpikeinBlue_ind-BlueOnset]; [distFromSoma' AUC]];
%         SpikeAUC{g}=[[NaN RheoBase_trace]; [distFromSoma' AUC]]; %RheoBase
%         SpikeAUC{g}(2:end,2:end)=SpikeAUC{g}(2:end,2:end)./mean(SpikeAUC{g}(perisomaInd+1,2:end),[1]);
%         %SpikeAUC{g}(2:end,:)=SpikeAUC{g}(dOrder+1,:);
%
%         g=g+1;
%     end
% end
%
% figure(14); clf; tiledlayout(1,2);
% nexttile([1 1]);
% [binnedZ binX binY]=show3Dbinning(SpikeAmp, 11, 8, 'image'); hold all
% allX=[]; allY=[]; allZ=[];
%   for n = 1:numel(SpikeAmp)
%         timeAxis = SpikeAmp{n}(1, 2:end);
%         xAxis = SpikeAmp{n}(2:end, 1);
%         zData = SpikeAmp{n}(2:end, 2:end);
%         [X, Y] = meshgrid(timeAxis, xAxis);
%         allX = [allX; X(:)];
%         allY = [allY; Y(:)];
%         allZ = [allZ; zData(:)];
%     end
% %scatter3(allX,allY,allZ,20,[0.7 0.7 0.7],'filled')
% set(gca, 'YDir', 'reverse');
% colormap("turbo")
% %zlim([0.3 1.6])
% xlabel('Optical Rheobase')
% ylabel('Distance from soma (\mum)')
% title('Normalized Spike Amplitude')
% %view(0,90)
%
% nexttile([1 1]);
% [binnedZ binX binY]=show3Dbinning(SpikeAUC, 11, 8, 'image'); hold all
% allX=[]; allY=[]; allZ=[];
%   for n = 1:numel(SpikeAUC)
%         timeAxis = SpikeAUC{n}(1, 2:end);
%         xAxis = SpikeAUC{n}(2:end, 1);
%         zData = SpikeAUC{n}(2:end, 2:end);
%         [X, Y] = meshgrid(timeAxis, xAxis);
%         allX = [allX; X(:)];
%         allY = [allY; Y(:)];
%         allZ = [allZ; zData(:)];
%     end
% %scatter3(allX,allY,allZ,20,[0.8 0.8 0.8],'filled')
% set(gca, 'YDir', 'reverse');
% colormap("turbo")
% %zlim([0.3 1.4])
% xlabel('Optical Rheobase')
% ylabel('Distance from soma (\mum)')
% title('Normalized Spike AUC')
% %view(0,90)

%% Dendrite targeted stimulation



