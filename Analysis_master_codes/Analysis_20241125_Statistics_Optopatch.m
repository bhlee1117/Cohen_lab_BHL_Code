
clear
clc;
cd '/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/Statistics_Optopatch_Prism';
[~, ~, raw] = xlsread(['/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/' ...
    'Prism_OptopatchData_Arrangement.xlsx'], 'Sheet1', 'B5:P175');

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


%%
% for f=[168]
% load(fullfile(fpath{f},'OP_Result.mat'),'Result');
%     interDendDist=[];
%     SkelDend = Skeletonize_dendrite(Result.ref_im,10,0.02,10);
%     nROI=size(Result.ftprnt,3);
%     for i=1%:nROI
%         i
%         for j=1:nROI
%             [interDendDist(i,j), path]=geodesic_distance(SkelDend,get_coord(Result.ftprnt(:,:,i)),get_coord(Result.ftprnt(:,:,j)));
%         end
%     end
%     Result.interDendDist=interDendDist;
%     save(fullfile(fpath{f},'OP_Result.mat'),'Result','-v7.3');
% end
% 
% %%
% [~, unqInd] = unique([Mouse NeuronInd] ,'row');
% for i=unqInd([37])'
% i
%     cat_trace=[];
%     cat_spike=[];
%     SameCellInd=find(Mouse==Mouse(i) & NeuronInd==NeuronInd(i) & isGoodCell==1);
%     for j=SameCellInd'
% load(fullfile(fpath{j},'OP_Result.mat'),'Result');
%     cat_trace=[cat_trace Result.normTraces];
%     cat_spike=[cat_spike Result.spike];
%     end
% 
%     cat_trace=cat_trace-movmedian(cat_trace,300,2);
%     [V D eigTrace]=get_eigvector(cat_trace',20);
%     ind2use=find(cumsum(D)/sum(D)>0.80,1)
%     F0_PCA=sqrt(sum((V(:,[1:ind2use]).*sqrt(D([1:ind2use]))').^2,2));
% 
%     for j=SameCellInd'
%         j
%         load(fullfile(fpath{j},'OP_Result.mat'),'Result');
%         Result.F0_PCA=F0_PCA;
%         save(fullfile(fpath{j},'OP_Result.mat'),'Result','-v7.3');
%     end
% end
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
        CatResult{wfn,g(wfn),g2}.pixelsize=PixelSize(f2read);
        CatResult{wfn,g(wfn),g2}.maintrunkROI=maintrunkROI{f2read};
        g(wfn)=g(wfn)+1;
    end
    g2=g2+1;
end

%% Short pulses Soma vs Dend
nTau=[-40:40];
g=ones(1,4); SP_STA=[]; distFromSoma=[]; 
dendriteaxis_bin=[-280:80:600];
perisomadist=[-50 50]; cmap=[0 0 0; 1 0 0];  cmap_light=[0.7 0.7 0.7; 1 0.7 0.7];

figure(10); clf;
tiledlayout(9,2)
for n=1:size(CatResult,3)
    cax=[];
    for wvf=[2 4]
        for rep=1:size(CatResult,2)            
            if ~isempty(CatResult{wvf,rep,n})
            normTr=CatResult{wvf,rep,n}.normTraces./CatResult{wvf,rep,n}.F0_PCA;
            %noi=CatResult{wvf,rep,n}.maintrunkROI; % only main trunks
            noi=[1:size(CatResult{wvf,rep,n}.normTraces,1)]; % all ROIs
            
            nROI=length(noi);
            nTime=size(CatResult{wvf,rep,n}.normTraces,2);
            noi_dist=find(ismember(CatResult{wvf,rep,n}.dist_order,noi));
            dOrder=CatResult{wvf,rep,n}.dist_order(noi_dist);
            %[~, dOrder]=sort(CatResult{wvf,rep,n}.dist_order(noi_dist),'ascend');
            som_spike=max(CatResult{wvf,rep,n}.SpClass([1 2],:),[],1);
            som_1st_spike=[]; som_2ndst_spike=[];
            for b=1:max(bwlabel(CatResult{wvf,rep,n}.Blue))
                pulseon=(bwlabel(CatResult{wvf,rep,n}.Blue)==b);
                spind=find((som_spike.*pulseon)>0);
                if length(spind)>0
                som_1st_spike=[som_1st_spike spind(1)];
                end
                if length(spind)>1
                som_2ndst_spike=[som_2ndst_spike spind(2)];
                end
            end
            sp_na_1=sum((som_1st_spike'+nTau)<0 | (som_1st_spike'+nTau)>nTime,2)==0;
            som_1st_spike=som_1st_spike(sp_na_1);
            
            SP_STA{wvf,g(wvf),1}=reshape(normTr(dOrder,som_1st_spike'+nTau),nROI,[],length(nTau));
            SP_STA{wvf,g(wvf),1}=squeeze(mean(SP_STA{wvf,g(wvf),1},2,'omitnan'));
            SP_STA{wvf,g(wvf),1}=SP_STA{wvf,g(wvf),1}-mean(SP_STA{wvf,g(wvf),1}(:,1:10),2);

            if ~isempty(som_2ndst_spike)
            sp_na_2=sum((som_2ndst_spike'+nTau)<0 | (som_2ndst_spike'+nTau)>nTime,2)==0;
            som_2ndst_spike=som_2ndst_spike(sp_na_2);
            SP_STA{wvf,g(wvf),2}=reshape(normTr(dOrder,som_2ndst_spike'+nTau),nROI,[],length(nTau));
            SP_STA{wvf,g(wvf),2}=squeeze(mean(SP_STA{wvf,g(wvf),2},2,'omitnan'));
            SP_STA{wvf,g(wvf),2}=SP_STA{wvf,g(wvf),2}-mean(SP_STA{wvf,g(wvf),2}(:,1:10),2);
            end

            if isempty(cax)
                cax=[prctile(tovec(SP_STA{wvf,g(wvf),1}),5) prctile(tovec(SP_STA{wvf,g(wvf),1}),99.9)];
            end

            Dsign=ones(1,size(CatResult{wvf,rep,n}.interDendDist,2));
            Dsign(CatResult{wvf,rep,n}.dist_order(1:find(CatResult{wvf,rep,n}.dist_order==1)-1))=-1;
            dendaxis=CatResult{wvf,rep,n}.interDendDist(1,:).*Dsign;
            distFromSoma{wvf,g(wvf)}=dendaxis*CatResult{wvf,rep,n}.pixelsize;
            distFromSoma{wvf,g(wvf)}=distFromSoma{wvf,g(wvf)}(dOrder);

            nexttile(2*n-2+wvf/2,[1 1])
            imagesc(SP_STA{wvf,g(wvf),1}(:,:),cax)
            title([num2str(n) '#,' num2str(wvf) ',' num2str(rep)])

            if find(dOrder==1)~=1
            set(gca,'YTick',[1 find(dOrder==1) length(dOrder)],'YTickLabel',num2str([min(distFromSoma{wvf,g(wvf)}) 0 max((distFromSoma{wvf,g(wvf)}))]',3))
        else
            set(gca,'YTick',[1 length(dOrder)],'YTickLabel',num2str([0 max((distFromSoma{wvf,g(wvf)}))]',3))
            end

            g(wvf)=g(wvf)+1;
            end
        end
    end
end
colormap(turbo);

figure(11); clf; l=[]; ls=[];
for wvf=[2 4]
    frstAmp_kink=[]; frstAmp=[];
    for n=1:sum(cellfun(@isempty,distFromSoma(wvf,:))==0,2) %neuron
        perisomaInd=find(distFromSoma{wvf,n}'>perisomadist(1) & distFromSoma{wvf,n}'<perisomadist(2));
        frstAmp{n}=[distFromSoma{wvf,n}' max(SP_STA{wvf,n,1}(:,-nTau(1)+[0:3]),[],2)];
        frstAmp{n}(:,2)=frstAmp{n}(:,2)./mean(frstAmp{n}(perisomaInd,2));
        %frstAmp_kink{n}=[distFromSoma{wvf,n}' max(SP_STA{wvf,n}(:,-nTau(1)+[0:3]),[],2)-mean(SP_STA{wvf,n}(:,-nTau(1)+[-4:-1]),2,'omitnan')];
        [~, maxfrm]=max(SP_STA{wvf,n,1}(:,1:-nTau(1)+4),[],2);
        AUC_kink=[]; snd_dip=[];
        for r=1:size(SP_STA{wvf,n,1},1)
            [~, snd_dip(r)]=min(SP_STA{wvf,n,1}(r,maxfrm(r)+[1:5]),[],2);
            [~, pre_dip(r)]=min(SP_STA{wvf,n,1}(r,maxfrm(r)+[-3:0]),[],2);
            postAmp=SP_STA{wvf,n,1}(r,maxfrm(r)+snd_dip(r)); preAmp=SP_STA{wvf,n,1}(r,maxfrm(r)+pre_dip(r)-4);
            if postAmp>preAmp
                AUC_kink(r,1)=sum(SP_STA{wvf,n,1}(r,maxfrm(r)+pre_dip(r)-4:maxfrm(r)+snd_dip(r)),2)-(postAmp-preAmp)*(snd_dip(r)-pre_dip(r)+5)/2-preAmp*(snd_dip(r)-pre_dip(r)+5);
            else
                subTr=SP_STA{wvf,n,1}(r,maxfrm(r)+pre_dip(r)-4:maxfrm(r)+snd_dip(r))-preAmp;
                AUC_kink(r,1)=sum(subTr(subTr>0));
            end
        end
        frstAmp_kink{n}=[distFromSoma{wvf,n}' AUC_kink];
        frstAmp_kink{n}(:,2)=frstAmp_kink{n}(:,2)./mean(frstAmp_kink{n}(perisomaInd,2));
    end
    nexttile(1,[1 1])
    [mean_amplitudes std_amplitudes dendBin_center ind]=binning_data(frstAmp,dendriteaxis_bin);
    N_neuron=sum((cellfun(@sum,ind)),2)';
    for n=1:length(frstAmp)
        scatter(frstAmp{n}(:,1),frstAmp{n}(:,2),20,cmap_light(wvf/2,:),'filled'); hold all
    end
    l(wvf/2)=errorbar(dendBin_center, mean_amplitudes, std_amplitudes./sqrt(N_neuron), 'o-', 'CapSize', 5,'color',cmap(wvf/2,:),'LineWidth',2); hold all
    xlabel('Distance from Soma (\mum)')
    ylabel('Normalized bAP amplitude')


    nexttile(2,[1 1])
    [mean_amplitudes std_amplitudes dendBin_center ind]=binning_data(frstAmp_kink,dendriteaxis_bin);
    N_neuron=sum((cellfun(@sum,ind)),2)';
    for n=1:length(frstAmp)
        scatter(frstAmp_kink{n}(:,1),frstAmp_kink{n}(:,2),20,cmap_light(wvf/2,:),'filled'); hold all
    end
    ls(wvf/2)=errorbar(dendBin_center, mean_amplitudes, std_amplitudes./sqrt(N_neuron), 'o-', 'CapSize', 5,'color',cmap(wvf/2,:),'LineWidth',2); hold all
    xlabel('Distance from Soma (\mum)')
    ylabel('Normalized bAP AUC')

end
legend(l,{'Soma Stim.','Distal dend Stim.'})
legend(ls,{'Soma Stim.','Distal dend Stim.'})

figure(12); clf; wvf=2; l=[]; lk=[];
for s=[1:2]
    frstAmp_kink=[]; frstAmp=[];
    for n=1:sum(cellfun(@isempty,distFromSoma(wvf,:))==0,2)
        perisomaInd=find(distFromSoma{wvf,n}'>perisomadist(1) & distFromSoma{wvf,n}'<perisomadist(2));
        if ~isempty(SP_STA{wvf,n,s})
        frstAmp{n}=[distFromSoma{wvf,n}' max(SP_STA{wvf,n,s}(:,-nTau(1)+[0:3]),[],2)];
        frstAmp{n}(:,2)=frstAmp{n}(:,2)./mean(frstAmp{n}(perisomaInd,2));
        [~, maxfrm]=max(SP_STA{wvf,n,s}(:,1:-nTau(1)+4),[],2);        
        AUC_kink=[]; snd_dip=[]; pre_dip=[];
        for r=1:size(SP_STA{wvf,n,s},1)
            [~, snd_dip(r)]=min(SP_STA{wvf,n,s}(r,maxfrm(r)+[1:7]),[],2);
            [~, pre_dip(r)]=min(SP_STA{wvf,n,s}(r,maxfrm(r)+[-4:0]),[],2);
            postAmp=SP_STA{wvf,n,s}(r,maxfrm(r)+snd_dip(r)); preAmp=SP_STA{wvf,n,s}(r,maxfrm(r)+pre_dip(r)-5);
            if postAmp>preAmp
        AUC_kink(r,1)=sum(SP_STA{wvf,n,s}(r,maxfrm(r)+pre_dip(r)-5:maxfrm(r)+snd_dip(r)),2)-(postAmp-preAmp)*(snd_dip(r)-pre_dip(r)+6)/2-preAmp*(snd_dip(r)-pre_dip(r)+6);
            else
                subTr=SP_STA{wvf,n,s}(r,maxfrm(r)+pre_dip(r)-5:maxfrm(r)+snd_dip(r))-preAmp;
        AUC_kink(r,1)=sum(subTr(subTr>0));
            end
        end
        frstAmp_kink{n}=[distFromSoma{wvf,n}' AUC_kink];
        frstAmp_kink{n}(:,2)=frstAmp_kink{n}(:,2)./mean(frstAmp_kink{n}(perisomaInd,2));
        else
        frstAmp{n}=[distFromSoma{wvf,n}' NaN(length(distFromSoma{wvf,n}),1)];
        frstAmp_kink{n}=[distFromSoma{wvf,n}' NaN(length(distFromSoma{wvf,n}),1)];
        end
    end
nexttile(1,[1 1])
[mean_amplitudes std_amplitudes dendBin_center ind]=binning_data(frstAmp,dendriteaxis_bin);
N_neuron=sum((cellfun(@sum,ind)),2)';
for n=1:length(frstAmp) 
    scatter(frstAmp{n}(:,1),frstAmp{n}(:,2),20,cmap_light(s,:),'filled'); hold all
end
l(s)=errorbar(dendBin_center, mean_amplitudes, std_amplitudes./sqrt(N_neuron), 'o-', 'CapSize', 5,'color',cmap(s,:),'LineWidth',2); hold all
xlabel('Distance from Soma (\mum)')
ylabel('Normalized bAP amplitude')

nexttile(2,[1 1])
for n=1:length(frstAmp_kink) 
    scatter(frstAmp_kink{n}(:,1),frstAmp_kink{n}(:,2),20,cmap_light(s,:),'filled'); hold all
end
[mean_amplitudes std_amplitudes dendBin_center ind]=binning_data(frstAmp_kink,dendriteaxis_bin);
N_neuron=sum((cellfun(@sum,ind)),2)';
lk(s)=errorbar(dendBin_center, mean_amplitudes, std_amplitudes./sqrt(N_neuron), 'o-', 'CapSize', 5,'color',cmap(s,:),'LineWidth',2); hold all
xlabel('Distance from Soma (\mum)')
ylabel('Normalized bAP AUC')
end
legend(l,{'1st spike','2nd spike'})
legend(lk,{'1st spike','2nd spike'})


%% Pulse (Soma Stim)
nTau=[-70:70]; nTauPeriSp=[-3:7]; nTauPeriSp2=[-3:3]; conductionspeed=170; %um/ms
g=ones(1,4); SP_STA=[]; distFromSoma=[];
dendriteaxis_bin=[-80:80:600];
perisomadist=[-50 50]; cmap=[0 0 0; 1 0 0];  cmap_light=[0.7 0.7 0.7; 1 0.7 0.7];
BlueStimN=[2 5 1;3 3 1;3 4 1;4 2 1;4 3 1;4 4 1;5 3 1;5 5 1;6 4 1;6 5 1;6 3 2;6 4 2;6 5 2;7 1 1;7 2 1;7 3 1;]; % Neuron number, Blue Pulse # , N th session
%BlueStimN=[2 6 1;3 6 1;4 6 1;5 6 1;6 6 1;6 6 2;7 6 1];

SpikeAmp=[]; SpikeAUC=[]; g=1;
for n=1:size(BlueStimN,1) % neuron
    cax=[];
    Neuron=BlueStimN(n,1);
    rep=BlueStimN(n,3);
    BluePulseN=BlueStimN(n,2);
    wvf=1;
    if ~isempty(CatResult{wvf,rep,Neuron})

        normTr=CatResult{wvf,rep,Neuron}.normTraces./CatResult{wvf,rep,Neuron}.F0_PCA;
        %noi=[1:size(CatResult{wvf,rep,Neuron}.normTraces,1)];
        noi=CatResult{wvf,rep,Neuron}.maintrunkROI;
        normTr=normTr(noi,:);

        nROI=length(noi);
        nTime=size(CatResult{wvf,rep,Neuron}.normTraces,2);
        [~, dOrder]=sort(CatResult{wvf,rep,Neuron}.dist_order(noi),'ascend');
        som_spike=max(CatResult{wvf,rep,Neuron}.SpClass([1 2],:),[],1);
        bwBlue=bwlabel(CatResult{wvf,rep,Neuron}.Blue);
        bwBlue((bwBlue-(max(bwBlue)-5))<0)=0;
        bwBlue=bwlabel(bwBlue>0);
        SpikeinBlue=som_spike.*(bwBlue==BluePulseN);
        SpikeinBlue_ind=find(SpikeinBlue);
        BlueOnset=find(bwBlue==BluePulseN,1);
        normTr=normTr-prctile(normTr(:,BlueOnset+[-500:-200]),30,2);

        Dsign=ones(1,size(CatResult{wvf,rep,Neuron}.interDendDist,2));
        Dsign(CatResult{wvf,rep,Neuron}.dist_order(1:find(CatResult{wvf,rep,Neuron}.dist_order==1)-1))=-1;
        dendaxis=CatResult{wvf,rep,Neuron}.interDendDist(1,:).*Dsign;
        dendaxis=dendaxis(noi);
        distFromSoma=dendaxis*CatResult{wvf,rep,Neuron}.pixelsize;
        perisomaInd=find(distFromSoma'>perisomadist(1) & distFromSoma'<perisomadist(2));

        SpikeMat=permute(reshape(normTr(:,SpikeinBlue_ind'+nTauPeriSp),nROI,[],length(nTauPeriSp)),[1 3 2]); %1: ROI, 2:time, 3:event
        SpikeAmp{g}=[[NaN SpikeinBlue_ind-BlueOnset]; [distFromSoma' squeeze(max(SpikeMat,[],2))]];
        SpikeAmp{g}(2:end,2:end)=SpikeAmp{g}(2:end,2:end)./mean(SpikeAmp{g}(perisomaInd+1,2:end),[1]);
        SpikeAmp{g}(2:end,:)=SpikeAmp{g}(dOrder+1,:);

        SPdelay=round(abs(distFromSoma)/conductionspeed);
        sub=[];
        for r=1:size(normTr,1)
            sp_vec=ind2vec(size(normTr,2),find(som_spike)+SPdelay(r),1,0);
        [~, sub(r,:)]=get_subthreshold(normTr(r,:),sp_vec,2*(2+SPdelay(r))+1,15);
        end
        normTr_sub=normTr-sub;

        SpikeMat_sub=permute(reshape(normTr(:,SpikeinBlue_ind'+nTauPeriSp),nROI,[],length(nTauPeriSp)),[1 3 2]); %1: ROI, 2:time, 3:event
        AUC=squeeze(sum(SpikeMat_sub,2,'omitnan'));
        % 
        % for r=1:size(SpikeMat,1)
        %     for s=1:size(SpikeMat,3)
        %     AUC(r,s)=get_AUC(squeeze(SpikeMat(r,:,s)),-nTauPeriSp(1)+1+SPdelay(r),3,3);
        %     end
        % end

        SpikeAUC{g}=[[NaN SpikeinBlue_ind-BlueOnset]; [distFromSoma' AUC]];
        SpikeAUC{g}(2:end,2:end)=SpikeAUC{g}(2:end,2:end)./mean(SpikeAUC{g}(perisomaInd+1,2:end),[1]);
        SpikeAUC{g}(2:end,:)=SpikeAUC{g}(dOrder+1,:);

        g=g+1;
    end
end

figure(13); clf; tiledlayout(1,2);
nexttile([1 1]);
[binnedZ binX binY]=show3Dbinning(SpikeAmp, 11, 9, 'image'); hold all
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
xlabel('Time after blue onset (ms)')
ylabel('Distance from soma (\mum)')
title('Normalized Spike Amplitude')
%view(0,90)

nexttile([1 1]);
[binnedZ binX binY]=show3Dbinning(SpikeAUC, 11, 9, 'image'); hold all
allX=[]; allY=[]; allZ=[];
shading flat
  for n = 1:numel(SpikeAUC)
        timeAxis = SpikeAUC{n}(1, 2:end);
        xAxis = SpikeAUC{n}(2:end, 1);
        zData = SpikeAUC{n}(2:end, 2:end);
        [X, Y] = meshgrid(timeAxis, xAxis);
        allX = [allX; X(:)];
        allY = [allY; Y(:)];
        allZ = [allZ; zData(:)];
    end
%scatter3(allX,allY,allZ,20,[0.8 0.8 0.8],'filled')
set(gca, 'YDir', 'reverse');
colormap("turbo")
%zlim([0.3 1.4])
xlabel('Time after blue onset (ms)')
ylabel('Distance from soma (\mum)')
title('Normalized Spike AUC')
%view(0,90)

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
    if ~isempty(CatResult{wvf,rep,Neuron})

        normTr=CatResult{wvf,rep,Neuron}.normTraces./CatResult{wvf,rep,Neuron}.F0_PCA;
        %noi=[1:size(CatResult{wvf,rep,Neuron}.normTraces,1)];
        noi=CatResult{wvf,rep,Neuron}.maintrunkROI;
        normTr=normTr(noi,:);

        nROI=length(noi);
        nTime=size(CatResult{wvf,rep,Neuron}.normTraces,2);
        [~, dOrder]=sort(CatResult{wvf,rep,Neuron}.dist_order(noi),'ascend');
        som_spike=max(CatResult{wvf,rep,Neuron}.SpClass([1 2],:),[],1);
        bwBlue=bwlabel(CatResult{wvf,rep,Neuron}.Blue);
        bwBlue((bwBlue-(max(bwBlue)-5))<0)=0;
        bwBlue=bwlabel(bwBlue>0);
        SpikeinBlue=som_spike.*(bwBlue==BluePulseN);
        SpikeinBlue_ind=find(SpikeinBlue);
        BlueOnset=find(bwBlue==BluePulseN,1);
        Rheobase=mean(CatResult{wvf,rep,Neuron}.Blue(SpikeinBlue_ind(1:2)));

        
        normTr=normTr-prctile(normTr(:,BlueOnset+[-500:-200]),30,2);

        Dsign=ones(1,size(CatResult{wvf,rep,Neuron}.interDendDist,2));
        Dsign(CatResult{wvf,rep,Neuron}.dist_order(1:find(CatResult{wvf,rep,Neuron}.dist_order==1)-1))=-1;
        dendaxis=CatResult{wvf,rep,Neuron}.interDendDist(1,:).*Dsign;
        dendaxis=dendaxis(noi);
        distFromSoma=dendaxis*CatResult{wvf,rep,Neuron}.pixelsize;
        perisomaInd=find(distFromSoma'>perisomadist(1) & distFromSoma'<perisomadist(2));

        SpikeMat=permute(reshape(normTr(:,SpikeinBlue_ind'+nTauPeriSp),nROI,[],length(nTauPeriSp)),[1 3 2]); %1: ROI, 2:time, 3:event
        %SpikeAmp{g}=[[NaN SpikeinBlue_ind-BlueOnset]; [distFromSoma' squeeze(max(SpikeMat,[],2))]]; 
        RheoBase_trace=CatResult{wvf,rep,Neuron}.Blue(SpikeinBlue_ind)/Rheobase;
        SpikeAmp{g}=[[NaN RheoBase_trace]; [distFromSoma' squeeze(max(SpikeMat,[],2))]]; %RheoBase
        SpikeAmp{g}(2:end,2:end)=SpikeAmp{g}(2:end,2:end)./mean(SpikeAmp{g}(perisomaInd+1,2:end),[1]);
        SpikeAmp{g}(2:end,:)=SpikeAmp{g}(dOrder+1,:);

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

figure(14); clf; tiledlayout(1,2);
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

nexttile([1 1]);
[binnedZ binX binY]=show3Dbinning(SpikeAUC, 11, 8, 'image'); hold all
shading flat
allX=[]; allY=[]; allZ=[];
  for n = 1:numel(SpikeAUC)
        timeAxis = SpikeAUC{n}(1, 2:end);
        xAxis = SpikeAUC{n}(2:end, 1);
        zData = SpikeAUC{n}(2:end, 2:end);
        [X, Y] = meshgrid(timeAxis, xAxis);
        allX = [allX; X(:)];
        allY = [allY; Y(:)];
        allZ = [allZ; zData(:)];
    end
%scatter3(allX,allY,allZ,20,[0.8 0.8 0.8],'filled')
set(gca, 'YDir', 'reverse');
colormap("turbo")
%zlim([0.3 1.4])
xlabel('Optical Rheobase')
ylabel('Distance from soma (\mum)')
title('Normalized Spike AUC')
%view(0,90)

%% Pulse (WF Stim)
nTau=[-70:70]; nTauPeriSp=[-3:7]; nTauPeriSp2=[-3:3]; conductionspeed=170; %um/ms
g=ones(1,4); SP_STA=[]; distFromSoma=[];
dendriteaxis_bin=[-80:80:600];
perisomadist=[-50 50]; cmap=[0 0 0; 1 0 0];  cmap_light=[0.7 0.7 0.7; 1 0.7 0.7];
BlueStimN=[10 3 1;10 4 1;10 5 1;11 4 1;11 5 1;12 3 1;12 4 1;12 5 1;13 4 1;13 5 1; 14 5 1;15 4 1;15 5 1]; % Neuron number, Blue Pulse # , N th session
%BlueStimN=[2 6 1;3 6 1;4 6 1;5 6 1;6 6 1;6 6 2;7 6 1]; 

SpikeAmp=[]; SpikeAUC=[]; g=1;
for n=1:size(BlueStimN,1) % neuron
    cax=[];
    Neuron=BlueStimN(n,1);
    rep=BlueStimN(n,3);
    BluePulseN=BlueStimN(n,2);
    wvf=5;
    if ~isempty(CatResult{wvf,rep,Neuron})

        normTr=CatResult{wvf,rep,Neuron}.normTraces./CatResult{wvf,rep,Neuron}.F0_PCA;
        %noi=[1:size(CatResult{wvf,rep,Neuron}.normTraces,1)];
        noi=CatResult{wvf,rep,Neuron}.maintrunkROI;
        normTr=normTr(noi,:);

        nROI=length(noi);
        nTime=size(CatResult{wvf,rep,Neuron}.normTraces,2);
        [~, dOrder]=sort(CatResult{wvf,rep,Neuron}.dist_order(noi),'ascend');
        som_spike=max(CatResult{wvf,rep,Neuron}.SpClass([1 2],:),[],1);
        bwBlue=bwlabel(CatResult{wvf,rep,Neuron}.Blue);
        bwBlue((bwBlue-(max(bwBlue)-5))<0)=0;
        bwBlue=bwlabel(bwBlue>0);
        SpikeinBlue=som_spike.*(bwBlue==BluePulseN);
        SpikeinBlue_ind=find(SpikeinBlue);
        BlueOnset=find(bwBlue==BluePulseN,1);
        normTr=normTr-prctile(normTr(:,BlueOnset+[-500:-200]),30,2);

        Dsign=ones(1,size(CatResult{wvf,rep,Neuron}.interDendDist,2));
        Dsign(CatResult{wvf,rep,Neuron}.dist_order(1:find(CatResult{wvf,rep,Neuron}.dist_order==1)-1))=-1;
        dendaxis=CatResult{wvf,rep,Neuron}.interDendDist(1,:).*Dsign;
        dendaxis=dendaxis(noi);
        distFromSoma=dendaxis*CatResult{wvf,rep,Neuron}.pixelsize;
        perisomaInd=find(distFromSoma'>perisomadist(1) & distFromSoma'<perisomadist(2));

        SpikeMat=permute(reshape(normTr(:,SpikeinBlue_ind'+nTauPeriSp),nROI,[],length(nTauPeriSp)),[1 3 2]); %1: ROI, 2:time, 3:event
        SpikeAmp{g}=[[NaN SpikeinBlue_ind-BlueOnset]; [distFromSoma' squeeze(max(SpikeMat,[],2))]];
        SpikeAmp{g}(2:end,2:end)=SpikeAmp{g}(2:end,2:end)./mean(SpikeAmp{g}(perisomaInd+1,2:end),[1]);
        SpikeAmp{g}(2:end,:)=SpikeAmp{g}(dOrder+1,:);

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

        SpikeAUC{g}=[[NaN SpikeinBlue_ind-BlueOnset]; [distFromSoma' AUC]];
        SpikeAUC{g}(2:end,2:end)=SpikeAUC{g}(2:end,2:end)./mean(SpikeAUC{g}(perisomaInd+1,2:end),[1]);
        SpikeAUC{g}(2:end,:)=SpikeAUC{g}(dOrder+1,:);

        g=g+1;
    end
end

figure(13); clf; tiledlayout(1,2);
nexttile([1 1]);
[binnedZ binX binY]=show3Dbinning(SpikeAmp, 9, 6, 'image'); hold all
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
scatter3(allX,allY,allZ,20,[0.7 0.7 0.7],'filled')
set(gca, 'YDir', 'reverse');
colormap("turbo")
%zlim([0.3 1.6])
xlabel('Time after blue onset (ms)')
ylabel('Distance from soma (\mum)')
title('Normalized Spike Amplitude')
%view(0,90)

nexttile([1 1]);
[binnedZ binX binY]=show3Dbinning(SpikeAUC, 9, 7, 'image'); hold all
allX=[]; allY=[]; allZ=[];
  for n = 1:numel(SpikeAUC)
        timeAxis = SpikeAUC{n}(1, 2:end);
        xAxis = SpikeAUC{n}(2:end, 1);
        zData = SpikeAUC{n}(2:end, 2:end);
        [X, Y] = meshgrid(timeAxis, xAxis);
        allX = [allX; X(:)];
        allY = [allY; Y(:)];
        allZ = [allZ; zData(:)];
    end
%scatter3(allX,allY,allZ,20,[0.8 0.8 0.8],'filled')
set(gca, 'YDir', 'reverse');
colormap("turbo")
%zlim([0.3 1.4])
xlabel('Time after blue onset (ms)')
ylabel('Distance from soma (\mum)')
title('Normalized Spike AUC')
%view(0,90)

%% Ramp (WF Stim)
nTau=[-70:70]; nTauPeriSp=[-2:3]; nTauPeriSp2=[-3:3]; conductionspeed=170; %um/ms
g=ones(1,4); SP_STA=[]; distFromSoma=[];
dendriteaxis_bin=[-80:80:600];
perisomadist=[-10 10]; cmap=[0 0 0; 1 0 0];  cmap_light=[0.7 0.7 0.7; 1 0.7 0.7];
BlueStimN=[10 6 1;12 6 1;13 6 1;14 6 1];

SpikeAmp=[]; SpikeAUC=[]; g=1;
for n=1:size(BlueStimN,1) % neuron
    cax=[];
    Neuron=BlueStimN(n,1);
    rep=BlueStimN(n,3);
    BluePulseN=BlueStimN(n,2);
    wvf=5;
    if ~isempty(CatResult{wvf,rep,Neuron})

        normTr=CatResult{wvf,rep,Neuron}.normTraces./CatResult{wvf,rep,Neuron}.F0_PCA;
        %noi=[1:size(CatResult{wvf,rep,Neuron}.normTraces,1)];
        noi=CatResult{wvf,rep,Neuron}.maintrunkROI;
        normTr=normTr(noi,:);

        nROI=length(noi);
        nTime=size(CatResult{wvf,rep,Neuron}.normTraces,2);
        [~, dOrder]=sort(CatResult{wvf,rep,Neuron}.dist_order(noi),'ascend');
        som_spike=max(CatResult{wvf,rep,Neuron}.SpClass([1 2],:),[],1);
        bwBlue=bwlabel(CatResult{wvf,rep,Neuron}.Blue);
        bwBlue((bwBlue-(max(bwBlue)-5))<0)=0;
        bwBlue=bwlabel(bwBlue>0);
        SpikeinBlue=som_spike.*(bwBlue==BluePulseN);
        SpikeinBlue_ind=find(SpikeinBlue);
        BlueOnset=find(bwBlue==BluePulseN,1);
        Rheobase=mean(CatResult{wvf,rep,Neuron}.Blue(SpikeinBlue_ind(1:2)));

        normTr=normTr-prctile(normTr(:,BlueOnset+[-500:-200]),30,2);

        Dsign=ones(1,size(CatResult{wvf,rep,Neuron}.interDendDist,2));
        Dsign(CatResult{wvf,rep,Neuron}.dist_order(1:find(CatResult{wvf,rep,Neuron}.dist_order==1)-1))=-1;
        dendaxis=CatResult{wvf,rep,Neuron}.interDendDist(1,:).*Dsign;
        dendaxis=dendaxis(noi);
        distFromSoma=dendaxis*CatResult{wvf,rep,Neuron}.pixelsize;
        perisomaInd=find(distFromSoma'>perisomadist(1) & distFromSoma'<perisomadist(2));

        SpikeMat=permute(reshape(normTr(:,SpikeinBlue_ind'+nTauPeriSp),nROI,[],length(nTauPeriSp)),[1 3 2]); %1: ROI, 2:time, 3:event
        %SpikeAmp{g}=[[NaN SpikeinBlue_ind-BlueOnset]; [distFromSoma' squeeze(max(SpikeMat,[],2))]]; 
        RheoBase_trace=CatResult{wvf,rep,Neuron}.Blue(SpikeinBlue_ind)/Rheobase;
        SpikeAmp{g}=[[NaN RheoBase_trace]; [distFromSoma' squeeze(max(SpikeMat,[],2))]]; %RheoBase
        SpikeAmp{g}(2:end,2:end)=SpikeAmp{g}(2:end,2:end)./mean(SpikeAmp{g}(perisomaInd+1,2:end),[1]);
        %SpikeAmp{g}(2:end,:)=SpikeAmp{g}(dOrder+1,:);

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
        %SpikeAUC{g}(2:end,:)=SpikeAUC{g}(dOrder+1,:);

        g=g+1;
    end
end

figure(14); clf; tiledlayout(1,2);
nexttile([1 1]);
[binnedZ binX binY]=show3Dbinning(SpikeAmp, 11, 8, 'image'); hold all
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

nexttile([1 1]);
[binnedZ binX binY]=show3Dbinning(SpikeAUC, 11, 8, 'image'); hold all
allX=[]; allY=[]; allZ=[];
  for n = 1:numel(SpikeAUC)
        timeAxis = SpikeAUC{n}(1, 2:end);
        xAxis = SpikeAUC{n}(2:end, 1);
        zData = SpikeAUC{n}(2:end, 2:end);
        [X, Y] = meshgrid(timeAxis, xAxis);
        allX = [allX; X(:)];
        allY = [allY; Y(:)];
        allZ = [allZ; zData(:)];
    end
%scatter3(allX,allY,allZ,20,[0.8 0.8 0.8],'filled')
set(gca, 'YDir', 'reverse');
colormap("turbo")
%zlim([0.3 1.4])
xlabel('Optical Rheobase')
ylabel('Distance from soma (\mum)')
title('Normalized Spike AUC')
%view(0,90)

%% Dendrite targeted stimulation

