
clear
clc;
cd '/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/Statistics_Optopatch_Prism';
[~, ~, raw] = xlsread(['/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/' ...
    'Prism_OptopatchData_Arrangement.xlsx'], 'Sheet1', 'B5:O175');

fpath=raw(:,1);
Mouse=cell2mat(raw(:,2));
NeuronInd=cell2mat(raw(:,5));
CamType=raw(:,3);
StructureData=raw(:,10);
StimROI=raw(:,6);
StimWfn=raw(:,7);
PixelSize=cell2mat(cellfun(@(x) (str2num(num2str(x))),raw(:,12),'UniformOutput',false));
refROI=cellfun(@(x) (str2num(num2str(x))),raw(:,14),'UniformOutput',false);

place_bin=150; time_segment=15000; overlap=200;
alignedMovFN = {'STA_Mat_SS','STA_Mat_CS','STA_Mat_dSP'};
bound=6;
title_str={'Basal','Apical','Peri-Soma'};
set(0,'DefaultFigureWindowStyle','docked')


%%
% for f=1:5
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

% %%
[~, unqInd] = unique([Mouse NeuronInd] ,'row');
for i=unqInd(32)'
i
    cat_trace=[];
    cat_spike=[];
    SameCellInd=find(Mouse==Mouse(i) & NeuronInd==NeuronInd(i));
    for j=SameCellInd'
load(fullfile(fpath{j},'OP_Result.mat'),'Result');
    cat_trace=[cat_trace Result.normTraces];
    cat_spike=[cat_spike Result.spike];
    end

    cat_trace=cat_trace-movmedian(cat_trace,300,2);
    [V D eigTrace]=get_eigvector(cat_trace',10);
    ind2use=find(cumsum(D)/sum(D)>0.90,1);
    F0_PCA=sqrt(sum((V(:,[1:ind2use]).*sqrt(D([1:ind2use]))').^2,2));

    for j=SameCellInd'
        load(fullfile(fpath{j},'OP_Result.mat'),'Result');
        Result.F0_PCA=F0_PCA;
        save(fullfile(fpath{j},'OP_Result.mat'),'Result','-v7.3');
    end
end
%% Concatenate data
StimROI_Ind={'Soma','Distal Dend'};
StimWfn_Ind={'Ramp Stim','Short Pulse'};
CatResult=[];
foi_somDend=[1 22 14 18 10 26 25 32 46];
g2=1;
[~, unqInd] = unique([Mouse NeuronInd] ,'row');
for i=unqInd(foi_somDend)'

    SameCellInd=find(Mouse==Mouse(i) & NeuronInd==NeuronInd(i));
    isSoma = ~cellfun(@isempty,strfind(StimROI(SameCellInd), StimROI_Ind{1}));
    isDend = ~cellfun(@isempty,strfind(StimROI(SameCellInd), StimROI_Ind{2}));
    isRamp = ~cellfun(@isempty,strfind(StimWfn(SameCellInd), StimWfn_Ind{1}));
    isSP   = ~cellfun(@isempty,strfind(StimWfn(SameCellInd), StimWfn_Ind{2}));

    ROIwvf_ind=[isSoma isDend isRamp isSP];
    validind=find(sum(ROIwvf_ind,2)>=2);
    patterns = [1 0 1 0; 1 0 0 1; 0 1 1 0; 0 1 0 1];
    values = [1; 2; 3; 4]; % 1=soma, ramp; 2=soma, sp; 3=dend, ramp; 4=dend, sp;

    g=ones(1,4);
    for j=1:length(validind)
        f2read=SameCellInd(validind(j));
        load(fullfile(fpath{f2read},'OP_Result.mat'),'Result');
        wfn = values(find(ismember(patterns, ROIwvf_ind(validind(j),:), 'rows')));
        CatResult{wfn,g(wfn),g2}=Result;
        CatResult{wfn,g(wfn),g2}.pixelsize=PixelSize(f2read);
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
tiledlayout(size(CatResult,3),2)
for n=1:size(CatResult,3)
    cax=[];
    for wvf=[2 4]
        for rep=1:size(CatResult,2)            
            if ~isempty(CatResult{wvf,rep,n})
            normTr=CatResult{wvf,rep,n}.normTraces./CatResult{wvf,rep,n}.F0_PCA;
            nROI=size(CatResult{wvf,rep,n}.normTraces,1);
            nTime=size(CatResult{wvf,rep,n}.normTraces,2);
            dOrder=CatResult{wvf,rep,n}.dist_order;
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
            

            SP_STA{wvf,g(wvf),1}=reshape(normTr(:,som_1st_spike'+nTau),nROI,[],length(nTau));
            SP_STA{wvf,g(wvf),1}=squeeze(mean(SP_STA{wvf,g(wvf),1},2,'omitnan'));
            SP_STA{wvf,g(wvf),1}=SP_STA{wvf,g(wvf),1}-mean(SP_STA{wvf,g(wvf),1}(:,1:10),2);

            if ~isempty(som_2ndst_spike)
            sp_na_2=sum((som_2ndst_spike'+nTau)<0 | (som_2ndst_spike'+nTau)>nTime,2)==0;
            som_2ndst_spike=som_2ndst_spike(sp_na_2);
            SP_STA{wvf,g(wvf),2}=reshape(normTr(:,som_2ndst_spike'+nTau),nROI,[],length(nTau));
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

            nexttile(2*n-2+wvf/2,[1 1])
            imagesc(SP_STA{wvf,g(wvf),1}(dOrder,:),cax)
            title([num2str(n) '#,' num2str(wvf) ',' num2str(rep)])
            g(wvf)=g(wvf)+1;
            end
        end
    end
end
colormap(turbo);

figure(11); clf; l=[]; ls=[];
for wvf=[2 4]
    frstAmp_kink=[]; frstAmp=[];
    for n=1:sum(cellfun(@isempty,distFromSoma(wvf,:))==0,2)
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

figure(12); clf; wvf=4; l=[]; lk=[];
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


%%
nTau=[-70:70]; nTauPeriSp=[-3:7]; nTauPeriSp2=[-3:3]; conductionspeed=170; %um/ms
g=ones(1,4); SP_STA=[]; distFromSoma=[];
dendriteaxis_bin=[-280:80:600];
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
        nROI=size(CatResult{wvf,rep,Neuron}.normTraces,1);
        nTime=size(CatResult{wvf,rep,Neuron}.normTraces,2);
        dOrder=CatResult{wvf,rep,Neuron}.dist_order;
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
[binnedZ binX binY]=show3Dbinning(SpikeAmp, 8, 6, 'surf'); hold all
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
zlim([0.3 1.6])
xlabel('Time after blue onset (ms)')
ylabel('Distance from soma (\mum)')
zlabel('Normalized Spike Amplitude')
view(36,38)

nexttile([1 1]);
[binnedZ binX binY]=show3Dbinning(SpikeAUC, 8, 6, 'surf'); hold all
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
zlim([0.3 1.4])
xlabel('Time after blue onset (ms)')
ylabel('Distance from soma (\mum)')
zlabel('Normalized Spike AUC')
view(36,38)

%%
i=89; bound=6;
load(fullfile(fpath{i},'OP_Result.mat'),'Result');
load([fpath{i} '/mcTrace' num2str(1,'%02d') '.mat']);
mov_mc=readBinMov_BHL(fpath{i});
mean_F=squeeze(mean(mov_res(bound:end-bound,bound:end-bound,:),[1 2]));

mov_res= mov_mc-mean(mov_mc,3);
bkg = zeros(1, size(mov_mc,3));

[~, blueOff]=get_blueoffTrace(mean_F,[Result.Blue],70);
[y_fit]=expfitDM_2(find(blueOff)',mean_F(find(blueOff)),[1:size(mov_mc,3)]',1000);
bkg(1,:)=y_fit;
mov_res = SeeResiduals(mov_res,Result.mc);
mov_res = SeeResiduals(mov_res,Result.mc.^2);
mov_res = SeeResiduals(mov_res,Result.mc(:,1).*Result.mc(:,end));
mov_res= SeeResiduals(mov_res,bkg,1);
CS_s=find(Result.SpClass(2,:)==1);
CS_STAmov=reshape(mov_res(:,:,CS_s'+[-100:100]),size(mov_res,1),size(mov_res,2),[],201);
CS_STAmov=squeeze(mean(CS_STAmov,3));
CS_STAmov=-(CS_STAmov-mean(CS_STAmov(:,:,1:30),3));
CS_STAmov=CS_STAmov.*double(max(Result.bvMask,[],3)==0);
STAmovVec=tovec(CS_STAmov);

%%

% CS_avg=Result.normTraces(:,CS_s'+[-100:100]);
% CS_avg=squeeze(mean(reshape(Result.normTraces(:,CS_s'+[-100:100]),size(Result.ftprnt,3),[],201),2));
% CS_avg=CS_avg-mean(CS_avg(:,1:30),2);

coord_1d=dim_reduce(get_coord(CatResult{2}.ftprnt));
coord_1d=coord_1d-coord_1d(1);
[~, dsort]=sort(coord_1d,'ascend');

% nROI=size(Result.ftprnt,3);
% maxfrm=[]; F0slope=[]; Fslope=[]; Fslope_weight=[];
% F0=tovec(Result.ref_im)-100;
% F0_filter=tovec(imgaussfilt(Result.ref_im,2))-100;
%F_ref=mean(STAmovVec(:,100+[20:30]),2,'omitnan')'*tovec(Result.ftprnt);
%F_ref=mean(CS_avg(:,100+[20:30]),2,'omitnan')';


figure(13); clf;
for n=1:nROI
    
px=find(tovec(Result.ftprnt(:,:,n)>0) & tovec(max(Result.bvMask,[],3)==0));
[~, maxfrm(n)]=max(mean(STAmovVec(px,1:105),1,'omitnan'));

dF=STAmovVec(:,maxfrm(n));
%dF=max(STAmovVec,[],2);
F0_weight=F0.*rescale(tovec(Result.ftprnt(:,:,n)));
dF_weight=dF.*rescale(tovec(Result.ftprnt(:,:,n)));

px2=px(find(F0(px)>prctile(F0(px),30) & F0(px)<prctile(F0(px),98)));
px2_weight=px(find(F0_weight(px)>prctile(F0_weight(px),30) & F0_weight(px)<prctile(F0_weight(px),98)));

nexttile([1 1])
plot(F0(px),dF(px),'.'); hold all
%plot(F0_weight(px),dF_weight(px),'.'); hold all
[p]=polyfit(F0(px2), dF(px2), 1); hold all
[p_weight]=polyfit(F0_weight(px2_weight), dF_weight(px2_weight), 1); hold all
% Get fitted values
y_fit = polyval(p, F0(px2));
plot(F0(px2),y_fit,'r')
SS_res = sum((dF(px2) - y_fit).^2);
SS_tot = sum((dF(px2) - mean(dF(px2))).^2);
R_squared = 1 - (SS_res / SS_tot);

y_fit_weight = polyval(p, F0_weight(px2_weight));
SS_res = sum((dF_weight(px2_weight) - y_fit_weight).^2);
SS_tot = sum((dF_weight(px2_weight) - mean(dF(px2_weight))).^2);
R_squared_weight = 1 - (SS_res / SS_tot);

title(['d(\DeltaF)/dF0 : ' num2str(p(1),2) ', R^2 : ' num2str(R_squared,2)])
Fslope(n)=p(1);
Fslope_weight(n,1)=p_weight(1);
Rsq(n)=R_squared;
F0slope(n,1)=CS_avg(n,maxfrm(n))/Fslope(n);
F0slope_weight(n,1)=CS_avg(n,maxfrm(n))/Fslope_weight(n);
end

figure(14); clf;
nexttile([1 1])
l=plot(CS_avg(dsort,:)'./F0slope_weight(dsort)');
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(turbo(nROI),2))
xlabel('Time (ms)'); ylabel('\DeltaF/F (robust)');

nexttile([1 1])
imagesc(CS_avg(dsort,:)./F0slope_weight(dsort));
colormap(turbo);
xlabel('Time (ms)'); ylabel('Basal to apical ROI');

nexttile([1 1])
Dsign=ones(1,size(Result.interDendDist,2));
Dsign(dsort(1:find(dsort==1)-1))=-1;
dendaxis=Result.interDendDist(1,:).*Dsign;

plot(dendaxis(dsort),max(CS_avg(dsort,100:104)./F0slope_weight(dsort),[],2),'.-'); hold all
plot(dendaxis([19 1:8]),max(CS_avg([19 1:8],100:104)./F0slope_weight([19 1:8]),[],2),'ro')
xlabel('Distance from soma (\mum)'); ylabel('\DeltaF/F (robust)');

nexttile([1 1])
l=plot(CS_avg(dsort,:)'./F_ref(dsort));
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(turbo(nROI),2))
xlabel('Time (ms)'); ylabel('Spike Amplitude (\DeltaF/F_{ref})');

nexttile([1 1])
imagesc(CS_avg(dsort,:)./F_ref(dsort)');
colormap(turbo);
xlabel('Time (ms)'); ylabel('Basal to apical ROI');

nexttile([1 1])
Dsign=ones(1,size(Result.interDendDist,2));
Dsign(dsort(1:find(dsort==1)-1))=-1;
dendaxis=Result.interDendDist(1,:).*Dsign;

plot(dendaxis(dsort),max(CS_avg(dsort,100:104)./F_ref(dsort)',[],2),'.-'); hold all
plot(dendaxis([19 1:8]),max(CS_avg([19 1:8],100:104)./F_ref([19 1:8])',[],2),'ro')
xlabel('Distance from soma (\mum)'); ylabel('Spike Amplitude (\DeltaF/F_{ref})');
%%

figure(15); clf; tiledlayout(5,1); ax1=[];
ax1=[ax1 nexttile([1 1])];
NormTrace=CatResult{4,1}.normTraces./F_ref';
NormTrace=NormTrace-prctile(NormTrace,10,2);
plot(NormTrace([1],:)')
title('Soma Stimulation')
ylabel('\DeltaF/F_{ref}')

ax1=[ax1 nexttile([1 1])];
imagesc(NormTrace(dsort,:),[-0.5 2.5]);
if find(dsort==1)~=1
    set(gca,'YTick',[1 find(dsort==1) nROI],'YTickLabel',num2str([min(interDendDist(1,:).*Dsign) 0 max((interDendDist(1,:).*Dsign))]',3))
else
    set(gca,'YTick',[1 nROI],'YTickLabel',num2str([0 max((interDendDist(1,:)))]',3))
end

ax1=[ax1 nexttile([1 1])];
NormTrace=CatResult{3,2}.normTraces./F_ref';
NormTrace=NormTrace-prctile(NormTrace,10,2);
plot(NormTrace([1],:)')
title('Distal Dendrite Stimulation')
ylabel('\DeltaF/F_{ref}')

ax1=[ax1 nexttile([1 1])];
imagesc(NormTrace(dsort,:),[-0.5 2.5]); colormap(turbo);
if find(dsort==1)~=1
    set(gca,'YTick',[1 find(dsort==1) nROI],'YTickLabel',num2str([min(interDendDist(1,:).*Dsign) 0 max((interDendDist(1,:).*Dsign))]',3))
else
    set(gca,'YTick',[1 nROI],'YTickLabel',num2str([0 max((interDendDist(1,:)))]',3))
end

ax1=[ax1 nexttile([1 1])];
plot(CatResult{4,1}.Blue)
xlabel('Time (ms)')
ylabel('Time (ms)')
linkaxes(ax1,'x')