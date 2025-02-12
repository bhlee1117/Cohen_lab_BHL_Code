clear
clc;
[~, ~, raw] = xlsread(['/Volumes/cohen_lab/Lab/Labmembers' ...
    '/Byung Hun Lee/Data/PrismPCdata_Arrangement.xlsx'], 'Sheet1', 'C5:Z31');

ref_ROI=cellfun(@(x) (str2num(num2str(x))),raw(:,9),'UniformOutput',false);

oblique_ROI=cellfun(@(x) (str2num(num2str(x))),raw(:,10),'UniformOutput',false);
PeriSoma_ROI=cellfun(@(x) (str2num(num2str(x))),raw(:,11),'UniformOutput',false);
basal_ROI=cellfun(@(x) (str2num(num2str(x))),raw(:,12),'UniformOutput',false);
apical_ROI=cellfun(@(x) (str2num(num2str(x))),raw(:,13),'UniformOutput',false);
distal_ROI=cellfun(@(x) (str2num(num2str(x))),raw(:,14),'UniformOutput',false);

fpath=raw(:,1)';
StructureData=raw(:,8);
BadROI=cellfun(@(x) (str2num(num2str(x))),raw(:,17),'UniformOutput',false);
EndFrame=cell2mat(raw(:,15));
ifmotionReject=cell2mat(raw(:,16));
ifdirtRemov=cell2mat(raw(:,18));
Pixelsize=cell2mat(raw(:,6));
save_figto='/Volumes/BHL_WD18TB/PP72_PlaceCellResults';
place_bin=150; time_segment=15000; overlap=200;
alignedMovFN = {'STA_Mat_SS','STA_Mat_CS','STA_Mat_dSP'};
bound=6;
title_str={'Basal','Apical','Peri-Soma'};
PlaceFieldList=cellfun(@(x) (str2num(num2str(x))),raw(:,21),'UniformOutput',false);
PlaceFieldBin=cellfun(@(x) (str2num(num2str(x))),raw(:,22),'UniformOutput',false);
set(0,'DefaultFigureWindowStyle','docked')
%foi=[1 4 5 6 8 10 11 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27];
foi=[1 4 5 6 8 10 11 15 16 17 18 19 20 21 22 23 24 25 26 27];
%foi=20;
%%
nTau={[-200:50],[-200:50],[-200:50]}; %SS, CS, Brst
clear SpikeInd MatSpike MatSTA MatBlue MatCStrace MatSub SpikeList NormalizedTrace_ch NormalizedTrace_dirt SpikeIndBlueOff Dist_order allSpikeMat noi interDendDist noi_dist derivSub LapSubSilent
clear Subthreshold dendaxis BrstOrder roisD roisD_order

for f=foi
    f
    load(fullfile(fpath{f},'PC_Result.mat'),'Result') %load traces
    if ~isempty(find(ismember(Result.spike(1,:),Result.SpClass(1,:))==0)) | ~isempty(find(ismember(Result.spike(1,:),Result.SpClass(2,:))==0))
        error('Spclass and spike mat is not matching');
    end

    Result.spike=Result.spike>0; Result.SpClass=Result.SpClass>0;
    Dist_order{f}=Result.BrancDist_order;
    interDendDist{f}=Result.interDendDist*Pixelsize(f);
    Dsign=ones(1,size(interDendDist{f},1));
    Dsign(Dist_order{f}(1:find(Dist_order{f}==1)-1))=-1;
    perisomaROI=setdiff(find(interDendDist{f}(1,:)<40),BadROI{f}); % ROIs < 40 um from soma
    noi=setdiff([1:size(Result.ftprnt,3)],BadROI{f});
    noi_dist{f}=ismember(Dist_order{f},noi);
    roisD(f,:)={basal_ROI{f},PeriSoma_ROI{f},apical_ROI{f},oblique_ROI{f},distal_ROI{f}}; %set the ROIs
    for dClass=1:size(roisD,2)
        g=1;
        if ~isnan(roisD{f,dClass})
      for d=roisD{f,dClass}
          dind=setdiff(find(Result.BranchLabel==d),BadROI{f});
          roisD_order{f}{dClass,g}=ismember(Dist_order{f}(noi_dist{f}),dind)*d;
          g=g+1;
      end
        end
    end

    dendaxis{f}=interDendDist{f}(1,:).*Dsign;
    dendaxis{f}=dendaxis{f}(Dist_order{f}(noi_dist{f}));

    NormalizedTrace=(Result.normTraces)./Result.F0_PCA;
    bAP_STA=get_STA(NormalizedTrace, Result.spike(1,:), 30, 20);
    bAP_STA=bAP_STA-prctile(bAP_STA,20,2);
    SpikeHeight=max(mean(bAP_STA(perisomaROI,:),1),[],2);

    NormalizedTrace=NormalizedTrace/SpikeHeight;
    NormalizedTrace_dirt{f,1}=NormalizedTrace;
    NormalizedTrace_dirt{f,1}(:,Result.motionReject>0)=NaN;
    NormalizedTrace_ch(f,:)=cellfun(@(x) x./Result.F_ref,Result.norm_trace_check,'UniformOutput',false);
    NormalizedTrace_ch{f,1}(:,Result.motionReject>0)=NaN; NormalizedTrace_ch{f,2}(:,Result.motionReject>0)=NaN;
    if ifdirtRemov(f)
        NormalizedTrace_dirt{f,1}(Result.dirtTrace>0)=NaN;
        NormalizedTrace_ch{f,1}(Result.dirtTrace>0)=NaN; NormalizedTrace_ch{f,2}(Result.dirtTrace>0)=NaN;
    end
    NormalizedTrace_dirt{f,1}=NormalizedTrace_dirt{f,1}(Dist_order{f}(noi_dist{f},:),:);
    NormalizedTrace_ch{f,1}=NormalizedTrace_ch{f,1}(Dist_order{f}(noi_dist{f},:),:);
    NormalizedTrace_ch{f,2}=NormalizedTrace_ch{f,2}(Dist_order{f}(noi_dist{f},:),:);

    Subthreshold{f}=get_subthreshold(NormalizedTrace_dirt{f},max(Result.spike(1,:),[],1)>0,7,17);
    Subthreshold{f}(isnan(NormalizedTrace_dirt{f}))=NaN;

    % spike class, SS:1, CS:2, dSP:3, BS:4
    ss_time=find(Result.SpClass(1,:)); % BS is subclass of SS
    brst=bwlabel((ss_time(2:end)-ss_time(1:end-1))<=20); % SSs that have an ISI shorter than 20 ms are BS.
    SpClass=Result.SpClass; BS_trace=zeros(1,size(Result.traces,2));
    for b=1:max(brst)
        bwn=find(brst==b);
        SpClass(1,ss_time([bwn bwn(end)+1]))=0;
        SpClass(4,ss_time([bwn(1)]))=1;
        BS_trace(1,[ss_time(bwn): ss_time(bwn(end)+1)])=b;
    end
    SpClass=SpClass([1 2 4],:);
    Classvec = get_Class2index(SpClass);
    SpikeClassvec=Classvec.*Result.spike(1,:);

    BrstOrder = get_BurstOrder(Result.spike(1,:), 20) - SpClass(1,:);
    BrstOrder(find(SpClass(3,:)))=1;
    ComplexSpikeOrder=get_spikeOrder(Result.CStrace,Result.spike(1,:));

    for stype=1:3
        s_list=find(SpClass(stype,:)>0);
        [~, MatSpClass{f,stype}]=get_STA(SpikeClassvec,SpClass(stype,:),-nTau{stype}(1),nTau{stype}(end));
        [~, MatSpike{f,stype}]=get_STA(Result.spike(1,:),SpClass(stype,:),-nTau{stype}(1),nTau{stype}(end));
        [~, MatBlue{f,stype}]=get_STA(Result.Blue,SpClass(stype,:),-nTau{stype}(1),nTau{stype}(end));
        [~, MatCStrace{f,stype}]=get_STA(Result.CStrace,SpClass(stype,:),-nTau{stype}(1),nTau{stype}(end));
        [~, MatSub{f,stype}]=get_STA(Subthreshold{f},SpClass(stype,:),-nTau{stype}(1),nTau{stype}(end));
        [~, MatOrder{f,stype}]=get_STA(BrstOrder,SpClass(stype,:),-nTau{stype}(1),nTau{stype}(end));

        switch stype
            case 1
                [~, MatSTA{f,stype,1} sptime]=get_STA(NormalizedTrace_dirt{f,1},SpClass(stype,:),-nTau{stype}(1),nTau{stype}(end));
                SpikeList{f,stype,1}=ind2vec(size(NormalizedTrace_dirt{f},2),sptime,1);
            case 2
                for ns=1:5
                    [~, MatSTA{f,stype,ns} sptime]=get_STA(NormalizedTrace_dirt{f,1},(ComplexSpikeOrder)==ns,-nTau{stype}(1),nTau{stype}(end));
                    SpikeList{f,stype,ns}=ind2vec(size(NormalizedTrace_dirt{f},2),sptime,1);
                end
            case 3
                for ns=1:5
                    [~, MatSTA{f,stype,ns} sptime]=get_STA(NormalizedTrace_dirt{f,1},(BrstOrder.*(1-Result.CStrace))==ns,-nTau{stype}(1),nTau{stype}(end));
                    SpikeList{f,stype,ns}=ind2vec(size(NormalizedTrace_dirt{f},2),sptime,1);
                end
        end
    end

    allSpikeMat{f}=double(Result.spike);
    allSpikeMat{f}(:,Result.motionReject>0)=NaN;
    if ifdirtRemov(f)
        allSpikeMat{f}(Result.dirtTrace>0)=NaN;
    end
    allSpikeClassMat{f}=SpClass;

    BlueStim{f}=Result.Blue;
    VRtrack{f}=Result.VR;
    CStrace{f}=Result.CStrace;
    Ftprnts{f}=Result.ftprnt;
    AvgImg{f}=Result.ref_im;

    % LapFR{f}=PlaceTrigger_average(double(allSpikeMat{f}(1,:)==1),150,VRtrack{f},0.002,115); %total trace
    % LapV{f}=PlaceTrigger_average(NormalizedTrace_dirt{f},150,VRtrack{f},0.002,115); %total trace
    % LapSub{f}=PlaceTrigger_average(Subthreshold{f},150,VRtrack{f},0.002,115); %total trace
    % 
    % subthreshold_silent=Subthreshold{f};
    % subthreshold_silent(:,unique(get_perispikeIndex(allSpikeMat{f}(1,:),[-8:30])))=NaN;
    % LapSubSilent{f}=PlaceTrigger_average(subthreshold_silent,150,VRtrack{f},0.002,115); %total trace
end

MatSTA=cellfun(@(x) permute(x,[1 3 2]),MatSTA,'UniformOutput',false);
MatSpClass=cellfun(@(x) permute(x,[1 3 2]),MatSpClass,'UniformOutput',false);
MatSpike=cellfun(@(x) permute(x,[1 3 2]),MatSpike,'UniformOutput',false);
MatBlue=cellfun(@(x) permute(x,[1 3 2]),MatBlue,'UniformOutput',false);
MatSub=cellfun(@(x) permute(x,[1 3 2]),MatSub,'UniformOutput',false);
MatOrder=cellfun(@(x) permute(x,[1 3 2]),MatOrder,'UniformOutput',false);
MatCStrace=cellfun(@(x) permute(x,[1 3 2]),MatCStrace,'UniformOutput',false);

%% show STA kymo
figure(101); clf; tiledlayout(9,9)
stype_str={'SS','CS','BS'};
sub_time=[-19:-5]; %from spike
ax1=[];
for f=23
    [~, BlueOffSpike]=find_wCondition(SpikeList{f,1,1},BlueStim{f}==0);
    refkymo=mean(MatSTA{f,1,1}(:,:,BlueOffSpike)-median(MatSTA{f,1,1}(:,-nTau{1}(1)+sub_time,BlueOffSpike),2,'omitnan'),3,'omitnan');
    cax=[prctile(refkymo(:),0.5) prctile(refkymo(:),99.5)];

    Dsign=ones(1,size(interDendDist{f},1));
    Dsign(Dist_order{f}(1:find(Dist_order{f}==1)-1))=-1;
    for stype=1:3
        [~, BlueOffSpike]=find_wCondition(SpikeList{f,stype,1},BlueStim{f}==0);
        if ~isempty(BlueOffSpike) % If there is no spike
            imshow=mean(MatSTA{f,stype,1}(:,:,BlueOffSpike)-median(MatSTA{f,stype,1}(:,-nTau{stype}(1)+sub_time,BlueOffSpike),2,'omitnan'),3,'omitnan');
        else
            imshow=zeros(sum(noi_dist{f}),length(nTau{stype}));
        end
        %imshow=mean(STAmat{f,stype}(:,:,SpikeIndBlueOff{f,sty  pe}),3,'omitnan');
        if isnan(cax); cax=[0 1]; end;
        ax1=[ax1 nexttile([1 1])];
        imagesc(imshow,cax)
        title(['Cell# ' num2str(f) ', ' stype_str{stype} ', N=' num2str(length(BlueOffSpike))])
        if find(Dist_order{f}(noi_dist{f})==1)~=1
            set(gca,'YTick',[1 find(Dist_order{f}(noi_dist{f})==1) sum(noi_dist{f})],'YTickLabel',num2str([min(interDendDist{f}(1,:).*Dsign) 0 max((interDendDist{f}(1,:).*Dsign))]',3))
        else
            set(gca,'YTick',[1 sum(noi_dist{f})],'YTickLabel',num2str([0 max((interDendDist{f}(1,:)))]',3))
        end
        set(gca,'XTick',[1 -nTau{stype}(1)+1 length(nTau{stype})],'XTickLabel',num2str([nTau{stype}(1) 0 nTau{stype}(end)]',3))
    end
end
linkaxes(ax1,'x')
colormap(turbo)

figure(100); clf;
for f=foi
    nexttile([1 1])
    imshow2(AvgImg{f},[])
    title(num2str(f))
end

%% Measure the Spike attenuation rate
clear SpikeAmp SpikeAmp_kink STAcrop STAmat_crop_kink SpikeDelay SpikeOrder AllSpikeAmp AllSpikeAmp_kink preSub AllSpikepreSub Allactivity preSubprox
pre_spike_window=[-5:-1];
pre_spike_window2=[-20:-1];
Amp_window=[0];
for f=foi
    nROI=size(NormalizedTrace_dirt{f},1);
    catMat=[];
    for stype=1:3
        if ~isempty(MatSTA{f,stype})
            [~, BlueOffSpike]=find_wCondition(SpikeList{f,stype,1},BlueStim{f}==0);
            catMat=cat(3,catMat,MatSTA{f,stype}(:,:,BlueOffSpike));
        end
    end

    AllSTA=mean(catMat,3,'omitnan');
    AllSTA=AllSTA-prctile(AllSTA,15,2);
    [~, maxAmptime]=max(AllSTA,[],2);
    maxAmptime=maxAmptime+nTau{1}(1);

    ISItr=get_BurstOrder(max([allSpikeMat{f}(1,:); CStrace{f}]),50);

    for stype=1:3
        for ns=1:size(MatSTA,3)
            if ~isempty(MatSTA{f,stype,ns})
                [~, SpInterest]=find_wCondition(SpikeList{f,stype,ns},BlueStim{f}==0);
                STAcrop{f,stype,ns}=mean(MatSTA{f,stype,ns}(:,-nTau{stype}(1)+[-3:5],SpInterest),3,'omitnan');
                SpikeDelay{f,stype}(:,ns)=get_delay(STAcrop{f,stype,ns},100);
                SpikeDelay{f,stype}=SpikeDelay{f,stype}-SpikeDelay{f,stype}(find(Dist_order{f}(noi_dist{f})==1),:);

                meanSTA=mean(MatSTA{f,stype,ns},3,'omitnan')-prctile(mean(MatSTA{f,stype,ns},3,'omitnan'),20,2);
                peripeakSTA=[];
                for n=1:size(meanSTA,1)
                    peripeakSTA(n,:)=meanSTA(n,maxAmptime(n)-nTau{stype}(1)+[-5:5]);
                    SpikeAUC{f,stype}(n,ns)=get_AUC(meanSTA(n,:),maxAmptime(n)-nTau{stype}(1),3,3);
                end
                SpikeAmp{f,stype}(:,ns)=max(peripeakSTA,[],2);
            end

            if ns==1
                [~, SpInterest]=find_wCondition(SpikeList{f,stype,ns},ISItr==1);
                meanSTA=mean(MatSTA{f,stype,ns}(:,:,SpInterest),3,'omitnan');
                preSub{f,stype}=mean(meanSTA(:,-nTau{stype}(1)+pre_spike_window2),2,'omitnan');
                preSubprox{f,stype}=mean(meanSTA(:,-nTau{stype}(1)+pre_spike_window),2,'omitnan');
            end
        end
    end
end

% Show triggering average traces
figure(102); clf;
nthspike=[1 5 5]; % ss cs bs
tiledlayout(sum(nthspike),26); g2=1;

for f=foi
    Dsign=ones(1,size(interDendDist{f},1));
    Dsign(Dist_order{f}(1:find(Dist_order{f}==1)-1))=-1;
    cax=[prctile(STAcrop{f,1,1}(:),1) prctile(STAcrop{f,1,1}(:),99)];
    g=1;
    for stype=1:3
        for ns=1:nthspike(stype)
            nexttile(g2+26*(g-1),[1 1]);
            imagesc(STAcrop{f,stype,ns},cax)
            if g==1;
                title(num2str(f))
            end
            g=g+1;

            if find(Dist_order{f}(noi_dist{f})==1)~=1
                set(gca,'YTick',[1 find(Dist_order{f}(noi_dist{f})==1) sum(noi_dist{f})],'YTickLabel',num2str([min(interDendDist{f}(1,:).*Dsign) 0 max((interDendDist{f}(1,:).*Dsign))]',3))
            else
                set(gca,'YTick',[1 sum(noi_dist{f})],'YTickLabel',num2str([0 max((interDendDist{f}(1,:)))]',3))
            end
        end
    end

    g2=g2+1;
end
colormap(turbo)

figure(103); clf;
numsp2show=[1 5 5]; %SS, CS, BS
conduction_Speed=NaN(max(foi),sum(numsp2show));
conduction_Speed_basal=NaN(max(foi),sum(numsp2show));

for f=foi
    g=1;
    for stype=1:3
        for s=1:numsp2show(stype)

            xaxis=dendaxis{f};
            ap=find(xaxis>=0);
            bs=find(xaxis<=0);

            if ~isempty(MatSTA{f,stype,s})

                [p S] = polyfit(xaxis(ap), SpikeDelay{f,stype}(ap,s), 1);
                y_fit = polyval(p, xaxis(ap));

                SS_res = sum((SpikeDelay{f,stype}(ap,s)' - y_fit).^2);
                SS_tot = sum((SpikeDelay{f,stype}(ap,s)' - mean(SpikeDelay{f,stype}(ap,s)')).^2);
                R_squared = 1 - (SS_res / SS_tot);

                if length(bs)>2
                    [p_bs S] = polyfit(-xaxis(bs), SpikeDelay{f,stype}(bs,s), 1);
                    y_fit_bs = polyval(p_bs, -xaxis(bs));
                    SS_res = sum((SpikeDelay{f,stype}(bs,s)' - y_fit_bs).^2);
                    SS_tot = sum((SpikeDelay{f,stype}(bs,s)' - mean(SpikeDelay{f,stype}(bs,s)')).^2);
                    R_squared_bs = 1 - (SS_res / SS_tot);
                else
                    p_bs=NaN;
                    R_squared_bs=0;
                end

                if R_squared>0.5
                    conduction_Speed(f,g)=1/p(1);
                else
                    conduction_Speed(f,g)=NaN;
                end

                if R_squared>0.7 && max(abs(xaxis(bs)))>80
                    conduction_Speed_basal(f,g)=1/p_bs(1);
                else
                    conduction_Speed_basal(f,g)=NaN;
                end

                g=g+1;
            else
                conduction_Speed(f,g)=NaN;
                conduction_Speed_basal(f,g)=NaN;
                g=g+1;
            end
        end
    end
end

cmap=hsv(size(conduction_Speed,2));
boxplot(conduction_Speed, 'PlotStyle', 'traditional', 'Colors', 'k'); hold all
for j = 1:size(conduction_Speed,2)
    % Scatter plot points for each category with a slight horizontal offset for visibility
    jitter = (rand(size(conduction_Speed,1), 1) - 0.5) * 0.1; % Small random jitter for horizontal spacing
    scatter(j + jitter, conduction_Speed(:, j), 35, cmap(j,:), 'filled', 'MarkerFaceAlpha', 0.6);
end
set(gca,'XTick',[1:size(conduction_Speed,2)*2],'XTickLabel',[{'SS'} counting_string([1:5],'CS ') counting_string([1:5],'BS ')])
ylabel('Conduction speed (\mum/ms)')

% Show attenuation of spike
dendriteaxis_bin=[-360 -200 -100 -50 50 100 200 300 500];
ylimval=[0 1.4]; ylimval2=[0 2.5];
figure(104); clf; tiledlayout(2,3)
cmap=distinguishable_colors(6); cmap=cmap([1 2 3 4 6 5],:);
nthspike=[1 4 4];
for stype=1:3
    nexttile(stype,[1 1])
    l=[];
    for s=1:nthspike(stype) % show Amplitude
        NormAmp_X=[];
        for f=foi
            if ~isempty(SpikeAmp{f,stype}) & s<=size(SpikeAmp{f,stype},2)
                NormAmp_X{f}=[dendaxis{f}' SpikeAmp{f,stype}(:,s)];
            else
                NormAmp_X{f}=[NaN NaN];
            end
            if stype==1
                scatter(NormAmp_X{f}(:,1),NormAmp_X{f}(:,2),20,[0.6 0.6 0.6],'filled','o','MarkerFaceAlpha',0.7); hold all
            end
        end

        [mean_amplitudes std_amplitudes dendBin_center ind]=binning_data(NormAmp_X(foi),dendriteaxis_bin);
        N_neuron=sum((cellfun(@sum,ind)),2)';
        l(s)=errorbar(dendBin_center, mean_amplitudes, std_amplitudes./sqrt(N_neuron), 'o-', 'CapSize', 5,'color',cmap(s,:),'LineWidth',2); hold all
    end
    xlabel('Distance from Soma (\mum)')
    ylabel(['Normalized bAP amplitude of ' stype_str{stype}])
    legend(l,counting_string([1:nthspike(stype)],[stype_str{stype} ' ']))
    ylim(ylimval);

    l=[];
    for s=1:nthspike(stype) % show AUC
        NormAmp_X=[];
        for f=foi
            nexttile(stype+3,[1 1])
            if ~isempty(SpikeAmp{f,stype}) & s<=size(SpikeAmp{f,stype},2)
                NormAmp_X{f}=[dendaxis{f}' SpikeAUC{f,stype}(:,s)];
            else
                NormAmp_X{f}=[NaN NaN];
            end
            if stype==1
                scatter(NormAmp_X{f}(:,1),NormAmp_X{f}(:,2),20,[0.6 0.6 0.6],'filled','o','MarkerFaceAlpha',0.7); hold all
            end
        end
        [mean_amplitudes std_amplitudes dendBin_center ind]=binning_data(NormAmp_X(foi),dendriteaxis_bin);
        N_neuron=sum((cellfun(@sum,ind)),2)';
        l(s)=errorbar(dendBin_center, mean_amplitudes, std_amplitudes./sqrt(N_neuron), 'o-', 'CapSize', 5,'color',cmap(s,:),'LineWidth',2); hold all
    end
    xlabel('Distance from Soma (\mum)')
    ylabel(['Normalized bAP AUC of ' stype_str{stype}])
    legend(l,counting_string([1:nthspike(stype)],[stype_str{stype} ' ']))
    ylim(ylimval2);
end

% Comparing amplitude of SS, CS, BS
figure(105); clf;
cmap=distinguishable_colors(6);
nexttile([1 1])
for stype=1:3
    NormAmp_X=[];
    for f=foi
        if ~isempty(SpikeAUC{f,stype})
            NormAmp=mean(SpikeAmp{f,stype}(:,1),2,'omitnan');
            NormAmp_X{f}=[dendaxis{f}' NormAmp];
        else
            NormAmp_X{f}=[NaN NaN];
        end
    end
    [mean_amplitudes std_amplitudes dendBin_center ind]=binning_data(NormAmp_X(foi),dendriteaxis_bin);
    N_neuron=sum((cellfun(@sum,ind)),2)';
    errorbar(dendBin_center, mean_amplitudes, std_amplitudes./sqrt(N_neuron), 'o-', 'CapSize', 5,'color',cmap(stype,:),'LineWidth',2); hold all
end
xlabel('Distance from Soma (\mum)')
ylabel('Normalized bAP amplitude')
legend({'SS','CS 1st','BS 1st'})
ylim(ylimval);

nexttile([1 1])
for stype=1:3
    NormAmp_X=[];
    for f=foi
        if ~isempty(SpikeAUC{f,stype})
            NormAmp=mean(SpikeAUC{f,stype}(:,1),2,'omitnan');
            NormAmp_X{f}=[dendaxis{f}' NormAmp];
        else
            NormAmp_X{f}=[NaN NaN];
        end
    end
    [mean_amplitudes std_amplitudes dendBin_center ind]=binning_data(NormAmp_X(foi),dendriteaxis_bin);
    N_neuron=sum((cellfun(@sum,ind)),2)';
    errorbar(dendBin_center, mean_amplitudes, std_amplitudes./sqrt(N_neuron), 'o-', 'CapSize', 5,'color',cmap(stype,:),'LineWidth',2); hold all
end
xlabel('Distance from Soma (\mum)')
ylabel('bAP AUC')
legend({'SS','CS 1st','BS 1st'})
ylim(ylimval2);

% Pre-subthreshold
figure(106); clf;
cmap=distinguishable_colors(6);
for stype=1:3
    PresubAmp=[];
    for f=foi
        if ~isempty(SpikeAmp{f,stype,1})

            PresubAmp{f}=[dendaxis{f}' preSub{f,stype}];
            PresubAmpProx{f}=[dendaxis{f}' preSubprox{f,stype}];

        else
            PresubAmp{f}=[NaN NaN];
            PresubAmpProx{f}=[NaN NaN];
        end
    end
    [mean_amplitudes std_amplitudes dendBin_center ind]=binning_data(PresubAmp(foi),dendriteaxis_bin);
    N_neuron=sum((cellfun(@sum,ind)),2)';
    nexttile(1,[1 1])
    errorbar(dendBin_center, mean_amplitudes, std_amplitudes./sqrt(N_neuron), 'o-', 'CapSize', 5,'color',cmap(stype,:),'LineWidth',2); hold all

    [mean_amplitudesProx std_amplitudesProx dendBin_center ind]=binning_data(PresubAmpProx(foi),dendriteaxis_bin);
    N_neuron=sum((cellfun(@sum,ind)),2)';
    nexttile(2,[1 1])
    errorbar(dendBin_center, mean_amplitudesProx, std_amplitudesProx./sqrt(N_neuron), 'o-', 'CapSize', 5,'color',cmap(stype,:),'LineWidth',2); hold all
end
nexttile(1,[1 1])
xlabel('Distance from Soma (\mum)')
ylabel('Mean Pre-spike amplitude (-20 ms to -1 ms)')
legend({'SS','CS','BS'})

nexttile(2,[1 1])
xlabel('Distance from Soma (\mum)')
ylabel('Mean Pre-spike amplitude (-5 ms to -1 ms)')
legend({'SS','CS','BS'})

%% Pre-spike subthreshold STA and Pre-spike spike rate STA 
figure(107); clf; 
tiledlayout(11,6);
clear STAspike STAspclass STAsub
prespike_timebin=[nTau{1}(1):20:20];
for f=foi
    cax=[prctile(MatSub{f,1,1}(:),5) prctile(MatSub{f,1,1}(:),95)];
    for stype=1:3
        ns=1;
        if ~isempty(MatSTA{f,stype,ns})
            [~, SpInterest]=find_wCondition(SpikeList{f,stype,ns},BlueStim{f}==0);
            if length(SpInterest)>10
            STAsub{f,stype}=mean(MatSub{f,stype,ns}(:,1:-nTau{1}(1),SpInterest),3,'omitnan');
            STAspike{f,stype}=[[nTau{1}(1):-1]' mean(MatSpike{f,stype,ns}(:,1:-nTau{1}(1),SpInterest),3,'omitnan')'];            
            for show_stype=1:3
            STAspclass{f,stype,show_stype}=[[nTau{1}(1):-1]' mean(MatSpClass{f,stype,ns}(show_stype,1:-nTau{1}(1),SpInterest),3,'omitnan')'];
            end
            end
            nexttile([1 1])
            imagesc(STAsub{f,stype},cax)
            title(['Neuron #' num2str(f) ',' stype_str{stype}])
            if find(Dist_order{f}(noi_dist{f})==1)~=1
                set(gca,'YTick',[1 find(Dist_order{f}(noi_dist{f})==1) sum(noi_dist{f})],'YTickLabel',num2str([min(dendaxis{f}) 0 max((dendaxis{f}))]',3))
            else
                set(gca,'YTick',[1 sum(noi_dist{f})],'YTickLabel',num2str([0 max((interDendDist{f}(1,:)))]',3))
            end
            set(gca,'XTick',[1 round(-nTau{stype}(1)/2) -nTau{stype}(1)],'XTickLabel',num2str([nTau{stype}(1) round(nTau{stype}(1)/2) -1]',3))
        end
    end
end
colormap(turbo)

figure(109); clf; 
tiledlayout(1,3);
edgesT=[nTau{stype}(1):5:0];
edgesX=[-150 -100 -50 50 100 200 300 500];
for stype=1:3
    showSTAsub=[]; g=1;
for f=foi
    if ~isempty(STAsub{f,stype})
    showSTAsub{g}=[[NaN [nTau{stype}(1):-1]]; [dendaxis{f}' STAsub{f,stype}]];
    g=g+1;
    end
end
nexttile([1 1])
[binnedZ binX binY]=show3Dbinning(showSTAsub, edgesT, edgesX, 'image'); hold all
caxis([0 0.3]);
xlabel('Pre-spike time (ms)')
ylabel('Distance from soma (\mum)')
title(['Mean subthreshold before ' stype_str{stype}])
set(gca, 'YDir', 'reverse');
colormap(turbo);
end

figure(108); clf;
tiledlayout(1,4);
ax1=nexttile(1,[1 1]);
for stype=1:3
 [mean_amplitudes std_amplitudes dendBin_center ind]=binning_data(STAspike(foi,stype),prespike_timebin);
    N_neuron=sum((cellfun(@sum,ind)),2)';
    errorbar(dendBin_center, 1000*mean_amplitudes, std_amplitudes./sqrt(N_neuron), 'o-', 'CapSize', 5,'color',cmap(stype,:),'LineWidth',2); hold all
end 
xlabel('Pre-spike time (ms)')
ylabel('Mean firing rate (Hz)')
legend({'SS','CS','BS'},'Location','northwest')
title('Mean firing rate before SS/CS/BS')

ax2=[];
for show_stype=1:3
ax2=[ax2 nexttile(1+show_stype,[1 1])];
for stype=1:3
 [mean_amplitudes std_amplitudes dendBin_center ind]=binning_data(STAspclass(foi,stype,show_stype),prespike_timebin);
    N_neuron=sum((cellfun(@sum,ind)),2)';
    errorbar(dendBin_center, 1000*mean_amplitudes, std_amplitudes./sqrt(N_neuron), 'o-', 'CapSize', 5,'color',cmap(stype,:),'LineWidth',2); hold all
end 
xlabel('Pre-spike time (ms)')
ylabel('Mean SS firing rate (Hz)')
legend({'SS','CS','BS'},'Location','northwest')
title(['Mean ' stype_str{show_stype}  ' firing rate before SS/CS/BS'])
end
linkaxes([ax2],'xy')

%% Subthreshold dynamics: Hierachy analysis and PCA analysis

area_threshold=50; perispike_time=[-2:20]; prespike_time=[-20:-2]; postspike_time=[0:20];
figure(112); clf; tiledlayout(11,8);
figure(113); clf; tiledlayout(11,8);
figure(114); clf; tiledlayout(11,8);
figure(115); clf; tiledlayout(11,8);
figure(116); clf; tiledlayout(11,8);
D_exc=NaN(max(foi),5); D_sub=NaN(max(foi),5); D_preSp=NaN(max(foi),5);
for f=foi
 Subth_trace=Subthreshold{f};
 nTime=size(Subthreshold{f},2);
 %Subth_trace=Subth_trace-movprc(Subth_trace,300,30,2);
    EXC_kymo=max(Subth_trace > prctile(Subth_trace(:,BlueStim{f}==0),90,2),[],1);
    EXC_kymobw=bwlabel(EXC_kymo);
    excbwprops=struct2array(regionprops(EXC_kymobw,'Area'));
    
    EXC_kymobw(ismember(EXC_kymobw,find(excbwprops<area_threshold)))=0;

    perispike_frame=unique([tovec(find(double(allSpikeMat{f}(1,:)==1))'+perispike_time); find(CStrace{f})']);
    perispike_frame(perispike_frame<=0 | perispike_frame>nTime)=[];
    nonvalid_frame=find(sum(isnan(Subthreshold{f}),1)>0);
    prespike_frame=unique([tovec(find(double(max(allSpikeClassMat{f})>0))'+prespike_time)]);
    postspike_frame=unique([tovec(find(double(allSpikeMat{f}(1,:)==1))'+postspike_time); find(CStrace{f})']);

    Blue_on_frame=find(imdilate(BlueStim{f}>0, [ones(1, 1), 1, ones(1, 200)]));
    Blue_off_frame=find(imdilate(BlueStim{f}==0, [1, ones(1, 50)]));

    Badframe=unique([perispike_frame' Blue_on_frame nonvalid_frame]);
    Goodframe=setdiff([1:nTime],Badframe);
    
    Badframe_preSp=unique([postspike_frame' Blue_on_frame nonvalid_frame]);
    Goodframe_preSp=setdiff(prespike_frame,Badframe);
    Goodframe_preSp(Goodframe_preSp<1)=[];

    EXC_kymobw(Badframe)=NaN;

    ExcCandFrame = bwlabel(imdilate(EXC_kymobw>0, [ones(1,5) 1 ones(1,5)]));
    ExcSub=Subthreshold{f}.*(ExcCandFrame>0);
    preSpFrame=bwlabel(ind2vec(nTime,Goodframe_preSp,1)>0);

    ExcPatt=[];
    for b=1:max(ExcCandFrame)
        bseg=find(ExcCandFrame==b);
        ExcPatt=[ExcPatt mean(Subthreshold{f}(:,bseg),2,'omitnan')];
    end
    % ExcCandFrame(ismember(ExcCandFrame,find(sum(isnan(ExcPatt),1)>0)))=0;
    % ExcCandFrame=bwlabel(ExcCandFrame>0);
    ExcPatt(:,sum(isnan(ExcPatt),1)>0)=[];

    ExcPatt_preSp=[];
    for b=1:max(preSpFrame)
        bseg=find(preSpFrame==b);
        ExcPatt_preSp=[ExcPatt_preSp mean(Subthreshold{f}(:,bseg),2,'omitnan')];
    end
    ExcPatt_preSp(:,sum(isnan(ExcPatt_preSp),1)>0)=[];

    [ExcPatt_reduce NclusterExc Excindmat]= hierachyCluster_BH(ExcPatt,4,0);
    [preSpPatt_reduce NclusterpreSp preSpindmat]= hierachyCluster_BH(ExcPatt_preSp,4,0);

    figure(112);
    for n=1:NclusterExc %Silent during blue off
        nexttile([1 1])
        showim=max(double(Ftprnts{f}(:,:,Dist_order{f}(noi_dist{f}))>0).*reshape((rescale(ExcPatt_reduce(:,n))+0.1),1,1,[]),[],3);
        imshow2(showim,[0 1])
        title(['Exc. Patt. #' num2str(n), ', Fraction :', num2str(sum(Excindmat(:,n)/size(Excindmat,1)),2)])
    end
    colormap(turbo)

    figure(113);
    for n=1:NclusterpreSp %Silent during blue off
        nexttile([1 1])
        showim=max(double(Ftprnts{f}(:,:,Dist_order{f}(noi_dist{f}))>0).*reshape((rescale(preSpPatt_reduce(:,n))+0.1),1,1,[]),[],3);
        imshow2(showim,[0 1])
        title(['Pre-Sp. Patt. #' num2str(n), ', Fraction :', num2str(sum(preSpindmat(:,n)/size(preSpindmat,1)),2)])
    end
    colormap(turbo)

    [ExcV ExcD ExcTrace]=get_eigvector(ExcSub(:,Goodframe)',size(ExcSub,1));
    [PreSpV PreSpD PreSpTrace]=get_eigvector(Subthreshold{f}(:,preSpFrame>0)',size(ExcSub,1));
    [SubV SubD subTrace]=get_eigvector(Subthreshold{f}(:,sum(isnan(Subthreshold{f}),1)==0)',size(ExcSub,1));
    
    cumD=cumsum([ExcD'; SubD'; PreSpD'],2)./sum([ExcD SubD PreSpD]',2);
    D_exc(f,:)=cumD(1,1:5);
    D_sub(f,:)=cumD(2,1:5);
    D_preSp(f,:)=cumD(3,1:5);

    figure(114);
    for pc=1:4
        nexttile([1 1])
        showScaleImage(Ftprnts{f}(:,:,Dist_order{f}(noi_dist{f}))>0, (ExcV(:,pc)),'turbo',[prctile(ExcV(:),1) prctile(ExcV(:),99)]);
        axis equal tight off
        if pc==1
        title(['N# ' num2str(f) ', Exc., ' num2str(ExcD(pc)/sum(ExcD)*100,2),' %'])
        else
        title([num2str(ExcD(pc)/sum(ExcD)*100,2),' %'])
        end
    end
    colormap(Aurora)

    figure(115);
    for pc=1:4
        nexttile([1 1])
        showScaleImage(Ftprnts{f}(:,:,Dist_order{f}(noi_dist{f}))>0, (PreSpV(:,pc)),'turbo',[prctile(PreSpV(:),1) prctile(PreSpV(:),99)]);
        axis equal tight off
        if pc==1
        title(['N# ' num2str(f) ', PreSp., ' num2str(PreSpD(pc)/sum(PreSpD)*100,2),' %'])
        else
        title([num2str(PreSpD(pc)/sum(PreSpD)*100,2),' %'])
        end
    end
    colormap(Aurora)

    figure(116);
    for pc=1:4
        nexttile([1 1])
        showScaleImage(Ftprnts{f}(:,:,Dist_order{f}(noi_dist{f}))>0, (SubV(:,pc)),'turbo',[prctile(SubV(:),1) prctile(SubV(:),99)]);
        axis equal tight off
        if pc==1
        title(['N# ' num2str(f) ', Sub., ' num2str(SubD(pc)/sum(SubD)*100,2),' %'])
        else
        title([num2str(SubD(pc)/sum(SubD)*100,2),' %'])
        end
    end
    colormap(Aurora)
end

figure(117); tiledlayout(1,3);
edgeVc=[5:5:105]/100;
for n=1:5
    nexttile(1,[1 1])
histogram(D_exc(:,n),edgeVc); hold all
nexttile(2,[1 1])
histogram(D_sub(:,n),edgeVc); hold all
nexttile(3,[1 1])
histogram(D_preSp(:,n),edgeVc); hold all
end
title('Centroid rank');
legend({'All subthreshold','Excitatory subthreshold','pre-Spike'})


%% Correlation matrix between ROIs
figure(118); clf
perispike_time=[-3:5];
Corrcoeff=[];
for f=foi
    nTime=size(Subthreshold{f},2);
    sub_ch=[]; sub_ch_dend=[]; 
    for ch=1:2
        sub_ch_dend{ch}=[];
        sub_ch{ch}=get_subthreshold(NormalizedTrace_ch{f,ch},max(allSpikeMat{f}(1,:),[],1)>0,7,17);
    end
    perispike_frame=unique([tovec(find(double(allSpikeMat{f}(1,:)==1))'+perispike_time); find(CStrace{f})']);
    perispike_frame(perispike_frame<=0 | perispike_frame>nTime)=[];
    nonvalid_frame=find(sum(isnan(Subthreshold{f}),1)>0);

    Blue_on_frame=find(imdilate(BlueStim{f}>0, [ones(1, 1), 1, ones(1, 200)]));
    Badframe=unique([perispike_frame' Blue_on_frame nonvalid_frame]);
    Goodframe=setdiff([1:nTime],Badframe);

    roisD_order_ind=cellfun(@find,roisD_order{f},'UniformOutput',false);
    labelvec=NaN(1,size(NormalizedTrace_dirt{f},1));
    labelclass=[];
    for dClass=1:5
        labelvec(cell2mat(roisD_order_ind(dClass,:)'))=dClass;
        for dend=1:size(roisD_order_ind,2)
            if ~isempty(roisD_order_ind{dClass,dend})
        labelclass=[labelclass dClass];
        sub_ch_dend{1}=[sub_ch_dend{1}; mean(sub_ch{1}(roisD_order_ind{dClass,dend},:),1,'omitnan')];
        sub_ch_dend{2}=[sub_ch_dend{2}; mean(sub_ch{2}(roisD_order_ind{dClass,dend},:),1,'omitnan')];
            end
        end
    end

    nexttile([1 1]);
    corrMat=get_corrMat(sub_ch{1},sub_ch{2},Goodframe);
    corrMat_dend=get_corrMat(sub_ch_dend{1},sub_ch_dend{2},Goodframe);
    Corrcoeff{f}=NaN(5,5);

    for dClassI=1:5
        for dClassJ=dClassI:5
            if dClassI==dClassJ
                if length(find(labelclass==dClassI))>1
                corrmat_tmp=corrMat_dend(find(labelclass==dClassI),find(labelclass==dClassJ));
                [a]=triu(ones(size(corrMat_dend(find(labelclass==dClassI),find(labelclass==dClassI)))),1);
                Corrcoeff{f}(dClassI,dClassJ)=mean(tovec(corrmat_tmp(a>0)),'omitnan');                
                else
                    Corrcoeff{f}(dClassI,dClassJ)=NaN;
                end
            else
                Corrcoeff{f}(dClassI,dClassJ)=mean(tovec(corrMat_dend(find(labelclass==dClassI),find(labelclass==dClassJ))),'omitnan');
            end
        end
    end

    cmap_label=hsv(5);
    cmap_label=cmap_label([4 2 3 5 1],:);
distMatrix = sqrt(2*(1 - corrMat));
distMatrix(logical(eye(size(distMatrix)))) = 0;
Z_ref = linkage(squareform(distMatrix), 'average');
leafOrder = optimalleaforder(Z_ref,squareform(distMatrix));
Cluster_ref= switchlabel(cluster(Z_ref, 'maxclust', 5));
[~, order_cluster]=sort(Cluster_ref,'ascend');

    l=imshow_label(corrMat(leafOrder,leafOrder),labelvec(leafOrder),cmap_label,{'Basal','Soma','Trunk','Oblique','Apical (distal)'});
    axis equal tight off
    if f==foi(end)
    else
        legend off
    end
    caxis([-0.3 1])
    title(num2str(f))
    drawnow;
end
colormap(turbo(256))

figure(119); clf;
label_str={'B','S','T','O','D'}; pair_str=[];
Corr_to_plot=[1 1; 3 3; 4 4; 5 5; 1 2; 1 3; 1 4; 1 5; 2 3; 2 4; 2 5; 3 4; 3 5; 4 5;]; g=1; show_c=[];
for f=foi
for p=1:size(Corr_to_plot,1)
    show_c(g,p)=Corrcoeff{f}(Corr_to_plot(p,1),Corr_to_plot(p,2));
    pair_str{p}=[label_str{Corr_to_plot(p,1)} '-' label_str{Corr_to_plot(p,2)}];
end
g=g+1;
end

%bar([1:size(Corr_to_plot,1)],mean(show_c,1,'omitnan'),'FaceColor',[0.5 .5 .5]); hold all
plot(1:size(Corr_to_plot,1),show_c','marker','o','color',[0.7 0.7 0.7]); hold all
x = repelem(1:size(show_c, 2), size(show_c, 1))'; % Group identifiers
y = show_c(:); % Flatten data for boxchart
b=boxchart(x, y,'BoxFaceColor', [0.2 0.5 0.9]); hold all
b.MarkerStyle = 'none'; b.BoxFaceAlpha=0.6;
b.WhiskerLineColor= [0.1 0.3 0.5]; b.LineWidth= 1.5;
%errorbar(1:size(Corr_to_plot,1),mean(show_c,1,'omitnan'),std(show_c,0,1,'omitnan')./sqrt(sum(~isnan(show_c))),'color','k','linewidth',2,'LineStyle','none','Marker','+')

set(gca,'XTick',[1:size(Corr_to_plot,1)],'XTickLabel',pair_str);
xlim([0.5 size(Corr_to_plot,1)+0.5]);
grid on
ylabel('Correlation coefficient')

%% Basal/Apical subthreshold dynamics

figure(124); clf;
figure(123); clf;
figure(121); clf;
figure(122); clf;
run_threshold=0.002;
PSD_sub=[];
preCS_SubBA=[];
preSS_SubBA=[];
BAxcorr_silent=NaN(max(foi),4001);
BAxcorr=NaN(max(foi),4001);
thetaPowerthres=0.4; phaseTarget=0;
FR_thetapeak=NaN(max(foi),3,2);
Run_thetapower=NaN(max(foi),2,2);
Run_BAunbal=NaN(max(foi),4);
Run_BA=NaN(max(foi),4);
PlaceField_BAunbal=NaN(max(foi),4);
PlaceField_BA=NaN(max(foi),4);
PlaceField_thetapower=NaN(max(foi),2);

FilterFreq=[5 12];
ax4=[];

for f=foi
    f
    nTime=size(Subthreshold{f},2);

    roisD_order_ind=cellfun(@find,roisD_order{f},'UniformOutput',false);

    if isempty(cell2mat(roisD_order_ind(1,:)'))
        basalind=cell2mat(roisD_order_ind(2,:)'); %if there is no basal, use soma
    else
        basalind=cell2mat(roisD_order_ind(1,:)');
    end
    apicalind=cell2mat(roisD_order_ind(5,:)'); %apical = distal dend

    filteredSubthreshold = pcafilterTrace(Subthreshold{f}([basalind; apicalind],:), 3); %filter subthreshold
    %filteredSubthreshold = NormalizedTrace_dirt{f};

    BasalSub = mean(filteredSubthreshold([1:length(basalind)],:),1,'omitnan');
    BasalSub(BlueStim{f}>0)=NaN;
    ApicalSub = mean(filteredSubthreshold([length(basalind)+1:end],:),1,'omitnan');
    ApicalSub(BlueStim{f}>0)=NaN;
    Basaltr = mean(NormalizedTrace_dirt{f}(basalind,:),1,'omitnan');
    Basaltr(BlueStim{f}>0)=NaN;
    Apicaltr = mean(NormalizedTrace_dirt{f}(apicalind,:),1,'omitnan');
    Apicaltr(BlueStim{f}>0)=NaN;

    [~, BApreCS]=get_STA([BasalSub; ApicalSub], allSpikeClassMat{f}(2,:).*(BlueStim{f}==0),15,-2); %pre complex spike
    BApreCS=mean(BApreCS,3,'omitnan');

    [~, BApreSS]=get_STA([BasalSub; ApicalSub], allSpikeClassMat{f}(1,:).*(BlueStim{f}==0),15,-2); %pre simple spike
    BApreSS=mean(BApreSS,3,'omitnan');

    [~, BApostSS]=get_STA([Basaltr; Apicaltr], allSpikeClassMat{f}(1,:).*(BlueStim{f}==0),0,3); %pre simple spike
    BApostSS=mean(BApostSS,3,'omitnan');

    [~, BApostCS1st]=get_STA([Basaltr; Apicaltr], allSpikeClassMat{f}(2,:).*(BlueStim{f}==0),0,3);
    BApostCS1st=mean(BApostCS1st,3,'omitnan');

    [~, BA_periSS]=get_STA([Basaltr; Apicaltr], allSpikeClassMat{f}(1,:).*(BlueStim{f}==0),1,3);

    Cstr=CStrace{f}.*(BlueStim{f}==0);  %CS trace during blue-off
    CSendpoint=ind2vec(nTime,find((Cstr(2:end)-Cstr(1:end-1))==-1),1);
    [~, BApostCS]=get_STA([BasalSub; ApicalSub], CSendpoint,5,20);
    BApostCS=mean(BApostCS,3,'omitnan');

    prespiketime=unique(find(max(allSpikeClassMat{f}).*(BlueStim{f}==0))'+[-40:-1]);
    prespiketime_vec=ind2vec(nTime,prespiketime(prespiketime>0),1);
    silenttime=setdiff([1:nTime],unique(find(max([allSpikeMat{f}(1,:); allSpikeClassMat{f}; CStrace{f}; BlueStim{f}])>0)'+[-5:30]));
    silenttime_vec=ind2vec(nTime,silenttime,1);
    runtime= double(imdilate(movmean((VRtrack{f}(end,:)),200)>run_threshold,[ones(1,2001)]));
    runtime(BlueStim{f}>0)=NaN;
    BASub_silent=[BasalSub; ApicalSub];
    BASub_silent(:,~silenttime_vec)=NaN;

    %Show Scatter heatmap of subthreshold dynamics
    % [V D]=get_eigvector(BASub_silent(:,sum(isnan(BASub_silent))==0),2);
    % clf; tiledlayout(2,2);
    % ax1=nexttile([1 1]); l=[];
    % scatter_heatmap2(BASub_silent(1,:),BASub_silent(2,:),linspace(-2,3,100),linspace(-2,3,100)); hold all
    % l=plot(tovec(squeeze(BA_periSS(1,:,:))'),tovec(squeeze(BA_periSS(2,:,:))'),'w');
    % plot([-2 3],[0 0],'color',[0.7 0.7 0.7]); plot([0 0],[-2 3],'color',[0.7 0.7 0.7]);
    % quiver(mean(BASub_silent(1,:),'omitnan'), mean(BASub_silent(2,:),'omitnan'), V(1,1)*sqrt(D(1)), V(2,1)*sqrt(D(1)), 'r', 'LineWidth', 1.5, 'DisplayName', 'Eigenvector 1');
    % quiver(mean(BASub_silent(1,:),'omitnan'), mean(BASub_silent(2,:),'omitnan'), V(1,2)*sqrt(D(2)), V(2,2)*sqrt(D(2)), 'g', 'LineWidth', 1.5, 'DisplayName', 'Eigenvector 2');
    % colormap(turbo); xlabel('Basal'); ylabel('Apical'); legend(l,{'SS'})
    % ax2=nexttile([1 1]);
    % scatter_heatmap2(BASub_silent(1,:),BASub_silent(2,:),linspace(-2,3,100),linspace(-2,3,100)); hold all
    % l=plot(Basaltr(Cstr>0),Apicaltr(Cstr>0),'w');
    % plot([-2 3],[0 0],'color',[0.7 0.7 0.7]); plot([0 0],[-2 3],'color',[0.7 0.7 0.7]);
    % quiver(mean(BASub_silent(1,:),'omitnan'), mean(BASub_silent(2,:),'omitnan'), V(1,1)*sqrt(D(1)), V(2,1)*sqrt(D(1)), 'r', 'LineWidth', 1.5, 'DisplayName', 'Eigenvector 1');
    % quiver(mean(BASub_silent(1,:),'omitnan'), mean(BASub_silent(2,:),'omitnan'), V(1,2)*sqrt(D(2)), V(2,2)*sqrt(D(2)), 'g', 'LineWidth', 1.5, 'DisplayName', 'Eigenvector 2');
    % colormap(turbo); xlabel('Basal'); ylabel('Apical'); legend(l,{'CS'})
    % linkaxes([ax1 ax2],'xy')
    % ax3=nexttile([1 1]); l=[];
    % scatter_heatmap2(BASub_silent(1,:),BASub_silent(2,:),linspace(-2,3,100),linspace(-2,3,100)); hold all
    % plot([-2 3],[0 0],'color',[0.7 0.7 0.7]); plot([0 0],[-2 3],'color',[0.7 0.7 0.7]);
    % l(1)=scatter(BApreCS(1,:),BApreCS(2,:),40,'filled','Marker','>','MarkerFaceColor',[1 0 0],'MarkerFaceAlpha',0.7); hold all
    % l(2)=scatter(BApreSS(1,:),BApreSS(2,:),40,'filled','Marker','o','MarkerFaceColor',[0 0 0],'MarkerFaceAlpha',0.7);
    % quiver(mean(BASub_silent(1,:),'omitnan'), mean(BASub_silent(2,:),'omitnan'), V(1,1)*sqrt(D(1)), V(2,1)*sqrt(D(1)), 'r', 'LineWidth', 1.5, 'DisplayName', 'Eigenvector 1');
    % quiver(mean(BASub_silent(1,:),'omitnan'), mean(BASub_silent(2,:),'omitnan'), V(1,2)*sqrt(D(2)), V(2,2)*sqrt(D(2)), 'g', 'LineWidth', 1.5, 'DisplayName', 'Eigenvector 2');
    % colormap(turbo); xlabel('Basal'); ylabel('Apical'); legend(l,{'pre CS','pre SS'})
    % ax4=nexttile([1 1]); l=[];
    % scatter_heatmap2(BASub_silent(1,:),BASub_silent(2,:),linspace(-2,3,100),linspace(-2,3,100)); hold all
    % plot([-2 3],[0 0],'color',[0.7 0.7 0.7]); plot([0 0],[-2 3],'color',[0.7 0.7 0.7]);
    % l(1)=scatter(BApostCS1st(1,:),BApostCS1st(2,:),40,'filled','Marker','>','MarkerFaceColor',[1 0 0],'MarkerFaceAlpha',0.7); hold all
    % l(2)=scatter(BApostSS(1,:),BApostSS(2,:),40,'filled','Marker','o','MarkerFaceColor',[0 0 0],'MarkerFaceAlpha',0.7);
    % quiver(mean(BASub_silent(1,:),'omitnan'), mean(BASub_silent(2,:),'omitnan'), V(1,1)*sqrt(D(1)), V(2,1)*sqrt(D(1)), 'r', 'LineWidth', 1.5, 'DisplayName', 'Eigenvector 1');
    % quiver(mean(BASub_silent(1,:),'omitnan'), mean(BASub_silent(2,:),'omitnan'), V(1,2)*sqrt(D(2)), V(2,2)*sqrt(D(2)), 'g', 'LineWidth', 1.5, 'DisplayName', 'Eigenvector 2');
    % colormap(turbo); xlabel('Basal'); ylabel('Apical'); legend(l,{'3-8 ms after 1st CS','3-8 ms after 1st SS'})
    % linkaxes([ax3 ax4],'xy')
    % 
    %BAxcorr_silent(f,:)=nanXCorr(BASub_silent(1,:),BASub_silent(2,:),2000,1);
    %BAxcorr(f,:)=nanXCorr(BasalSub,ApicalSub,2000,1);

    [freqPSD ~, PSD_bs]=nanPSD(BasalSub,1000,10000);
    [freqPSD ~, PSD_ap]=nanPSD(ApicalSub,1000,10000);
    PSD_sub{f,1}=[freqPSD PSD_bs];
    PSD_sub{f,2}=[freqPSD PSD_ap];

    BasalSub_intNan=mean(interpolateNaN(Subthreshold{f}(basalind,:)),1,'omitnan');
    ApicalSub_intNan=mean(interpolateNaN(Subthreshold{f}(apicalind,:)),1,'omitnan');

    % [Basalwvt freq]=get_waveletTransform(BasalSub_intNan,1000,[0 50]);
    % [Apicalwvt freq]=get_waveletTransform(ApicalSub_intNan,1000,[0 50]);
    % thetaInd=find(freq>5 & freq<12);
    % thetaPowerBA=[mean(abs(Basalwvt(thetaInd,:)),1,'omitnan'); mean(abs(Apicalwvt(thetaInd,:)),1,'omitnan')];

    [BasalPhase BasalSubFilt BasalthetaPower] = get_phase(BasalSub_intNan, 1000, FilterFreq);
    [ApicalPhase ApicalSubFilt ApicalthetaPower] = get_phase(ApicalSub_intNan, 1000, FilterFreq);
    Run_thetapower(f,:,1)=[mean(BasalthetaPower(runtime>0)) mean(BasalthetaPower(runtime==0))];
    Run_thetapower(f,:,2)=[mean(ApicalthetaPower(runtime>0)) mean(ApicalthetaPower(runtime==0))];

    crossingVecApical = [((ApicalPhase(1:end-1) < phaseTarget) & (ApicalPhase(2:end) > phaseTarget))' 0];
    crossingVecBasal = [((BasalPhase(1:end-1) < phaseTarget) & (BasalPhase(2:end) > phaseTarget))' 0];
    ApicalphaseVec=(crossingVecApical & zscore(ApicalthetaPower)>1.5);
    ApicalphaseVec_silent=(ApicalphaseVec & silenttime_vec);
    BasalphaseVec=(crossingVecBasal & zscore(BasalthetaPower)>1.5);
    BasalphaseVec_silent=(BasalphaseVec & silenttime_vec);

    PhaseSTAApical=get_STA(Subthreshold{f},ApicalphaseVec_silent,200,200);
    PhaseSTABasal=get_STA(Subthreshold{f},BasalphaseVec_silent,200,200);
    PhaseSTASpclassBA=get_STA(allSpikeClassMat{f},BasalphaseVec_silent,1000,1000);
    %
    FR_thetapeak(f,[1:3],1)=mean(PhaseSTASpclassBA(:,150+[-80:60]),2,'omitnan')'*1000;
    PhaseSTASpclass=get_STA(allSpikeClassMat{f},ApicalphaseVec_silent,1000,1000);

    [PhaseSTASpclass_bin binaxis]= Bin_Vector(PhaseSTASpclass, [-1000:1000], [-1000:25:1000]);
    [PhaseSTASpclassBA_bin binaxis]= Bin_Vector(PhaseSTASpclassBA, [-1000:1000], [-1000:25:1000]);
    FR_thetapeak(f,[1:3],2)=mean(PhaseSTASpclass(:,150+[-80:60]),2,'omitnan')'*1000;

    [~, sortind]=sort(dendaxis{f},'ascend');
    [Xq, Yq] = meshgrid([-200:200], min(dendaxis{f}):30:max(dendaxis{f}));
    normTrace_showqApical = interp2([-200:200], dendaxis{f}(sortind), PhaseSTAApical(sortind,:), Xq, Yq, 'linear');
    normTrace_showqBasal = interp2([-200:200], dendaxis{f}(sortind), PhaseSTABasal(sortind,:), Xq, Yq, 'linear');

    BAunbalance=BASub_silent(1,:)-BASub_silent(2,:);
    Run_BAunbal(f,:)=[std(BAunbalance(runtime>0),'omitnan') std(BAunbalance(runtime==0),'omitnan') mean(BAunbalance(runtime>0),'omitnan') mean(BAunbalance(runtime==0),'omitnan')];
    Run_BA(f,:)=[mean(BASub_silent(:,runtime>0),2,'omitnan')' mean(BASub_silent(:,runtime==0),2,'omitnan')'];

    if ~isempty(PlaceFieldList{f}) % in place field
        binTrack=(ceil(VRtrack{f}(5,:)/((115)/150)));
        PFvec=zeros(1,nTime);
        for p=1:length(PlaceFieldBin{f})/2
            if PlaceFieldBin{f}(2*(p-1)+1)>PlaceFieldBin{f}(2*(p-1)+2)
                Pvec=~(binTrack<(PlaceFieldBin{f}(2*(p-1)+1)) & binTrack>(PlaceFieldBin{f}(2*(p-1)+2)));
            else
                Pvec=(binTrack>(PlaceFieldBin{f}(2*(p-1)+1)) & binTrack<(PlaceFieldBin{f}(2*(p-1)+2)));
            end
            Lapvec=(VRtrack{f}(8,:)>PlaceFieldList{f}(2*(p-1)+1) & VRtrack{f}(8,:)<PlaceFieldList{f}(2*(p-1)+2));
            PFvec=PFvec| (Lapvec & Pvec);
        end
        PlaceField_BAunbal(f,:)=[std(BAunbalance(PFvec>0),'omitnan') std(BAunbalance(PFvec==0),'omitnan') mean(BAunbalance(PFvec>0),'omitnan') mean(BAunbalance(PFvec==0),'omitnan')];
    end

    PlaceField_thetapower(f,:,1)=[mean(BasalthetaPower(PFvec>0)) mean(BasalthetaPower(PFvec==0))];
    PlaceField_thetapower(f,:,2)=[mean(ApicalthetaPower(PFvec>0)) mean(ApicalthetaPower(PFvec==0))];
    PlaceField_BA(f,:)=[mean(BASub_silent(:,PFvec>0),2,'omitnan')' mean(BASub_silent(:,PFvec==0),2,'omitnan')'];

    figure(121); clf;
    nexttile([1 1])
    imagesc([-200:200],Yq(:,1),normTrace_showqBasal); colormap(turbo);
    xlabel('Time (ms)')
    ylabel('Distance from soma (\mum)')
    title(['Basal ,' num2str(f)]);
    drawnow;
    
     nexttile([1 1])
    imagesc([-200:200],Yq(:,1),normTrace_showqApical); colormap(turbo);
    xlabel('Time (ms)')
    ylabel('Distance from soma (\mum)')
    title(['Apical ,' num2str(f)]);
    drawnow;

    figure(122);
    nexttile([1 1])
    plot(binaxis,movmean(PhaseSTASpclass_bin*1000,3,2)');
    xlabel('Time (ms)')
    ylabel('Spike rate (s^{-1})')
    title(f);
    drawnow;

    figure(123);
    nexttile([1 1])
    plot(binaxis,movmean(PhaseSTASpclassBA_bin*1000,3,2)');
    xlabel('Time (ms)')
    ylabel('Spike rate (s^{-1})')
    title(f);
    drawnow;

    figure(124);
    ax4=[ax4 nexttile([1 1])];
    for spclass=1:3
        BAunbal_STA=get_STA(BAunbalance,allSpikeClassMat{f}(spclass,:).*(BlueStim{f}==0),500,400);
        plot(BAunbal_STA); hold all
    end
    title(f);
    drawnow;

end
linkaxes(ax4,'xy');

%plot temporal variance

% spike rate around theta peak 
figure(126); clf; cmap=lines(3);
for spclass=1:3
    plot([1 2]+2*(spclass-1),squeeze(FR_thetapeak(foi,spclass,:)),'marker','o','color',[0.6 0.6 0.6]); hold all
    FRtmp=squeeze(FR_thetapeak(foi,spclass,:));
    M=mean(FRtmp,1,'omitnan'); S=std(FRtmp,0,1,'omitnan')./sqrt(length(foi));
    errorbar([1 2]+2*(spclass-1),M,S,'color',cmap(spclass,:))
    [~, p]=ttest(FRtmp(:,1),FRtmp(:,2));
    line([1 2]+2*(spclass-1),[max(FRtmp(:))*1.1 max(FRtmp(:))*1.1],'color',[0 0 0]);
    text(1.5+2*(spclass-1),max(FRtmp(:))*1.1+0.5,['p = ' num2str(p,3)],'HorizontalAlignment','center')
end
ylabel('Spike rate (s^{-1})')
xlim([0.5 6.5])
set(gca,'XTick',[1:6],'XTickLabel',{'Basal peak SS','Apical peak SS','Basal peak CS','Apical peak CS','Basal peak BS','Apical peak BS'})


figure(125); clf; cmap=lines(3); tiledlayout(2,3);

% Basal & apical during run or rest
nexttile([1 1]);
for r=1:2
    plot([1 2]+(r-1)*2,Run_BA(foi,[1 3]+r-1),'marker','o','color',[0.6 0.6 0.6]); hold all
    M=mean(Run_BA(foi,[1 3]+r-1),1,'omitnan'); S=std(Run_BA(foi,[1 3]+r-1),0,1,'omitnan')./sqrt(length(foi));
    errorbar([1 2]+(r-1)*2,M,S,'color',cmap(2,:));
    [p]=signrank(Run_BA(foi,[1]+r-1),Run_BA(foi,[3]+r-1));
    mm=max(tovec(Run_BA(foi,[1 3]+r-1)));
    line([1 2]+(r-1)*2,[mm*1.1 mm*1.1],'color',[0 0 0]);
    text(1.5+(r-1)*2,mm*1.1+0.01,['p = ' num2str(p,2)],'HorizontalAlignment','center')
end
xlim([0.5 4.5]);
ylabel('Voltage (spike height')
set(gca,'XTick',[1:4],'XTickLabel',{'Basal, Run','Basal, Rest','Apical, Run','Apical, Rest'})

% Basal & apical in PF or out PF
nexttile([1 1]);
for r=1:2
    plot([1 2]+(r-1)*2,PlaceField_BA(foi,[1 3]+r-1),'marker','o','color',[0.6 0.6 0.6]); hold all
    M=mean(PlaceField_BA(foi,[1 3]+r-1),1,'omitnan'); S=std(PlaceField_BA(foi,[1 3]+r-1),0,1,'omitnan')./sqrt(length(foi));
    errorbar([1 2]+(r-1)*2,M,S,'color',cmap(2,:));
    [p]=signrank(PlaceField_BA(foi,[1]+r-1),PlaceField_BA(foi,[3]+r-1));
    mm=max(tovec(PlaceField_BA(foi,[1 3]+r-1)));
    line([1 2]+(r-1)*2,[mm*1.1 mm*1.1],'color',[0 0 0]);
    text(1.5+(r-1)*2,mm*1.1+0.01,['p = ' num2str(p,2)],'HorizontalAlignment','center')
end
xlim([0.5 4.5]);
ylabel('Voltage (spike height')
set(gca,'XTick',[1:4],'XTickLabel',{'Basal, in PF','Basal, out PF','Apical, in PF','Apical, out PF'})

% Basal - apical during run or rest
nexttile([1 1]);
plot([1 2],Run_BAunbal(foi,3:4),'marker','o','color',[0.6 0.6 0.6]); hold all
M=mean(Run_BAunbal(foi,3:4),1,'omitnan'); S=std(Run_BAunbal(foi,3:4),0,1,'omitnan')./sqrt(length(foi));
errorbar([1 2],M,S,'color',cmap(1,:));
[p]=signrank(Run_BAunbal(foi,3),Run_BAunbal(foi,4));
xlim([0.5 2.5])
set(gca,'XTick',[1 2],'XTickLabel',{'Run','Rest'})
title(['p = ' num2str(p,2)])
ylabel('Basal- Apical (Spike height)')

% Basal - Apical in PF or Out PF
nexttile([1 1]);
plot([1 2],PlaceField_BAunbal(foi,3:4),'marker','o','color',[0.6 0.6 0.6]); hold all
M=mean(PlaceField_BAunbal(foi,3:4),1,'omitnan'); S=std(PlaceField_BAunbal(foi,3:4),0,1,'omitnan')./sqrt(length(foi));
errorbar([1 2],M,S,'color',cmap(3,:));
[p]=signrank(PlaceField_BAunbal(foi,3),PlaceField_BAunbal(foi,4));
xlim([0.5 2.5])
ylabel('Basal- Apical (Spike height)')
set(gca,'XTick',[1 2],'XTickLabel',{'In PF','Out PF'})
title(['p = ' num2str(p,2)])

% Theta power during run or rest
nexttile([1 1]);
for r=1:2
    plot([1 2]+(r-1)*2,Run_thetapower(foi,:,r),'marker','o','color',[0.6 0.6 0.6]); hold all
    M=mean(Run_thetapower(foi,:,r),1,'omitnan'); S=std(Run_thetapower(foi,:,r),0,1,'omitnan')./sqrt(length(foi));
    errorbar([1 2]+(r-1)*2,M,S,'color',cmap(2,:));
    [p]=signrank(Run_thetapower(foi,1,r),Run_thetapower(foi,2,r));
    mm=max(tovec(Run_thetapower(foi,:,r)));
    line([1 2]+(r-1)*2,[mm*1.1 mm*1.1],'color',[0 0 0]);
    text(1.5+(r-1)*2,mm*1.1+0.01,['p = ' num2str(p,2)],'HorizontalAlignment','center')
end
xlim([0.5 4.5]);
ylabel('Theta-band power')
set(gca,'XTick',[1:4],'XTickLabel',{'Basal, Run','Basal, Rest','Apical, Run','Apical, Rest'})

% Theta power in PF or out PF
nexttile([1 1]);
for r=1:2
    plot([1 2]+(r-1)*2,PlaceField_thetapower(foi,:,r),'marker','o','color',[0.6 0.6 0.6]); hold all
    M=mean(PlaceField_thetapower(foi,:,r),1,'omitnan'); S=std(PlaceField_thetapower(foi,:,r),0,1,'omitnan')./sqrt(length(foi));
    errorbar([1 2]+(r-1)*2,M,S,'color',cmap(2,:));
    [p]=signrank(PlaceField_thetapower(foi,1,r),PlaceField_thetapower(foi,2,r));
    mm=max(tovec(PlaceField_thetapower(foi,:,r)));
    line([1 2]+(r-1)*2,[mm*1.1 mm*1.1],'color',[0 0 0]);
    text(1.5+(r-1)*2,mm*1.1+0.01,['p = ' num2str(p,2)],'HorizontalAlignment','center')
end
xlim([0.5 4.5]);
ylabel('Theta-band power')
set(gca,'XTick',[1:4],'XTickLabel',{'Basal, in PF','Basal, out PF','Apical, in PF','Apical, out PF'})






figure(127); clf;
plot([-2000:2000],BAxcorr_silent(foi,:)','color',[0.8 0.8 0.8]); hold all
errorbar_shade([-2000:2000],mean(BAxcorr_silent(foi,:),1,'omitnan'),std(BAxcorr_silent(foi,:),0,1,'omitnan'));
xlabel('\tau (ms)')
ylabel('R(\tau)')

figure(128); clf; cmap=distinguishable_colors(2);
[mean_PSDbs std_PSDbs b_cent]=binning_data(PSD_sub(foi,1),[0:1:20]);
[mean_PSDap std_PSDap b_cent]=binning_data(PSD_sub(foi,2),[0:1:20]);
errorbar_shade(b_cent,mean_PSDbs,std_PSDbs,cmap(1,:)); hold all
errorbar_shade(b_cent,mean_PSDap,std_PSDap,cmap(2,:));

% 
% sample_gridSS=linspace(-1,2,50);
% sample_gridCS=linspace(-1,2,30);
% interp_grid=linspace(-1,2,100);
% [Xq,Yq] = meshgrid(interp_grid,interp_grid);
% 
% clf;
% nexttile([1 1])
% [Prob_preCSsubBA xc yc]=scatter_heatmap2(preCS_SubBA(1,:),preCS_SubBA(2,:),sample_gridCS,sample_gridCS);
% shading faceted
% [Xmesh, Ymesh] = meshgrid(xc,yc);
% Prob_preCSsubBA=interp2(Xmesh,Ymesh,Prob_preCSsubBA,Xq,Yq);
% 
% nexttile([1 1])
% [Prob_preSSsubBA xc yc]=scatter_heatmap2(preSS_SubBA(1,:),preSS_SubBA(2,:),sample_gridSS,sample_gridSS);
% [Xmesh, Ymesh] = meshgrid(xc,yc);
% Prob_preSSsubBA=interp2(Xmesh,Ymesh,Prob_preSSsubBA,Xq,Yq);
% colormap(turbo)
% ax = gca; % Get current axes
% ax.GridColor = 'w'; % Set grid color (e.g., black)
% ax.GridAlpha = 0.5; % Set transparency of the grid lines























