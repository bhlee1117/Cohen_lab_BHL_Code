% Analysis on AAV expression sample and plot, YQ601
% 2025/07/20, Byung Hun Lee
clear
clc;
[~, ~, raw] = xlsread(['/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/NaVInactivationData_Arrangement.xlsx'], ...
    'Sheet1', 'C5:Q415');
fpath=raw(:,1);
Mouse=cell2mat(raw(:,2));
NeuronInd=cell2mat(raw(:,3));
Drug=raw(:,4);
pixelsize=cell2mat(raw(:,5));
set(0,'DefaultFigureWindowStyle','docked')
%% Load data
NeuronMat=[]; CatTraces=[]; g=1;
BlueWvf_str={'tri','pulse3','RP','IncT'};
for f=1:length(fpath)
    f
    load([fpath{f} '/SomOP_Result.mat'])
    for n=find(Result.isGood>0)
        Subth=get_subthreshold(Result.normtrace(n,:),Result.spike(n,:),9,20);
        %F0PCA=get_F0PCA(Subth,1);
        F0PCA=sqrt(mean(Subth(2:end).*Subth(1:end-1),'omitnan'));
        NormTr=Result.normtrace(n,:)./F0PCA;
        CatTraces{g,1}=NormTr;
        CatTraces{g,2}=Result.Blue;
        CatTraces{g,3}=Result.spike(n,:);
        CatTraces{g,4}=Result.ref_im;
        CatTraces{g,5}=Result.SpClass(:,:,n);
        CatTraces{g,6}=Result.CStrace(n,:);
        try
        CatTraces{g,7}=Result.ftprnt(:,:,n);
        catch
        CatTraces{g,7}=Result.c_ftprnt(:,:,n);
        end
        switch max(bwlabel(Result.Blue))
            case 2
                BlueWvf=1; %Triangle
            case 3
                BlueWvf=2;
            case 29
                BlueWvf=4; %increasing blue pulse width
            case 24
                BlueWvf=3; %RP
            case 20
                BlueWvf=4; %increasing blue pulse width
        end
        NeuronMat=[NeuronMat; [f Mouse(f) NeuronInd(f) max(bwlabel(Result.Blue)) BlueWvf Result.NtypeInd contains(Drug{f},'Ket') pixelsize(f)]];
        g=g+1;
    end
end
fieldName={'FileID','MouseID','NeuronID','N_blueStim','BlueWvf','Ntype','IsKet','pixelsize'};
NeuronMat=array2table(NeuronMat,'VariableNames',fieldName);

fieldName={'TraceV','Blue','Spike','ref_im','SpClass','CStrace','ftprnt'};
CatTraces=array2table(CatTraces,'VariableNames',fieldName);
CatTraces_double=table2array(CatTraces);
%% Show traces
lengthV=cellfun(@length,CatTraces_double);
ShowInd=find(NeuronMat.BlueWvf==3 & lengthV(:,1)<15000 & NeuronMat.Ntype==1); %RP
FindBlueOn=cellfun(@(x,y) find(bwlabel(x>0)==(y-4),1)-400,CatTraces_double(:,2),num2cell(NeuronMat.N_blueStim',1)','UniformOutput',false);

figure(22); clf; hold all;
l=cellfun(@(x,y,z) plot(rescale(x(z:end))+y),CatTraces_double(ShowInd,1),num2cell([1:length(ShowInd)], 1)',FindBlueOn(ShowInd));
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(turbo(length(l)),2));
plot(rescale(CatTraces_double{ShowInd(1),2}(FindBlueOn{1}:end))*2-1.5,'color',[0 0.6 1]);
axis off tight;  xlim([0 8000]); ylim([-2 length(l)+1.5]);

ShowInd=find(NeuronMat.BlueWvf==1 & NeuronMat.Ntype==1); %Tri
lengthV=cellfun(@length,CatTraces_double);
FindBlueOn=cellfun(@(x,y) find(bwlabel(x>0)==(y-4),1)-500,CatTraces_double(:,2),num2cell(NeuronMat.N_blueStim',1)','UniformOutput',false);
figure(23); clf; hold all;
l=cellfun(@(x,y) plot(rescale(x)+y),CatTraces_double(ShowInd,1),num2cell([1:length(ShowInd)], 1)');
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(turbo(length(l)),2));
plot(rescale(CatTraces_double{ShowInd(1),2})*2-1.5,'color',[0 0.6 1]);
axis off tight; xlim([0 14000]); ylim([-2 length(l)+1.5]);

ShowInd=find(NeuronMat.BlueWvf==4 & NeuronMat.Ntype==1); %IncT
lengthV=cellfun(@length,CatTraces_double);
FindBlueOn=cellfun(@(x,y) find(bwlabel(x>0)==(y-4),1)-500,CatTraces_double(:,2),num2cell(NeuronMat.N_blueStim',1)','UniformOutput',false);
figure(24); clf; hold all;
l=cellfun(@(x,y) plot(rescale(x)+y),CatTraces_double(ShowInd,1),num2cell([1:length(ShowInd)], 1)');
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(turbo(length(l)),2));
plot(rescale(CatTraces_double{ShowInd(1),2})*2-1.5,'color',[0 0.6 1]);
axis off tight; xlim([0 14000]); ylim([-2 length(l)+1.5]);

ShowInd=find(NeuronMat.BlueWvf==4 & NeuronMat.Ntype==2); %IncT
lengthV=cellfun(@length,CatTraces_double);
FindBlueOn=cellfun(@(x,y) find(bwlabel(x>0)==(y-4),1)-500,CatTraces_double(:,2),num2cell(NeuronMat.N_blueStim',1)','UniformOutput',false);
figure(25); clf; hold all;
l=cellfun(@(x,y) plot(rescale(x)+y),CatTraces_double(ShowInd,1),num2cell([1:length(ShowInd)], 1)');
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(turbo(length(l)),2));
plot(rescale(CatTraces_double{ShowInd(1),2})*2-1.5,'color',[0 0.6 1]);
axis off tight; xlim([0 14000]); ylim([-2 length(l)+1.5]);
%% Statistics Cell type
nTau_celltype=[20 20];
ADP_nTau=[4 6]; ISI_edge=[0:2:40]; timeinterp=[-2:0.2:2];
CelltypeMat=[]; Celltype_STAs=[]; Celltype_refim=[];
for n=1:size(CatTraces,1)
    if ~all(isnan(CatTraces.TraceV{n}))
    STAtr=get_STA(CatTraces.TraceV{n},CatTraces.Spike{n},nTau_celltype(1),nTau_celltype(2));
    STAtr_norm=STAtr./(STAtr(21)-mean(STAtr(1:nTau_celltype(1))));
    STAtr_norminterp=interp1([-2:2],STAtr_norm(nTau_celltype(1)+1+[-2:2]),timeinterp,'spline');
    Spikewidth=sum(STAtr_norminterp>max(STAtr_norminterp)/2)*0.2;
    ADP=mean(STAtr_norm(nTau_celltype(1)+[ADP_nTau(1):ADP_nTau(2)]))-mean(STAtr_norm(1:nTau_celltype(1)));
    sptimes=find(CatTraces.Spike{n});
    [counts] = histcounts(diff(sptimes), ISI_edge); 
    [~, max_idx] = max(counts); 
    PeakISI = (ISI_edge(max_idx) + ISI_edge(max_idx + 1)) / 2;
    Celltype_STAs{n}=STAtr_norm';

    bound=round(25/NeuronMat.pixelsize(n));
    if bound<size(CatTraces.ref_im{n},1)/3    
    ref_im=CatTraces.ref_im{n}-medfilt2(CatTraces.ref_im{n},[bound bound]);
    ref_im=ref_im(round(bound/2):end-round(bound/2),round(bound/2):end-round(bound/2));
    else
    ref_im=CatTraces.ref_im{n};
    end
    [imcounts, imedge]=histcounts(ref_im,100);
    level = otsuthresh(imcounts);
    im_threshold=imedge(round(level * (length(imcounts) - 1)));
    Radius=struct2array(regionprops(ref_im>im_threshold,'EquivDiameter')).*NeuronMat.pixelsize(n);
    Radius=max(Radius);
    CelltypeMat=[CelltypeMat; [ADP PeakISI im_threshold Radius Spikewidth NeuronMat.Ntype(n)]];
    Celltype_refim{n}=ref_im>im_threshold;
    else
    CelltypeMat=[CelltypeMat; NaN(1,6)];
    Celltype_STAs{n}=NaN;
    Celltype_refim{n}=NaN;
    end
end

fieldName={'ADP','PeakISI','im_threshold','Radius','Spikewidth','ExpectedNtype'};
CelltypeMat=array2table(CelltypeMat,'VariableNames',fieldName);

CelltypeMat.ExpectedNtype([83 88 89 95 100:103],1) = 1;
%% Clustering PYR and INTs
figure(44); clf;
showind=~isnan(CelltypeMat.ADP) & CelltypeMat.PeakISI>3;
cmap=hsv(2);
nexttile([1 1]);
scatter(CelltypeMat.ADP(showind),CelltypeMat.Radius(showind),20,cmap(CelltypeMat.ExpectedNtype(showind),:))
xlabel('ADP'); ylabel('Diameter (\mum)')
nexttile([1 1]);
scatter(CelltypeMat.ADP(showind),CelltypeMat.PeakISI(showind),20,cmap(CelltypeMat.ExpectedNtype(showind),:))
xlabel('ADP'); ylabel('ISI peak (ms)');
nexttile([1 1]);
scatter(CelltypeMat.PeakISI(showind),CelltypeMat.Spikewidth(showind),20,cmap(CelltypeMat.ExpectedNtype(showind),:))
xlabel('ISI peak'); ylabel('Spike width (ms)');
nexttile([1 1]);
scatter(CelltypeMat.Radius(showind),CelltypeMat.PeakISI(showind),20,cmap(CelltypeMat.ExpectedNtype(showind),:))
xlabel('Diameter (\mum)'); ylabel('ISI peak');

% K-means clustering and 3D plot
n2use=~isnan(CelltypeMat.ADP) & CelltypeMat.PeakISI>3;
n2useINT=find(CelltypeMat.ExpectedNtype(n2use)==2);
n2usePYR=find(CelltypeMat.ExpectedNtype(n2use)==1);
data = [CelltypeMat.ADP(n2use) CelltypeMat.PeakISI(n2use) CelltypeMat.Radius(n2use)];  % Replace
figure(45); clf; ax3=[];
for k = 2:3  % Adjust
[idx, C] = kmeans(data, k, 'Replicates', 10);  % More stable
gm = fitgmdist(data, k);
idxgmm = cluster(gm, data);
ax3=[ax3 nexttile([1 1])];
scatter3(data(:,1), data(:,2), data(:,3), 10, idx, 'filled'); hold all
scatter3(data(n2useINT,1), data(n2useINT,2), data(n2useINT,3), 40,'ro'); 
xlabel('ADP'); ylabel('Peak ISI'); zlabel('Diameter (\mum)');
end
linkprop([ax3], {'CameraPosition', 'CameraUpVector'});

Celltype_stas_sub=cell2mat(Celltype_STAs(n2use)');
PYRstas=Celltype_stas_sub(n2usePYR,:);
PYRstas=PYRstas-median(PYRstas(:,1:15),2);
INTstas=Celltype_stas_sub(n2useINT,:);
INTstas=INTstas-median(INTstas(:,1:15),2);

figure(46); clf; tiledlayout(1,3,'TileSpacing','tight'); ax4=[];
ax4=[ax4 nexttile([1 1])];
plot(PYRstas','color',[0.6 0.6 0.6]); hold all
plot(mean(PYRstas,1,'omitnan'),'color',[0 0 0]); axis off;
ax4=[ax4 nexttile([1 1])];
plot(INTstas','color',[0.6 0.6 0.6]); hold all
plot(mean(INTstas,1,'omitnan'),'color',[0 0 0]); axis off;
linkaxes(ax4,'xy');
nexttile([1 1]);
scatter3(data(n2usePYR,1), data(n2usePYR,2), data(n2usePYR,3), 30, [0 0 0], 'filled'); hold all;
scatter3(data(n2useINT,1), data(n2useINT,2), data(n2useINT,3), 30, [1 0 0], 'filled'); 
legend({'PYR','INT'});
xlabel('ADP'); ylabel('Peak ISI'); zlabel('Diameter (\mum)');
%% Statistics IncT 
t2integrate=[0:100]; IncTmat=[];
for Ntype=1:2
    NOIs=find(NeuronMat.BlueWvf==4 & NeuronMat.Ntype==Ntype); %IncT

    for n=NOIs'
        SpikeTr=CatTraces.Spike{n};
        m=NeuronMat.MouseID(n);
        SpClassTr=CatTraces.SpClass{n};
        blueTr=CatTraces.Blue{n};
        VoltTr=CatTraces.TraceV{n};
        SpikeHeight=mean(VoltTr(SpikeTr>0));
        VoltTr=VoltTr/SpikeHeight;
        CSTr=CatTraces.CStrace{n};
        bwCS=bwlabel(CSTr);
        bwCS(bwCS==0)=NaN;
        bwblue=bwlabel(blueTr);

        for b=1:max(bwblue)
            bframe=(bwblue==b);
            bonset=find(bframe,1);
            spikefindwindow=bonset+[0:sum(bframe)+20];
            if all(ismember(bonset+t2integrate,[1:size(VoltTr,2)])) %within the frame
                AUC=sum(VoltTr(bonset+t2integrate));
                Nspike=sum(SpikeTr(spikefindwindow));
                isCS=0;
                CSinBlue=unique(bwCS(bframe)); CSinBlue(isnan(CSinBlue))=[];
                FirstCSframe=[];
                for cs=CSinBlue
                    FirstCSframe=find(bwCS==cs,1);
                end
                CSinBlue(FirstCSframe<bonset)=[];
                if ~isempty(CSinBlue)
                    isCS=1;
                    intTime=find(ismember(bwCS,bwCS(bframe)));
                else
                    intTime=unique(find((SpikeTr.*ind2vec(size(VoltTr,2),spikefindwindow,1))>0)'+[-9:10]);
                end
                AUCcs=sum(VoltTr(intTime));
                IncTmat=[IncTmat; [n b bonset sum(bframe) AUC AUCcs Nspike isCS Ntype length(intTime) m]];
            end
        end
    end
end

fieldName={'Neuron','BlueSeq','BlueOnset','BlueWidth','AUC','AUCspike','Nspike','isCS','Celltype','AUCintTime','MouseID'};
IncTmat=array2table(IncTmat,'VariableNames',fieldName);

figure(27); clf;
Stats=[]; cmap=[1 0.5 0.5; 0.5 0.8 1]; cmap_Mean=[1 0 0; 0 0.8 1]; l=[];
Ntype_str={'PRY','INT'}; usedNeuron=[];
for Ntype=1:2
    neurons=unique(IncTmat.Neuron(IncTmat.Celltype==Ntype));
    usedNeuron{Ntype}=[];
    g=1;
    for n=neurons'
        n2show=find(IncTmat.Neuron==n & IncTmat.Nspike>=1);
        BlueWidth=IncTmat.BlueWidth(n2show);
        AUCspike=IncTmat.AUCspike(n2show);
        Nspike=IncTmat.Nspike(n2show);
        isCS=IncTmat.isCS(n2show);
        if (Ntype==1)*(mean(AUCspike(BlueWidth<=30))<45)*(mean(AUCspike(BlueWidth>=60))>27)+(Ntype==2)
            usedNeuron{Ntype}(g)=n;
            nexttile(4*(Ntype-1)+1,[1 1]);
            scatter(BlueWidth,AUCspike,15,cmap(Ntype,:),'filled','MarkerFaceAlpha',0.5); hold all
            %scatter(IncTmat.BlueWidth(n2show),IncTmat.isCS(n2show))
            %title(g);
            xlabel('Pulse width (ms)'); ylabel('AUC (A.U.)');
            
            nexttile(4*(Ntype-1)+2,[1 1]);
            scatter(BlueWidth,Nspike,15,cmap(Ntype,:),'filled','MarkerFaceAlpha',0.5); hold all
            xlabel('Pulse width (ms)'); ylabel('Spike #');

            % nexttile(4*(Ntype-1)+3,[1 1]);
            % scatter(Nspike,isCS,15,cmap(Ntype,:),'filled','MarkerFaceAlpha',0.3); hold all
            
            nexttile(4*(Ntype-1)+4,[1 1]);
            scatter(Nspike,AUCspike,15,cmap(Ntype,:),'filled','MarkerFaceAlpha',0.3); hold all
            xlabel('Spike #'); ylabel('AUC (A.U.)');

            Stats{Ntype}{g,1}=[IncTmat.BlueWidth(n2show) IncTmat.AUCspike(n2show)];
            Stats{Ntype}{g,2}=[IncTmat.BlueWidth(n2show) IncTmat.Nspike(n2show)];
            Stats{Ntype}{g,3}=[IncTmat.Nspike(n2show) IncTmat.isCS(n2show)];
            Stats{Ntype}{g,4}=[IncTmat.Nspike(n2show) IncTmat.AUCspike(n2show)];
            g=g+1;
        end
    end
    nexttile(4*(Ntype-1)+1,[1 1]);
    [M S xc N]=binning_data_median(Stats{Ntype}(:,1),[5:10:160]);
    Sall=S./sqrt(g);
    l(Ntype)=errorbar(xc,mean(M,1,'omitnan'),Sall,'color',cmap_Mean(Ntype,:),'LineWidth',1.5);
    title(Ntype_str(Ntype));
    xlim([0 105]); ylim([0 60]);

    nexttile(4*(Ntype-1)+2,[1 1]);
    [M S xc N]=binning_data(Stats{Ntype}(:,2),[5:10:160]);
    Sall=S./sqrt(g);
    l(Ntype)=errorbar(xc,mean(M,1,'omitnan'),Sall,'color',cmap_Mean(Ntype,:),'LineWidth',1.5);
    title(Ntype_str(Ntype));
    xlim([0 150]); ylim([0 10]);

    nexttile(4*(Ntype-1)+3,[1 1]);
    %[M S xc N]=binning_data(Stats{Ntype}(:,3),[0.5:1:8.5]);
    [~, CSProb, xc, N]=zscore_binning(Stats{Ntype}(:,3),[0.5:1:8.5]);
    CSProb=cell2mat(cellfun(@(x) x(:,2),CSProb,'UniformOutput',false)')';
    CSProbMean=mean(CSProb,1,'omitnan'); CSProbStd=std(CSProb,0,1,'omitnan')./sqrt(sum(~isnan(CSProb),1));
    plot(xc,CSProb','color',[0.6 0.6 0.6]); hold all
    l(Ntype)=errorbar(xc,CSProbMean,CSProbStd,'color',cmap_Mean(Ntype,:),'LineWidth',1.5);
    title(Ntype_str(Ntype));
    xlim([0 8.5]); ylim([0 1]);
    xlabel('Spike #'); ylabel('CS probability');

    nexttile(4*(Ntype-1)+4,[1 1]);
    [M S xc N]=binning_data(Stats{Ntype}(:,4),[0.5:1:8.5]);
    % [~, MeanAUC_spikebin, xc, N]=zscore_binning(Stats{Ntype}(:,4),[0.5:1:8.5]);
    % MeanAUC_spikebin=cell2mat(cellfun(@(x) x(:,2),MeanAUC_spikebin,'UniformOutput',false)')';
    % MeanSSAUC=cellfun(@(x,y) mean(x(y,2),'omitnan'),Stats{Ntype}(:,4),N(1,:)');
    % MeanSSAUC=MeanSSAUC*xc;
    % SupraRatio=MeanAUC_spikebin./MeanSSAUC;

    % DatMean=mean(MeanAUC_spikebin,1,'omitnan'); Datstd=std(MeanAUC_spikebin,0,1,'omitnan')./sqrt(sum(~isnan(MeanAUC_spikebin),1));
    % SSAUCMean=mean(MeanSSAUC,1,'omitnan'); SSAUCstd=std(MeanSSAUC,0,1,'omitnan')./sqrt(sum(~isnan(MeanSSAUC),1));
    % l(Ntype)=errorbar(xc,DatMean,Datstd,'color',cmap_Mean(Ntype,:),'LineWidth',1.5); hold all
    % l2(Ntype)=errorbar(xc,SSAUCMean,SSAUCstd,'color',[0.3 0.3 0.3],'LineWidth',1.5); hold all
    l(Ntype)=errorbar(xc,M,S./sqrt(sum(cellfun(@sum,N),2))','color',cmap_Mean(Ntype,:),'LineWidth',1.5); hold all
    lbase(Ntype)=plot(xc,M(1)*xc,'color',[0.3 0.3 0.3],'LineWidth',1.5,'LineStyle','--'); hold all
    title(Ntype_str(Ntype));
    xlim([0 8.5]); ylim([0 65]);
    % pval=[];
    % for sp=1:size(MeanAUC_spikebin,2)
    %     %[pval(sp)]=ranksum(MeanAUC_spikebin(:,sp),MeanSSAUC(:,sp));
    %     [~, pval(sp)] = ttest(SupraRatio(:,sp),1);
    % end
    % pval
end
%legend(l,{'PYR','INT'})
box off;

%% Show example IncT
NOIs=find(NeuronMat.BlueWvf==4 & NeuronMat.Ntype==1); %IncT
g=1; offset=1.5;
figure(41); clf; showHalf=[1 1 1 1 2];
tiledlayout(5,5);
for n=NOIs([13 15 16 17 18])'
    nexttile(2,[5 4])
    avgImg=CatTraces.ref_im{n};
    SpikeTr=CatTraces.Spike{n};
    SpClassTr=CatTraces.SpClass{n};
    blueTr=CatTraces.Blue{n};
    VoltTr=CatTraces.TraceV{n};
    SpikeHeight=mean(VoltTr(SpikeTr>0));
    VoltTr=VoltTr/SpikeHeight;
    CStr=VoltTr; CStr(CatTraces.CStrace{n}==0)=NaN;
    t_show=[(showHalf(g)-1)*7580+600:(showHalf(g))*7580-200];
    %t_show=[1:size(VoltTr,2)];
    plot(VoltTr(t_show)-g*offset,'k'); hold all
    plot(CStr(t_show)-g*offset,'r'); hold all

    nexttile(1+5*(g-1),[1 1])
    imshow2(avgImg,[])
    g=g+1;
end
nexttile(16,[1 1]);
drawScaleBar(10/(1.17),'horizontal')
nexttile(2,[5 4])
plot(blueTr(t_show)/2,'color',[0 0.6 1]); axis off;

NOIs=find(NeuronMat.BlueWvf==4 & NeuronMat.Ntype==2); %IncT
g=1; offset=1.5;
figure(43); clf; showHalf=[1 1 1 1 1 1];
tiledlayout(7,5);
for n=NOIs([3 6 7 8 11 12])'
    nexttile(2,[6 4])
    avgImg=CatTraces.ref_im{n};
    SpikeTr=CatTraces.Spike{n};
    SpClassTr=CatTraces.SpClass{n};
    blueTr=CatTraces.Blue{n};
    VoltTr=CatTraces.TraceV{n};
    SpikeHeight=mean(VoltTr(SpikeTr>0));
    VoltTr=VoltTr/SpikeHeight;
    CStr=VoltTr; CStr(CatTraces.CStrace{n}==0)=NaN;
    t_show=[(showHalf(g)-1)*7580+600:(showHalf(g))*7580-200];
    %t_show=[1:size(VoltTr,2)];
    plot(VoltTr(t_show)-g*offset,'k'); hold all
    plot(CStr(t_show)-g*offset,'r'); hold all

    nexttile(1+5*(g-1),[1 1])
    imshow2(avgImg,[])
    g=g+1;
end
nexttile(2,[6 4])
plot(blueTr(t_show)/2,'color',[0 0.6 1]); axis off;

%% Show example RP + Ramp
NOIs=find(NeuronMat.BlueWvf==3 & NeuronMat.Ntype==1); %IncT
g=1; offset=1.5;
figure(42); clf; showHalf=[1 1 1 1];
tiledlayout(5,5);
for n=NOIs([2 15 18 25])'
    nexttile(2,[5 4])
    avgImg=CatTraces.ref_im{n};
    SpikeTr=CatTraces.Spike{n};
    SpClassTr=CatTraces.SpClass{n};
    blueTr=CatTraces.Blue{n};
    VoltTr=CatTraces.TraceV{n};
    SpikeHeight=mean(VoltTr(SpikeTr>0));
    VoltTr=VoltTr/SpikeHeight;
    CStr=VoltTr; CStr(CatTraces.CStrace{n}==0)=NaN;
    t_show=[3300:length(VoltTr)-1200];
    plot(VoltTr(t_show)+g*offset,'color',[0.4 0.4 0.4]); hold all
    plot(CStr(t_show)+g*offset,'r'); hold all

    nexttile(1+5*(g-1),[1 1])
    imshow2(avgImg,[])
    g=g+1;
end
nexttile(16,[1 1]);
drawScaleBar(10/(6.5/25),'horizontal')
nexttile(2,[5 4])
plot(blueTr(t_show)/2,'color',[0 0.6 1]); axis off;

%% Show example Ramp
NOIs=find(NeuronMat.BlueWvf==3 & NeuronMat.Ntype==1); %IncT
g=1; offset=1.5;
figure(42); clf; showHalf=[1 1 1 1];
tiledlayout(5,5);
for n=NOIs([2 15 18 25])'
    nexttile(2,[5 4])
    avgImg=CatTraces.ref_im{n};
    SpikeTr=CatTraces.Spike{n};
    SpClassTr=CatTraces.SpClass{n};
    blueTr=CatTraces.Blue{n};
    VoltTr=CatTraces.TraceV{n};
    SpikeHeight=mean(VoltTr(SpikeTr>0));
    VoltTr=VoltTr/SpikeHeight;
    CStr=VoltTr; CStr(CatTraces.CStrace{n}==0)=NaN;
    t_show=[8500:length(VoltTr)-1000];
    plot(VoltTr(t_show)+g*offset,'k'); hold all
    plot(CStr(t_show)+g*offset,'r'); hold all

    nexttile(1+5*(g-1),[1 1])
    imshow2(avgImg,[])
    g=g+1;
end
nexttile(2,[5 4])
plot(blueTr(t_show)/2,'color',[0 0.6 1]); axis off;

%% Show example tri-Ramp
NOIs=find(NeuronMat.BlueWvf==1 & NeuronMat.Ntype==1 & NeuronMat.IsKet==0); %IncT
g=1; offset=1.5;
figure(42); clf; showHalf=[1 1 1 2 1];
%tiledlayout(5,5);
for n=NOIs([4 6 10 11 15])'
    nexttile(2,[5 4])
    avgImg=CatTraces.ref_im{n};
    SpikeTr=CatTraces.Spike{n};
    SpClassTr=CatTraces.SpClass{n};
    blueTr=CatTraces.Blue{n};
    VoltTr=CatTraces.TraceV{n};
    SpikeHeight=mean(VoltTr(SpikeTr>0));
    VoltTr=VoltTr/SpikeHeight;
    CStr=VoltTr; CStr(CatTraces.CStrace{n}==0)=NaN;
    t_show=[(showHalf(g)-1)*6300+50:(showHalf(g))*6300+50];
    %t_show=[1:length(VoltTr)];
    plot(VoltTr(t_show)+g*offset,'k'); hold all
    plot(CStr(t_show)+g*offset,'r'); hold all

    % nexttile(1+5*(g-1),[1 1])
    % imshow2(avgImg,[])
    g=g+1;
end
nexttile(2,[5 4])
plot(blueTr(t_show)/2,'color',[0 0.6 1]); axis off;

%% Statistics IncT
NOIs=find(NeuronMat.BlueWvf==1 & NeuronMat.Ntype==1 & NeuronMat.IsKet==0); %IncT
CS_mat=[]; SS_mat=[];
Rheobin=[0 1.5 3 5 7 9]; timebin=[0:1000:6000];
for n=NOIs'
    SpikeTr=CatTraces.Spike{n};
    bwCS=bwlabel(CatTraces.CStrace{n});
    CSspikeTr=CatTraces.Spike{n}.*bwCS;
    SpClassTr=CatTraces.SpClass{n};
    blueTr=CatTraces.Blue{n};
    VoltTr=CatTraces.TraceV{n};
    SpikeHeight=mean(VoltTr(SpikeTr>0));
    bwblue=bwlabel(blueTr);
    BlueSpike=blueTr(find(SpikeTr)); Rheoconst=mean(BlueSpike(1:5));
    RheoTr=blueTr/Rheoconst;

    for b=1:max(bwblue)
        toi=(bwblue==b);
       sptime=SpikeTr(toi);
       RheoTr_sub=RheoTr(toi);
       bwCS_sub=bwCS(toi);
       %CStime=find(SpClassTr(2,toi));
       CStime=find(CSspikeTr(toi)>0);
       SStime=find(SpClassTr(1,toi));
       Ncs=[];
       for cs=1:length(CStime)
       Ncs(cs)=sum(CSspikeTr==bwCS_sub(CStime(cs)));
       end
       CS_mat=[CS_mat; [repmat([n, b],length(CStime),1), RheoTr_sub(CStime)' CStime']];
       SS_mat=[SS_mat; [repmat([n, b],length(SStime),1), RheoTr_sub(SStime)' SStime']];
    end
end
RatioBin=[]; g=1; CSrheobin_count=[]; SSrheobin_count=[]; CStimebin_count=[];
for n=unique([CS_mat(:,1); SS_mat(:,1)])'
    % CSrheobin_count=histcounts2(CS_mat(CS_mat(:,1)==n,4),CS_mat(CS_mat(:,1)==n,3),Rheobin);
    % SSrheobin_count=histcounts2(SS_mat(SS_mat(:,1)==n,4),SS_mat(SS_mat(:,1)==n,3),Rheobin);
    % RatioBin(:,:,g)=CSrheobin_count./(CSrheobin_count+SSrheobin_count);
    CSrheobin_count(g,:,1)=histcounts(CS_mat(CS_mat(:,1)==n & CS_mat(:,4)<7000,3),Rheobin);
    SSrheobin_count(g,:,1)=histcounts(SS_mat(SS_mat(:,1)==n & SS_mat(:,4)<7000,3),Rheobin);
    %CSrheobin_count(g,:,2)=histcounts(CS_mat(CS_mat(:,1)==n & CS_mat(:,4)>=3000,3),Rheobin);
    %SSrheobin_count(g,:,2)=histcounts(SS_mat(SS_mat(:,1)==n & SS_mat(:,4)>=3000,3),Rheobin);
    CStimebin_count(g,:,1)=histcounts(CS_mat(CS_mat(:,1)==n,4),timebin);
    SStimebin_count(g,:,1)=histcounts(SS_mat(SS_mat(:,1)==n,4),timebin);
    g=g+1;
end
RatioBin=CSrheobin_count./(CSrheobin_count+SSrheobin_count);
RatioBin(RatioBin==inf)=NaN;
Rheobin_c=mean([Rheobin(2:end); Rheobin(1:end-1)]);
timebin_c=mean([timebin(2:end); timebin(1:end-1)])/1000;

figure(32); clf; %tiledlayout(1,5);
%nexttile(4,[1 2]);
M=mean(RatioBin,1,'omitnan'); S=std(RatioBin,0,1,'omitnan'); N=sum(~isnan(RatioBin),1);
plot(Rheobin_c,RatioBin,'color',[0.7 0.7 0.7]); hold all;
errorbar(Rheobin_c,M(:,:,1),S(:,:,1)./sqrt(N(:,:,1)),'color',[0 0 0],'LineWidth',1.5);
%errorbar(Rheobin_c,M(:,:,2),S(:,:,2)./sqrt(N(:,:,2)),'color',[0.2 0.8 0.3],'LineWidth',1.5);
xlim([0 9]); ylim([0 1]);
xlabel('Optical rheobase'); ylabel('Fraction of complex spike'); box off;
% 
% % nexttile(1,[1 3]); l=[];
% MC=mean(CStimebin_count,1,'omitnan'); SC=std(CStimebin_count,0,1,'omitnan'); N=size(SStimebin_count,1);
% MS=mean(SStimebin_count,1,'omitnan'); SS=std(SStimebin_count,0,1,'omitnan');
% MR=mean(CStimebin_count./(CStimebin_count+SStimebin_count),1,'omitnan'); SR=std(CStimebin_count./(CStimebin_count+SStimebin_count),0,1,'omitnan');
% l(1)=errorbar(timebin_c,MC,SC./sqrt(N),'color',[1 0 0],'LineWidth',1.5); hold all
% ylabel('Spike count');
% l(2)=errorbar(timebin_c,MS,SS./sqrt(N),'color',[0 0 0],'LineWidth',1.5); yyaxis right;
% l(3)=errorbar(timebin_c,MR,SR./sqrt(N),'color',[0.93 0.69 0.13],'LineWidth',1.5);
% xlabel('Time (s)'); ylabel('Fraction of complex spike'); box off;
% legend(l,{'CS','SS'})
% clf;
% scatter(CS_mat(:,4),CS_mat(:,3),20,[1 0 0],'filled')
% hold all
% scatter(SS_mat(:,4),SS_mat(:,3),20,[0 0 0],'filled')
% pfsz = size(LapMat_sub);
% pfcentroid = (LapMat_sub * [1:pfsz(2)]') ./ sum(LapMat_sub,2);
% pfsize = sqrt(sum(LapMat_sub .* (repmat([1:pfsz(2)],pfsz(1),1) - pfcentroid).^2,2)./ sum(LapMat_sub,2));

%% Statistics Ramp
omitN = [57 78 79 81 82 90 94 96 97 98 99 104 105 106 107 110 111];
NOIs=find(NeuronMat.BlueWvf==3 & NeuronMat.Ntype==1 & NeuronMat.IsKet==0); %RP
NOIs=setdiff(NOIs,omitN);
SomStimSpikeMat=[];
figure(333); clf;
for n=NOIs'
    SpikeTr=CatTraces.Spike{n};
    bwCS=bwlabel(CatTraces.CStrace{n});
    CSspikeTr=CatTraces.Spike{n}.*bwCS;
    SpClassTr=CatTraces.SpClass{n};
    blueTr=CatTraces.Blue{n};
    VoltTr=CatTraces.TraceV{n};
    SpikeHeight=mean(VoltTr(SpikeTr>0));
    bwblue=bwlabel(blueTr>0);
    BlueSpike=blueTr(find(SpikeTr.*(bwblue==max(bwblue)))); 
    %BlueSpike=blueTr(find(SpikeTr.*(bwblue>0))); 
    Rheoconst=mean(BlueSpike(1));
    RheoTr=blueTr/Rheoconst;

    normtr=rescale(VoltTr);
    plot(normtr'+n); hold all
    plot(find(bwCS>0), normtr(:,bwCS>0)'+n,'r');

    for b=1:max(bwblue)
        toi=(bwblue==b);
       sptime=SpikeTr(toi);
       RheoTr_sub=RheoTr(toi);
       bwCS_sub=bwCS(toi);
       %CStime=find(SpClassTr(2,toi));
       CStime=find(CSspikeTr(toi)>0);
       SStime=find(SpClassTr(1,toi));
       Ncs=[];
       for cs=1:length(CStime)
       Ncs(cs)=sum(CSspikeTr==bwCS_sub(CStime(cs)));
       end
       SomStimSpikeMat=[SomStimSpikeMat; [repmat([n, b-max(bwblue)],length(CStime),1), RheoTr_sub(CStime)' CStime' ones(length(CStime),1)]];
       SomStimSpikeMat=[SomStimSpikeMat; [repmat([n, b-max(bwblue)],length(SStime),1), RheoTr_sub(SStime)' SStime' zeros(length(SStime),1)]];
    end
end

fieldName={'NeuronID','StimOrder','RheoBase','SpikeTime','IsCS'};
SomStimSpikeMat=array2table(SomStimSpikeMat,'VariableNames',fieldName);

Rheobin=[0.7:0.35:4.5]; timebin=[0:1000:6000]; minSpikeN=20;
Rheobin2=[0 1.4 9]; 
RatioBin=[]; RatioBin2=[]; g=1; 
CSrheobin_count=NaN(length(unique(SomStimSpikeMat.NeuronID)),length(Rheobin)-1); CStimebin_count=[];
SSrheobin_count=NaN(length(unique(SomStimSpikeMat.NeuronID)),length(Rheobin)-1); 
CStimebin_count=NaN(length(unique(SomStimSpikeMat.NeuronID)),length(timebin)-1);
SStimebin_count=NaN(length(unique(SomStimSpikeMat.NeuronID)),length(timebin)-1);
CSrheobin_count2=NaN(length(unique(SomStimSpikeMat.NeuronID)),length(Rheobin2)-1); 
SSrheobin_count2=NaN(length(unique(SomStimSpikeMat.NeuronID)),length(Rheobin2)-1); 
for n=unique(SomStimSpikeMat.NeuronID)'
    allspN=SomStimSpikeMat.NeuronID==n & SomStimSpikeMat.StimOrder==0;
    countindCS=SomStimSpikeMat.NeuronID==n & SomStimSpikeMat.StimOrder==0 & SomStimSpikeMat.IsCS==1; % Only Ramp, CS
    countindSS=SomStimSpikeMat.NeuronID==n & SomStimSpikeMat.StimOrder==0 & SomStimSpikeMat.IsCS==0; % Only Ramp, SS
    if ~isempty(countindCS) & ~isempty(countindSS) & sum(allspN)>=minSpikeN
    CSrheobin_count(g,:,1)=histcounts(SomStimSpikeMat.RheoBase(countindCS),Rheobin);
    CStimebin_count(g,:,1)=histcounts(SomStimSpikeMat.SpikeTime(countindCS),timebin);
    SSrheobin_count(g,:,1)=histcounts(SomStimSpikeMat.RheoBase(countindSS),Rheobin);
    SStimebin_count(g,:,1)=histcounts(SomStimSpikeMat.SpikeTime(countindSS),timebin);

    CSrheobin_count2(g,:,1)=histcounts(SomStimSpikeMat.RheoBase(countindCS),Rheobin2);
    SSrheobin_count2(g,:,1)=histcounts(SomStimSpikeMat.RheoBase(countindSS),Rheobin2);
    g=g+1;
    end
end
RatioBin=CSrheobin_count./(CSrheobin_count+SSrheobin_count);
RatioBin(RatioBin==inf)=NaN;
Rheobin_c=mean([Rheobin(2:end); Rheobin(1:end-1)]);
timebin_c=mean([timebin(2:end); timebin(1:end-1)])/1000;

RatioBin2=CSrheobin_count2./(CSrheobin_count2+SSrheobin_count2);
RatioBin2(RatioBin2==inf)=NaN;

MeanCSweakvsstrong_2ndhalf=mean(RatioBin2,1,'omitnan');
StdCSweakvsstrong_2ndhalf=std(RatioBin2,0,1,'omitnan')./sqrt(sum(~isnan(RatioBin2),1));
fprintf(['Rheobase <1.4 : %2.2f ± %2.2f, > 1.4: %2.2f ± %2.2f  \n'],...
    MeanCSweakvsstrong_2ndhalf(1),StdCSweakvsstrong_2ndhalf(1),MeanCSweakvsstrong_2ndhalf(2),StdCSweakvsstrong_2ndhalf(2));

figure(32); clf; %tiledlayout(1,5);
%nexttile(4,[1 2]);
ShowIndCS=SomStimSpikeMat.StimOrder==0 & SomStimSpikeMat.IsCS==1;
ShowIndSS=SomStimSpikeMat.StimOrder==0 & SomStimSpikeMat.IsCS==0;
h(1)=histogram([SomStimSpikeMat.RheoBase(ShowIndSS);SomStimSpikeMat.RheoBase(ShowIndCS)],Rheobin,'FaceColor',[0.4 0.4 0.4]); hold all
h(2)=histogram(SomStimSpikeMat.RheoBase(ShowIndCS),Rheobin,'FaceColor',[1 0 0]); 
ylabel('Spike Counts');
yyaxis right;
M=mean(RatioBin,1,'omitnan'); S=std(RatioBin,0,1,'omitnan'); N=sum(~isnan(RatioBin),1);
errorbar(Rheobin_c,M(:,:,1),S(:,:,1)./sqrt(N(:,:,1)),'color',[0 0 1],'LineWidth',1.5);
xlim([0.6 4.5]); ylim([0 0.9]);
xlabel('Optical rheobase'); ylabel('Fraction of complex spike'); box off;
legend(h,{'Total spikes','Complex spikes'});
%% Statistics Pulse
NOIs=find(NeuronMat.BlueWvf==3 & NeuronMat.Ntype==1 & NeuronMat.IsKet==0); %RP
NOIs=setdiff(NOIs,omitN);
RheoThreshold=1.45;

Rheobin=[1:0.25:9.5]; timebin=[0:60:120 260:120:500]; minSpikeN=1;
timebin2=[250 500];
Rheobin2=[0 1.5 2 3 4 5 7 9]; 
RatioBin=[]; g=1; 
CStimebin_count=NaN(length(unique(SomStimSpikeMat.NeuronID)),length(timebin)-1,2);
SStimebin_count=NaN(length(unique(SomStimSpikeMat.NeuronID)),length(timebin)-1,2);
CStimebin_count2=NaN(length(unique(SomStimSpikeMat.NeuronID)),length(timebin2)-1,2);
SStimebin_count2=NaN(length(unique(SomStimSpikeMat.NeuronID)),length(timebin2)-1,2);
for n=unique(SomStimSpikeMat.NeuronID)'
    allspN=SomStimSpikeMat.NeuronID==n & ismember(SomStimSpikeMat.StimOrder,-[1:5]);
    countindCSweak=SomStimSpikeMat.NeuronID==n & ismember(SomStimSpikeMat.StimOrder,-[1:5]) & SomStimSpikeMat.IsCS==1 & SomStimSpikeMat.RheoBase<=RheoThreshold; % Only Ramp, CS
    countindSSweak=SomStimSpikeMat.NeuronID==n & ismember(SomStimSpikeMat.StimOrder,-[1:5]) & SomStimSpikeMat.IsCS==0 & SomStimSpikeMat.RheoBase<=RheoThreshold; % Only Ramp, SS
    countindCSstrong=SomStimSpikeMat.NeuronID==n & ismember(SomStimSpikeMat.StimOrder,-[1:5]) & SomStimSpikeMat.IsCS==1 & SomStimSpikeMat.RheoBase>RheoThreshold; % Only Ramp, CS
    countindSSstrong=SomStimSpikeMat.NeuronID==n & ismember(SomStimSpikeMat.StimOrder,-[1:5]) & SomStimSpikeMat.IsCS==0 & SomStimSpikeMat.RheoBase>RheoThreshold; % Only Ramp, SS

    if sum(allspN)>=minSpikeN
    CStimebin_count(g,:,1)=histcounts(SomStimSpikeMat.SpikeTime(countindCSweak),timebin);
    SStimebin_count(g,:,1)=histcounts(SomStimSpikeMat.SpikeTime(countindSSweak),timebin);
    CStimebin_count(g,:,2)=histcounts(SomStimSpikeMat.SpikeTime(countindCSstrong),timebin);
    SStimebin_count(g,:,2)=histcounts(SomStimSpikeMat.SpikeTime(countindSSstrong),timebin);

    CStimebin_count2(g,:,1)=histcounts(SomStimSpikeMat.SpikeTime(countindCSweak),timebin2);
    SStimebin_count2(g,:,1)=histcounts(SomStimSpikeMat.SpikeTime(countindSSweak),timebin2);
    CStimebin_count2(g,:,2)=histcounts(SomStimSpikeMat.SpikeTime(countindCSstrong),timebin2);
    SStimebin_count2(g,:,2)=histcounts(SomStimSpikeMat.SpikeTime(countindSSstrong),timebin2);
    g=g+1;
    end
end

RatioBin=CStimebin_count./(CStimebin_count+SStimebin_count);
RatioBin(RatioBin==inf)=NaN;
timebin_c=mean([timebin(2:end); timebin(1:end-1)]);
M=mean(RatioBin,1,'omitnan'); S=std(RatioBin,0,1,'omitnan'); N=sum(~isnan(RatioBin),1); p=[];
for b=1:size(RatioBin,2)
%p_tmp=get_pValue({RatioBin(:,b,1),RatioBin(:,b,2)},0);
p_tmp=get_pValue(squeeze(RatioBin(:,b,:)),1);
p(b)=p_tmp(1,2);
end
p

RatioBin2=CStimebin_count2./(CStimebin_count2+SStimebin_count2);
RatioBin2(RatioBin2==inf)=NaN;
Ratestmp=squeeze(RatioBin2(:,:,:));
MeanCSweakvsstrong_2ndhalf=mean(Ratestmp,1,'omitnan');
StdCSweakvsstrong_2ndhalf=std(Ratestmp,0,1,'omitnan')./sqrt(sum(~isnan(Ratestmp),1));
fprintf(['CS fraction from 2nd half of stimulation, Weak: %2.2f ± %2.2f, Strong: %2.2f ± %2.2f  \n'],...
    MeanCSweakvsstrong_2ndhalf(1),StdCSweakvsstrong_2ndhalf(1),MeanCSweakvsstrong_2ndhalf(2),StdCSweakvsstrong_2ndhalf(2));


figure(33); clf; %tiledlayout(1,5);
%nexttile(4,[1 2]);
nexttile([1 1]);
% ShowIndCS=ismember(SomStimSpikeMat.StimOrder,-[1:5]) & SomStimSpikeMat.IsCS==1 & SomStimSpikeMat.RheoBase<=RheoThreshold;
% ShowIndSS=ismember(SomStimSpikeMat.StimOrder,-[1:5]) & SomStimSpikeMat.IsCS==0 & SomStimSpikeMat.RheoBase<=RheoThreshold;
% histogram([SomStimSpikeMat.SpikeTime(ShowIndSS);SomStimSpikeMat.SpikeTime(ShowIndCS)],timebin,'FaceColor',[0.4 0.4 0.4]); hold all
% histogram(SomStimSpikeMat.SpikeTime(ShowIndCS),timebin,'FaceColor',[1 0 0]); yyaxis right;
errorbar(timebin_c,M(:,:,1),S(:,:,1)./sqrt(N(:,:,1)),'color',[0.5 0.5 0.5],'LineWidth',1.5); hold all
% nexttile([1 1]);
% ShowIndCS=ismember(SomStimSpikeMat.StimOrder,-[1:5]) & SomStimSpikeMat.IsCS==1 & SomStimSpikeMat.RheoBase>RheoThreshold;
% ShowIndSS=ismember(SomStimSpikeMat.StimOrder,-[1:5]) & SomStimSpikeMat.IsCS==0 & SomStimSpikeMat.RheoBase>RheoThreshold;
% histogram([SomStimSpikeMat.SpikeTime(ShowIndSS);SomStimSpikeMat.SpikeTime(ShowIndCS)],timebin,'FaceColor',[0.4 0.4 0.4]); hold all
% histogram(SomStimSpikeMat.SpikeTime(ShowIndCS),timebin,'FaceColor',[1 0 0]); yyaxis right;
errorbar(timebin_c,M(:,:,2),S(:,:,2)./sqrt(N(:,:,2)),'color',[0 0 0],'LineWidth',1.5);
xlim([0 500]); ylim([0 1]);
xlabel('Time after stimulation onset (ms)'); ylabel('Fraction of complex spike'); box off;
legend({'Weak stimulation','Strong stimulation'})

figure(34); clf;
[~, ii]=unique([SomStimSpikeMat.NeuronID SomStimSpikeMat.StimOrder],'row');
stim_Rheo_FR_CSratio=[];
for stimInd=ii'
    if ismember(SomStimSpikeMat.StimOrder(stimInd),-[1:5])
    Ind2see=find([SomStimSpikeMat.NeuronID==SomStimSpikeMat.NeuronID(stimInd) & SomStimSpikeMat.StimOrder==SomStimSpikeMat.StimOrder(stimInd)]);
    stim_Rheo_FR_CSratio=[stim_Rheo_FR_CSratio; [SomStimSpikeMat.RheoBase(Ind2see(1)) length(Ind2see)*2 sum(SomStimSpikeMat.IsCS(Ind2see))]];
    end
end

%% Ketamien administration
NOIs=find(NeuronMat.IsKet==1); 
figure(55); clf; g=1; scale=10;
for n=NOIs([1 6])'
   nexttile([1 1])
   imshow2(CatTraces.ref_im{n},[]); drawScaleBar(10/0.26,'horizontal','color',[1 1 1])
   NeuronIDs=find(NeuronMat.MouseID==NeuronMat.MouseID(n) & NeuronMat.NeuronID==NeuronMat.NeuronID(n));
   [~, nuniqe]=unique([NeuronMat.BlueWvf(NeuronIDs) NeuronMat.IsKet(NeuronIDs)],'row');
   nexttile([1 1])
   g2=1;
   for n2=NeuronIDs(nuniqe)'
       if NeuronMat.BlueWvf(n2)==3
           time2show=[3000:13500];
       else
           time2show=[1:size(CatTraces.TraceV{n2},2)];
       end
   TrwCS=CatTraces.TraceV{n2};
   TrwCS(CatTraces.CStrace{n2}==0)=NaN;
   plot(CatTraces.TraceV{n2}(time2show)+scale*g2,'k'); hold all
   plot(TrwCS(time2show)+scale*g2,'r');
   plot(rescale(CatTraces.Blue{n2}(time2show))*2+scale*(g2-1/5),'color',[0 0.6 1]); axis off;
   drawScaleBar(1000,'horizontal');
   g2=g2+1;
   end
   g=g+1;
end
CSfraction=[]; g=1;
for n=NOIs'
   
   NeuronIDs=find(NeuronMat.MouseID==NeuronMat.MouseID(n) & NeuronMat.NeuronID==NeuronMat.NeuronID(n));
   experimentlist=[NeuronMat.BlueWvf(NeuronIDs) NeuronMat.IsKet(NeuronIDs)];
   [experimentlistuniq, nuniqe]=unique(experimentlist,'row');
   KetamineNeuronID=NeuronIDs(nuniqe(experimentlistuniq(:,2)>0));
   IsoNeuronID=NeuronIDs(nuniqe(experimentlistuniq(:,2)==0));
   CSspcat_ket=[]; CSspcat_iso=[];
   for k=KetamineNeuronID'
       CSspcat_ket=[CSspcat_ket [CatTraces.Spike{k}.*CatTraces.CStrace{k}; CatTraces.Spike{k}]];
   end
   for k=IsoNeuronID'
       CSspcat_iso=[CSspcat_iso [CatTraces.Spike{k}.*CatTraces.CStrace{k}; CatTraces.Spike{k}]];
   end
   CSfraction(g,:)=[sum(CSspcat_ket(1,:),2)/sum(CSspcat_ket(2,:),2) sum(CSspcat_iso(1,:),2)/sum(CSspcat_iso(2,:),2)];
   g=g+1;
end
figure(56); clf;
p=Boxplot_wPoints2(CSfraction([1 3 4 5 6],[2 1]),[0 0 0; 0.5 0.5 0.5]);
drawPValueLines(p,0,'TextYOffset',0.05); box off;
ylabel('Fraction of complex spike')
set(gca,'xtick',[1 2],'XTickLabel',{'Iso','Iso/ketamine'})
set_fontsize(13);

%% Classfy complex spike (Figure S)

ExampleTrace=CatTraces.TraceV{11};
sp_soma=CatTraces.Spike{11};
CStrace=CatTraces.CStrace{11};
STA=get_STA(ExampleTrace,sp_soma,1,1);
spikeheight=STA(2);
SS_frm=4889+[-30:60];
BS_frm=4691+[-40:150];
CS_frm=7585+[-40:150];

% Default parameters
CS_thres = [4 1];
N_Spike = 3;
N_Spike2ISI = 3;
Max_AUC = 150;

% Process subthreshold trace
tr_hi=ExampleTrace-movprc(ExampleTrace,400,30,2);

tr_sub= get_subthreshold(tr_hi,sp_soma,7,17);
CStraceTr=tr_hi;
CStraceTr(~CStrace)=NaN;
nTime = length(tr_sub);

% Create figure and main plot
f = figure(30); clf;
nexttile([1 1])
plot(ExampleTrace,'color',[0 0 0]); hold all;
drawScaleBar(1000,'horizontal','color',[0 0 0],'Linewidth',3,'Position',[14500 -2]);
drawScaleBar(4,'vertical','color',[0 0 0],'Linewidth',3,'Position',[14500 -2]);
axis tight off;

nexttile([1 1]);
tr2show=tr_hi(SS_frm);
plot(tr2show, 'k','linewidth',1); hold on; grid on;
plot(tr_sub(SS_frm),'color',[1 0.5 0],'linewidth',1.5);
plot(find(sp_soma(SS_frm)>0), tr2show(find(sp_soma(SS_frm)>0)),'ro','linewidth',1);
box off;
ylabel('Voltage (Z score)'); xlabel('Time (ms)');

nexttile([1 1]);
tr2show=tr_hi(BS_frm);
plot(tr2show, 'k','linewidth',1); hold on; grid on;
plot(tr_sub(BS_frm),'color',[1 0.5 0],'linewidth',1.5);
plot(find(sp_soma(BS_frm)>0), tr2show(find(sp_soma(BS_frm)>0)),'ro','linewidth',1);
box off;
ylabel('Voltage (Z score)'); xlabel('Time (ms)');

nexttile([1 1]);
tr2show=tr_hi(CS_frm);
plot(tr2show, 'k','linewidth',1); hold on; grid on;
plot(tr_sub(CS_frm),'color',[1 0.5 0],'linewidth',1.5);
plot(find(sp_soma(CS_frm)>0), tr2show(find(sp_soma(CS_frm)>0)),'ro','linewidth',1);
box off;
ylabel('Voltage (Z score)'); xlabel('Time (ms)');
set_fontsize(13);

%% CS amplitude distritubion
CSamp_ints=[];
for n=1:157
tr_hi=CatTraces.TraceV{n}-movprc(CatTraces.TraceV{n},300,30,2);
tr_sub = movmean(tr_hi, 20, 2);
tr_sub = tr_sub - movmedian(tr_sub, 300, 2);
[trans, tr_trace] = detect_transient2(tr_sub, [1.5 1], CatTraces.Spike{n}, 10);
[~, csfrstfrm]=unique(bwlabel(CatTraces.CStrace{n})); csfrstfrm=csfrstfrm(2:end);
cscandints=tr_trace(csfrstfrm);
STA=get_STA(CatTraces.TraceV{n},CatTraces.Spike{n},1,1);
spikeheight=STA(2);
CSamp_ints{n}=[trans.amp(cscandints(cscandints>0))' trans.int(cscandints(cscandints>0))' repmat(spikeheight,length(cscandints(cscandints>0)),1)];
end
nonemptycell=cell2mat(cellfun(@(x) ~isempty(x),CSamp_ints,'UniformOutput',false));
minCSamps=cellfun(@(x) [prctile(x(:,1),3) prctile(x(:,2),3) prctile(x(:,1),3)/x(1,3) prctile(x(:,2),3)/x(1,3)],CSamp_ints(nonemptycell),'UniformOutput',false);
minCSamps=cell2mat(minCSamps');
%%
g=1; offset=1.5;
figure(52); clf; showHalf=[1 1 1 1];
%tiledlayout(5,5);
g=1;
for n=NOIs'
    nexttile(2,[5 4])
    avgImg=CatTraces.ref_im{n};
    SpikeTr=CatTraces.Spike{n};
    SpClassTr=CatTraces.SpClass{n};
    blueTr=CatTraces.Blue{n};
    VoltTr=CatTraces.TraceV{n};
    SpikeHeight=mean(VoltTr(SpikeTr>0));
    VoltTr=VoltTr/SpikeHeight;
    CStr=VoltTr; CStr(CatTraces.CStrace{n}==0)=NaN;
    t_show=[3300:length(VoltTr)-1200];
    plot(VoltTr(t_show)+g*offset,'color',[0.4 0.4 0.4]); hold all
    plot(CStr(t_show)+g*offset,'r'); hold all

    nexttile(1+5*(g-1),[1 1])
    imshow2(avgImg,[])
    g=g+1;
end
nexttile(16,[1 1]);
drawScaleBar(10/(6.5/25),'horizontal')
nexttile(2,[5 4])
plot(blueTr(t_show)/2,'color',[0 0.6 1]); axis off;



