clear
clc;
cd '/Volumes/BHL18TB_D2/YQ601_PlaceCellResults';
[~, fpath, raw] = xlsread(['/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/' ...
    'PlaceCellData_Arrangement.xlsx'], 'Sheet1', 'C8:N60');

save_figto='/Volumes/BHL18TB_D2/YQ601_PlaceCellResults';
place_bin=150;
velocity_threshold=0.002;
VRsession=cell2mat(raw(:,12));
MouseInd=cell2mat(raw(:,2));
NeuroInd=cell2mat(raw(:,3));
time_size=15000; %segment size
overlap=200;
set(0,'DefaultFigureWindowStyle','docked')
foi=setdiff([1:length(fpath)],[36]);

[~, ~, PFs] = xlsread(['/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/' ...
    'SomV_PlaceField_Arrangement.xlsx'], 'Sheet1', 'A5:G52');
PFs=cell2mat(PFs);
PFfieldName={'FileInd','Nind','LapStart','LapEnd','World','PFstart','PFend'};
PFs=array2table(PFs,'VariableNames',PFfieldName);

%% Load data
SomPCresult=[];
for f=foi
    f
    load([fpath{f} '/SomPC_Result.mat'])
    SomPCresult{f}=Result;
end

%% Show all PFs
figure(100); clf;
PF_list=[];
for pf=1:size(PFs,1)
    f=PFs.FileInd(pf);
    n=PFs.Nind(pf);
    [~, Lap_firstFrame]=unique(SomPCresult{f}.VR(8,:));
    try
        Lap_world=SomPCresult{f}.VR(2,Lap_firstFrame+500);
    catch
        Lap_world=SomPCresult{f}.VR(2,Lap_firstFrame(1:end-1)+500);
        Lap_world(end)=SomPCresult{f}.VR(2,Lap_firstFrame(end)+60);
    end
    FRmap_world2=NaN(size(SomPCresult{f}.Lap_FR(:,:,1,n)));
    FRmap_world3=NaN(size(SomPCresult{f}.Lap_FR(:,:,1,n)));
    FRmap=ringmovMean(SomPCresult{f}.Lap_FR(:,:,1,n),5);
    % FRmap_world2(Lap_world==2,:)=FRmap(Lap_world==2,:);
    % FRmap_world3(Lap_world~=2,:)=FRmap(Lap_world~=2,:);
    % nexttile([1 1])
    % imagesc(FRmap);
    % nexttile([1 1])
    % imagesc(FRmap_world3);
    CSSSoverlap=cat(3,mat2gray(SomPCresult{f}.Lap_CS(2:end-1,:,1,n)),mat2gray(SomPCresult{f}.Lap_SS(2:end-1,:,1,n)));
    CSSSoverlap=ringmovMean(CSSSoverlap,5)*3;
    CSSS_world2=NaN([size(SomPCresult{f}.Lap_FR(:,:,1,n)) 2]);
    CSSS_world3=NaN([size(SomPCresult{f}.Lap_FR(:,:,1,n)) 2]);
    CSSS_world2(Lap_world(2:end-1)==2,:,:)=CSSSoverlap(Lap_world(2:end-1)==2,:,:);
    CSSS_world3(Lap_world(2:end-1)~=2,:,:)=CSSSoverlap(Lap_world(2:end-1)~=2,:,:);
    nexttile([1 1])
    imagesc(im_merge(CSSSoverlap,[1 0 0; 1 1 1]));
    title(['Neuron ID#', num2str(pf)])
    % nexttile([1 1])
    % imagesc(im_merge(CSSS_world3,[1 0 0; 1 1 1]));
    % title(['Neuron ID#', num2str(pf) ', World 3'])
    st=1;
    % while ~isempty(st)
    %     st=input(['Neuron : ' num2str(n) ',Start lap?\n']);
    %     if isempty(st)
    %         break;
    %     end
    %     ed=input(['Neuron : ' num2str(n) ',End lap? \n']);
    %     wd=input(['Neuron : ' num2str(n) ',World? \n']);
    %     PF_list=[PF_list; [f n [[st ed]+1] wd]];
    % end
end

%% Representative image

N2show=table2array(PFs(10,:));
if isfield(SomPCresult{N2show(1)},'Blue')
    blueLaps=unique(double(SomPCresult{N2show(1)}.Blue>0).*SomPCresult{N2show(1)}.VR(8,:));
    blueLaps(blueLaps==0)=[];
else
    blueLaps=[];
end
Subth=get_subthreshold(SomPCresult{N2show(1)}.normtrace(N2show(2),:),SomPCresult{N2show(1)}.spike(N2show(2),:),7,40);
normConst=get_F0PCA(Subth);
%normConst=sqrt(mean(Subth(1,1:end-1).*Subth(1,2:end)));
%normTr=zscore(SomPCresult{N2show(1)}.normtrace(N2show(2),:));
normTr=SomPCresult{N2show(1)}.normtrace(N2show(2),:)./normConst;
Subth=Subth./normConst;
nTime=size(normTr,2);
CSvolt=NaN(1,nTime);
CSvolt(SomPCresult{N2show(1)}.CStrace(N2show(2),:)>0)=normTr(1,SomPCresult{N2show(1)}.CStrace(N2show(2),:)>0);
spt=find(SomPCresult{N2show(1)}.spike(N2show(2),:));
figure(111); clf;
plot(normTr,'color',[0 0 0]); hold all
%plot(Subth,'color',[0 0 0]);
plot(CSvolt,'r');
%plot(SomPCresult{N2show(1)}.Blue,'b');
plot(rescale(SomPCresult{N2show(1)}.VR(5,:))*3-7,'color',[0 1 0.2],'LineWidth',1.5);
%plot([spt; spt]',[9 10]','color',[0 0 0],'linewidth',1.5)
%drawScaleBar(60000,'horizontal'); axis off;
axis off;

figure(113); clf;
lap2show=[14:24 26 27];
PF2show_center=66.7;
scaleoffset=show_traces_align_Position_wCS(normTr, SomPCresult{N2show(1)}.CStrace(N2show(2),:)>0, Subth, PF2show_center, [-4000:3000], SomPCresult{N2show(1)}.VR,lap2show,10); box off;
%show_traces_align_Position_wCS(SomPCresult{N2show(1)}.subThresholdSlow(N2show(2),:)/normConst, [], [], PF2show_center, [-4000:3000], SomPCresult{N2show(1)}.VR,lap2show,scaleoffset,1.5);
% figure(114);
% plot(normTr,'color',[0.5 0.5 0.5]); hold all
% plot(Subth,'color',[0 0 0]);
% plot(CSvolt,'r');
% %plot(SomPCresult{N2show(1)}.Blue,'b');
% %plot(rescale(SomPCresult{N2show(1)}.VR(5,:))*3-6,'color',[0 1 0.2]);
% plot([spt; spt]',[6 7]','color',[0 0 0],'linewidth',1.5)
% xlim([37000 152500]);
% drawScaleBar(10000,'horizontal'); axis off;

figure(112); clf; tiledlayout(4,1,'padding','tight');
cmap_rateMap=[1 0 0; 0.5 0.5 0.5; 0 0 0]; show_lap=setdiff([N2show(3):N2show(4)],blueLaps)-1;
LapFR2show=ringmovMean(SomPCresult{N2show(1)}.Lap_FR(2:end-1,:,1,N2show(2)),5);
Lapsub2show=ringmovMean(SomPCresult{N2show(1)}.Lap_SubSlow(2:end-1,:,1,N2show(2)),5);
LapCS2show=ringmovMean(SomPCresult{N2show(1)}.Lap_CS(2:end-1,:,1,N2show(2)),5)*1000;
LapSS2show=ringmovMean(SomPCresult{N2show(1)}.Lap_SS(2:end-1,:,1,N2show(2)),5)*1000;
CatRateMap=cat(3,mat2gray(SomPCresult{N2show(1)}.Lap_CS(2:end-1,:,1,N2show(2))),mat2gray(SomPCresult{N2show(1)}.Lap_SS(2:end-1,:,1,N2show(2))),Lapsub2show,LapFR2show);
CatRateMap=ringmovMean(CatRateMap,5)*3;
MeanRate=permute(mean(CatRateMap(show_lap,:,:),1,'omitnan'),[2 3 1]);
posaxis=[1:150]/150*2; % in m
lapaxis=[show_lap];
%nexttile([1 1])
% imagesc(posaxis,lapaxis,LapFR2show(lapaxis,:)); title('AP');
% colorbar; box off;
csax=nexttile([1 1]);
imagesc(posaxis,lapaxis+1,LapCS2show(lapaxis,:,1)); title('CS');
cb=colorbar; box off; cb.Label.String = 'Firing rate (Hz)'; ylabel('VR trials');
nexttile([1 1])
imagesc(posaxis,lapaxis+1,LapSS2show(lapaxis,:,1)); title('SS');
cb=colorbar; box off; cb.Label.String = 'Firing rate (Hz)'; ylabel('VR trials');
nexttile([1 1])
imagesc(posaxis,lapaxis+1,im_merge(CatRateMap(lapaxis,:,1:2),[1 0 0; 1 1 1])); title('Merged rate map of CS and SS');
ylabel('VR trials');
%colorbar; box off;
% subax=nexttile([1 1]);
% imagesc(posaxis,lapaxis+1,Lapsub2show(lapaxis,:)); title('Subthreshold');
% colorbar; box off;
nexttile([1 1])
colormap(gray);
l=plot(posaxis,rescale2(MeanRate(:,[1 2]),1),'linewidth',1.5); legend({'CS','SS'});
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(cmap_rateMap([1 2],:),2));
xlabel('Position (m)'); ylabel('Tuning curve'); box off; %xlim([0.7 1.5]);
%colormap(subax,gen_colormap([0 0.2 1;1 1 1;1 0 0]));
colormap(csax,gen_colormap([0 0 0;1 0 0]));

%figure(114); clf; cmap_rateMap=[1 0 0; 0.5 0.5 0.5; 0 0 0]; show_lap=setdiff([N2show(3):N2show(4)],blueLaps)-1;
% LapFR2show=ringmovMean(SomPCresult{N2show(1)}.Lap_FR(2:end-1,:,1,N2show(2)),5);
% Lapsub2show=ringmovMean(SomPCresult{N2show(1)}.Lap_SubSlow(2:end-1,:,1,N2show(2)),5);
% LapCS2show=ringmovMean(SomPCresult{N2show(1)}.Lap_CS(2:end-1,:,1,N2show(2)),5);
% LapSS2show=ringmovMean(SomPCresult{N2show(1)}.Lap_SS(2:end-1,:,1,N2show(2)),5);
% CatRateMap=cat(3,mat2gray(SomPCresult{N2show(1)}.Lap_CS(2:end-1,:,1,N2show(2))),mat2gray(SomPCresult{N2show(1)}.Lap_SS(2:end-1,:,1,N2show(2))),Lapsub2show,LapFR2show);
% CatRateMap=ringmovMean(CatRateMap,5)*3;
% MeanRate=permute(mean(CatRateMap(show_lap,:,:),1,'omitnan'),[2 3 1]);
% posaxis=[1:150]/150*2; % in m
% lapaxis=[show_lap];
% pfsz = size(LapMat_sub);
%     pfcentroid = (LapMat_sub * [1:pfsz(2)]') ./ sum(LapMat_sub,2);
%     pfsize = sqrt(sum(LapMat_sub .* (repmat([1:pfsz(2)],pfsz(1),1) - pfcentroid).^2,2)./ sum(LapMat_sub,2));

%% Show voltage traces
figure(113); clf; g=1;
tr2showind=unique([PFs.FileInd PFs.Nind],'rows');
for pf=1:size(tr2showind,1)
    F2show=tr2showind(pf,1);
    Neuron2show=tr2showind(pf,2);
    nTime=size(SomPCresult{F2show}.normtrace,2);
    tr2show=rescale(SomPCresult{F2show}.normtrace(Neuron2show,1:end-100));
    CSvolt=NaN(1,nTime);
    CSvolt(SomPCresult{F2show}.CStrace(Neuron2show,1:end-100)>0)=tr2show(1,SomPCresult{F2show}.CStrace(Neuron2show,1:end-100)>0);
    plot(tr2show+g,'k'); hold all
    plot(CSvolt+g,'r'); hold all
    drawnow;
    g=g+1;

    LapFR2show=ringmovMean(SomPCresult{N2show(1)}.Lap_FR(2:end-1,:,1,N2show(2)),5);
    Lapsub2show=ringmovMean(SomPCresult{N2show(1)}.Lap_SubSlow(2:end-1,:,1,N2show(2)),5);
    LapCS2show=ringmovMean(SomPCresult{N2show(1)}.Lap_CS(2:end-1,:,1,N2show(2)),5);
    LapSS2show=ringmovMean(SomPCresult{N2show(1)}.Lap_SS(2:end-1,:,1,N2show(2)),5);
end
axis off;

%% Measure PFs
CSinPFmat=[]; SSinPFmat=[]; SSpositionNorm=[]; CSpositionNorm=[]; Lap_normalizedPF=[];
g=1; 
subplacefield_limit=4; % S.d.
NormPF_edge=linspace(-subplacefield_limit,subplacefield_limit,40);
figure(25); clf;
for pf=1:size(PFs,1)
    SSpositionNorm{pf}=[]; CSpositionNorm{pf}=[];
    F2show=PFs.FileInd(pf);
    Neuron2show=PFs.Nind(pf);
    nTime=size(SomPCresult{F2show}.normtrace,2);

    if isfield(SomPCresult{F2show},'Blue')
        blueLaps=unique(double(SomPCresult{F2show}.Blue>0).*SomPCresult{F2show}.VR(8,:));
        blueLaps(blueLaps==0)=[];
    else
        blueLaps=[];
    end

    Lap2show=[PFs.LapStart(pf):PFs.LapEnd(pf)];
    if PFs.PFstart(pf)>PFs.PFend(pf)
    Bin2show=[PFs.PFstart(pf):PFs.PFend(pf)+150];        
    else
    Bin2show=[PFs.PFstart(pf):PFs.PFend(pf)];
    end
    Lap2show=setdiff(Lap2show,blueLaps);

    LapMat=cat(3,SomPCresult{F2show}.Lap_CS(:,:,1,Neuron2show),SomPCresult{F2show}.Lap_SS(:,:,1,Neuron2show) ...
        ,SomPCresult{F2show}.Lap_FR(:,:,1,Neuron2show));
    LapMat2=repmat(LapMat,1,3);
    LapMat=repmat(LapMat,1,2); 
    LapMat_subFR=LapMat(Lap2show,Bin2show,3);
    LapMat_suball=LapMat2(Lap2show,:,:);

    bin_dist=ceil(SomPCresult{F2show}.VR(5,:)/((115)/150));
    lap_dist=SomPCresult{F2show}.VR(8,:);
    isRun=movmean(SomPCresult{F2show}.VR(end,:),50)>0.002;
    Lap_staytime=NaN(max(lap_dist),max(bin_dist));
    for l=1:max(lap_dist)
        for b=1:max(bin_dist)
            Lap_staytime(l,b)=sum(bin_dist==b & lap_dist==l);
        end
    end

    bin_vec=ismember(bin_dist,Bin2show);
    lap_vec=SomPCresult{F2show}.VR(8,:).*ismember(SomPCresult{F2show}.VR(8,:),Lap2show);
    %PFinCS_vec=lap_vec(1:nTime).*bin_vec(1:nTime).*SomPCresult{F2show}.SpClass(2,:,Neuron2show);
    PFinCS_vec=lap_vec(1:nTime).*bin_vec(1:nTime).*SomPCresult{F2show}.spike(Neuron2show,:).*SomPCresult{F2show}.CStrace(Neuron2show,:);
    PFinCSfrst_vec=lap_vec(1:nTime).*bin_vec(1:nTime).*SomPCresult{F2show}.SpClass(2,:,Neuron2show);
    PFinSS_vec=lap_vec(1:nTime).*bin_vec(1:nTime).*SomPCresult{F2show}.SpClass(1,:,Neuron2show);
    CSinds=unique(PFinCS_vec); CSinds(CSinds==0)=[];

    LapMat_subFR(isnan(LapMat_subFR))=0;
    pfsz = size(LapMat_subFR);
    pfcentroid = (LapMat_subFR * [1:pfsz(2)]') ./ sum(LapMat_subFR,2);
    pfsize = sqrt(sum(LapMat_subFR .* (repmat([1:pfsz(2)],pfsz(1),1) - pfcentroid).^2,2)./ sum(LapMat_subFR,2));
    
    for lapsub=1:size(LapMat_subFR,1)
        if pfsize(lapsub)>0.01
        pfrange=[-pfsize(lapsub)*subplacefield_limit:pfsize(lapsub)*subplacefield_limit];
        Lap_normpftmp=LapMat_suball(lapsub,round(pfcentroid(lapsub)+Bin2show(1)+150+pfrange),:);
        for sptype=1:3
        Lap_normalizedPF{pf}(lapsub,:,sptype)=interp1(pfrange/pfsize(lapsub),Lap_normpftmp(:,:,sptype),NormPF_edge);
        end
        else
        Lap_normalizedPF{pf}(lapsub,:,sptype)=NaN;    
        end
    end

    nexttile([1 1]);
    imagesc(ringmovMean(LapMat_subFR,5)); hold all
    scatter(pfcentroid,[1:pfsz(1)],'k','filled');
    scatter(pfcentroid-pfsize,[1:pfsz(1)],'k','filled');
    scatter(pfcentroid+pfsize,[1:pfsz(1)],'k','filled');
    title(num2str(pf))

    LapwCS=ismember(Lap2show,CSinds);
    for Laps=1:length(Lap2show)
        Distfromcentroid=bin_dist(PFinCS_vec==Lap2show(Laps))-pfcentroid(Laps)-Bin2show(1);
        Distfromcentroid_SS=bin_dist(PFinSS_vec==Lap2show(Laps))-pfcentroid(Laps)-Bin2show(1);
        SSpositionNorm{pf}=[SSpositionNorm{pf} (bin_dist(PFinSS_vec==Lap2show(Laps))-pfcentroid(Laps)-Bin2show(1))/pfsize(Laps)];
        CSpositionNorm{pf}=[CSpositionNorm{pf} (bin_dist(PFinCS_vec==Lap2show(Laps))-pfcentroid(Laps)-Bin2show(1))/pfsize(Laps)];
        SS_sptime=find(PFinSS_vec==Lap2show(Laps)); CS_sptime=find(PFinCS_vec==Lap2show(Laps));
        
        scatter(bin_dist(PFinCS_vec==Lap2show(Laps))-Bin2show(1),Laps,'r','filled')
        CSinPFmat=[CSinPFmat; [repmat([pf F2show Neuron2show pfsize(Laps) Lap2show(Laps)],length(Distfromcentroid),1) Distfromcentroid' CS_sptime' isRun(CS_sptime)']];        
        SSinPFmat=[SSinPFmat; [repmat([pf F2show Neuron2show pfsize(Laps) Lap2show(Laps)],length(Distfromcentroid_SS),1) Distfromcentroid_SS' SS_sptime' isRun(SS_sptime)']];        
    end
    g=g+1;
    drawnow;
end

FieldName={'PFind','Fileind','Nind','PFsize','LapNumber','DistfromCent','SpikeTime','IsRun'};
CSinPFmat=array2table(CSinPFmat,'VariableNames',FieldName);
FieldName={'PFind','Fileind','Nind','PFsize','LapNumber','DistfromCent','SpikeTime','IsRun'};
SSinPFmat=array2table(SSinPFmat,'VariableNames',FieldName);
%% Prevalence of CS in the periphery of place field, Figure 5K
%omitPFs=[4 12 14 15 42 43 46];
figure(26); clf; tiledlayout(2,1); ax1=[];
pf2use=setdiff([1:size(PFs,1)],[]); PFsizeRange=[7.75 12];
meanPFsize=NaN(size(PFs,1),1);
for pf=unique(SSinPFmat.PFind)'; meanPFsize(pf,:)=mean(SSinPFmat.PFsize(SSinPFmat.PFind==pf)); end;
pf2use2=setdiff(find(meanPFsize>= PFsizeRange(1) & meanPFsize< PFsizeRange(2)),[]);

sp2useSS=ismember(SSinPFmat.PFind,pf2use) & SSinPFmat.IsRun;
sp2useCS=ismember(CSinPFmat.PFind,pf2use) & CSinPFmat.IsRun;

SSpositionNorm_noInf = cellfun(@(x) subsasgn(x, substruct('()', {find(abs(x)>100)}), NaN(1,length(find(abs(x)>100)))), SSpositionNorm, 'UniformOutput', false);
distancebin=linspace(-0.3,0.3,30);
fieldsizebin=linspace(0.05,0.2,30);
[X, Y] = ndgrid(distancebin,fieldsizebin);
[DensityAP,~] = ksdensity([[SSinPFmat.DistfromCent(sp2useSS);CSinPFmat.DistfromCent(sp2useCS)] [SSinPFmat.PFsize(sp2useSS); CSinPFmat.PFsize(sp2useCS)]]/150*2, [X(:), Y(:)]);
[DensitySS,~] = ksdensity([SSinPFmat.DistfromCent(sp2useSS)/150*2 SSinPFmat.PFsize(sp2useSS)/150*2], [X(:), Y(:)]);
[DensityCS,~] = ksdensity([CSinPFmat.DistfromCent(sp2useCS)/150*2 CSinPFmat.PFsize(sp2useCS)/150*2], [X(:), Y(:)]);
DensityAP=reshape(DensityAP,size(X,1),size(X,2)); DensitySS=reshape(DensitySS,size(X,1),size(X,2)); DensityCS=reshape(DensityCS,size(X,1),size(X,2));
FractionCShist=DensityCS./(DensityAP).*(DensitySS>0.0);

ax1=[ax1 nexttile([1 1])];
% [~, cont]=contour(X,Y,DensityAP,5); hold all;
% cont.EdgeColor=[1 0 0];
scatter_density(CSinPFmat.DistfromCent(sp2useCS)/150*2,CSinPFmat.PFsize(sp2useCS)/150*2,35,[],[20 40],0.15); hold all;
[~, cont]=contour(X,Y,DensityAP,5); hold all;
cont.EdgeColor=[1 0 0]; box off; colorbar off;
%CS_count=scatter_heatmap2(CSinPFmat.DistfromCent/150*2,CSinPFmat.PFsize/150*2,distancebin,fieldsizebin);

%densityScatterChart(CSinPFmat.DistfromCent,CSinPFmat.PFsize/150*2);
xlabel('Distance from place field centroid (m)'); ylabel('Place field size (m)'); 
% cb=colorbar; caxis([20 40]);
% cb.Label.String='Probability density';

ax1=[ax1 nexttile([1 1])];
scatter_density(SSinPFmat.DistfromCent(sp2useSS)/150*2,SSinPFmat.PFsize(sp2useSS)/150*2,35,[],[15 35],0.07); hold all;
[~, cont]=contour(X,Y,DensityAP,5); hold all;
cont.EdgeColor=[1 0 0]; box off;
%SS_count=scatter_heatmap2(SSinPFmat.DistfromCent/150*2,SSinPFmat.PFsize/150*2,distancebin,fieldsizebin);
%[~, cont]=contour(X,Y,DensityCS,5);
%densityScatterChart(CSinPFmat.DistfromCent,CSinPFmat.PFsize/150*2);
xlabel('Distance from place field centroid (m)'); ylabel('Place field size (m)'); 
cb=colorbar; caxis([20 40]);
cb.Label.String='Probability density';
linkaxes(ax1,'xy')
xlim([-0.3 0.3]); ylim([0.05 0.2]);
%%
figure(27); clf; tiledlayout(2,1); ax1=[];
nexttile([1 1]);
meanFiringRate_normPF=cell2mat(cellfun(@(x) permute(mean(x,1,'omitnan'),[2 1 3]),Lap_normalizedPF,'UniformOutput',false));
meanFiringRate_normPF=1000*permute(meanFiringRate_normPF,[2 1 3]);
meanFiringRate_normPF=mean(cat(4,meanFiringRate_normPF,fliplr(meanFiringRate_normPF)),4,'omitnan');
meanFiringRate_normPF=fliplr(meanFiringRate_normPF(:,1:length(NormPF_edge)/2,:));
NormPF_edge_abs=fliplr(-NormPF_edge(1:length(NormPF_edge)/2));
Fration_normPF=meanFiringRate_normPF(:,:,1)./meanFiringRate_normPF(:,:,3);
M_fr_normPF=mean(meanFiringRate_normPF(pf2use2,:,:),1,'omitnan'); Sem_fr_normPF=std(meanFiringRate_normPF(pf2use2,:,:),0,1,'omitnan')./sqrt(sum(~isnan(meanFiringRate_normPF(pf2use2,:,:)),1));
M_ratio_normPF=mean(Fration_normPF(pf2use2,:),1,'omitnan'); Sem_ratio_normPF=std(Fration_normPF(pf2use2,:),0,1,'omitnan')./sqrt(sum(~isnan(meanFiringRate_normPF(pf2use2,:,:)),1));

errorbar_shade(NormPF_edge_abs,M_fr_normPF(:,:,3)',Sem_fr_normPF(:,:,3)',[0 0 0]); hold all
errorbar_shade(NormPF_edge_abs,M_fr_normPF(:,:,1)',Sem_fr_normPF(:,:,1)',[1 0 0]); 
ylabel('Firing rate (Hz)'); yyaxis right;
errorbar_shade(NormPF_edge_abs,M_ratio_normPF(:,:,1)',Sem_ratio_normPF(:,:,1)',[0 0 1]);
p=get_pValue([Fration_normPF(pf2use2,1) mean(Fration_normPF(pf2use2,4:6),2)],1)
xlabel('Normalized place field (s.d.)'); ylabel('Fraction of complex spikes'); xlim([0 4]); box off;

nexttile([1 1]);
meanFiringRate_normPF=cell2mat(cellfun(@(x) permute(mean(x,1,'omitnan'),[2 1 3]),Lap_normalizedPF,'UniformOutput',false));
meanFiringRate_normPF=1000*permute(meanFiringRate_normPF,[2 1 3]);
meanFiringRate_normPF=movmean(meanFiringRate_normPF,5,2);
allFR_normPF=meanFiringRate_normPF(:,:,3);
%allFR_normPF=movmean(meanFiringRate_normPF(:,:,3),5,2,'omitnan');
Fration_normPF=meanFiringRate_normPF(:,:,1)./allFR_normPF;
M_fr_normPF=mean(meanFiringRate_normPF(pf2use2,:,:),1,'omitnan'); Sem_fr_normPF=std(meanFiringRate_normPF(pf2use2,:,:),0,1,'omitnan')./sqrt(sum(~isnan(meanFiringRate_normPF(pf2use2,:,:)),1));
M_ratio_normPF=mean(Fration_normPF(pf2use2,:),1,'omitnan'); Sem_ratio_normPF=std(Fration_normPF(pf2use2,:),0,1,'omitnan')./sqrt(sum(~isnan(meanFiringRate_normPF(pf2use2,:,:)),1));

errorbar_shade(NormPF_edge,M_fr_normPF(:,:,3)',Sem_fr_normPF(:,:,3)',[0 0 0]); hold all
errorbar_shade(NormPF_edge,M_fr_normPF(:,:,1)',Sem_fr_normPF(:,:,1)',[1 0 0]); 
ylabel('Firing rate (Hz)'); yyaxis right;
errorbar_shade(NormPF_edge,M_ratio_normPF(:,:,1)',Sem_ratio_normPF(:,:,1)',[0 0 1]);
p=get_pValue([Fration_normPF(pf2use2,20) mean(Fration_normPF(pf2use2,16:17),2) mean(Fration_normPF(pf2use2,25:27),2)],1)
xlabel('Normalized place field (s.d.)'); ylabel('Fraction of complex spikes'); xlim([-2.5 2.5]); box off;
set_fontsize(13);
% plot(NormPF_edge_abs,0.2*(p(1,:)<0.05),'c')
% plot(NormPF_edge_abs,0.3*(p(1,:)<0.01),'g')
% ax1=[ax1 nexttile([1 1])];
% imagesc(distancebin,fieldsizebin,im_merge(cat(3,mat2gray(SS_count),mat2gray(CS_count)),[1 1 1;1 0 0]));
% 
% ax1=[ax1 nexttile([1 1])];
% scatter_density(SSinPFmat.DistfromCent/150*2,SSinPFmat.PFsize/150*2,35,[]);
% hold all
% scatter(CSinPFmat.DistfromCent/150*2,CSinPFmat.PFsize/150*2,15,[1 0 0],'filled')
% %densityScatterChart(CSinPFmat.DistfromCent,CSinPFmat.PFsize/150*2);
% xlabel('Distance from place field centroid (m)'); ylabel('Place field size (m)'); 
% cb=colorbar;
% cb.Label.String='Relative density';

%SShist=histcounts2(SSinPFmat.DistfromCent/150*2,SSinPFmat.PFsize/150*2,distancebin,fieldsizebin);
%[SShist Xc Yc]=histcounts2(SSinPFmat.DistfromCent/150*2,SSinPFmat.PFsize/150*2,'BinMethod','auto');
%[CShist Xc Yc]=histcounts2(CSinPFmat.DistfromCent/150*2,CSinPFmat.PFsize/150*2,distancebin,fieldsizebin);
%CShist=histcounts2(CSinPFmat.DistfromCent/150*2,CSinPFmat.PFsize/150*2,Xc,Yc);

% [X, Y] = meshgrid(mean([distancebin(1:end-1); distancebin(2:end)]), mean([fieldsizebin(1:end-1); fieldsizebin(2:end)]));
% FractionCShist=CShist./(SShist+CShist);
%%
HistMat=[]; normpf_bin=[-3:0.2:-2 -1.8 -1.2 -0.8 -0.5 -0.2 0 0.2 0.5 0.8 1.2 1.8 2:0.2:3]; g=1;
for pf=setdiff([1:size(PFs,1)],omitPFs)
HistMat{g,1}=histcounts(([SSpositionNorm_noInf{pf} CSpositionNorm{pf}]),normpf_bin,'Normalization','countdensity');
HistMat{g,2}=histcounts(([SSpositionNorm{pf}]),normpf_bin,'Normalization','countdensity');
HistMat{g,3}=histcounts(([CSpositionNorm{pf}]),normpf_bin,'Normalization','countdensity');

%ksdensity([CSinPFmat.DistfromCent/150*2 CSinPFmat.PFsize/150*2], [X(:), Y(:)]);
g=g+1;
end
APcountmat=cell2mat(HistMat(:,1));
CScountmat=cell2mat(HistMat(:,3));
Fractionmat=CScountmat./APcountmat;
MeanAP=mean(APcountmat,1,'omitnan'); SemAP=std(APcountmat,0,1,'omitnan')./sqrt(sum(~isnan(APcountmat),1));
MeanCS=mean(CScountmat,1,'omitnan'); SemCS=std(CScountmat,0,1,'omitnan')./sqrt(sum(~isnan(CScountmat),1));
MeanFrac=mean(Fractionmat,1,'omitnan'); SemFrac=std(Fractionmat,0,1,'omitnan')./sqrt(sum(~isnan(Fractionmat),1));

ax1=[ax1 nexttile([1 1])];
errorbar_shade(mean([normpf_bin(2:end); normpf_bin(1:end-1)],1),MeanAP,SemAP,[0 0 0]); hold all; 
errorbar_shade(mean([normpf_bin(2:end); normpf_bin(1:end-1)],1),MeanCS,SemCS,[1 0 0]); yyaxis right;
errorbar_shade(mean([normpf_bin(2:end); normpf_bin(1:end-1)],1),MeanFrac,SemFrac,[1 0 1]);

linkaxes(ax1,'xy')
%%
PFsizeRange=[0.05 0.2]/2*150;
cs2showInd=CSinPFmat.PFsize>PFsizeRange(1) & CSinPFmat.PFsize<PFsizeRange(2);
for pf=unique(CSinPFmat.PFind(cs2showInd))'
    LapInRange=CSinPFmat.LapNumber(cs2showInd);
    unique(LapInRange(CSinPFmat.PFind(cs2showInd)==pf))
end


%%
figure(27); clf; PFbinEdge=[-3:0.3:-1 -0.9:0.1:0.9 1:0.3:3]; %tiledlayout(1,2);
CScountsbinned=[]; SScountsbinned=[]; APcountsbinned=[]; CScountsbinned_pdf=[]; g=1;
PFsizeRange=[0.05 0.15]/2*150;
pf2showInd=CSinPFmat.PFsize>PFsizeRange(1) & CSinPFmat.PFsize<PFsizeRange(2);
for pf=unique(CSinPFmat.PFind(pf2showInd))'
%[CScountsbinned(pf,:)]=histcounts(CSinPFmat.DistfromCent(CSinPFmat.PFind==pf)./CSinPFmat.PFsize(CSinPFmat.PFind==pf),PFbinEdge);
APcountsbinned(g,:)=histcounts([SSpositionNorm{pf} CSpositionNorm{pf}],PFbinEdge);
SScountsbinned(g,:)=histcounts([SSpositionNorm{pf}],PFbinEdge);
CScountsbinned(g,:)=histcounts(CSpositionNorm{pf},PFbinEdge);
%CScountsbinned_pdf(g,:)=histcounts(CSpositionNorm{pf},PFbinEdge,'Normalization','probability');
g=g+1;
end
PFbin_c=mean([PFbinEdge(1:end-1); PFbinEdge(2:end)]);
Ratiobinned=CScountsbinned./(SScountsbinned+CScountsbinned);
Ratiobinned(Ratiobinned==Inf)=NaN;
%plot(PFbin_c,CScountsbinned','color',[0.7 0.7 0.7]); hold all
MC=mean(CScountsbinned,1,'omitnan'); SC=std(CScountsbinned,0,1,'omitnan')./sqrt(size(CScountsbinned,1));
MS=mean(SScountsbinned,1,'omitnan'); SS=std(SScountsbinned,0,1,'omitnan')./sqrt(size(SScountsbinned,1));
MR=mean(Ratiobinned,1,'omitnan'); SR=std(Ratiobinned,0,1,'omitnan')./sqrt(size(SScountsbinned,1));

nexttile([1 1])
errorbar(PFbin_c,MC,SC,'r','linewidth',1.5); hold all; box off;
errorbar(PFbin_c,MS,SS,'k','linewidth',1.5); hold all; yyaxis right;
%plot(PFbin_c,Ratiobinned','color',[0.7 0.7 0.7]); hold all
errorbar(PFbin_c,MR,SR,'b','linewidth',1.5); hold all; box off;
xlabel('Normalized place field'); ylabel('Fraction of CS'); axis tight; ylim([0 0.5]);

% nexttile([1 1])
% plot(PFbin_c,CScountsbinned_pdf','color',[0.7 0.7 0.7]); hold all
% errorbar(PFbin_c,M,S,'r','linewidth',1.5); hold all; box off;
% xlabel('Normalized place field'); ylabel('CS probability');  axis tight; ylim([0 0.4]);

% nexttile([1 1])
% plot(PFbin_c,CScountsbinned','color',[0.7 0.7 0.7]); hold all
% errorbar(PFbin_c,MC,SC,'r'); hold all; box off;
% xlabel('Normalized place field'); ylabel('CS count'); axis tight;