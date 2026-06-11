%% Set the path
clear
clc;
cd '/Users/bhlee1117/Documents/GitHub/Cohen_lab_BHL_Code/Analysis_master_codes';
[~, ~, raw] = xlsread(['/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/' ...
    'NaVInactivationData_Arrangement.xlsx'], 'Sheet1', 'C5:T442');

save_to='/Volumes/BHL18TB_D2/Arranged_Data/Prism_OptopatchResult';
fpath=raw(:,1);
set(0,'DefaultFigureWindowStyle','docked')
time_segment=25000;
DrugInd=raw(:,4);
%% load data
[DrugLegend DrugLabels]=unique(DrugInd([412:438]));

DrugInjectionResult=[]; g=1;
for i=[412:438]
    VoltResult=importdata(fullfile(fpath{i},'OP_Result.mat'));
    fpsplit=split(fpath{i},'/');
    [~, labelIdx] = ismember(DrugInd(i), DrugLegend);
    DrugInjectionResult{g}=VoltResult;
    DrugInjectionResult{g}.Drug=DrugLegend(labelIdx);
    DrugInjectionResult{g}.DrugInd=labelIdx;
    DrugInjectionResult{g}.filename=fpsplit(end);
    DrugInjectionResult{g}.CSfraction=sum(VoltResult.spike.*(VoltResult.CStrace>0),2)./sum(VoltResult.spike,2);
    g=g+1;
end
disp('Data loaded');
%% Show representative trace and CS fraction
DrugList=cellfun(@(x) x.Drug,DrugInjectionResult)';
FileList=cellfun(@(x) x.filename,DrugInjectionResult)';
IsoF=[1 16; 1 15; 3 9; 4 6; 12 12;];
MK801F=[2 15;2 34;5 6;6 2;9 9];

figure(22); clf; tiledlayout(1,2,'Padding','tight'); offsetScale=1.5;
nexttile([1 1])
for cell2show=1:length(IsoF)
tr2show=DrugInjectionResult{IsoF(cell2show,1)}.traces(IsoF(cell2show,2),:);
nTime=size(tr2show,2);
Bluetr=DrugInjectionResult{IsoF(cell2show,1)}.Blue;
Spiketr=DrugInjectionResult{IsoF(cell2show,1)}.spike(IsoF(cell2show,2),:);
CStr=DrugInjectionResult{IsoF(cell2show,1)}.CStrace(IsoF(cell2show,2),:);

tr2show_nonblue=get_blueoffTrace(movmean(tr2show,20),Bluetr,300,50)';
tr2show=tr2show-movmedian(tr2show_nonblue,1000);
subtr2show=get_subthreshold(tr2show,Spiketr,7,17);
SpikeHeight=get_STA(tr2show,Spiketr,0,0);
%tr2showZscored=tr2show./std(subtr2show,0,2);
tr2showZscored=tr2show./SpikeHeight;
CS2showZscored=tr2showZscored;
CS2showZscored(CStr==0)=NaN;

plot([1:nTime]*0.001,tr2showZscored+offsetScale*cell2show,'k'); hold all
plot([1:nTime]*0.001,CS2showZscored+offsetScale*cell2show,'r'); hold all
end
plot([1:nTime]*0.001,rescale(Bluetr)*2/3+0.2,'color',[0 0.6 1]);
axis off;
drawScaleBar(1,'horizontal','color','k','position',[10.2 1.2]);
drawScaleBar(1,'vertical','color','k','position',[10.2 1.2]);
xlim([0.2 10.25]);
title('1% Iso.');

nexttile([1 1]);
for cell2show=1:length(MK801F)
tr2show=DrugInjectionResult{MK801F(cell2show,1)}.traces(MK801F(cell2show,2),:);
nTime=size(tr2show,2);
Bluetr=DrugInjectionResult{MK801F(cell2show,1)}.Blue;
Spiketr=DrugInjectionResult{MK801F(cell2show,1)}.spike(MK801F(cell2show,2),:);
CStr=DrugInjectionResult{MK801F(cell2show,1)}.CStrace(MK801F(cell2show,2),:);

tr2show_nonblue=get_blueoffTrace(movmean(tr2show,20),Bluetr,300,50)';
tr2show=tr2show-movmedian(tr2show_nonblue,1000);
subtr2show=get_subthreshold(tr2show,Spiketr,7,17);
SpikeHeight=get_STA(tr2show,Spiketr,0,0);
%tr2showZscored=tr2show./std(subtr2show,0,2);
tr2showZscored=tr2show./SpikeHeight;
CS2showZscored=tr2showZscored;
CS2showZscored(CStr==0)=NaN;

plot([1:nTime]*0.001,tr2showZscored+offsetScale*cell2show,'k'); hold all
plot([1:nTime]*0.001,CS2showZscored+offsetScale*cell2show,'r'); hold all;
end
plot([1:nTime]*0.001,rescale(Bluetr)*2/3+0.2,'color',[0 0.6 1]);
axis off;
drawScaleBar(1,'horizontal','color','k','position',[10.2 1.2]);
drawScaleBar(1,'vertical','color','k','position',[10.2 1.2]);
xlim([0.2 10.25]);
title('1% Iso. + MK-801 0.5 µg');
set_font('Arial'); set_fontsize(13);
set_figsize(300,150);
copygraphics(gcf, 'ContentType', 'vector');
%%
CompareDrug=[1 2];
fprintf('Comparing drug label %s and %s.\n',DrugLegend{CompareDrug(1)},DrugLegend{CompareDrug(2)});
DrugIndList=cellfun(@(x) x.DrugInd,DrugInjectionResult)';
CSfraction=cellfun(@(x) x.CSfraction(sum(x.spike,2)>0),DrugInjectionResult,'UniformOutput',false)';
CSfractionCat={cell2mat(CSfraction(DrugIndList==CompareDrug(1))),cell2mat(CSfraction(DrugIndList==CompareDrug(2)))}; %iso & iso+MK801;
figure(23); clf;
p=Violin_wPoints(CSfractionCat,[0 0 0; 0.5 0.5 0.5]);
ylim([0 0.95]);
drawPValueLines(p,0,'TextYOffset',0.05)
set(gca,'XTick',[1 2],'XTickLabel',{'Iso.',sprintf('Iso. + MK-801')});
ylabel('Fraction of spikes in CS');
set_font('Arial'); set_fontsize(14);
fprintf('p = %3.0d, N = 1st cond: %3.0f, 2nd cond: %3.0f, 2 mice\n',p(1,2),length(CSfractionCat{1}),length(CSfractionCat{2}))
set_figsize(60,100); ylim([0 1]);
copygraphics(gcf, 'ContentType', 'vector','Resolution',300);


% show_footprnt_contour(DrugInjectionResult{IsoF}.ftprnt,DrugInjectionResult{IsoF}.ref_im);
% figure(2); clf;
% show_footprnt_contour(DrugInjectionResult{MK801F}.ftprnt,DrugInjectionResult{MK801F}.ref_im);

