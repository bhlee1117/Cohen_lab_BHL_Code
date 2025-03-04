clear
clc;
[~, ~, raw] = xlsread(['/Volumes/cohen_lab/Lab/Labmembers' ...
    '/Byung Hun Lee/Data/PrismPCdata_Arrangement.xlsx'], 'Sheet1', 'C5:W31');

ref_ROI=cellfun(@(x) (str2num(num2str(x))),raw(:,10),'UniformOutput',false);
basal_ROI=cellfun(@(x) (str2num(num2str(x))),raw(:,11),'UniformOutput',false);
apical_ROI=cellfun(@(x) (str2num(num2str(x))),raw(:,12),'UniformOutput',false);
fpath=raw(:,1)';
StructureData=raw(:,10);
BadROI=cellfun(@(x) (str2num(num2str(x))),raw(:,15),'UniformOutput',false);
EndFrame=cell2mat(raw(:,13));
ifmotionReject=cell2mat(raw(:,14));
ifdirtRemov=cell2mat(raw(:,16));
Pixelsize=cell2mat(raw(:,6));
save_figto='/Volumes/BHL_WD18TB/PP72_PlaceCellResults';
place_bin=150; time_segment=15000; overlap=200;
alignedMovFN = {'STA_Mat_SS','STA_Mat_CS','STA_Mat_dSP'};
bound=6;
title_str={'Basal','Apical','Peri-Soma'};
PlaceFieldList=cellfun(@(x) (str2num(num2str(x))),raw(:,19),'UniformOutput',false);
PlaceFieldBin=cellfun(@(x) (str2num(num2str(x))),raw(:,20),'UniformOutput',false);
set(0,'DefaultFigureWindowStyle','docked')
foi=[1 4 5 6 8 10 11 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27];

%%

for f=27
    load(fullfile(fpath{f},'PC_Result.mat'),'Result');
    noi=setdiff([1:size(Result.ftprnt,3)],BadROI{f});
    %noi=[1:size(Result.ftprnt,3)];
    noi_dist=ismember(Result.dist_order,noi);
    %validFrame=find(Result.motionReject==0 & sum(Result.dirtTrace,1)==0);
    validFrame=find(Result.motionReject==0);
    [V D eigTrace]=get_eigvector(movmean(Result.normTraces(noi,validFrame),10,2)',10);
    %ind2use=find(cumsum(D)/sum(D)>0.80,1);
    ind2use=[2];
    F0_PCA=ones(1,size(Result.ftprnt,3));
    F0_PCA(noi)=sqrt(sum((V(:,[1:ind2use]).*sqrt(D([1:ind2use]))').^2,2));
    %normTr=Result.normTraces./F0_PCA';
    normTr=Result.normTraces./Result.F_ref;
    normTr=normTr(Result.dist_order(noi_dist),:);
    SpMat=Result.SpClass(1:2,:);
    SpMat_modified=interactiveFrameCheck(normTr,SpMat,[200]);
    
    Result.SpClass(1:2,:)=SpMat_modified;
    CS_trace_add=tovec(find(SpMat_modified(2,:))+[-2:50]');
    Result.CStrace(CS_trace_add)=1;
    unmarkspike=find((SpMat(1,:)-SpMat_modified(1,:))==1 & Result.CStrace==0);
    markspike=find((SpMat_modified(1,:)-SpMat(1,:))==1);

    Result.spike(ref_ROI{f},unmarkspike)=0;
    Result.spike(ref_ROI{f},markspike)=1;
    Result.spike(ref_ROI{f},find(max(Result.SpClass(1:2,:),[],1)>0))=1;

    save(fullfile(fpath{f},'PC_Result.mat'),'Result','-v7.3');
    load(fullfile(fpath{f},'PC_Result.mat'),'Result');
    backupServer(fpath{f},'BHL18TB_D2','ByungHun2TB/Prism_PlaceCell_backup','PC_Result.mat')
end

%%
%foi=[1 4 5 6 8 10 11 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27];

for f=11
load(fullfile(fpath{f},'PC_Result.mat'),'Result');
noi=setdiff([1:size(Result.ftprnt,3)],BadROI{f});
    %noi=[1:size(Result.ftprnt,3)];
    noi_dist=ismember(Result.dist_order,noi);
    %validFrame=find(Result.motionReject==0 & sum(Result.dirtTrace,1)==0);
    %sptime=ind2vec(size(Result.SpClass,2),unique(find(max(Result.SpClass(1:2,:),[],1))'+[-10:30]),1);
    validFrame=find(Result.motionReject==0);
    %validFrame=[1:size(Result.SpClass,2)];
    [V D eigTrace]=get_eigvector(movmean(Result.normTraces(noi,validFrame),10,2)',25);
    ind2use=find(cumsum(D)/sum(D)>0.95,1)
    %cumsum(D)/sum(D)
    ind2use=[6];
    F0_PCA=ones(1,size(Result.ftprnt,3));
    F0_PCA(noi)=sqrt(sum((V(:,[1:ind2use]).*sqrt(D([1:ind2use]))').^2,2));
    normTr=Result.normTraces./F0_PCA';    
    normTr=normTr(Result.dist_order(noi_dist),:);
    normTr_ref=Result.normTraces./Result.F_ref;
    normTr_ref=normTr_ref(Result.dist_order(noi_dist),:);

    Result.F0_PCA=F0_PCA;
    
bwCS=bwlabel(Result.CStrace);
ValidCS=unique(bwCS.*Result.SpClass(2,:));
ValidCS=ValidCS(2:end);
Result.CStrace=ismember(bwCS,ValidCS);

bwCS=bwlabel(Result.CStrace);
newCS=zeros(1,size(bwCS,2));
for b=1:max(bwCS)
frmb=find(bwCS==b);
CS1st=find(Result.SpClass(2,:)==1 & bwCS==b,1);
    newCS(CS1st:frmb(end))=1;
end
Result.CStrace=newCS;

tr_ref=Result.normTraces(ref_ROI{f}(1),:);
for sp_type=1:2
sp=max(Result.SpClass(sp_type,:),[],1);
[~, shift]=max(reshape(tr_ref(find(sp>0)+[-1:0]'),2,[]),[],1);
shift=shift-2;
sp_time = find(sp>0)+shift;
sp_new=zeros(1,length(tr_ref));
sp_new(sp_time)=1;
Result.SpClass(sp_type,:)=sp_new;
end
notdSP=find(Result.SpClass(3,:)>0 & max(Result.SpClass([1:2],:),[],1)>0);
Result.SpClass(3,notdSP)=0;

sp=max(Result.spike(ref_ROI{f},:),[],1);
[~, shift]=max(reshape(tr_ref(find(sp>0)+[-1:0]'),2,[]),[],1);
shift=shift-2;
sp_time = find(sp>0)+shift;
sp_new=zeros(1,length(tr_ref));
sp_new(sp_time)=1;
Result.spike(ref_ROI{f},:)=repmat(sp_new,length(ref_ROI{f}),1);

clf; tiledlayout(5,1)
ax1=nexttile([2 1]);
imagesc(normTr,[-0.003 0.01]);
colormap(turbo);
ax3=nexttile([2 1]);
imagesc(normTr_ref,[-3 8]);
colormap(turbo);
ax2=nexttile([1 1]);
plot(Result.SpClass'); hold all
plot(Result.CStrace);
linkaxes([ax1 ax2 ax3],'x')
legend({'SS','CS','dS','CStr'})

save(fullfile(fpath{f},'PC_Result.mat'),'Result','-v7.3');
    load(fullfile(fpath{f},'PC_Result.mat'),'Result');
    backupServer(fpath{f},'BHL18TB_D2','cohen_lab/Lab/Labmembers/Byung Hun Lee/Data','PC_Result.mat')
end


%%

for f=27
    load(fullfile(fpath{f},'PC_Result.mat'),'Result');
    noi=setdiff([1:size(Result.ftprnt,3)],BadROI{f});
    %noi=[1:size(Result.ftprnt,3)];
    noi_dist=ismember(Result.dist_order,noi);
    normTr=Result.normTraces./Result.F0_PCA';
    normTr=normTr(Result.dist_order(noi_dist),:);
    SpMat=Result.spike(1,:);
    SpMat_modified=interactiveFrameCheck(normTr,SpMat,[100]);
    
tr_ref=Result.normTraces(ref_ROI{f}(1),:);
    sp=SpMat_modified;
[~, shift]=max(reshape(tr_ref(find(sp>0)+[-1:0]'),2,[]),[],1);
shift=shift-2;
sp_time = find(sp>0)+shift;
sp_new=zeros(1,length(tr_ref));
sp_new(sp_time)=1;
Result.spike(ref_ROI{f},:)=repmat(sp_new,length(ref_ROI{f}),1);
end
%%
for f=[5 6 8 10 11 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27]
    f
       load(fullfile(fpath{f},'PC_Result.mat'),'Result');
Result.spike(ref_ROI{f},:)=repmat(shrinkToMax(Result.spike(1,:), Result.normTraces(ref_ROI{f}(1),:)),length(ref_ROI{f}),1);
save(fullfile(fpath{f},'PC_Result.mat'),'Result','-v7.3');
    load(fullfile(fpath{f},'PC_Result.mat'),'Result');
    backupServer(fpath{f},'BHL18TB_D2','cohen_lab/Lab/Labmembers/Byung Hun Lee/Data','PC_Result.mat')
end

%%

for f=foi
    f
load(fullfile(fpath{f},'PC_Result.mat'),'Result');
Result.spike(1,find(Result.SpClass(1,:)==1))=1;
sptime=find(Result.spike(1,:));
consSp=find((sptime(2:end)-sptime(1:end-1))==1);
for ss=consSp
[~, arg]=max(Result.normTraces(ref_ROI{f}(1),sptime(ss)+[0 1]));
if arg==1
Result.spike(ref_ROI{f},sptime(ss))=1;
Result.spike(ref_ROI{f},sptime(ss)+1)=0;
else
Result.spike(ref_ROI{f},sptime(ss))=0;
Result.spike(ref_ROI{f},sptime(ss)+1)=1;
end
end
Result.SpClass(1,:)=Result.spike(1,:);
Result.SpClass(1,:)=Result.SpClass(1,:).*(Result.CStrace==0);
save(fullfile(fpath{f},'PC_Result.mat'),'Result','-v7.3');
load(fullfile(fpath{f},'PC_Result.mat'),'Result');
backupServer(fpath{f},'BHL18TB_D2','cohen_lab/Lab/Labmembers/Byung Hun Lee/Data','PC_Result.mat')
end

    %%

    
