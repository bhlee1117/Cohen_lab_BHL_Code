clear
clc;
[~, ~, raw] = xlsread(['/Volumes/cohen_lab/Lab/Labmembers' ...
    '/Byung Hun Lee/Data/PrismPCdata_Arrangement.xlsx'], 'Sheet1', 'C5:S31');

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
set(0,'DefaultFigureWindowStyle','docked')
foi=[1 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27];

%%
nTau_bAP=[-25:15];
nTau=[-70:50];
for f=[27]
f
    load(fullfile(fpath{f},'PC_Result.mat'),'Result')
 load([fpath{f} '/output_data.mat'])
    sz=double(Device_Data{1, 3}.ROI([2 4]));
    
alignmovlist=dir(fullfile(fpath{f},[alignedMovFN{1} '*.bin']));
AlignMov=[];
nReadmov=min([length(alignmovlist) 4]);
for l=1:nReadmov
AlignMov=cat(3,AlignMov,readBinMov(fullfile(fpath{f},alignmovlist(l).name),sz(2),sz(1)));
end

if length(alignmovlist)>4
AlignMovrsh=double(reshape(AlignMov,size(AlignMov,1),size(AlignMov,2),51,[]));
else
    if isfield(Result,'StackedSpike')
AlignMovrsh=double(reshape(AlignMov,size(AlignMov,1),size(AlignMov,2),[],size(Result.StackedSpike{1},2)));
    else
AlignMovrsh=double(reshape(AlignMov,size(AlignMov,1),size(AlignMov,2),[],sum(Result.SpClass(1,:))));        
    end
end

som_spike=find(Result.spike(1,:)); 
bAP_ref=[];
for s=som_spike
    isnearby=sum(ismember(s+nTau_bAP,som_spike))>1;
    isnearbyCS=sum(ismember(s+nTau_bAP,find(Result.CStrace)))>1;
    if ~isnearby & ~isnan(s) & ~isnearbyCS
        bAP_ref=[bAP_ref s];
    end
end
sp_na=sum((bAP_ref'+nTau_bAP)<0 | (bAP_ref'+nTau_bAP)>size(Result.traces,2),2)==0;
bAP_ref=bAP_ref(sp_na);

if isfield(Result,'StackedSpike')
Refindex=ismember(Result.StackedSpike{1,1}(2,1:size(AlignMovrsh,4)),bAP_ref);
else
Refindex=ismember(find(Result.SpClass(1,:)) ,bAP_ref);    
Refindex=Refindex(1:size(AlignMovrsh,4));
end

STAmov=mean(-AlignMovrsh(:,:,:,find(Refindex)),4);
STAmov=STAmov-mean(STAmov(:,:,1:15),3);
STAmovVec=tovec(STAmov);
end
%%

Fslope=[]; Rsq=[]; dff=[]; dff_weight=[]; dff_Ref=[]; dff_Ref_weight=[];
Fslope_weight=[]; 
bound=6;
figure(f); clf;
F0=tovec(Result.ref_im)-100;
F0_filter=tovec(imgaussfilt(Result.ref_im,2))-100;
F_ref=mean(STAmovVec(:,-nTau_bAP(1)+[8:11]),2,'omitnan');
for n=1:size(Result.ftprnt,3)
    
px=find(tovec(Result.ftprnt(:,:,n)>0));
[~, maxfrm]=max(mean(STAmovVec(px,:),1,'omitnan'));

dF=STAmovVec(:,maxfrm);
%dF=max(STAmovVec,[],2);
F0_weight=F0.*rescale(tovec(Result.ftprnt(:,:,n)));
dF_weight=dF.*rescale(tovec(Result.ftprnt(:,:,n)));

px2=px(find(F0(px)>prctile(F0(px),30) & F0(px)<prctile(F0(px),98)));
px2_weight=px(find(F0_weight(px)>prctile(F0_weight(px),30) & F0_weight(px)<prctile(F0_weight(px),98)));

nexttile([1 1])
plot(F0(px),dF(px),'.'); hold all
plot(F0_weight(px),dF_weight(px),'.'); hold all
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
Fslope_weight(n)=p_weight(1);
Rsq(n)=R_squared;

dff(n)=sum((dF).*(tovec(Result.ftprnt(:,:,n))>0),1,'omitnan')./sum((F0_filter).*(tovec(Result.ftprnt(:,:,n))>0),1,'omitnan');
dff_weight(n)=sum((dF).*rescale(tovec(Result.ftprnt(:,:,n))),1,'omitnan')./sum((F0_filter).*rescale(tovec(Result.ftprnt(:,:,n))),1,'omitnan');

dff_Ref(n)=sum((dF).*(tovec(Result.ftprnt(:,:,n))>0),1,'omitnan')./sum((F_ref).*(tovec(Result.ftprnt(:,:,n))>0),1,'omitnan');
dff_Ref_weight(n)=sum((dF).*rescale(tovec(Result.ftprnt(:,:,n))),1,'omitnan')./sum((F_ref).*rescale(tovec(Result.ftprnt(:,:,n))),1,'omitnan');
end

figure(3); clf;

Dsign=ones(1,size(Result.interDendDist,1));
Dsign(Result.dist_order(1:find(Result.dist_order==1)-1))=-1;
dendaxis=Result.interDendDist(1,:).*Dsign;

nexttile([1 1])
plot(dendaxis(Result.dist_order),Fslope(Result.dist_order)/mean(Fslope(ref_ROI{f}))); hold all
plot(dendaxis(Result.dist_order),Fslope_weight(Result.dist_order)/mean(Fslope_weight(ref_ROI{f})))
plot(dendaxis(Result.dist_order([2 7 12 15 18 22 25 31 36 38 40 42 44 46 47])),ones(1,15)*1.5,'ro')
title('Robust dFF')
legend({'Equal','Weighted'})
nexttile([1 1])
plot(dendaxis(Result.dist_order),dff(Result.dist_order)/mean(dff(ref_ROI{f}))); hold all
plot(dendaxis(Result.dist_order),dff_weight(Result.dist_order)/mean(dff_weight(ref_ROI{f})))
plot(dendaxis(Result.dist_order([2 7 12 15 18 22 25 31 36 38 40 42 44 46 47])),ones(1,15)*1.5,'ro')
legend({'Equal','Weighted'})
title('dFF')
nexttile([1 1])
plot(dendaxis(Result.dist_order),dff_Ref(Result.dist_order)/mean(dff_Ref(ref_ROI{f}))); hold all
plot(dendaxis(Result.dist_order),dff_Ref_weight(Result.dist_order)/mean(dff_Ref_weight(ref_ROI{f})))
plot(dendaxis(Result.dist_order([2 7 12 15 18 22 25 31 36 38 40 42 44 46 47])),ones(1,15)*1.5,'ro')
legend({'Equal','Weighted'})
title('F_{ref}')


%%
clf;
plot(F_ref,Fslope,'.')

%%
bound=7;
nTau={[-70:50],[-70:50],[-70:50]}; %SS, CS, Brst

for f=foi(end)
load(fullfile(fpath{f},'PC_Result.mat'),'Result')
% load(fullfile(fpath{f},"output_data.mat"))
% sz=double(Device_Data{1, 3}.ROI([2 4]));
% mov_mc=double(readBinMov([fpath{f} '/mc_ShutterReg' num2str(2,'%02d') '.bin'],sz(2),sz(1)));
% mov_res= mov_mc-mean(mov_mc,3);
% load([fpath{f} '/mcTrace' num2str(2,'%02d') '.mat']);
% mov_res= mov_mc-mean(mov_mc,3);
% mov_res = SeeResiduals(mov_res,mcTrace.xymean);
% mov_res = SeeResiduals(mov_res,mcTrace.xymean.^2);
% mov_res = SeeResiduals(mov_res,mcTrace.xymean(:,1).*mcTrace.xymean(:,end));
% mov_res = mov_res.*double(max(Result.bvMask,[],3)==0);
% mov_res_crop=mov_res(bound:end-bound,bound:end-bound,:);
% 
% [V D eigTrace]=get_eigvector(tovec(mov_res_crop(:,:,1:7780)),10);

%F_ref=Result.F_ref;
NormalizedTrace=(Result.normTraces);

NormalizedTrace_dirt=NormalizedTrace;
NormalizedTrace_dirt(:,Result.motionReject>0)=NaN;
if isfield(Result,'dirtTrace')
NormalizedTrace_dirt(Result.dirtTrace>0)=NaN;
end

Subthreshold=get_subthreshold(NormalizedTrace_dirt,max(Result.spike(1,:),[],1)>0,7,17);
Subthreshold(isnan(NormalizedTrace_dirt))=NaN;

noi=setdiff([1:size(NormalizedTrace_dirt,1)],BadROI{f});
noi_dist=ismember(Result.dist_order,noi);

validtime=find(sum(isnan(Subthreshold),1)==0);
% [V D eigTrace]=get_eigvector(Subthreshold(:,validtime)',10);
% imshow_patch(toimg(squeeze(max(tovec(Result.ftprnt>0).*reshape(rescale2(V,1)+0.1,1,size(Result.ftprnt,3),10),[],2)),size(Result.ref_im,1),size(Result.ref_im,2)))
% 
% F0_PCA=sqrt(sum(((V(:,[1:3]).^2).*D([1:3])'),2));

F0_PCA=get_F0PCA(Subthreshold(:,validtime));

NormalizedTrace_dirt=NormalizedTrace_dirt./Result.F_ref;
 % Isolated Somatic spike
    SS_s=[];
    som_spike=find(Result.spike(1,:));
    for s=som_spike
        if sum((s+nTau{1})<0 | (s+nTau{1})>size(Result.traces,2),2)==0
            isnearby=sum(ismember(s+nTau{1},som_spike))>1;
            isnearbyCS=sum(ismember(s+nTau{1},find(Result.CStrace)))>1;
            ispartCS=sum(Result.CStrace(s+nTau{1}))>0;
            if ~isnearby & ~isnan(s) & ~isnearbyCS & ~ispartCS
                SS_s=[SS_s s];
            end
        end
    end
    sp_na=sum((SS_s'+nTau{1})<0 | (SS_s'+nTau{1})>size(Result.traces,2),2)==0;
    SS_s=SS_s(sp_na);
  STAmat=reshape(NormalizedTrace_dirt(:,SS_s'+nTau{1}),size(Result.traces,1),[],length(nTau{1}));
  SS_STA=squeeze(mean(STAmat(Result.dist_order(noi_dist),:,:),2,'omitnan'));
  SS_STA=SS_STA-mean(SS_STA(:,1:15),2);
  figure(15); clf;
  nexttile([1 1])
  imagesc(SS_STA)
  imshow_patch(toimg(squeeze(max(tovec(Result.ftprnt(:,:,Result.dist_order(noi_dist))>0).*reshape(rescale(F0_PCA(Result.dist_order(noi_dist)))+0.1,1,sum(noi_dist),1),[],2)),size(Result.ref_im,1),size(Result.ref_im,2)))
  Result.F0_PCA=F0_PCA;
  save(fullfile(fpath{f},'PC_Result.mat'),'Result','-v7.3')
  load(fullfile(fpath{f},'PC_Result.mat'),'Result');
backupServer(fpath{f},'BHL18TB_D2','cohen_lab/Lab/Labmembers/Byung Hun Lee/Data','PC_Result.mat')
end