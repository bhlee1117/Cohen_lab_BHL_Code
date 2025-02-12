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

for f=foi(21)
load(fullfile(fpath{f},'PC_Result.mat'),'Result');
SkelDend = Skeletonize_dendrite(Result.ref_im,10,0.02,10);
    Result.BranchLabel = interactive_ROIlabel(Result.ftprnt,Result.ref_im);
    Result.BrachDistfromSoma=[];
    Result.BrancDist_order=[];
    dendmask=[];
    for d=1:max(Result.BranchLabel)
        dendmask(:,:,d)=max(Result.ftprnt(:,:,Result.BranchLabel==d),[],3)>0;
[Result.BrachDistfromSoma(d), path]=geodesic_distance(SkelDend,get_coord(Result.ftprnt(:,:,1)),get_coord(dendmask(:,:,d)));       
    end

    coord_1d=dim_reduce(get_coord(dendmask));
if coord_1d(1)<median(coord_1d)
    [~, denddistOrder]=sort(coord_1d,'ascend');

    for d2=denddistOrder'
        d2ind=find(Result.BranchLabel==d2);

        [~, d2ind_sort]=sort(Result.interDendDist(1,d2ind),'ascend');
        Result.BrancDist_order=[Result.BrancDist_order; d2ind(d2ind_sort)'];
    end
else
    [~, denddistOrder]=sort(coord_1d,'descend');
    for d2=denddistOrder'
        d2ind=find(Result.BranchLabel==d2);

        [~, d2ind_sort]=sort(Result.interDendDist(1,d2ind),'ascend');
        Result.BrancDist_order=[Result.BrancDist_order; d2ind(d2ind_sort)'];
    end
end

figure(f); clf;
show_footprnt_contour(Result.ftprnt(:,:,Result.BrancDist_order),Result.ref_im);
drawnow;

save(fullfile(fpath{f},'PC_Result.mat'),'Result','-v7.3');
pause(3);
load(fullfile(fpath{f},'PC_Result.mat'),'Result');
backupServer(fpath{f},'BHL18TB_D2','cohen_lab/Lab/Labmembers/Byung Hun Lee/Data','PC_Result.mat')
end

%%
figure(10); clf;
for f=foi(21:22)
    f
load(fullfile(fpath{f},'PC_Result.mat'),'Result');

perisomaROI=setdiff(find(Result.interDendDist(1,:)<40),BadROI{f}); % ROIs < 40 um from soma
NormalizedTrace=(Result.normTraces)./Result.F_ref;
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

nTime=size(Result.normTraces,2);
perispike_frame=unique([tovec(find(double(Result.spike(1,:)==1))'+[-5:20]); find(Result.CStrace)']);
perispike_frame(perispike_frame<=0 | perispike_frame>nTime)=[];
nonvalid_frame=find(sum(isnan(NormalizedTrace_dirt{f}),1)>0);
Blue_on_frame=find(imdilate(Result.Blue>0, [ones(1, 1), 1, ones(1, 200)]));
Badframe=unique([perispike_frame' Blue_on_frame nonvalid_frame]);
Goodframe=setdiff([1:nTime],Badframe);
sub_ch=[];
for ch=1:2
sub=get_subthreshold( NormalizedTrace_ch{f,ch},max(Result.spike(1,:),[],1)>0,7,17);
for d=1:max(Result.BranchLabel)
    dind=find(Result.BranchLabel==d);
    dind=setdiff(dind,BadROI{f});
   sub_ch{ch}(d,:)=mean(sub(dind,:),1,'omitnan');
end
end

corrMat=get_corrMat(sub_ch{1},sub_ch{2},Goodframe);

dendmask=[];

for d=1:max(Result.BranchLabel)
    dendmask(:,:,d)=max(Result.ftprnt(:,:,find(Result.BranchLabel==d)),[],3);
end
l=nexttile([1 1]);
show_footprnt_contour(dendmask,Result.ref_im);
title(f);
%colormap(l,gray)

l2=nexttile([1 1]);
imagesc(corrMat); axis equal tight;
title(f);
%colormap(l2,turbo);
colormap(turbo(256))

 drawnow;
end
