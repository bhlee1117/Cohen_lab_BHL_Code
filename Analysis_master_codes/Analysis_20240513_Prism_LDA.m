clear
clc;
cd '/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/FromBackup/PP72_PlaceCellResults';
[~, ~, raw] = xlsread(['/Volumes/cohen_lab/Lab/Labmembers' ...
    '/Byung Hun Lee/Data/PrismPCdata_Arrangement.xlsx'], 'Sheet1', 'C5:P23');

% [~, ~, NeuronsToUse]=xlsread(['/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/' ...
%     'PlaceCellData_Arrangement.xlsx'], 'Sheet1', 'L8:M46');
%
% NeuronsToUse=cellfun(@(x) (str2num(num2str(x))),NeuronsToUse,'UniformOutput',false);
ref_ROI=cellfun(@(x) (str2num(num2str(x))),raw(:,10),'UniformOutput',false);
basal_ROI=cellfun(@(x) (str2num(num2str(x))),raw(:,11),'UniformOutput',false);
apical_ROI=cellfun(@(x) (str2num(num2str(x))),raw(:,12),'UniformOutput',false);
fpath=raw(:,1)';
StructureData=raw(:,10);
EndFrame=cell2mat(raw(:,13));
ifmotionReject=cell2mat(raw(:,14));
save_figto='/Volumes/BHL_WD18TB/PP72_PlaceCellResults';
place_bin=150; time_segment=15000; overlap=200;
alignedMovFN = {'STA_Mat_SS','STA_Mat_CS','STA_Mat_dSP'};
bound=6;
save_at='/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Updates/2024/20240408_Movs_Figs/';
%%
f=18; load(fullfile(fpath{f},'PC_Result.mat'),'Result')
nROI=size(Result.normTraces,1);
nTau={[-55:15],[-60:150],[-20:20]}; %SS, CS, dSP
nTau_bAP=[-20:20];
% Isolated Somatic spike
som_spike=find(Result.spike(1,:));
tr_ref=Result.normTraces(ref_ROI{f},:);
tr_sub=mean(tr_ref,1)-movprc(mean(tr_ref,1),200,20,2);
tr_sub=get_subthreshold(tr_sub,Result.spike(1,:),5,10);
[trans tr_trace]=detect_transient2(tr_sub,[5 1.5],Result.spike(1,:),15);

% Calculate distance order
SkelDend = Skeletonize_dendrite(Result.ref_im,4,0.02,25);
interDendDist=[];
for i=1:size(Result.normTraces,1)
    for j=1:size(Result.normTraces,1)
        [interDendDist(i,j), ~]=geodesic_distance(SkelDend,get_coord(Result.ftprnt(:,:,i)),get_coord(Result.ftprnt(:,:,j)));
    end
end
coord_1d=dim_reduce(get_coord(Result.ftprnt));
[~, dist_order]=sort(coord_1d,'ascend');
som_roi=find(dist_order==1);
geodist=interDendDist(1,:)'.*sign(coord_1d-coord_1d(1));
show_footprnt_contour(Result.ftprnt(:,:,dist_order),Result.ref_im)

% Isolated Somatic spike
bAP_s=[];
for s=som_spike
    isnearby=sum(ismember(s+nTau{1},som_spike))>1;
    isnearbyCS=sum(ismember(s+nTau{1},find(Result.CStrace)))>1;
    ispartCS=tr_trace(s)>0;
    if ~isnearby & ~isnan(s) & ~isnearbyCS & ~ispartCS
        bAP_s=[bAP_s s];
    end
end

% Reference bAPs for F_ref
bAP_ref=[];
for s=som_spike
    isnearby=sum(ismember(s+nTau_bAP,som_spike))>1;
    isnearbyCS=sum(ismember(s+nTau_bAP,find(Result.CStrace)))>1;
    ispartCS=tr_trace(s)>0;
    if ~isnearby & ~isnan(s) & ~isnearbyCS & ~ispartCS
        bAP_ref=[bAP_ref s];
    end
end

% Isolated Complex spike
CS_s=[];
C_list=find(Result.SpClass(2,:));
CS_label=max(Result.spike,[],1).*bwlabel(Result.CStrace);
CS_list=[];
for g=1:length(C_list)
    s=C_list(g);
    s_tmp=Result.spike(1,:);
    s_tmp(find(CS_label==g))=0;
    isnearbyCS=max(bwlabel(Result.CStrace(s+nTau{2})))>1;
    isnearbyS=sum(ismember(s+nTau{2},find(s_tmp)))>0;
    if ~isnearbyCS & ~isnearbyS
        CS_s=[CS_s s];
        CS_list=[CS_list g];
    end
end

prc_normTr=Result.normTraces;
%prc_normTr=Result.normTraces-movprc(Result.normTraces,500,35,2);
STA_SSmat=reshape(prc_normTr(:,bAP_s'+nTau{1}),nROI,[],length(nTau{1}));
STA_SS=squeeze(mean(reshape(prc_normTr(:,bAP_ref'+nTau_bAP),nROI,[],length(nTau_bAP)),2));
F_ref=mean(STA_SS(:,-nTau_bAP(1)+[10:14]),2);
STA_SSmat=STA_SSmat./F_ref;
prc_normTrCS=Result.normTraces;
%prc_normTrCS=Result.normTraces-movprc(Result.normTraces,500,35,2);
STA_CSmat=reshape(prc_normTrCS(:,CS_s'+nTau{2}),nROI,[],length(nTau{2}));
STA_CSmat=STA_CSmat./F_ref;
CSpike=Result.spike.*Result.CStrace;
STA_CSpikemat=reshape(CSpike(:,CS_s'+nTau{2}),nROI,[],length(nTau{2}));

figure(2); clf;
cax=[-0.5 1.5];
tiledlayout(2,2)
nexttile([1 1])
imagesc(squeeze(mean(STA_SSmat(dist_order,:,:),2)),cax); hold all
plot(1,som_roi,'marker','>','MarkerFaceColor','r')
ylabel('ROI (distance order)')
title('STA of SS')
nexttile([1 1])
imagesc(squeeze(mean(STA_CSmat(dist_order,:,:),2)),cax); hold all
plot(1,som_roi,'marker','>','MarkerFaceColor','r')
ylabel('ROI (distance order)')
title('STA of CS')
colormap(turbo)
nexttile([1 1])
l=plot(nTau{1},squeeze(mean(STA_SSmat(dist_order,:,:),2))');
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(turbo(nROI),2))
xlabel('Peri-spike time (ms)')
nexttile([1 1])
l=plot(nTau{2},squeeze(mean(STA_CSmat(dist_order,:,:),2))');
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(turbo(nROI),2))
xlabel('Peri-spike time (ms)')
for n=1:nROI
    roi_l{n}=['ROI#' num2str(n)];
end
legend(roi_l,'FontSize',8)

figure(3); clf; cmap=distinguishable_colors(6); ax1=[];
rois={basal_ROI{f},apical_ROI{f}};
tiledlayout(2,2)
ax1=[ax1 nexttile([1 1])];
errorbar_shade(nTau{1},squeeze(mean(mean(STA_SSmat(rois{1},:,:),1),2))',squeeze(std(mean(STA_SSmat(rois{1},:,:),1),0,2))',cmap(1,:)); hold all
errorbar_shade(nTau{1},squeeze(mean(mean(STA_SSmat(rois{2},:,:),1),2))',squeeze(std(mean(STA_SSmat(rois{2},:,:),1),0,2))',cmap(2,:));
xlabel('Peri-spike time (ms)')
ax1=[ax1 nexttile([1 1])];
errorbar_shade(nTau{2},squeeze(mean(mean(STA_CSmat(rois{1},:,:),1),2))',squeeze(std(mean(STA_CSmat(rois{1},:,:),1),0,2))',cmap(1,:)); hold all
errorbar_shade(nTau{2},squeeze(mean(mean(STA_CSmat(rois{2},:,:),1),2))',squeeze(std(mean(STA_CSmat(rois{2},:,:),1),0,2))',cmap(2,:));
xlabel('Peri-spike time (ms)')
ax1=[ax1 nexttile([1 1])];
errorbar_shade(nTau{1},squeeze(mean(mean(STA_SSmat(rois{1},:,:),1),2))',squeeze(std(mean(STA_SSmat(rois{1},:,:),1),0,2))',cmap(3,:)); hold all
errorbar_shade(nTau{2},squeeze(mean(mean(STA_CSmat(rois{1},:,:),1),2))',squeeze(std(mean(STA_CSmat(rois{1},:,:),1),0,2))',cmap(6,:));
xlabel('Peri-spike time (ms)')
ax1=[ax1 nexttile([1 1])];
errorbar_shade(nTau{1},squeeze(mean(mean(STA_SSmat(rois{2},:,:),1),2))',squeeze(std(mean(STA_SSmat(rois{2},:,:),1),0,2))',cmap(3,:)); hold all
errorbar_shade(nTau{2},squeeze(mean(mean(STA_CSmat(rois{2},:,:),1),2))',squeeze(std(mean(STA_CSmat(rois{2},:,:),1),0,2))',cmap(6,:));
xlabel('Peri-spike time (ms)')
linkaxes(ax1,'xy')

%% 
cmap=distinguishable_colors(6);
prespike_time=50; %ms
sum_bin= 7; %ms
training_fraction=0.5; noi=setdiff([1:nROI],9);
LDA_CSmat=STA_CSmat(dist_order(noi),:,-nTau{2}(1)-prespike_time:-nTau{2}(1)-1);
LDA_SSmat=STA_SSmat(dist_order(noi),:,-nTau{1}(1)-prespike_time:-nTau{1}(1)-1);
LDA_CSmat=movsum(LDA_CSmat,sum_bin,3); LDA_SSmat=movsum(LDA_SSmat,sum_bin,3); 
LDA_SSmat=tovec(permute(LDA_SSmat(:,:,[round(sum_bin/2):sum_bin:end]),[1 3 2])); 
LDA_CSmat=tovec(permute(LDA_CSmat(:,:,[round(sum_bin/2):sum_bin:end]),[1 3 2]));

featureInput=[LDA_SSmat LDA_CSmat]';
featureClass=[zeros(1,size(STA_SSmat,2)) (zeros(1,size(STA_CSmat,2))+1)]';
trainingset=randi(size(featureClass,1),1,round(size(featureClass,1)*training_fraction));
predictset=setdiff([1:size(featureInput,1)],trainingset);
randomInput = median(featureInput) + 2*range(featureInput) .* (rand(5000,size(featureInput,2))-0.5);

mdl = fitcdiscr(featureInput(trainingset,:), featureClass(trainingset));
mlg = fitglm(featureInput(trainingset,:), featureClass(trainingset), 'Distribution', 'binomial', 'Link', 'logit');

covMat = featureInput'*featureInput;
[V, D] = eig(covMat);
D = diag(D);
D = D(end:-1:1);
V = V(:,end:-1:1);
vSign = sign(max(V) - max(-V));  % make the largest value always positive
V = V.*vSign;
figure(22); clf;
for i=1:20
nexttile([1 1]);
imagesc(reshape(V(:,i),length(noi),[]))
title(num2str(i));
end
colormap(turbo)

eigTraces=featureInput*V;
randTraces=randomInput*V;

coefficients = mdl.Coeffs(1,2).Linear;
%featureInput_reduced = featureInput * mdl.Coeffs(1,2).Linear;
featureInput_reduced = featureInput * mdl.Mu';

predictedLabels = predict(mlg, featureInput(predictset,:)) > 0.5;
predictedLabels2 = predict(mdl, featureInput) > 0.5;
randomLabels=predict(mdl, randomInput)>0.5;
accuracy = mean((predictedLabels) == featureClass(predictset))

figure(20); clf; cax=[-1 3];
nexttile([1 1])
imagesc(reshape(squeeze(mean(LDA_SSmat,2)),length(noi),[]),cax)
title('STA of SS')
nexttile([1 1])
imagesc(reshape(squeeze(mean(LDA_CSmat,2)),length(noi),[]),cax)
title('STA of CS')
nexttile([1 1])
imagesc(reshape(mdl.Mu(1,:),length(noi),[]),cax)
title('LDA mean for SS')
nexttile([1 1])
imagesc(reshape(mdl.Mu(2,:),length(noi),[]),cax)
title('LDA mean for CS')
colormap(turbo)

figure(19); clf;
ax3=nexttile([1 1]);
show_ind=[4 5 6]; cmap_random=[0.3 0.3 1;1 0.3 0.3];
scatter3(randTraces(randomLabels==0,show_ind(1)), randTraces(randomLabels==0,show_ind(2)),randTraces(randomLabels==0,show_ind(3)), 5, 'filled');hold all
scatter3(randTraces(randomLabels==1,show_ind(1)), randTraces(randomLabels==1,show_ind(2)),randTraces(randomLabels==1,show_ind(3)), 5, 'filled');hold all
xlabel(['PC#' num2str(show_ind(1))]);
ylabel(['PC#' num2str(show_ind(2))]);
zlabel(['PC#' num2str(show_ind(3))]);
title('Random Data'); grid on;
legend({'SS','CS'})
ax4=nexttile([1 1]);
scatter3(eigTraces(featureClass==0,show_ind(1)), eigTraces(featureClass==0,show_ind(2)), eigTraces(featureClass==0,show_ind(3)), 20, 'filled'); hold all
scatter3(eigTraces(featureClass==1,show_ind(1)), eigTraces(featureClass==1,show_ind(2)), eigTraces(featureClass==1,show_ind(3)), 20, 'filled'); hold all
xlabel(['PC#' num2str(show_ind(1))]);
ylabel(['PC#' num2str(show_ind(2))]);
zlabel(['PC#' num2str(show_ind(3))]);
title('LDA Data'); grid on;
legend({'SS','CS'})
linkprop([ax3, ax4], {'CameraPosition', 'CameraTarget', 'CameraUpVector', 'CameraViewAngle'});
%%
covMat = featureInput*featureInput';
[V, D] = eig(covMat);
D = diag(D);
D = D(end:-1:1);
V = V(:,end:-1:1);
vSign = sign(max(V) - max(-V));  % make the largest value always positive
V = V.*vSign;
eigImgs = (V'*featureInput);

figure(21); clf;
nexttile([1 1])
plot(cumsum(D)/sum(D))

figure(22); clf;
for i=1:20
nexttile([1 1]);
imagesc(reshape(eigImgs(i,:),length(noi),[]))
title(num2str(i));
end
colormap(turbo)

pc_ind=[1:8];
pairs=[1 2; 1 3; 1 4; 1 5; 1 6; 1 7; 2 3; 2 4; 2 5; 2 6; 2 7; 3 4; 3 5; 3 6;3 7;...
       4 5; 4 6; 4 7; 5 6; 5 7; 6 7; 7 8; 7 9; 7 10];
proj_val=V;
trainingset=randi(size(featureClass,1),1,round(size(featureClass,1)*training_fraction));
predictset=setdiff([1:size(featureInput,1)],trainingset);
pclda = fitcdiscr(proj_val(trainingset,pc_ind), featureClass(trainingset));

predictedLabels = predict(pclda, proj_val(predictset,pc_ind));
accuracy = mean(predictedLabels == featureClass(predictset))

pcldaImgs=pclda.Mu*eigImgs(pc_ind,:);

figure(23); clf; cax=[-1 3];
nexttile([1 1])
imagesc(reshape(squeeze(mean(LDA_SSmat,2)),length(noi),[]),cax)
title('STA of SS')
nexttile([1 1])
imagesc(reshape(squeeze(mean(LDA_CSmat,2)),length(noi),[]),cax)
title('STA of CS')
nexttile([1 1])
imagesc(reshape(pcldaImgs(1,:),length(noi),[]),cax)
title('LDA mean for SS')
nexttile([1 1])
imagesc(reshape(pcldaImgs(2,:),length(noi),[]),cax)
title('LDA mean for CS')
colormap(turbo)
% 
% 
% X_lda = proj_val(:,pc_ind) * pclda.Coeffs(1,2).Linear;
%     minX = min(X_lda);
%     maxX = max(X_lda);
% [x1Grid, x2Grid] = meshgrid(linspace(minX(1), maxX(1), 50), ...
%                            linspace(minX(2), maxX(2), 50));
% gridPoints = [x1Grid(:), x2Grid(:)];
% coeffInv = pinv(pclda.Coeffs(1,2).Linear');
% originalSpaceGridPoints = gridPoints * coeffInv;
% predictedLabels = predict(pclda, originalSpaceGridPoints);
% labelGrid = reshape(predictedLabels, size(x1Grid));
% 
% figure(23); clf;
% imagesc(linspace(minX(1), maxX(1), 50), linspace(minX(2), maxX(2), 50), labelGrid);
% set(gca, 'YDir', 'normal'); % Adjust the y-axis to display correctly
% colormap('jet'); hold all
% scatter(X_lda(:,1), X_lda(:,2), 10, cmap(featureClass,:), 'filled'); hold all
% 
