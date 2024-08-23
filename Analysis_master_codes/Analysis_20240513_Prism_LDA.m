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
sum_bin= 12; %ms
training_fraction=0.45; noi=setdiff([1:nROI],[8 9]);
LDA_CSmat=STA_CSmat(dist_order(noi),:,-nTau{2}(1)-prespike_time:-nTau{2}(1)-1);
LDA_SSmat=STA_SSmat(dist_order(noi),:,-nTau{1}(1)-prespike_time:-nTau{1}(1)-1);
LDA_CSmat=movsum(LDA_CSmat,sum_bin,3); LDA_SSmat=movsum(LDA_SSmat,sum_bin,3); 
LDA_SSmat=tovec(permute(LDA_SSmat(:,:,[round(sum_bin/2):sum_bin:end]),[1 3 2])); 
LDA_CSmat=tovec(permute(LDA_CSmat(:,:,[round(sum_bin/2):sum_bin:end]),[1 3 2]));
SS_nSet= size(LDA_SSmat,2); CS_nSet=size(LDA_CSmat,2);

featureInput=[LDA_SSmat LDA_CSmat]';
featureClass=[zeros(1,size(STA_SSmat,2)) (zeros(1,size(STA_CSmat,2))+1)]';
maxlda=0; maxlog=0;
for iter=1:1000
[nSet, minarg]=min([SS_nSet CS_nSet]);
randSS=randi(SS_nSet,1,nSet); randCS=randi(CS_nSet,1,nSet);
trainingset=[randSS(1:round(nSet*training_fraction)) randCS(1:round(nSet*training_fraction))+SS_nSet];
predictset=[randSS(round(nSet*training_fraction)+1:end) randCS(round(nSet*training_fraction)+1:end)+SS_nSet];
randomInput = median(featureInput) + 2*range(featureInput) .* (rand(5000,size(featureInput,2))-0.5);

mdl = fitcdiscr(featureInput(trainingset,:), featureClass(trainingset));
mlg = fitglm(featureInput(trainingset,:), featureClass(trainingset), 'Distribution', 'binomial', 'Link', 'logit');

coefficients = mdl.Coeffs(1,2).Linear;
featureInput_reduced = featureInput * mdl.Mu';

predictedLabels_lda = predict(mdl, featureInput(predictset,:)) > 0.5;
predictedLabels_logit = predict(mlg, featureInput(predictset,:)) > 0.5;
%randomLabels=predict(mdl, randomInput)>0.5;
accuracy_lda(iter) = mean((predictedLabels_lda) == featureClass(predictset));
accuracy_logit(iter) = mean((predictedLabels_logit) == featureClass(predictset));

if accuracy_lda(iter)>maxlda
max_mdl=mdl;
end

if accuracy_logit(iter)>maxlog
max_mlg=mlg;
end

maxlda=max(accuracy_lda);
maxlog=max(accuracy_logit);
end
%%

figure(20); clf; cax=[-1 5];
[pattcorr, p]=corr(max_mdl.Mu(1,:)',max_mdl.Mu(2,:)');
nexttile([1 1])
imagesc(reshape(squeeze(mean(LDA_SSmat,2)),length(noi),[]),cax)
title('STA of SS')
nexttile([1 1])
imagesc(reshape(squeeze(mean(LDA_CSmat,2)),length(noi),[]),cax)
title('STA of CS')
nexttile([1 1])
imagesc(reshape(max_mdl.Mu(1,:),length(noi),[]),cax)
title('LDA mean for SS')
nexttile([1 1])
imagesc(reshape(max_mdl.Mu(2,:),length(noi),[]),cax)
title('LDA mean for CS')
nexttile([1 1])
imagesc(reshape(max_mdl.Mu(1,:),length(noi),[])-reshape(squeeze(mean(LDA_SSmat,2)),length(noi),[]))
title('LDA mean - STA for SS')
colorbar
nexttile([1 1])
imagesc(reshape(max_mdl.Mu(2,:),length(noi),[])-reshape(squeeze(mean(LDA_CSmat,2)),length(noi),[]))
title('LDA mean - STA for CS')
colorbar
colormap(turbo)
nexttile([1 1])
imagesc(reshape(table2array(max_mlg.Coefficients(2:end,1)),length(noi),[]))
title('Weight vector for Logistic Regression')
colormap(turbo)
nexttile([1 1])
imagesc(reshape(max_mdl.Coeffs(1,2).Linear,length(noi),[]))
title('Weight vector from LDA')
colormap(turbo)

figure(21); clf; 
nexttile([1 1])
h=histogram(accuracy_lda,30);
xlabel('Accuracy')
ylabel('# of trial')
title("Linear discrimination analysis")
nexttile([1 1])
h=histogram(accuracy_logit,30);
xlabel('Accuracy')
ylabel('# of trial')
title("Logistic Regression")

mdl_project=featureInput*max_mdl.Mu';
mdl_Corr=corr(featureInput',max_mdl.Mu');
mdl_Corr_rand=corr(randomInput',max_mdl.Mu');
mdl_project_rand=randomInput*max_mdl.Mu';
predictedLabels_randlda = predict(max_mdl, randomInput) > 0.5;
figure(30); clf; ax1=[];
ax1=[ax1 nexttile([1 1])];
h=histogram(mdl_Corr(1:SS_nSet,1),30,'Normalization','probability'); hold all
histogram(mdl_Corr(SS_nSet+1:end,1),h.BinEdges,'Normalization','probability')
hold all
plot([pattcorr pattcorr],[0 0.2],'r')
xlabel('Corrlation coefficient')
title('Correlation btw pattern #1 and SS/CS')
legend({'SS','CS'})
ax1=[ax1 nexttile([1 1])];
histogram(mdl_Corr(1:SS_nSet,2),h.BinEdges,'Normalization','probability'); hold all
histogram(mdl_Corr(SS_nSet+1:end,2),h.BinEdges,'Normalization','probability')
plot([pattcorr pattcorr],[0 0.2],'r')
xlabel('Corrlation coefficient')
legend({'SS','CS'})
title('Correlation btw pattern #2 and SS/CS')
ax1=[ax1 nexttile([1 1])];
histogram(mdl_Corr(1:SS_nSet,1),h.BinEdges,'Normalization','probability'); hold all
histogram(mdl_Corr(1:SS_nSet,2),h.BinEdges,'Normalization','probability')
plot([pattcorr pattcorr],[0 0.2],'r')
xlabel('Corrlation coefficient')
legend({'Patt #1','Patt #2'})
title('Correlation btw SS and pattern')
ax1=[ax1 nexttile([1 1])];
histogram(mdl_Corr(SS_nSet+1:end,1),h.BinEdges,'Normalization','probability'); hold all
histogram(mdl_Corr(SS_nSet+1:end,2),h.BinEdges,'Normalization','probability')
plot([pattcorr pattcorr],[0 0.2],'r')
xlabel('Corrlation coefficient')
legend({'Patt #1','Patt #2'})
title('Correlation btw CS and pattern')
linkaxes(ax1,'x')

%%

cmap=distinguishable_colors(6);
prespike_time=50; %ms
sum_bin= 5; %ms
training_fraction=0.45; noi=setdiff([1:nROI],[8 9]);
LDA_CSmat=STA_CSmat(dist_order(noi),:,-nTau{2}(1)-prespike_time:-nTau{2}(1)-1);
LDA_SSmat=STA_SSmat(dist_order(noi),:,-nTau{1}(1)-prespike_time:-nTau{1}(1)-1);
LDA_CSmat=movsum(LDA_CSmat,sum_bin,3); LDA_SSmat=movsum(LDA_SSmat,sum_bin,3); 
LDA_SSmat=tovec(permute(LDA_SSmat(:,:,[round(sum_bin/2):sum_bin:end]),[1 3 2])); 
LDA_CSmat=tovec(permute(LDA_CSmat(:,:,[round(sum_bin/2):sum_bin:end]),[1 3 2]));
SS_nSet= size(LDA_SSmat,2); CS_nSet=size(LDA_CSmat,2);

featureInput=[LDA_SSmat LDA_CSmat]';
featureClass=[zeros(1,size(STA_SSmat,2)) (zeros(1,size(STA_CSmat,2))+1)]';

covMat = featureInput*featureInput';
[V, D] = eig(covMat);
D = diag(D);
D = D(end:-1:1);
V = V(:,end:-1:1);
vSign = sign(max(V) - max(-V));  % make the largest value always positive
V = V.*vSign;
eigImgs = (V'*featureInput);

% figure(21); clf;
% nexttile([1 1])
% plot(cumsum(D)/sum(D))

figure(22); clf;
for i=1:20
nexttile([1 1]);
imagesc(reshape(eigImgs(i,:),length(noi),[]))
title(['PC #' num2str(i) ' , Fraction: ' num2str(D(i)/sum(D),2)]);
end
colormap(turbo)

pc_ind=[1:20];
pairs=[1 2; 1 3; 1 4; 1 5; 1 6; 1 7; 2 3; 2 4; 2 5; 2 6; 2 7; 3 4; 3 5; 3 6;3 7;...
       4 5; 4 6; 4 7; 5 6; 5 7; 6 7; 7 8; 7 9; 7 10];
featureInputPCA=V(:,pc_ind);
SS_nSet= size(LDA_SSmat,2); CS_nSet=size(LDA_CSmat,2);

featureClass=[zeros(1,size(STA_SSmat,2)) (zeros(1,size(STA_CSmat,2))+1)]';
maxldaPCA=0; maxlogPCA=0;
for iter=1:1000

[nSet, minarg]=min([SS_nSet CS_nSet]);
randSS=randi(SS_nSet,1,nSet); randCS=randi(CS_nSet,1,nSet);
trainingset=[randSS(1:round(nSet*training_fraction)) randCS(1:round(nSet*training_fraction))+SS_nSet];
predictset=[randSS(round(nSet*training_fraction)+1:end) randCS(round(nSet*training_fraction)+1:end)+SS_nSet];
randomInput = median(featureInputPCA) + 2*range(featureInputPCA) .* (rand(5000,size(featureInputPCA,2))-0.5);

mdl = fitcdiscr(featureInputPCA(trainingset,:), featureClass(trainingset));
mlg = fitglm(featureInputPCA(trainingset,:), featureClass(trainingset), 'Distribution', 'binomial', 'Link', 'logit');

coefficients = mdl.Coeffs(1,2).Linear;
featureInput_reduced = featureInputPCA * mdl.Mu';

predictedLabels_lda = predict(mdl, featureInputPCA(predictset,:)) > 0.5;
predictedLabels_logit = predict(mlg, featureInputPCA(predictset,:)) > 0.5;
%randomLabels=predict(mdl, randomInput)>0.5;
accuracy_ldaPCA(iter) = mean((predictedLabels_lda) == featureClass(predictset));
accuracy_logitPCA(iter) = mean((predictedLabels_logit) == featureClass(predictset));

if accuracy_ldaPCA(iter)>maxldaPCA
max_mdlPCA=mdl;
end

if accuracy_logitPCA(iter)>maxlogPCA
max_mlgPCA=mlg;
end

maxldaPCA=max(accuracy_ldaPCA);
maxlogPCA=max(accuracy_logitPCA);
end


%%

figure(20); clf; cax=[-1 5];
LDAmu2img=max_mdlPCA.Mu*eigImgs(pc_ind,:);
nexttile([1 1])
imagesc(reshape(squeeze(mean(LDA_SSmat,2)),length(noi),[]),cax)
title('STA of SS')
nexttile([1 1])
imagesc(reshape(squeeze(mean(LDA_CSmat,2)),length(noi),[]),cax)
title('STA of CS')
nexttile([1 1])
imagesc(reshape(LDAmu2img(1,:),length(noi),[]),cax)
title('LDA mean for SS')
nexttile([1 1])
imagesc(reshape(LDAmu2img(2,:),length(noi),[]),cax)
title('LDA mean for CS')
nexttile([1 1])
imagesc(reshape(LDAmu2img(1,:),length(noi),[])-reshape(squeeze(mean(LDA_SSmat,2)),length(noi),[]))
title('LDA mean - STA for SS')
colorbar
nexttile([1 1])
imagesc(reshape(LDAmu2img(2,:),length(noi),[])-reshape(squeeze(mean(LDA_CSmat,2)),length(noi),[]))
title('LDA mean - STA for CS')
colorbar
colormap(turbo)
nexttile([1 1])
imagesc(reshape(table2array(max_mlg.Coefficients(2:end,1)),length(noi),[]))
title('Weight vector for Logistic Regression')
colormap(turbo)
nexttile([1 1])
imagesc(reshape(max_mdl.Coeffs(1,2).Linear,length(noi),[]))
title('Weight vector for Logistic Regression')
colormap(turbo)

figure(21); clf; 
nexttile([1 1])
h=histogram(accuracy_ldaPCA,30);
xlabel('Accuracy')
ylabel('# of trial')
title("Linear discrimination analysis")
nexttile([1 1])
h=histogram(accuracy_logitPCA,30);
xlabel('Accuracy')
ylabel('# of trial')
title("Logistic Regression")

mdl_projectPCA=featureInputPCA*max_mdlPCA.Mu';
mdl_CorrPCA=corr(featureInputPCA',max_mdlPCA.Mu');
[pattcorr, p]=corr(max_mdlPCA.Mu(1,:)',max_mdlPCA.Mu(2,:)');
figure(30); clf; ax1=[];
ax1=[ax1 nexttile([1 1])];
h=histogram(mdl_CorrPCA(1:SS_nSet,1),30,'Normalization','probability'); hold all
histogram(mdl_CorrPCA(SS_nSet+1:end,1),h.BinEdges,'Normalization','probability')
hold all
plot([pattcorr pattcorr],[0 0.2],'r')
xlabel('Corrlation coefficient')
title('Correlation btw pattern #1 and SS/CS')
legend({'SS','CS'})
ax1=[ax1 nexttile([1 1])];
histogram(mdl_CorrPCA(1:SS_nSet,2),h.BinEdges,'Normalization','probability'); hold all
histogram(mdl_CorrPCA(SS_nSet+1:end,2),h.BinEdges,'Normalization','probability')
plot([pattcorr pattcorr],[0 0.2],'r')
xlabel('Corrlation coefficient')
legend({'SS','CS'})
title('Correlation btw pattern #2 and SS/CS')
ax1=[ax1 nexttile([1 1])];
histogram(mdl_CorrPCA(1:SS_nSet,1),h.BinEdges,'Normalization','probability'); hold all
histogram(mdl_CorrPCA(1:SS_nSet,2),h.BinEdges,'Normalization','probability')
plot([pattcorr pattcorr],[0 0.2],'r')
xlabel('Corrlation coefficient')
legend({'Patt #1','Patt #2'})
title('Correlation btw SS and pattern')
ax1=[ax1 nexttile([1 1])];
histogram(mdl_CorrPCA(SS_nSet+1:end,1),h.BinEdges,'Normalization','probability'); hold all
histogram(mdl_CorrPCA(SS_nSet+1:end,2),h.BinEdges,'Normalization','probability')
plot([pattcorr pattcorr],[0 0.2],'r')
xlabel('Corrlation coefficient')
legend({'Patt #1','Patt #2'})
title('Correlation btw CS and pattern')
linkaxes(ax1,'x')

pc_pair=[1 5;1 7;1 13;1 15;1 20;5 7;5 13;5 15;5 20;7 13;7 15;7 20;13 15;13 20;15 20];
figure(31); clf;
for p=1:size(pc_pair,1)
    nexttile([1 1])
    plot(V(1:SS_nSet,pc_pair(p,1)),V(1:SS_nSet,pc_pair(p,2)),'.'); hold all
    plot(V(SS_nSet+1:end,pc_pair(p,1)),V(SS_nSet+1:end,pc_pair(p,2)),'.'); hold all
    xlabel(['PC #' num2str(pc_pair(p,1))])
    ylabel(['PC #' num2str(pc_pair(p,2))])
end
legend({'SS','CS'});