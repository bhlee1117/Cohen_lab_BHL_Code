
clear
clc;
[~, ~, raw] = xlsread(['/Volumes/cohen_lab/Lab/Labmembers' ...
    '/Byung Hun Lee/Data/PrismPCdata_Arrangement.xlsx'], 'Sheet1', 'C5:AA31');

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
NSeesawComponent=cell2mat(raw(:,25));
save_figto='/Volumes/BHL18TB_D2/PP72_PlaceCellResults';
place_bin=150; time_segment=15000; overlap=200;
alignedMovFN = {'STA_Mat_SS','STA_Mat_CS','STA_Mat_dSP'};
bound=6;
title_str={'Basal','Apical','Peri-Soma'};
PlaceFieldList=cellfun(@(x) (str2num(num2str(x))),raw(:,21),'UniformOutput',false);
PlaceFieldBin=cellfun(@(x) (str2num(num2str(x))),raw(:,22),'UniformOutput',false);
set(0,'DefaultFigureWindowStyle','docked')
%foi=[1 4 5 6 8 10 11 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27];
foi=[1 4 5 6 8 10 11 15 16 17 18 19 20 21 22 23 24 25 26 27];
%foi=23;

%% Load parameter
win_pre = 50; win_post = 150;
time_segment=15000; bound=5; overlap=200;
nTauPeak=[50 150];
for f = [foi]
    f
    load(fullfile(fpath{f}, 'PC_Result.mat'), 'Result')
    load(fullfile(fpath{f},'RobustdFFfit.mat'))

    if ~isempty(find(ismember(find(max(Result.SpClass([1 2],:),[],1)>0),find(Result.spike(1,:)>0))==0))
        diffind=find(ismember(find(max(Result.SpClass([1 2],:),[],1)>0),find(Result.spike(1,:)>0))==0);
        ss=find(max(Result.SpClass([1 2],:),[],1)>0);
        disp([num2str(ss(diffind)) ' th frame is not matching']);
        error('Spclass and spike mat is not matching');
    end

    Result.spike = Result.spike > 0;
    Result.SpClass = Result.SpClass > 0;
    PCresult{f}.Dist_order = Result.BrancDist_order;
    PCresult{f}.BranchLabel = Result.BranchLabel;
    PCresult{f}.interDendDist = Result.interDendDist * Pixelsize(f);
    PCresult{f}.pixelsize = Pixelsize(f);
    PCresult{f}.swc=Result.SWC;
    PCresult{f}.Blue=Result.Blue;
    PCresult{f}.traces_bvMask=Result.traces_bvMask;
    PCresult{f}.SpikeHeight_fit=Result.SpikeHeight_fit;
    PCresult{f}.ftprnt=Result.ftprnt;
    PCresult{f}.ref_im=Result.ref_im;
    PCresult{f}.bvMask=Result.bvMask;
    PCresult{f}.normTraces=Result.normTraces;
    PCresult{f}.motionReject=Result.motionReject;
    PCresult{f}.F0_PCA=Result.F0_PCA;
    PCresult{f}.Scalefactor=RobustdFFfit.ScaleFactor;

    Dsign = ones(1, size(PCresult{f}.interDendDist, 1));
    Dsign(PCresult{f}.Dist_order(1:find(PCresult{f}.Dist_order == 1) - 1)) = -1;
    perisomaROI = setdiff(find(PCresult{f}.interDendDist(1, :) < 40), BadROI{f});
    noi = setdiff(1:size(Result.ftprnt, 3), BadROI{f});
    PCresult{f}.noi_dist = ismember(PCresult{f}.Dist_order, noi);

    PCresult{f}.sortDist=PCresult{f}.Dist_order(PCresult{f}.noi_dist);
    dendaxis = PCresult{f}.interDendDist(1,:) .* Dsign;
    PCresult{f}.dendaxis = dendaxis;% dendaxis(PCresult{f}.sortDist);

    roisD(f,:)={basal_ROI{f},PeriSoma_ROI{f},apical_ROI{f},oblique_ROI{f},distal_ROI{f}}; %set the ROIs
    for dClass=1:size(roisD,2) % 1:basal, 2:soma, 3:trunk, 4:obliqu, 5:distal
        g=1;
        if ~isnan(roisD{f,dClass})
            for d=roisD{f,dClass} %branch label
                dind=setdiff(find(Result.BranchLabel==d),BadROI{f});
                PCresult{f}.roisD_order{dClass,g}=ismember(PCresult{f}.sortDist,dind)*d;
                PCresult{f}.roi_dClass(ismember(PCresult{f}.sortDist,dind))=dClass;
                PCresult{f}.branch_dClass(ismember(PCresult{f}.sortDist,dind))=dClass;
                g=g+1;
            end
        end
    end

    NormalizedTrace = Result.normTraces ./ Result.F0_PCA;
    [bAP_STA MatbAP]= get_STA(NormalizedTrace, Result.spike(1,:), win_pre, win_post);
    bAP_STA = bAP_STA - prctile(bAP_STA, 30, 2);
    PCresult{f}.MatbAP=permute(MatbAP(PCresult{f}.sortDist,:,:),[1 3 2]);
    [PCresult{f}.SpikeHeight] = max(mean(bAP_STA(perisomaROI,:), 1), [], 2);
    [~, PCresult{f}.SpikeDelay] = max(bAP_STA, [], 2);
    PCresult{f}.SpikeDelay = PCresult{f}.SpikeDelay - (win_pre + 1);

    NormalizedTrace_dirt = NormalizedTrace;
    NormalizedTrace_dirt(:, Result.motionReject > 0) = NaN;
    NormTrCheck = cellfun(@(x) x ./ Result.F0_PCA, Result.norm_trace_check, 'UniformOutput', false);
    NormTrCheck{1}(:, Result.motionReject > 0) = NaN;
    NormTrCheck{2}(:, Result.motionReject > 0) = NaN;
    if ifdirtRemov(f)
        PCresult{f}.dirtTrace=Result.dirtTrace;
        NormalizedTrace_dirt(Result.dirtTrace > 0) = NaN;
        NormTrCheck{1}(Result.dirtTrace > 0) = NaN;
        NormTrCheck{2}(Result.dirtTrace > 0) = NaN;
    end
    NormalizedTrace_dirt = NormalizedTrace_dirt(PCresult{f}.sortDist, :);
    NormTrCheck{1} = NormTrCheck{1}(PCresult{f}.sortDist, :);
    NormTrCheck{2} = NormTrCheck{2}(PCresult{f}.sortDist, :);
    [nROI nTime]=size(NormalizedTrace_dirt);

    PCresult{f}.NormalizedTrace_dirt = NormalizedTrace_dirt;
    PCresult{f}.NormalizedTrace_ch = NormTrCheck;
    PCresult{f}.Subthreshold = get_subthreshold(NormalizedTrace_dirt, max(Result.spike(1,:), [], 1) > 0, 7, 17);
    PCresult{f}.Subthreshold(isnan(NormalizedTrace_dirt)) = NaN;

    sptime = find(Result.spike(1,:));
    brst = bwlabel((sptime(2:end) - sptime(1:end-1)) < 20);
    SpClass = Result.SpClass;
    BS_trace = zeros(1, nTime);
    CS_trace= Result.CStrace;
    bwCStrace= bwlabel(CS_trace);
    CS_fistSp=find(SpClass(2,:)>0);
    if any(CS_trace(CS_fistSp)==0)
        CS_trace(CS_fistSp(find(CS_trace(CS_fistSp)==0))'+[0:50])=1;
    end

    SubthInterp=interpolateNaN(mean(PCresult{f}.Subthreshold(roisD{f,2},:),1,'omitnan'))./PCresult{f}.SpikeHeight; %soma subthreshold
    bwCScand=bwlabel(SubthInterp>0.45); % Subthreshold level higher than 0.3 spike height is included to complex spike interval
    CS_trace_new=CS_trace;
    for b=1:max(bwCStrace)
        candInd=unique(bwCScand(bwCStrace==b));
        candInd(candInd==0)=[];
        CS_trace_new(ismember(bwCScand,candInd))=1;
    end
    % brst_fistSp=[];
    % for b=1:max(brst)
    % brst_fistSp(b)=sptime(find(brst==b,1));
    % end
    CS_trace_new(CS_fistSp-1)=0;    %do not connect two CSs
    %CS_trace_new(brst_fistSp-1)=0;    %do not connect CSs and BSs

    bwCStrace= bwlabel(CS_trace_new);
    for b=1:max(bwCStrace)
        if any(SpClass(2,bwCStrace==b)>0)
        else
            CS_trace_new(bwCStrace==b)=0; %remove the CS traces without spikes
        end
    end
    CS_trace=CS_trace_new;
    bwCStrace= bwlabel(CS_trace);

    for b = 1:max(brst)
        bind=find(brst == b);
        bwn = sptime([bind]);
        if any(CS_trace(bwn)) % part of CS?
            CSind=unique(bwCStrace(bwn));
            CSind(CSind==0)=[];
            if length(CSind)>1
                error(['Spike train in more than one CS, check frame #' num2str(bwn)]);
            end
            if SpClass(2,bwn(1))==0 & CS_trace(bwn(1))==0 % the first spike is not assigned to CS?
                originalCSinterval=find(bwCStrace==CSind);
                SpClass(:,bwn(1):originalCSinterval(end))=0;
                SpClass(2,bwn(1))=1;
                CS_trace(bwn(1):originalCSinterval(end))=1;
            end
        else
            SpClass(1, bwn(1):sptime(bind(end)+1)) = 0;
            SpClass(4, bwn(1)) = 1; % assign first spike of BS
            BS_trace(1, [bwn(1):sptime(bind(end)+1)+7]) = b;
        end
    end

    BS_trace=BS_trace(1:nTime); CS_trace=CS_trace(1:nTime);
    SpClass = SpClass([1 2 4 3], :); %SS, CS, BS ,dSP
    Classvec = get_Class2index(SpClass);
    SpikeClassvec = Classvec .* Result.spike(1,:);

    PCresult{f}.SpikeMat = double(Result.spike>0);
    PCresult{f}.SpikeMat(:, Result.motionReject > 0) = NaN;
    if ifdirtRemov(f)
        PCresult{f}.SpikeMat(Result.dirtTrace > 0) = NaN;
    end
    PCresult{f}.SpikeClassMat = SpClass; % first spike of class
    PCresult{f}.SpikeClassVecMat = SpikeClassvec; % all the spikes in class

    BrstOrder = get_BurstOrder(Result.spike(1,:), 20) - SpClass(1,:);
    BrstOrder(SpClass(3,:) > 0) = 1;
    ComplexSpikeOrder = get_spikeOrder(CS_trace, Result.spike(1,:));
    PCresult{f}.SpikeOrder=max([BrstOrder; ComplexSpikeOrder; SpClass(1,:)],[],1);
end
%% Compare F0 and Fstd
figure(20); clf;
Corr_F0FPCA=NaN(max(foi),2);
for f=23%foi
    nexttile([1 1])
    F0img=PCResult{f}.ref_im-100;
    dMask=point2img(PCResult{f}.SWC(:,1:2),2,size(F0img));
    F0img_masked=maskBinary(F0img,dMask,NaN);
    F0img_masked=medfilt2nan(F0img_masked,[13 13]);
    F0img_sub=F0img-F0img_masked;
    F0img_sub(isnan(F0img_sub))=0;
    F0_sub=tovec(F0img_sub)'*tovec(PCResult{f}.ftprnt)./sum(tovec(PCResult{f}.ftprnt),1);
    F0=tovec(F0img)'*tovec(PCResult{f}.ftprnt)./sum(tovec(PCResult{f}.ftprnt),1);
    Dist_order=PCResult{f}.dist_order;
    scatter(F0(Dist_order),PCResult{f}.F0_PCA(Dist_order),30,turbo(length(F0)),'filled'); 
    set(gca,'XScale','log','YScale','log');
    xlabel('F_0');
    ylabel('F_{std}');
    title(f);
    Corr_F0FPCA(f,2)=corr(F0_sub',PCResult{f}.F0_PCA);
    Corr_F0FPCA(f,1)=corr(F0',PCResult{f}.F0_PCA);
end
figure(21); clf;
show_footprnt_contour(PCResult{f}.ftprnt(:,:,Dist_order),PCResult{f}.ref_im);

%% Show example of std image
f2show=23; gaussiansz=1;
frm2read=[52000:63000];
dMask2=point2img(PCResult{f2show}.SWC(:,1:2),4,size(PCResult{f}.ref_im));

if isfile(fullfile(fpath{f},'F_std_img.mat'))
   load(fullfile(fpath{f},'F_std_img.mat'));
else
mov_mc=readBinMov_BHL_multiple(fpath{f2show},3,frm2read, 15000);
meanF=squeeze(mean(mov_mc,[1 2]));
y_fit=expfitDM_2([1:size(mov_mc,3)]',meanF,[1:size(mov_mc,3)]',10000);
y_fit=y_fit./mean(y_fit);
mc=PCResult{f2show}.mcTrace(frm2read,:);
sz=size(mov_mc);

mov_mc=mov_mc./reshape(y_fit,1,1,[]);
mov_res=mov_mc-median(mov_mc,3);
%mov_res=SeeResiduals(mov_res,y_fit);
mov_res=SeeResiduals(mov_res,mc);
mov_res=SeeResiduals(mov_res,mc.^2);
mov_res=SeeResiduals(mov_res,mc(:,1).*mc(:,end));

mov_res_sub=get_subthreshold(tovec(mov_res),PCResult{f2show}.spike(1,frm2read),7,20);
mov_res_sub=toimg(mov_res_sub,size(mov_res,1),size(mov_res,2));
mov_res_sub=maskEdge(mov_res_sub,5,0);

DirtExist=sum(PCResult{f}.dirtTrace(:,frm2read),1);
[F0_img_PCA VarianceExplained]=get_F0img_PCA_mask(imgaussfilt(mov_res_sub,gaussiansz),dMask2,find(DirtExist==0));
F0_img_F0=imgaussfilt(F0img,gaussiansz).*dMask2;
F0_std_mask=std(imgaussfilt(mov_res_sub(:,:,find(DirtExist==0)),gaussiansz),0,3).*dMask2;

[V D E]=get_eigvector(tovec(imgaussfilt(mov_res_sub(:,:,find(DirtExist==0)),gaussiansz))');
E_sub=reshape(E(:,1:10),sz(1),sz(2),[]);
figure(24); clf;
imshow2_patch(E_sub);
clear mov_res_sub mov_res mov_mc
save(fullfile(fpath{f},'F_std_img.mat'),'F0_img_PCA','VarianceExplained','D','F0_img_F0','F0_std_mask','E_sub','-v7.3');
end

figure(23); clf;
plot(D(1:30)/D(1),'color',[0 0.1 1],'linewidth',1.5); box off; axis tight; ylim([0 1]);
set(gca,'linewidth',1);
xlabel('Eigenvalue number'); ylabel('Normalized variance');
set_fontsize(18); set_font('Arial');
set_figsize(200, 70);

%% Compare F0 vs F_std vs std image
figure(25); tiledlayout(3,1,'padding','tight');
nexttile([1 1])
imshow2(F0_img_F0,[min(F0_img_F0(dMask))*3 max(F0_img_F0(dMask))*0.7])
nexttile([1 1])
imshow2(F0_img_PCA,[min(F0_img_PCA(dMask))*3 max(F0_img_PCA(dMask))*0.7])
nexttile([1 1])
imshow2(F0_std_mask,[min(F0_std_mask(dMask))*3 max(F0_std_mask(dMask))*0.7])

figure(26); clf;
scatter(F0_img_F0(dMask),F0_img_PCA(dMask),10,[0 0 0],'filled');
[R p]=corr(F0_img_F0(dMask),F0_img_PCA(dMask));
title(sprintf('R = %0.2f',R))
xlabel('F_0'); ylabel('F_{std}');

%%
STAdecay_cat=[];
showf=[23]; boi=[3 13 31];
RobustdFF_decay=[];
for f=foi
    SS_list=find(max(PCresult{f}.SpikeClassMat([1 2],:).*(PCresult{f}.Blue==0),[],1)>0);
    STA_tr=get_STA(PCresult{f}.normTraces(PCresult{f}.sortDist,:),SS_list,win_pre,win_post);
    STA_tr=STA_tr-median(STA_tr(:,1:win_pre/2),2);
    F0=tovec(PCresult{f}.ref_im-100)'*tovec(PCresult{f}.ftprnt);
    F0=F0(PCresult{f}.sortDist)';
    NormConst=get_threshold(PCresult{f}.traces_bvMask(PCresult{f}.sortDist,:),1);
    STAdecay_cat{f}={PCresult{f}.dendaxis(PCresult{f}.sortDist)', STA_tr, PCresult{f}.F0_PCA(PCresult{f}.sortDist), PCresult{f}.Scalefactor(PCresult{f}.sortDist), F0};
    fieldName={'dendaxis','STA','F0_PCA','RobustF0','F0'};
    STAdecay_cat{f}=array2table(STAdecay_cat{f},'VariableNames',fieldName);

    maxSTA=max(STAdecay_cat{f}.STA{1}(:,win_pre+[-2:4]),[],2);
    PerisomROIs=abs(STAdecay_cat{f}.dendaxis{1})<50;
    normFactor_Rob=mean(maxSTA(PerisomROIs)./STAdecay_cat{f}.RobustF0{1}(PerisomROIs));
    %normFactor_Rob=max(maxSTA./STAdecay_cat{f}.RobustF0{1});
    normFactor_F0=mean(maxSTA(PerisomROIs)./STAdecay_cat{f}.F0{1}(PerisomROIs));
    normFactor_PCA=mean(maxSTA(PerisomROIs)./STAdecay_cat{f}.F0_PCA{1}(PerisomROIs));

    RobdFF_dcay=maxSTA./STAdecay_cat{f}.RobustF0{1}./normFactor_Rob;
    %idx=STAdecay_cat{f}.dendaxis{1}>=0;
    %RobdFF_dcay=RobdFF_dcay(idx);
    xfit=abs(STAdecay_cat{f}.dendaxis{1});
    % wt=1./RobdFF_dcay;
    % wt(RobdFF_dcay>1.1)=0;
    %[ft, c_fit, R2]=expfit_wBd(xfit,RobdFF_dcay,[0:500],[1 -200],[0.5 -500],[2 -10],wt);
    %[ft p]=expfitDM_2(xfit,RobdFF_dcay,[0:500]',[200]);
    %c_fit(2)=-p;

    if isfile(fullfile(save_figto,'bAPdecay_fit.mat'))
    else
        fitresult = interactive_expfit(xfit, RobdFF_dcay);
        RobustdFF_decay=[RobustdFF_decay; {xfit, RobdFF_dcay, fitresult, fitresult.p, fitresult.Rsq}];
    end

    if f==showf
        ax1=[]; ax2=[];
        noi=setdiff([1:length(STAdecay_cat{f}.dendaxis{1})],boi);
        [dax,dsort]=sort(STAdecay_cat{f}.dendaxis{1}(noi));

        figure(26); clf;
        tiledlayout(5,3); cax=[-0.1 1.2];
        ax2=[ax2 nexttile(3,[3 1])];
        imagesc([-win_pre:win_post],[1:length(noi)],STAdecay_cat{f}.STA{1}(noi(dsort),:)./STAdecay_cat{f}.RobustF0{1}(noi(dsort))./normFactor_Rob,cax)
        cb=colorbar;
        cb.Label.String='Normalized amplitude';
        set(gca,'Ytick',[]); title('Robust ∆F/F')
        ax2=[ax2 nexttile(2,[3 1])];
        imagesc([-win_pre:win_post],[1:length(noi)],STAdecay_cat{f}.STA{1}(noi(dsort),:)./STAdecay_cat{f}.F0_PCA{1}(noi(dsort))./normFactor_PCA,cax)
        xlabel('Time'); title('∆F/F_{std}'); 
        set(gca,'Ytick',[]);
        ax2=[ax2 nexttile(1,[3 1])];
        imagesc([-win_pre:win_post],[1:length(noi)],STAdecay_cat{f}.STA{1}(noi(dsort),:)./STAdecay_cat{f}.F0{1}(noi(dsort))./normFactor_F0,cax);
        set_kymoYtick(dax); ylabel('Distance from soma (µm)');
        title('∆F/F_0'); 
        colormap(turbo);
        linkaxes(ax2,'xy'); xlim([-50 30]);

        ax1=[ax1 nexttile(12,[2 1])];
        scatter(abs(STAdecay_cat{f}.dendaxis{1}(noi)),maxSTA(noi)./STAdecay_cat{f}.RobustF0{1}(noi)./normFactor_Rob,20,[0 0 0],'filled'); hold all;
        ax1=[ax1 nexttile(11,[2 1])];
        scatter(abs(STAdecay_cat{f}.dendaxis{1}(noi)),maxSTA(noi)./STAdecay_cat{f}.F0_PCA{1}(noi)./normFactor_PCA,20,[0 0 0],'filled'); hold all;
        xlabel('Distance from soma (µm)');
        ax1=[ax1 nexttile(10,[2 1])];
        scatter(abs(STAdecay_cat{f}.dendaxis{1}(noi)),maxSTA(noi)./STAdecay_cat{f}.F0{1}(noi)./normFactor_F0,20,[0 0 0],'filled');
        ylabel('Normalized amplitude');
        linkaxes(ax1,'xy'); 
        ylim([0 1.4]); xlim([0 500]);
        set_font('Arial'); set_fontsize(18);
        set_figsize(260,150);
    end
end
if isfile(fullfile(save_figto,'bAPdecay_fit.mat'))
    load(fullfile(save_figto,'bAPdecay_fit.mat'));
else
    fieldName={'dax','maxSTA','fit_result','c_fit','R2'};
    RobustdFF_decay=array2table(RobustdFF_decay,'VariableNames',fieldName);
    save(fullfile(save_figto,'bAPdecay_fit.mat'),'RobustdFF_decay','STAdecay_cat','-v7.3');
end
%%
figure(27); clf; tiledlayout(5,3);
nexttile(10,[2 2]);
shade2show=[];
for f=1:size(RobustdFF_decay,1)
    [dax, dsort]=sort(RobustdFF_decay.dax{f});
    norm2show=median(RobustdFF_decay.maxSTA{f}(RobustdFF_decay.dax{f}<50));
    sta2show=RobustdFF_decay.maxSTA{f}(dsort)./norm2show;
    %sta2show=sta2show/max(sta2show);
    shade2show{f}=[dax sta2show];
    %plot(dax,RobustdFF_decay.maxSTA{f}(dsort),'color',[0.6 0.6 0.6],'linewidth',1); hold all;
    scatter(dax,sta2show,15,[0.8 0.8 0.8],'filled'); hold all;
end
% dat2bin=cellfun(@(x,y) [x(:,1) y(:,1)],RobustdFF_decay.dax,RobustdFF_decay.maxSTA,'UniformOutput',0)';
% dat2exp=cell2mat(cellfun(@(x) [x(2)],RobustdFF_decay.c_fit,'UniformOutput',0))';
STAdecay_bin=binning_data(shade2show,[-20:40:200 260 320 400 460 600]);
errorbar_shade(STAdecay_bin.centers,STAdecay_bin.mean,STAdecay_bin.sem,[1 0 0]);
xlim([0 500]); ylim([0 1.3]);
xlabel('Distance from soma (µm)'); ylabel('Normalized amplitude');
set(gca,'xtick',[0 500]);
nexttile(12,[2 1]);
Boxplot_wPoints2(dat2exp',[0.4 0.4 0.4]);
ylim([0 500]); box off; ylabel('Decay length (µm)');
set(gca,'xtick',[],'linewidth',1);
set_font('Arial'); set_fontsize(18); set_figsize(128,112);

%% Plot and compare different methods
omitROI=[6 7 20];
f2show=23;
figure(6); clf; cax=[-0.1 1];
show_ROI=setdiff([1:nROI],omitROI);
dax=PCResult{f2show}.interDendDist(1,PCResult{f2show}.dist_order)*Pixelsize(f2show);
dROI2norm=find(abs(dax)<70);
dax=dax(show_ROI);

normSTAF0=STAs./F0'; normSTAF0=normSTAF0./mean(normSTAF0(dROI2norm,51),1);
normSTAF0_rbst=dFF_robust; normSTAF0_rbst=normSTAF0_rbst./mean(normSTAF0_rbst(dROI2norm,51),1);
normSTAF0_ref=STAs./F0_ref'; normSTAF0_ref=normSTAF0_ref./mean(normSTAF0_ref(dROI2norm,51),1);
normSTAF0_PCA=STAs./F0PCA; normSTAF0_PCA=normSTAF0_PCA./mean(normSTAF0_PCA(dROI2norm,51),1);

nexttile([1 1]);
imagesc([-50:50],[1:length(show_ROI)],normSTAF0(show_ROI,:),cax); title('∆F/F_0'); set_kymoYtick(dax);
nexttile([1 1]);
imagesc([-50:50],[1:length(show_ROI)],normSTAF0_rbst(show_ROI,:),cax); title(['Robust \DeltaF/F,' newline '(Wong-Campos et al., 2023)']); set_kymoYtick(dax);
nexttile([1 1]);
imagesc([-50:50],[1:length(show_ROI)],normSTAF0_ref(show_ROI,:),cax); title(['∆F/F_{ref},' newline '(Park et al., 2025)']); set_kymoYtick(dax);
nexttile([1 1]);
imagesc([-50:50],[1:length(show_ROI)],normSTAF0_PCA(show_ROI,:),cax); title('∆F/F_{std}'); set_kymoYtick(dax);
colorbar;
colormap(turbo);

nexttile([1 1]);
scatter(dax,max(normSTAF0(show_ROI,50+[-2:4]),[],2),50,turbo(length(dax)),'filled');
ylim([0 1.5]); xlim([0 500]); 
xlabel('Distance (µm)'); ylabel('Normalized amplitude');
nexttile([1 1]);
scatter(dax,max(normSTAF0_rbst(show_ROI,50+[-2:4]),[],2),50,turbo(length(dax)),'filled');
ylim([0 1.5]);  xlim([0 500]); 
xlabel('Distance (µm)'); ylabel('Normalized amplitude');
nexttile([1 1]);
scatter(dax,max(normSTAF0_ref(show_ROI,50+[-2:4]),[],2),50,turbo(length(dax)),'filled');
ylim([0 1.5]); xlim([0 500]); 
xlabel('Distance (µm)'); ylabel('Normalized amplitude');
nexttile([1 1]);
scatter(dax,max(normSTAF0_PCA(show_ROI,50+[-2:4]),[],2),50,turbo(length(dax)),'filled');
ylim([0 1.5]); xlim([0 500]); 
xlabel('Distance (µm)'); ylabel('Normalized amplitude');
set_font('Arial'); set_fontsize(20);

figure(27); clf;
nexttile([1 1])
scatter(max(normSTAF0_PCA(show_ROI,50+[-2:4]),[],2),max(normSTAF0_rbst(show_ROI,50+[-2:4]),[],2),50,turbo(length(dax)),'filled');
xlabel('PCA'); ylabel('rbst');
nexttile([1 1])
scatter(max(normSTAF0(show_ROI,50+[-2:4]),[],2),max(normSTAF0_rbst(show_ROI,50+[-2:4]),[],2),50,turbo(length(dax)),'filled');
xlabel('F0'); ylabel('rbst');
nexttile([1 1])
scatter(max(normSTAF0_PCA(show_ROI,50+[-2:4]),[],2),max(normSTAF0(show_ROI,50+[-2:4]),[],2),50,turbo(length(dax)),'filled');
xlabel('PCA'); ylabel('F0');

%%
ShowROI=10;
ROI2show=PCResult{f}.ftprnt(:,:,ShowROI)>0;
figure; clf;
scatter(AmpImg(ROI2show),F0_img_F0(ROI2show),'.'); hold all; yyaxis right;
scatter(AmpImg(ROI2show),F0_img_PCA(ROI2show),'.');

%%
figure(8); clf;
STA_IsolatedSS_cat=[];
for f=foi
    BadROI_distorder=find(ismember(PCResult{f}.dist_order,BadROI{f}));
    STA_IsolatedSS_cat{f}=[];
    nexttile([1 1]);
    TrRaw=PCResult{f}.normTraces;
    F0=tovec(PCResult{f}.ref_im)'*tovec(PCResult{f}.ftprnt);
    [~, STmat, spike2average]=get_STA(PCResult{f}.SpClass,PCResult{f}.SpClass(1,:),50,-1);
    Isisolated=squeeze(sum(STmat,[1 3]))==0;
    spike2average=spike2average(Isisolated);
    [STA_IsolatedSS]=get_STA(TrRaw,spike2average,50,50);
    STA_IsolatedSS=STA_IsolatedSS-prctile(STA_IsolatedSS(:,1:40),30,2);
    STA_IsolatedSS=STA_IsolatedSS./PCResult{f}.F0_PCA;
    %STA_IsolatedSS=STA_IsolatedSS./F0';
    STA_IsolatedSS_cat{f}=STA_IsolatedSS(PCResult{f}.dist_order,:);
    STA_IsolatedSS_cat{f}=STA_IsolatedSS_cat{f}(setdiff([1:size(STA_IsolatedSS,1)],BadROI_distorder),:);
    imagesc(STA_IsolatedSS_cat{f})
    title(f)
end
%%
f2show=[1 4 10 11 23 25 26];
%f2show=[1 6 10 23 24];

figure(9); clf;
SpikeAmp_cat=[]; g=1; bAPdecayfit=[];

ft = fittype('a*exp(x/b)');
opts = fitoptions(ft);
opts.StartPoint = [1.2, -150];
opts.Lower = [1.4, -300];
opts.Upper = [2, Inf];

for f=f2show
    SpikeAmp_cat{g}=[];
    BadROI_distorder=find(ismember(PCResult{f}.dist_order,BadROI{f}));
    SpikeAmp=max(STA_IsolatedSS_cat{f}(:,50+[-2:4]),[],2);
    dax=PCResult{f}.interDendDist(1,PCResult{f}.dist_order)*Pixelsize(f);
    dax=dax(setdiff([1:length(dax)],BadROI_distorder));
    dROI2norm=find(abs(dax)<70);
    SpikeAmp=SpikeAmp./mean(SpikeAmp(dROI2norm),1);
    [dax, dsort]=sort(dax,'ascend');
    %nexttile([1 1]);
    scatter(dax,SpikeAmp(dsort),20,ones(1,3)*0.7,'filled'); hold all;
    SpikeAmp_cat{g}=[dax' SpikeAmp(dsort)];
    %[s, bAPdecayfit(g)]=expfitDM_2(dax',SpikeAmp(dsort),[0:500]',[150]);
    %[cFitted s]=expfit_wBd(dax',,[0:500],[1 -150],[0.5 -400],[2 0]);

    %opts.Weights = 1./(1 + dax.^0.5);
    opts.Weights = (1 + dax.^2);
    fitResult = fit(dax', SpikeAmp(dsort), ft, opts);
    s=fitResult([0:500]);
    bAPdecayfit(g,:)=[fitResult.a fitResult.b];

    % hold all
    % plot([0:500],s,'r');
    g=g+1;
end
[binResult] = binning_data(SpikeAmp_cat, [-20:40:300 360 420 500]);
errorbar_shade(binResult.centers,binResult.mean,binResult.sem,[0 0 0]);
xlabel('Distance from soma (µm)'); ylabel('Normalized bAP amplitude');
set_font('Arial'); set_fontsize(15); box off; 
xlim([0 360]); ylim([0 1.5]);

% opts = fitoptions(ft);
% opts.StartPoint = [1, -200];
% opts.Lower = [0, -Inf];
% opts.Upper = [Inf, 0];      % use 0 if expecting decay
% opts.Weights = 1 ./ (binResult.sem.^2);
% 
% [f, gof] = fit(binResult.centers(:), binResult.mean(:), ft, opts);
% hold all;
% plot(f([0:500]),'r')

fprintf('Decay Const. mean: %3.2f, std: %3.2f\n',mean(bAPdecayfit(:,2)),std(bAPdecayfit(:,2)))












%% Show STA
STAmovieSS=importdata('/Volumes/BHL18TB_D2/20240808/173430BHLm143_N3_VR2/STAmovieIsolatedSS.mat');

TrRaw=PCResult{f2show}.normTraces;
FtprintSum=squeeze(sum(PCResult{f2show}.ftprnt,[1 2]));
Ftprnt2show=PCResult{f2show}.ftprnt(:,:,PCResult{f2show}.dist_order);
%TrRaw=TrRaw./FtprintSum;
[nROI nTime]=size(TrRaw);
sz=size(PCResult{f2show}.ref_im);
TrRaw=TrRaw(PCResult{f2show}.dist_order,:);
TrSub=get_subthreshold(TrRaw,PCResult{f2show}.spike(1,:),7,17);
figure(23); clf;
l=plot(TrRaw');
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(gen_colormap(Plasma,nROI),2));

sp=STAmovieSS.IsolatedSS_frame;
[STAs]=get_STA(TrRaw,sp,50,50);
STAs=STAs-median(STAs(:,1:40),2);

F0PCA=PCResult{f2show}.F0_PCA(PCResult{f2show}.dist_order);
F0=tovec(PCResult{f2show}.ref_im-100)'*tovec(PCResult{f2show}.ftprnt(:,:,PCResult{f2show}.dist_order));

STAmovie=-STAmovieSS.STAmovieIsolatedSS/8.5;
%STAs2=tovec(-STAmovie)'*tovec(PCResult{f2show}.ftprnt(:,:,PCResult{f2show}.dist_order));
dFF_robust=zeros(nROI,size(STAmovie,3));
msg = '';
for t=1:size(STAmovie,3)
    fprintf(repmat('\b', 1, numel(msg)));
    msg = sprintf('Calculating robust dF/F.. (%d/%d)', t, size(STAmovie,3));
    fprintf('%s', msg);

[~, dFF_robust(:,t)]=get_robustdFF(STAmovie(:,:,t),Ftprnt2show,PCResult{f2show}.ref_im-100,0);
end
fprintf('\n');
dFF_robust=dFF_robust-median(dFF_robust(:,1:40),2);
F0_ref=mean(STAs(:,50+[10:19]),2)';

%% Comparing F0 image and Fstd
[tr_df polyROI]=polyLineKymo3(STAmovie,10,10);
STA_F0=STAmovie./F0_img_F0; STA_Fstd=STAmovie./F0_img_PCA;
STA_F0(isinf(STA_F0))=NaN; STA_Fstd(isinf(STA_Fstd))=NaN;
[tr_dff0]=apply_clicky(polyROI,STA_F0);
[tr_dffstd]=apply_clicky(polyROI,STA_Fstd);

AmpImg=max(STAmovie(:,:,STAmovieSS.nTauPeak(1)+[-1:4]),[],3);
AmpImg_F0=AmpImg./F0_img_F0.*dMask2;
AmpImg_F0_vec=tovec(AmpImg_F0);
SomaAmp_F0=mean(AmpImg_F0_vec(tovec(PCResult{f}.ftprnt(:,:,1))>0,:),1,'omitnan');
AmpImg_F0=AmpImg_F0./SomaAmp_F0;

AmpImg_Fstd=AmpImg./F0_img_PCA.*dMask2;
AmpImg_Fstd_vec=tovec(AmpImg_Fstd);
SomaAmp_Fstd=mean(AmpImg_Fstd_vec(tovec(PCResult{f}.ftprnt(:,:,1))>0,:),1,'omitnan');
AmpImg_Fstd=AmpImg_Fstd./SomaAmp_Fstd;

figure(24); clf; ax1=[];
ax1=[ax1 nexttile([1 1])];
imshow2(AmpImg,[]);
title('∆F image');
ax1=[ax1 nexttile([1 1])];
imshow2(AmpImg_F0,[0 1.5]);
title('∆F/F_0 image');
ax1=[ax1 nexttile([1 1])];
imshow2(AmpImg_Fstd,[0 1.5]);
title('∆F/F_{std} image');
linkaxes(ax1,'xy');
colormap(turbo);


%% Plot and compare different methods
omitROI=[6 7 20];
f2show=23;
figure(6); clf; cax=[-0.1 1];
show_ROI=setdiff([1:nROI],omitROI);
dax=PCResult{f2show}.interDendDist(1,PCResult{f2show}.dist_order)*Pixelsize(f2show);
dROI2norm=find(abs(dax)<70);
dax=dax(show_ROI);

normSTAF0=STAs./F0'; normSTAF0=normSTAF0./mean(normSTAF0(dROI2norm,51),1);
normSTAF0_rbst=dFF_robust; normSTAF0_rbst=normSTAF0_rbst./mean(normSTAF0_rbst(dROI2norm,51),1);
normSTAF0_ref=STAs./F0_ref'; normSTAF0_ref=normSTAF0_ref./mean(normSTAF0_ref(dROI2norm,51),1);
normSTAF0_PCA=STAs./F0PCA; normSTAF0_PCA=normSTAF0_PCA./mean(normSTAF0_PCA(dROI2norm,51),1);

nexttile([1 1]);
imagesc([-50:50],[1:length(show_ROI)],normSTAF0(show_ROI,:),cax); title('∆F/F_0'); set_kymoYtick(dax);
nexttile([1 1]);
imagesc([-50:50],[1:length(show_ROI)],normSTAF0_rbst(show_ROI,:),cax); title(['Robust \DeltaF/F,' newline '(Wong-Campos et al., 2023)']); set_kymoYtick(dax);
nexttile([1 1]);
imagesc([-50:50],[1:length(show_ROI)],normSTAF0_ref(show_ROI,:),cax); title(['∆F/F_{ref},' newline '(Park et al., 2025)']); set_kymoYtick(dax);
nexttile([1 1]);
imagesc([-50:50],[1:length(show_ROI)],normSTAF0_PCA(show_ROI,:),cax); title('∆F/F_{std}'); set_kymoYtick(dax);
colorbar;
colormap(turbo);

nexttile([1 1]);
scatter(dax,max(normSTAF0(show_ROI,50+[-2:4]),[],2),50,turbo(length(dax)),'filled');
ylim([0 1.5]); xlim([0 500]); 
xlabel('Distance (µm)'); ylabel('Normalized amplitude');
nexttile([1 1]);
scatter(dax,max(normSTAF0_rbst(show_ROI,50+[-2:4]),[],2),50,turbo(length(dax)),'filled');
ylim([0 1.5]);  xlim([0 500]); 
xlabel('Distance (µm)'); ylabel('Normalized amplitude');
nexttile([1 1]);
scatter(dax,max(normSTAF0_ref(show_ROI,50+[-2:4]),[],2),50,turbo(length(dax)),'filled');
ylim([0 1.5]); xlim([0 500]); 
xlabel('Distance (µm)'); ylabel('Normalized amplitude');
nexttile([1 1]);
scatter(dax,max(normSTAF0_PCA(show_ROI,50+[-2:4]),[],2),50,turbo(length(dax)),'filled');
ylim([0 1.5]); xlim([0 500]); 
xlabel('Distance (µm)'); ylabel('Normalized amplitude');
set_font('Arial'); set_fontsize(20);

figure(27); clf;
nexttile([1 1])
scatter(max(normSTAF0_PCA(show_ROI,50+[-2:4]),[],2),max(normSTAF0_rbst(show_ROI,50+[-2:4]),[],2),50,turbo(length(dax)),'filled');
xlabel('PCA'); ylabel('rbst');
nexttile([1 1])
scatter(max(normSTAF0(show_ROI,50+[-2:4]),[],2),max(normSTAF0_rbst(show_ROI,50+[-2:4]),[],2),50,turbo(length(dax)),'filled');
xlabel('F0'); ylabel('rbst');
nexttile([1 1])
scatter(max(normSTAF0_PCA(show_ROI,50+[-2:4]),[],2),max(normSTAF0(show_ROI,50+[-2:4]),[],2),50,turbo(length(dax)),'filled');
xlabel('PCA'); ylabel('F0');

%%
ShowROI=10;
ROI2show=PCResult{f}.ftprnt(:,:,ShowROI)>0;
figure; clf;
scatter(AmpImg(ROI2show),F0_img_F0(ROI2show),'.'); hold all; yyaxis right;
scatter(AmpImg(ROI2show),F0_img_PCA(ROI2show),'.');

%%
figure(8); clf;
STA_IsolatedSS_cat=[];
for f=foi
    BadROI_distorder=find(ismember(PCResult{f}.dist_order,BadROI{f}));
    STA_IsolatedSS_cat{f}=[];
    nexttile([1 1]);
    TrRaw=PCResult{f}.normTraces;
    F0=tovec(PCResult{f}.ref_im)'*tovec(PCResult{f}.ftprnt);
    [~, STmat, spike2average]=get_STA(PCResult{f}.SpClass,PCResult{f}.SpClass(1,:),50,-1);
    Isisolated=squeeze(sum(STmat,[1 3]))==0;
    spike2average=spike2average(Isisolated);
    [STA_IsolatedSS]=get_STA(TrRaw,spike2average,50,50);
    STA_IsolatedSS=STA_IsolatedSS-prctile(STA_IsolatedSS(:,1:40),30,2);
    STA_IsolatedSS=STA_IsolatedSS./PCResult{f}.F0_PCA;
    %STA_IsolatedSS=STA_IsolatedSS./F0';
    STA_IsolatedSS_cat{f}=STA_IsolatedSS(PCResult{f}.dist_order,:);
    STA_IsolatedSS_cat{f}=STA_IsolatedSS_cat{f}(setdiff([1:size(STA_IsolatedSS,1)],BadROI_distorder),:);
    imagesc(STA_IsolatedSS_cat{f})
    title(f)
end
%%
f2show=[1 4 10 11 23 25 26];
%f2show=[1 6 10 23 24];

figure(9); clf;
SpikeAmp_cat=[]; g=1; bAPdecayfit=[];

ft = fittype('a*exp(x/b)');
opts = fitoptions(ft);
opts.StartPoint = [1.2, -150];
opts.Lower = [1.4, -300];
opts.Upper = [2, Inf];

for f=f2show
    SpikeAmp_cat{g}=[];
    BadROI_distorder=find(ismember(PCResult{f}.dist_order,BadROI{f}));
    SpikeAmp=max(STA_IsolatedSS_cat{f}(:,50+[-2:4]),[],2);
    dax=PCResult{f}.interDendDist(1,PCResult{f}.dist_order)*Pixelsize(f);
    dax=dax(setdiff([1:length(dax)],BadROI_distorder));
    dROI2norm=find(abs(dax)<70);
    SpikeAmp=SpikeAmp./mean(SpikeAmp(dROI2norm),1);
    [dax, dsort]=sort(dax,'ascend');
    %nexttile([1 1]);
    scatter(dax,SpikeAmp(dsort),20,ones(1,3)*0.7,'filled'); hold all;
    SpikeAmp_cat{g}=[dax' SpikeAmp(dsort)];
    %[s, bAPdecayfit(g)]=expfitDM_2(dax',SpikeAmp(dsort),[0:500]',[150]);
    %[cFitted s]=expfit_wBd(dax',,[0:500],[1 -150],[0.5 -400],[2 0]);

    %opts.Weights = 1./(1 + dax.^0.5);
    opts.Weights = (1 + dax.^2);
    fitResult = fit(dax', SpikeAmp(dsort), ft, opts);
    s=fitResult([0:500]);
    bAPdecayfit(g,:)=[fitResult.a fitResult.b];

    % hold all
    % plot([0:500],s,'r');
    g=g+1;
end
[binResult] = binning_data(SpikeAmp_cat, [-20:40:300 360 420 500]);
errorbar_shade(binResult.centers,binResult.mean,binResult.sem,[0 0 0]);
xlabel('Distance from soma (µm)'); ylabel('Normalized bAP amplitude');
set_font('Arial'); set_fontsize(15); box off; 
xlim([0 360]); ylim([0 1.5]);

% opts = fitoptions(ft);
% opts.StartPoint = [1, -200];
% opts.Lower = [0, -Inf];
% opts.Upper = [Inf, 0];      % use 0 if expecting decay
% opts.Weights = 1 ./ (binResult.sem.^2);
% 
% [f, gof] = fit(binResult.centers(:), binResult.mean(:), ft, opts);
% hold all;
% plot(f([0:500]),'r')

fprintf('Decay Const. mean: %3.2f, std: %3.2f\n',mean(bAPdecayfit(:,2)),std(bAPdecayfit(:,2)))


%%
f=23;

nTauPeak=[150 150];
load([fpath{f} '/output_data.mat'])
load(fullfile(fpath{f},'PC_Result.mat'),'Result');
sz=double(Device_Data{1, 3}.ROI([2 4]));
Mov_PeakTA=double(readBinMov([fpath{f} '/PeakTriggered_movie.bin'],sz(2)*sz(1),301));
Mov_PeakTA=Mov_PeakTA-mean(Mov_PeakTA,2);



%% Clean up and norm
exclude_frq=[241.7 242]; %monitor
%exclude_frq2=[483.5 484]; %monitor
exclude_frq2=[55.5 56]; %motion
time_bin=15000; Fs=1000; %2nd trunk is the reliable trace
RatioMat=[]; Lag1AutoCorr=[]; F0intensity=[]; SWCs=[]; Ftprnts=[]; dFtr=[];
%PSD=[]; 
figure(11); clf;
cmap=gray(100);
for f=[18 20 23]%:length(fpath)
    f
    load(fullfile(fpath{f},'PC_Result.mat'),'Result')

    ref_trace=ref_ROI{f}(1);
    [nROI nTime]=size(Result.traces);
    sz=size(Result.ref_im);
    t_silent=tovec(squeeze(find(Result.spike(1,:))'+[-10:40])); t_silent(t_silent<1)=[]; t_silent(t_silent>nTime)=[];

    tr_res=Result.traces_bvMask./squeeze(sum(Result.ftprnt,[1 2]));
    meanF=mean(tr_res,1,'omitnan');
    y_fit=-expfitDM_2(setdiff([1:nTime],t_silent)',-meanF(setdiff([1:nTime],t_silent))',[1:nTime]',[10000 1000]);
    tr_res=squeeze(SeeResiduals(permute(tr_res,[1 3 2]),y_fit));

    freq_lowhigh=exclude_frq/(Fs/2);   [b, a] = butter(4, freq_lowhigh, 'stop');
    freq_lowhigh2=exclude_frq2/(Fs/2); [b2, a2] = butter(4, freq_lowhigh2, 'stop');

    clear traces_res_filtered noise noise_intp norm_trace sp_height SpHeight_intp sp_time tr_mc_imcorr tr_mc tr norm_trace_check traces_res_filtered_ch
    tN=[1:time_bin:nTime]; tN=[tN nTime];
    lwpass_fit2=NaN(nROI,nTime);
    mcTrace=squeeze(Result.mcTrace)';
    

    for n=1:size(Result.traces,1)
        tr=tr_res(n,1:nTime);
        % regress out motion frequency
        for t=1:length(tN)-1
            tr(tN(t):tN(t+1))=squeeze(SeeResiduals(reshape(tr(tN(t):tN(t+1)),1,1,[]),mcTrace(:,(tN(t):tN(t+1)))));
            tr(tN(t):tN(t+1))=squeeze(SeeResiduals(reshape(tr(tN(t):tN(t+1)),1,1,[]),mcTrace(:,(tN(t):tN(t+1))).^2));
            tr(tN(t):tN(t+1))=squeeze(SeeResiduals(reshape(tr(tN(t):tN(t+1)),1,1,[]),mcTrace(1,(tN(t):tN(t+1))).*mcTrace(end,(tN(t):tN(t+1)))));

            tr_mc(tN(t):tN(t+1))=tr(tN(t):tN(t+1));
        end

        % regress out monitor frequency
        traces_res_filtered(n,:) = filtfilt(b, a, tr);
        traces_res_filtered(n,:) = filtfilt(b2, a2, traces_res_filtered(n,:));

        norm_trace(n,:)=traces_res_filtered(n,:);%-movmedian(traces_res_filtered(n,:),500,2);
        lwpass_fit2(n,t_silent)=norm_trace(n,t_silent); lwpass_fit2(n,:)=movmedian(lwpass_fit2(n,:),20000,2,'omitnan');
        norm_trace(n,:)=norm_trace(n,:)-lwpass_fit2(n,:);
    end

    Sptime=find(Result.spike(1,:));
    SpikeHeight=norm_trace(ref_trace,Sptime);
    y_fit=expfitDM_2(Sptime',SpikeHeight',[1:nTime]',10000);
    y_fit=y_fit/y_fit(1);

    norm_trace=norm_trace./y_fit';
    if isfield(Result,'dirtTrace')
    norm_trace(Result.dirtTrace>0)=NaN;
    end
    norm_trace(:,Result.motionReject>0)=NaN;

    %[PSD{f}.Power, PSD{f}.frequency, PSD{f}.normPower]=nanPSD(norm_trace,1000,5000); %use parallel computing

    dF=norm_trace;
    dFtr{f}=norm_trace;
    dFfilt=pcafilterTrace(movmean(norm_trace,15,2,'omitnan'),[1 2 3]);
    dFfilt=dFfilt-mean(dFfilt,2,'omitnan');
    Lag1AutoCorr{f}=sqrt(mean(dFfilt(:,2:end).*dFfilt(:,1:end-1),2,'omitnan'));

    [xx,yy] = meshgrid(1:sz(2),1:sz(1));
    mask = false(sz(1),sz(2));
    Result.SWC(:,3)=Result.SWC(:,3)+3;
    Result.SWC(1,3)=12;
    r=Result.SWC(:,3);
    for i = 1:size(Result.SWC,1)
        mask = mask | ((xx - Result.SWC(i,1)).^2 + (yy - Result.SWC(i,2)).^2 <= r(i)^2);
    end
    strImg=Result.ref_im-100;%./imgaussfilt(Result.ref_im,50);
    % strImgMasked=strImg;
    % strImgMasked(mask==1 | max(Result.ftprnt>0.5,[],3)>0)=NaN;
    % strImgbkg = medfilt2nan(strImgMasked, ones(1,2)*10);
    % [~, idx]=bwdist(~isnan(strImgbkg));
    % strImgbkg(isnan(strImgbkg)) = strImgbkg(idx(isnan(strImgbkg)));  % assign nearest values
    % strImgbkg = imgaussfiltnan(strImgbkg, 3);
    % strImg=strImg-strImgbkg;
    F0intensity{f}=tovec(strImg)'*tovec(Result.ftprnt)./squeeze(sum(Result.ftprnt,[1 2]))';

    % figure(f); clf; tiledlayout(2,6);
    % nexttile([1 2])
    % showScaleScatter(Lag1AutoCorr, Result.SWC, Result.ftprnt,'turbo');
    % nexttile([1 2])
    % showScaleScatter(F0intensity, Result.SWC, Result.ftprnt,'turbo');
    % nexttile([1 2])
    % showScaleScatter(Lag1AutoCorr./F0intensity', Result.SWC, Result.ftprnt,'turbo');
    R=Lag1AutoCorr{f}./F0intensity{f}';
    Result.SWC(:,3)=3; Result.SWC(1,3)=15;
    SWCs{f}=Result.SWC;
    Ftprnts{f}=Result.ftprnt;    

    Dsign=ones(1,size(Result.interDendDist,2));
    Dsign(Result.dist_order(1:find(Result.dist_order==1)-1))=-1;
    dendaxis=Result.interDendDist(1,:).*Dsign;

    RatioMat=[RatioMat; repmat(f,nROI,1) dendaxis' R];
end

fieldName={'NeuronID','Distance','Ratio'};
RatioMat=array2table(RatioMat,'VariableNames',fieldName);
%%
figure(5); clf; tiledlayout(2,6)
R_cat=[]; cmap=hsv(3); g=1; pxsz=[0.936 1.17 1.17];
for f=[18 20 23]
    showind=RatioMat.NeuronID==f;
    R=RatioMat.Ratio(showind);
    %normPSD=PSD{f}.Power./sum(PSD{f}.Power(:,find(PSD{f}.frequency<5)),2);
    normPSD=PSD{f}.Power./var(dFtr{f},0,2,'omitnan');
    normthetaPower=mean(normPSD(:,PSD{f}.frequency>5 & PSD{f}.frequency<10),2,'omitnan');
    nexttile(2*g-1,[1 1])
    showScaleScatter(R/R(1), SWCs{f}, Ftprnts{f} ,'turbo',[0.5 1.5]);
    drawScaleBar(100/pxsz(g),'vertical');
    title('Std/F0');
    nexttile(2*g,[1 1])
    showScaleScatter(normthetaPower, SWCs{f}, Ftprnts{f} ,'turbo',[0.005 0.04]);
    title('\theta-power');

    showindN=RatioMat.NeuronID==f;
    showind_S=RatioMat.NeuronID==f & RatioMat.Distance<50;
    nexttile(7,[1 2])
    scatter(RatioMat.Distance(showindN),RatioMat.Ratio(showindN)./mean(RatioMat.Ratio(showind_S)),10,cmap(g,:),'filled'); hold all
    title('Std/F0'); xlabel('Distance (\mum)');
    nexttile(9,[1 2])
    scatter(RatioMat.Distance(showindN),Lag1AutoCorr{f},10,cmap(g,:),'filled'); hold all
    title('Std'); xlabel('Distance (\mum)');
    nexttile(11,[1 2])
    scatter(RatioMat.Distance(showindN),normthetaPower,10,cmap(g,:),'filled'); hold all
    title('\theta-power'); xlabel('Distance (\mum)'); ylabel('Normalized power (A.U.)')
    g=g+1;
end

figure(6); clf;
f=23;
l=plot(PSD{f}.frequency,(PSD{f}.Power./var(dFtr{f},0,2,'omitnan'))');
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(vec2cmap(RatioMat.Distance(RatioMat.NeuronID==23),'jet',[0 400]),2));
xlim([0 20])
set(gca,'yscale','log');
ylabel('Normalized power (A.U.)')
xlabel('Frequency (Hz)')
colorbar; colormap('jet');
cb=colorbar; cb.Label.String='Distance (\mum)';
cb.Ticks=[0 1]; cb.TickLabels=[0 400];

figure(4); clf; cmap=gray(40);
R_cat=[];
for f=foi
showind=RatioMat.NeuronID==f;
showind_S=RatioMat.NeuronID==f & RatioMat.Distance<50;
scatter(RatioMat.Distance(showind),RatioMat.Ratio(showind)./mean(RatioMat.Ratio(showind_S)),10,cmap(f,:),'filled'); hold all
R_cat=[R_cat; [RatioMat.Distance(showind) RatioMat.Ratio(showind)./mean(RatioMat.Ratio(showind_S))]];
end
hold all;
[M S Xc Nc]=binning_data({R_cat},[-250:50:450]);
errorbar(Xc,M,S./sqrt(cellfun(@sum,Nc))','r')
xlabel('Distance from soma (\mum)')
ylabel('Normalized std(dF)/F_0')


%%
mov_mc=[];
mov_mc=cat(3,mov_mc,double(readBinMov([fpath{f} '/mc_ShutterReg' num2str(1,'%02d') '.bin'],sz(2),sz(1))));
load([fpath{f} '/mcTrace' num2str(1,'%02d') '.mat']);
mc=mcTrace.xymean;
t_fit=[1:size(mov_mc,3)];
[bleaching_fit] = expfitDM_2(t_fit',squeeze(mean(mov_mc,[1 2])),t_fit',[100000 10000]);

mov_res= mov_mc-median(mov_mc,3);
mov_res = SeeResiduals(mov_res,mc);
mov_res = SeeResiduals(mov_res,mc.^2);
mov_res = SeeResiduals(mov_res,mc(:,1).*mc(:,end));
mov_res = SeeResiduals(mov_res,bleaching_fit,1);
mov_res=tovec(mov_res);
%F0img=get_F0img(toimg(mov_res,sz(2),sz(1)));
mov_res_sub=movmedian(mov_res,20,2);
F0img=get_F0img_PCA(imresize(toimg(mov_res_sub,sz(2),sz(1)),0.7),[3000:8000]);
F0img2=get_F0img(imresize(toimg(mov_res_sub,sz(2),sz(1)),0.7),[3000:8000]);

%%
[xx,yy] = meshgrid(1:sz(1),1:sz(2));
mask = false(sz(2),sz(1));
r=Result.SWC(:,3)+3;
r(1)=12;
for i = 1:size(Result.SWC,1)
    mask = mask | ((xx - Result.SWC(i,1)).^2 + (yy - Result.SWC(i,2)).^2 <= r(i)^2);
end
strImg=Result.ref_im-100;%./imgaussfilt(Result.ref_im,50);
strImgMasked=strImg;
strImgMasked(mask==1 | max(Result.ftprnt>0.5,[],3)>0)=NaN;
strImgbkg = medfilt2nan(strImgMasked, ones(1,2)*10);
[~, idx]=bwdist(~isnan(strImgbkg));
strImgbkg(isnan(strImgbkg)) = strImgbkg(idx(isnan(strImgbkg)));  % assign nearest values
strImgbkg = imgaussfiltnan(strImgbkg, 3);

F0imgMasked=F0img;
F0imgMasked(mask==1 | max(Result.ftprnt>0.5,[],3)>0)=NaN;
F0imgbkg = medfilt2nan(F0imgMasked, ones(1,2)*10);
[~, idx]=bwdist(~isnan(F0imgbkg));
F0imgbkg(isnan(F0imgbkg)) = F0imgbkg(idx(isnan(F0imgbkg)));  % assign nearest values
F0imgbkg = imgaussfiltnan(F0imgbkg, 3);

figure; clf; ax1=[];
ax1=[ax1 nexttile([1 1])];
imshow2(Result.ref_im,[])
ax1=[ax1 nexttile([1 1])];
imshow2(strImgbkg,[])
ax1=[ax1 nexttile([1 1])];
strImg=strImg-strImgbkg;
imshow2(strImg,[])
%imagesc(strImg./strImgbkg); colormap(turbo);

DR=[prctile(strImg(mask==0),70),prctile(strImg(mask==1),99)];
strImg_bin=grs2rgb(strImg,colormap('gray'),DR(1),DR(2));
strImg_bin=strImg_bin(:,:,1);
strImg_bin(strImg_bin<0.02)=0;
ax1=[ax1 nexttile([1 1])];
imshow2(strImg_bin,[])
linkaxes(ax1)
%%
figure(15); clf; tiledlayout(2,2);
nexttile([1 1]);
F0img_crop=F0img(5:end-5,5:end-5);
imshow2(F0img_crop,[]); title('Subthreshold s.d.')
nexttile([1 1]);
strImg_crop=strImg(5:end-5,5:end-5);
imshow2(strImg_crop,[0 2000]); title('Background subtracted F_0');
nexttile([1 1]);
%imshow2(strImg_bin(5:end-5,5:end-5)>0.03,[]); title('Mask')
imshow2(strImg_crop./F0img_crop,[-10 150]); title('F0/s.d.');
nexttile([1 1]);
scatter_density(tovec(strImg_crop(strImg_bin(5:end-5,5:end-5)>0.03)),tovec(F0img_crop(strImg_bin(5:end-5,5:end-5)>0.03)),10,[],[0 0.0001])
xlabel('Subthreshold s.d.'); ylabel('F_0'); title('Pixels in mask');

%%
noi=setdiff(1:size(Result.ftprnt,3),BadROI{f});
Dist_order{f}=Result.BrancDist_order;
noi_dist{f}=ismember(Dist_order{f},noi);
Ftprnt_order=Result.ftprnt(:,:,Dist_order{f}(noi_dist{f}));

dFF_range=[-1.8 1.8];
ROI2show=[3 37];
dFF_movie=[]; dFF_tmp=[]; colored_dFFmov=[]; g=1;
cmap=gen_colormap([0 0.2 1; 0 0 0; 1 0 0],256);
t2show=[1:301];
for n=ROI2show
    dFF_movie{g}=toimg(Mov_PeakTA(:,:,n),sz(2),sz(1));
    dFF_movie{g}=imgaussfilt(dFF_movie{g},1);
    dFF_movie{g}=movmean(dFF_movie{g},5,3);

    for t=t2show;
        dFF_tmp{g}=dFF_movie{g}(:,:,t);
        dFF_tmp{g}(max(Result.bvMask,[],3)>0)=NaN;
        dFF_tmp{g}=medfilt2nan(dFF_tmp{g},[8 8]);
        dFF_movie{g}(:,:,t) = imgaussfiltnan(dFF_tmp{g}, 2).*strImg_bin(:,:,1);
    end

    colored_dFFmov{g}=[];
    for t=t2show
        colored_dFFmov{g}(:,:,:,t) = grs2rgb(double(dFF_movie{g}(:,:,t)), cmap ,dFF_range(1),dFF_range(2)).*strImg_bin(:,:,1);
        colored_dFFmov{g}(:,:,:,t) = grs2rgb(double(strImg), colormap("gray"),1,3)/1.5 + colored_dFFmov{g}(:,:,:,t)*6;
        %colored_dFFmov(:,:,:,t) = colored_dFFmov(:,:,:,t).*mat2gray(strImg);
    end
    g=g+1;
end
% Generate movie;
figure(20); clf;
v = VideoWriter(fullfile(fpath{f},['PTAmovie_' num2str(ROI2show)]),'MPEG-4');
%v = VideoWriter([fpath2read '/SNAPT_movie'],'Uncompressed AVI');
v.FrameRate = 25;  %can adjust this, 5 - 10 works well for me
v.Quality= 100;
open(v);

for j = t2show
    clf;
    set(gca,'units','pixels','position',[200 0 1000 800])
    ax1=[]; tiledlayout(length(ROI2show),1,'padding','compact');
    for g=1:length(ROI2show)
        ax1=[ax1 nexttile([1 1])];
        imshow2(colored_dFFmov{g}(:,:,:,j),[0 1]);
        pbaspect([size(double(colored_dFFmov{g}(:,:,:,j)),2) size(double(colored_dFFmov{g}(:,:,:,j)),1) 1]),colormap(gray)
        ftprntboundary=bwboundaries(Ftprnt_order(:,:,ROI2show(g))); hold all
        ftprntboundary_center=mean(ftprntboundary{1});
        scatter(ftprntboundary_center(2),ftprntboundary_center(1)+15,50,[1 0 0],'filled','^');
        colormap(gen_colormap([0 0.2 1; 1 1 1; 1 0 0],256))
        cb=colorbar;
        cb.Ticks=[0 1]; cb.TickLabels=dFF_range;
    end
    drawScaleBar(100/1.17,'horizontal','color',[1 1 1],'Linewidth',3);
    text(50,150,'100 \mum','color','w','Fontsize',17);
    set_fontsize(13);
    axis off
    nexttile(1,[1 1]);
    text(7,12,[num2str((j)-150) ' ms'], 'FontSize', 20, 'color', [0.99 0.99 0.99])% the value 1. is to adjust timing by eyes
    pause(0.1)
    set(gcf,'color','w')    % Sets background to white
    frame = getframe(gcf);
    writeVideo(v,frame);
    pause(0.1);
end;
close(v);

%%
[xx,yy] = meshgrid(1:sz(1),1:sz(2));
mask = false(sz(2),sz(1));
r=Result.SWC(:,3)+3;
r(1)=10;
for i = 1:size(Result.SWC,1)
    mask = mask | ((xx - Result.SWC(i,1)).^2 + (yy - Result.SWC(i,2)).^2 <= r(i)^2);
end
strImg=Result.ref_im-100;%./imgaussfilt(Result.ref_im,50);
strImgMasked=strImg;
strImgMasked(mask==1 | max(Result.ftprnt>0.5,[],3)>0)=NaN;
strImgbkg = medfilt2nan(strImgMasked, ones(1,2)*35);
[~, idx]=bwdist(~isnan(strImgbkg));
strImgbkg(isnan(strImgbkg)) = strImgbkg(idx(isnan(strImgbkg)));  % assign nearest values
strImgbkg = imgaussfiltnan(strImgbkg, 5);

figure; clf; ax1=[];
ax1=[ax1 nexttile([1 1])];
imshow2(Result.ref_im,[])
ax1=[ax1 nexttile([1 1])];
imshow2(strImgbkg,[])
ax1=[ax1 nexttile([1 1])];
imshow2(strImg./strImgbkg,[])
%imagesc(strImg./strImgbkg); colormap(turbo);
strImg=strImg./strImgbkg;
DR=[prctile(strImg(mask==0),70),prctile(strImg(mask==1),95)];
strImg_bin=grs2rgb(strImg,colormap('gray'),DR(1),DR(2));
strImg_bin=strImg_bin(:,:,1);
strImg_bin(strImg_bin<0.05)=0;
ax1=[ax1 nexttile([1 1])];
imshow2(strImg_bin,[])
linkaxes(ax1)



%%




%%
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


%%