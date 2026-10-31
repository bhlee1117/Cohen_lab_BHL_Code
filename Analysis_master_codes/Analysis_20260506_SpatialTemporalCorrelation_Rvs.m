%% Load file path
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
%foi=23;
%% Load Result files

f=20;
f
load(fullfile(fpath{f}, 'PC_Result.mat'), 'Result');
load(fullfile(fpath{f},'RobustdFFfit.mat'))
cd(fpath{f});
sz=size(Result.ref_im);

if isfile(fullfile(fpath{f},'PC_Result_smallerROI.mat'))
else
    SkelDend = point2img(Result.SWC,5,size(Result.ref_im));
    SkelDend2 = point2img(Result.SWC,2.5,size(Result.ref_im));
    [~, removeMask]=get_ROI(SkelDend2);
    SkelDend2(max(removeMask,[],3)>0)=0;
    SkelDend2=imfill(SkelDend2,'holes');
    imshow2(SkelDend2,[])
    roi_out = divide_ROIs(Result.ftprnt>0, 300, SkelDend2);
    label_img = sum(roi_out .* permute(1:size(roi_out,3), [1 3 2]), 3);
    figure(3); clf; imagesc(label_img); axis image; colormap(lines); colorbar;

    [x, y] = meshgrid(1:sz(2), 1:sz(1));
    CheckerBoardPatt=mod(floor(x/1) + floor(y/1), 2);
    ftprint_smaller=roi_out.*max(Result.ftprnt,[],3);
    ftprint_smaller_even=ftprint_smaller.*CheckerBoardPatt;
    ftprint_smaller_odd=ftprint_smaller.*(1-CheckerBoardPatt);
    save(fullfile(fpath{f},'PC_Result_smallerROI.mat'),'ftprint_smaller_odd','ftprint_smaller_even','ftprint_smaller','SkelDend','SkelDend2','-v7.3');
    disp(['Smaller ROI ftprnt is saved to ' fullfile(fpath{f},'PC_Result_smallerROI.mat')])
end

%% Signal extraction from multiple movie files, in streaming mode
bound = 6;
load(fullfile(fpath{f},"output_data.mat"))
CamCounter=Device_Data{1, 2}.Counter_Inputs(1, 1).data;
CamTrigger=find(CamCounter(2:end)-CamCounter(1:end-1));

frm_end=EndFrame(f);
f_seg=[[1:time_segment:frm_end] frm_end+1]; f_seg(2:end)=f_seg(2:end)-1;
take_window=repmat([1 time_segment],length(f_seg)-1,1);
take_window(2:end,1)=take_window(2:end,1)+overlap; take_window(1:end-1,2)=take_window(1:end-1,2)+overlap;
take_window(end)=mod(f_seg(end),time_segment);
take_window(take_window==0)=time_segment;

take_window_real=take_window+(f_seg(1:end-1))';
take_window_real(1,:)=take_window_real(1,:)-1;

Resultsm=importdata(fullfile(fpath{f},'PC_Result_smallerROI.mat'));
Resultsm.traces=[];
Resultsm.traces_even=[];
Resultsm.traces_odd=[];
Resultsm.mc=[];
nROI=size(Resultsm.ftprint_smaller,3);

g=1; intDD=[]; pth=[]; pair=[];
for i=1:nROI
    fprintf('Calculating inter-ROI distance %2.0f / %2.0f \n',i,nROI);
    for k=i:nROI
        [intDD(i,k), pth(:,:,g)]=geodesic_distance(Resultsm.SkelDend2>0,get_coord(Resultsm.ftprint_smaller(:,:,i)),get_coord(Resultsm.ftprint_smaller(:,:,k)));
        pth(:,:,g)=pth(:,:,g);
        pair(g,:)=[i k];
        g=g+1;
    end
end

figure(4); clf;
for r=randi(nROI,1,6)
    nexttile([1 1])
    imshow2(im_merge(cat(3,pth(:,:,r),mat2gray(Result.ref_im),mat2gray(Resultsm.SkelDend2)),[1 0 0; 1 1 1; 0 0 1]),[]);
    title(r);
end
intDD=max(cat(3,intDD,intDD'),[],3);
Resultsm.interDendDist=intDD;

for j=1:length(f_seg)-1

    fprintf('Extracting signal from movie %2.0f/ %2.0f\n',j,length(f_seg)-1);
    mov_mc=double(readBinMov([fpath{f} '/mc_ShutterReg' num2str(j,'%02d') '.bin'],sz(1),sz(2)));
    load([fpath{f} '/mcTrace' num2str(j,'%02d') '.mat']);

    mov_mc=mov_mc(:,:,[take_window(j,1):take_window(j,2)]);
    motionTrace=movmean(mcTrace.xymean([take_window(j,1):take_window(j,2)],:),3,1);
    Resultsm.mc=[Resultsm.mc; motionTrace];

    bkg = zeros(1, size(mov_mc,3));
    mean_F=squeeze(mean(mov_mc,[1 2]));
    [~, blueOff]=get_blueoffTrace(mean_F,[Result.Blue(take_window_real(j,1):take_window_real(j,2))],70);
    [y_fit]=expfitDM_2(find(blueOff)',mean_F(find(blueOff)),[1:size(mov_mc,3)]',1000);
    bkg(1,:)=y_fit;
    bv_trace=tovec(mov_mc)'*tovec(Result.bvMask);
    V_trace=tovec(mov_mc)'*tovec(Result.ftprnt);
    bv_trace=SeeResiduals(permute(bv_trace,[2 3 1]),V_trace);
    bv_trace=squeeze(bv_trace)';

    mov_mc = mov_mc./reshape(bkg,1,1,[]);
    mov_res = mov_mc-mean(mov_mc,3);
    mov_res = SeeResiduals(mov_res,motionTrace);
    mov_res = SeeResiduals(mov_res,motionTrace.^2);
    mov_res = SeeResiduals(mov_res,motionTrace(:,1).*motionTrace(:,end));
    mov_res= SeeResiduals(mov_res,bv_trace,1);

    Resultsm.traces_odd=[Resultsm.traces_odd, -(tovec(mov_res)'*tovec(Resultsm.ftprint_smaller_odd))'];
    Resultsm.traces_even=[Resultsm.traces_even, -(tovec(mov_res)'*tovec(Resultsm.ftprint_smaller_even))'];
    Resultsm.traces=[Resultsm.traces, -(tovec(mov_res)'*tovec(Resultsm.ftprint_smaller))'];

    clear mov_mc mov_res
end

SpikeTime=find(Result.spike(1,:)>0);
SpikeHeight=mean(Resultsm.traces([1:5],SpikeTime),1);
[spheight_intime]=expfitDM_2(SpikeTime',SpikeHeight',[1:frm_end]',[20000]);
Resultsm.SpikeHeight_fit=spheight_intime';

Resultsm.traces=Resultsm.traces./spheight_intime';
Resultsm.traces_even=Resultsm.traces_even./spheight_intime';
Resultsm.traces_odd=Resultsm.traces_odd./spheight_intime';

figure(2*f); clf; tiledlayout(2,2)
nexttile([1 1])
show_footprnt_contour(Resultsm.ftprint_smaller,Result.ref_im)
nexttile([1 1])
show_footprnt_contour(Result.bvMask,Result.ref_im)
nexttile([1 2])
plot(rescale2(Resultsm.traces,2)'+[1:size(Resultsm.traces,1)]);
drawnow;

save([fpath{f} '/PC_Result_smallerROI.mat'],'Resultsm','-v7.3')
disp(['Smaller ftprnt projected to movie is saved to ' fullfile(fpath{f},'PC_Result_smallerROI.mat')])
%% Load traces
f=20;
maxLag=1500;
Resultsm=importdata(fullfile(fpath{f},'PC_Result_smallerROI.mat'));
Result=importdata(fullfile(fpath{f}, 'PC_Result.mat'));

%% Calculate Xcorrelation

[nROI nTime]= size(Resultsm.traces);
PeriSpike_Frame=tovec(find(max(Result.SpClass,[],1)>0)'+[-3:5])';
BlueOne_Frame=find(Result.Blue>0);
DirtContam_Frame=find(max(Result.dirtTrace,[],1)>0);
MotionReject_Frame=find(Result.motionReject>0);
Time2Reject=unique([MotionReject_Frame DirtContam_Frame BlueOne_Frame PeriSpike_Frame]);
Time2Reject2=unique([MotionReject_Frame DirtContam_Frame BlueOne_Frame]);

if ~isfile(fullfile(fpath{f},'xCorr_subthreshold.mat'))
    Traces_even=Resultsm.traces_even;
    Traces_odd=Resultsm.traces_odd;
    Traces=Resultsm.traces;
    Traces(:,Time2Reject)=NaN;
    Traces_even(:,Time2Reject)=NaN;
    Traces_odd(:,Time2Reject)=NaN;

    Traces2_even=Resultsm.traces_even;
    Traces2_odd=Resultsm.traces_odd;
    Traces2=Resultsm.traces;
    Traces2(:,Time2Reject2)=NaN;
    Traces2_even(:,Time2Reject2)=NaN;
    Traces2_odd(:,Time2Reject2)=NaN;

    crossCorr=[]; crossCorr_wSpike=[]; g=1;
    corrROI2see=[1:nROI];
    for i=[1:nROI]
        for j=i:nROI
            if i==j
                [crossCorr(:,i,j), lags, nPairs] = nanXCorr(Traces_even(i,:), Traces_odd(j,:), maxLag, true);
                [crossCorr_wSpike(:,i,j), lags, nPairs] = nanXCorr(Traces2_even(i,:), Traces2_odd(j,:), maxLag, true);
            else
                [crossCorr(:,i,j), lags, nPairs] = nanXCorr(Traces(i,:), Traces(j,:), maxLag, true);
                [crossCorr_wSpike(:,i,j), lags, nPairs] = nanXCorr(Traces2(i,:), Traces2(j,:), maxLag, true);
            end
            g=g+1;
            if mod(g,10)==0
                fprintf('calculating Xcorr of subthreshold dynamics (%3.0f/%3.0f)\n',g,nROI*(length(corrROI2see)+1)/2)
            end
        end
    end
    save(fullfile(fpath{f},'xCorr_subthreshold.mat'),'crossCorr','crossCorr_wSpike','lags','maxLag','corrROI2see','-v7.3');
    disp('Xcorrdata saved');
end

xCorrResult=importdata(fullfile(fpath{f},'xCorr_subthreshold.mat'));

% Mirror the correlation value
N = size(xCorrResult.crossCorr, 2);
for i = 1:N
    for j = 1:i-1
        xCorrResult.crossCorr(:, i, j) = xCorrResult.crossCorr(:, j, i);
        xCorrResult.crossCorr_wSpike(:, i, j) = xCorrResult.crossCorr_wSpike(:, j, i);
    end
end

%% Plot Xcorrelation
CorrCat=[];
ind2rmv=[1 2 6 8 17 18 20 21 22 25 39 40 45 55 41 42 44 48 50 51 54 60 63 68 69 73 83 93 95 96 97];
ind2show=setdiff([1:length(xCorrResult.corrROI2see)],ind2rmv);
%ind2show=[23 24];
for corrROI2see=ind2show

    dax=Resultsm.interDendDist(corrROI2see,:)*Pixelsize(f);
    %dax(57:end)=-dax(57:end);
    Xcorrimg2show=squeeze(xCorrResult.crossCorr(:,corrROI2see,:))';
    Xcorrimg2showwSp=squeeze(xCorrResult.crossCorr_wSpike(:,corrROI2see,:))';
    [y_fit_spatial, Lconst, ~, R2]=expfitDM_2(dax',Xcorrimg2show(:,xCorrResult.maxLag+1),[0:400]', [50]);
    [y_fit_spatialwSp LconstwSp]=expfitDM_2(dax',Xcorrimg2showwSp(:,xCorrResult.maxLag+1),[0:400]', [50]);
    idx2use_SP=[1:length(dax)];
    if R2<0.7 | Lconst<20 | Lconst>300
        fprintf('Fitting spatial correlation, ROI #%2.0f\n',corrROI2see)
        opt.weights=Xcorrimg2show(:,xCorrResult.maxLag+1);
        fitResult_Sp=interactive_expfit(dax,Xcorrimg2show(:,xCorrResult.maxLag+1),opt);
        y_fit_spatial=fitResult_Sp.yFit;
        Lconst=fitResult_Sp.b;
        idx2use_SP=fitResult_Sp.idxSel;

        if  fitResult_Sp.aborted
            break;
        end

        opt2.predefined_idx=fitResult_Sp.idxSel;
        opt2.weights=Xcorrimg2showwSp(:,xCorrResult.maxLag+1);
        fitResult_SP_wSp=interactive_expfit(dax,Xcorrimg2showwSp(:,xCorrResult.maxLag+1),opt2);
        y_fit_spatialwSp=fitResult_SP_wSp.yFit;
        LconstwSp=fitResult_SP_wSp.b;
    end

    %ind2fit=[xCorrResult.maxLag+1:length(xCorrResult.lags)];
    ind2fit=[xCorrResult.maxLag+[1:100]];
    tau2fit=xCorrResult.lags(ind2fit);
    [y_fit, tconst, ~, R2t]=expfitDM_2(tau2fit,xCorrResult.crossCorr(ind2fit,corrROI2see,corrROI2see),tau2fit,[100]);
    [y_fit_wSp tconst_wSp]=expfitDM_2(tau2fit,xCorrResult.crossCorr_wSpike(ind2fit,corrROI2see,corrROI2see),tau2fit,[100]);
    idx2use=[1:length(tau2fit)];
    if R2t<0.7 | tconst<10
        fprintf('Fitting temporal correlation, ROI #%2.0f\n',corrROI2see)
        opt.weights=xCorrResult.crossCorr(ind2fit,corrROI2see,corrROI2see);
        fitResult=interactive_expfit(tau2fit,xCorrResult.crossCorr(ind2fit,corrROI2see,corrROI2see),opt);
        y_fit=fitResult.yFit;
        tconst=fitResult.b;
        idx2use=fitResult.idxSel;
        if  fitResult.aborted
            break;
        end

        opt2.predefined_idx=fitResult.idxSel;
        opt2.weights=xCorrResult.crossCorr_wSpike(ind2fit,corrROI2see,corrROI2see);
        fitResultwSp=interactive_expfit(tau2fit,xCorrResult.crossCorr_wSpike(ind2fit,corrROI2see,corrROI2see),opt2);
        y_fit_wSp=fitResultwSp.yFit;
        tconst_wSp=fitResultwSp.b;
    end


    CorrCat=[CorrCat; {corrROI2see ,tau2fit,xCorrResult.crossCorr(ind2fit,corrROI2see,corrROI2see), xCorrResult.crossCorr_wSpike(ind2fit,corrROI2see,corrROI2see)...
        y_fit, tconst, y_fit_wSp, tconst_wSp, idx2use, ...
        dax,Xcorrimg2show(:,xCorrResult.maxLag+1), Xcorrimg2showwSp(:,xCorrResult.maxLag+1),...
        y_fit_spatial,Lconst, y_fit_spatialwSp, LconstwSp, idx2use_SP}];
end
fieldName={'ROIidx','tax','TCorr','TCorr_wSp','TCorr_fit','tconst','TCorr_fit_wSp','tconst_wSp','taxfitidx',...
    'dax','SCorr','SCorr_wSp','SCorr_fit','lconst','SCorr_fit_wSP','lconst_wSp','daxfitidx'};
CorrCat=array2table(CorrCat,'VariableNames',fieldName);

% save(fullfile(fpath{f},'xCorrSubOnly_fitted.mat'),'CorrCat','ind2show','-v7.3');
% disp('Xcorr fit data saved');

load(fullfile(fpath{f},'xCorrSubOnly_fitted.mat'));
% tt=cell2mat(CorrCat.tconst);
% figure;
% histogram(tt(:,1),100); hold all;
% %histogram(tt(:,2),100); hold all;

%% Plot Xcorrelation
%load(fullfile(fpath{f},'xCorrSubOnly_fitted.mat'))
figure(222); clf;
tcorr2show=cellfun(@(x,y,z) [x(z) y(z)],CorrCat.tax, CorrCat.TCorr, CorrCat.taxfitidx,'UniformOutput',0);
tcorrbinned=binning_data(tcorr2show',[0:1500]-0.5);
scorr2show=cellfun(@(x,y,z) [x(z)' y(z)],CorrCat.dax, CorrCat.SCorr,CorrCat.daxfitidx,'UniformOutput',0);
scorrbinned=binning_data(scorr2show',[0:30:500]-15);
tcorr2show_wSp=cellfun(@(x,y,z) [x(z) y(z)],CorrCat.tax, CorrCat.TCorr_wSp, CorrCat.taxfitidx,'UniformOutput',0);
tcorrbinned_wSp=binning_data(tcorr2show_wSp',[0:1500]-0.5);
scorr2show_wSp=cellfun(@(x,y,z) [x(z)' y(z)],CorrCat.dax, CorrCat.SCorr_wSp,CorrCat.daxfitidx,'UniformOutput',0);
scorrbinned_wSp=binning_data(scorr2show_wSp',[0:30:500]-15);

tiledlayout(1,6,'TileSpacing','tight','Padding','compact');
nexttile([1 2]);
ROI2show=58;
dax=Resultsm.interDendDist(ROI2show,ind2show)*Pixelsize(f);
[daxsorted, dsort2show]=sort(dax);
corrmap2show=squeeze(xCorrResult.crossCorr(:,ROI2show,ind2show(dsort2show)))';
imagesc(xCorrResult.lags,[1:length(daxsorted)],corrmap2show,[0 .6]);
set_kymoYtick(daxsorted);
cb=colorbar; cb.Label.String='Correlation'; colormap(turbo);
xlabel('Time lag (ms)'); ylabel('Distance (µm)'); xlim([-1000 1000]);

nexttile([1 2]); hold all;
scatterT=cellfun(@(x,y,z) plot(x(z),y(z),'color',[0.7 0.7 0.7]),CorrCat.tax, CorrCat.TCorr, CorrCat.taxfitidx ,'UniformOutput',0);
errorbar_shade(tcorrbinned.centers,tcorrbinned.mean,tcorrbinned.sem,[1 0 0]);
title(sprintf('T_c: %2.0f ± %2.0f ms \n',mean(cell2mat(CorrCat.tconst)),std(cell2mat(CorrCat.tconst))));
xlim([0 100]); ylim([0 1]);
xlabel('Time lag (ms)'); ylabel('Correlation');

nexttile([1 2]); hold all;
[~, dsort]=cellfun(@(x,y) sort(x(y)),CorrCat.dax, CorrCat.daxfitidx,'UniformOutput',0);
scatterT=cellfun(@(x,y,zz,z) plot(x(zz(z)),y(zz(z)),'color',[0.7 0.7 0.7]),CorrCat.dax, CorrCat.SCorr,CorrCat.daxfitidx ,dsort, 'UniformOutput',0);
errorbar_shade(scorrbinned.centers,scorrbinned.mean,scorrbinned.sem,[0.5 0 1]);
title(sprintf('L_c: %2.0f ± %2.0f µm \n',mean(cell2mat(CorrCat.lconst)),std(cell2mat(CorrCat.lconst))));
xlim([0 300]); ylim([0 1]);
xlabel('Distance (µm)'); ylabel('Correlation');
set_font('Arial'); set_fontsize(15);
set_figsize(295, 100);
% 
% nexttile([1 3]);
% corrmap2show_wSp=squeeze(xCorrResult.crossCorr_wSpike(:,ROI2show,ind2show(dsort2show)))';
% imagesc(xCorrResult.lags,[1:length(daxsorted)],corrmap2show_wSp,[0 .6]);
% set_kymoYtick(daxsorted);
% cb=colorbar; cb.Label.String='Correlation'; colormap(turbo);
% xlabel('Time lag (ms)'); ylabel('Distance (µm)'); xlim([-1000 1000]);
% 
% nexttile([1 2]); hold all;
% scatterT=cellfun(@(x,y,z) plot(x(z),y(z),'color',[0.7 0.7 0.7]),CorrCat.tax, CorrCat.TCorr_wSp, CorrCat.taxfitidx ,'UniformOutput',0);
% errorbar_shade(tcorrbinned_wSp.centers,tcorrbinned_wSp.mean,tcorrbinned_wSp.sem,[1 0 0]);
% title(sprintf('T_c: %2.0f ± %2.0f ms \n',mean(cell2mat(CorrCat.tconst_wSp)),std(cell2mat(CorrCat.tconst_wSp))));
% xlim([0 100]); ylim([0 1]);
% xlabel('Time lag (ms)'); ylabel('Correlation');
% 
% nexttile([1 2]); hold all;
% scatterT=cellfun(@(x,y,zz,z) plot(x(zz(z)),y(zz(z)),'color',[0.7 0.7 0.7]),CorrCat.dax, CorrCat.SCorr_wSp,CorrCat.daxfitidx ,dsort, 'UniformOutput',0);
% errorbar_shade(scorrbinned_wSp.centers,scorrbinned_wSp.mean,scorrbinned_wSp.sem,[0.5 0 1]);
% title(sprintf('L_c: %2.0f ± %2.0f µm \n',mean(cell2mat(CorrCat.lconst_wSp)),std(cell2mat(CorrCat.lconst_wSp))));
% xlim([0 300]); ylim([0 1]);
% xlabel('Distance (µm)'); ylabel('Correlation');
% set_font('Arial'); set_fontsize(18);
% set_figsize(340, 115);

%% Read movie segment
frm2read=[65000:75000];
mov_mc=readBinMov_BHL_multiple(fpath{f},3,frm2read,15000,'mc_ShutterReg');
fprintf('Extracting signal from movie %2.0f – %2.0f\n',frm2read(1),frm2read(end));
motionTrace=Result.mcTrace(frm2read,:);
SpikeTrace=find(max(Result.spike(:,frm2read),[],1)>0)'+[-5:10];

bkg = zeros(1, size(mov_mc,3));
mean_F=squeeze(mean(mov_mc,[1 2]));
t2fit=setdiff([1:length(frm2read)],SpikeTrace);
[y_fit]=expfitDM_2(t2fit',mean_F(t2fit),[1:size(mov_mc,3)]',1000);
bkg(1,:)=y_fit/mean(y_fit);

mov_mc = mov_mc./reshape(bkg,1,1,[]);
mov_res1 = mov_mc-mean(mov_mc,3);
mov_res = SeeResiduals(mov_res1 ,motionTrace);
mov_res = SeeResiduals(mov_res,motionTrace.^2);
mov_res = SeeResiduals(mov_res,motionTrace(:,1).*motionTrace(:,end));
F0=mean(mov_mc,3)-100;
%% Calculate Variance image
mov_res_zm = mov_res - mean(mov_res, 3, 'omitnan');
var_total = mean(mov_res_zm.^2, 3, 'omitnan');
cov_lag_t = mean(mov_res_zm(:,:,1:end-1) .* mov_res_zm(:,:,2:end), 3, 'omitnan');

% Shot noise estimate:
% Var(total residual) - temporally correlated component.
shotVar_t = var_total - cov_lag_t;

%estimate shot noise using neighboring-pixel covariance too
kernel = ones(3,3);
kernel(2,2) = 0;

% Count valid neighbors for edge pixels.
neighborCount = conv2(ones(size(F0)), kernel, 'same');

% Compute local neighbor average movie.
neighborAvg = NaN(size(mov_res_zm), 'like', mov_res_zm);

for t = 1:size(mov_res_zm,3)
    frame = mov_res_zm(:,:,t);

    neighborSum = conv2(frame, kernel, 'same');
    neighborAvg(:,:,t) = neighborSum ./ neighborCount;
end
cov_spatial_3x3 = mean(mov_res_zm .* neighborAvg, 3, 'omitnan');
shotVar_spatial_3x3 = var_total - cov_spatial_3x3;

shotVar_spatial_3x3(~isfinite(shotVar_spatial_3x3)) = NaN;
shotVar_spatial_3x3 = max(shotVar_spatial_3x3, 0);

var_before = mean(mov_res1.^2, 3, 'omitnan');
var_after  = mean(mov_res.^2,  3, 'omitnan');

motionVar = var_before - var_after;
motionVar(~isfinite(motionVar)) = NaN;
motionVar = max(motionVar, 0);

%% Calculate low bound threshold ∆F/F
dcrit = 5;              % detection criterion in d' units. Try 3 or 5.
maxPixels = 5000;       % max number of pixels to average
maxFrames = 100;         % max number of frames to average

M_list = 1:maxPixels;   % x-axis: number of pixels averaged
K_list = 1:maxFrames;   % y-axis: number of frames averaged

roiMask = max(Result.ftprnt(:,:,1:7),[],3)>0;

F0= mean(mov_mc-100, 3, 'omitnan');
validMask = roiMask & isfinite(F0) & F0 > 0;

% Choose which empirical shot variance to use.
% Temporal version is usually the safest starting point.
shotVar = shotVar_spatial_3x3;

% Alternative robust combination:
% shotVar = median(cat(3, shotVar_t, shotVar_x, shotVar_y), 3, 'omitnan');

validMask = validMask & isfinite(shotVar) & shotVar > 0;

% Flatten pixels inside ROI

F0_vec = F0(validMask);
shotVar_vec = shotVar(validMask);

nPixAvailable = numel(F0_vec);
maxPixels = min(maxPixels, nPixAvailable);
M_list = M_list(M_list <= maxPixels);

% Choose order in which pixels are averaged

score = F0_vec ./ sqrt(shotVar_vec);
[~, sortIdx] = sort(score, 'descend');

F0_sorted = F0_vec(sortIdx);
shotVar_sorted = shotVar_vec(sortIdx);

% Cumulative sums over pixels

cumF0 = cumsum(F0_sorted);
cumShotVar = cumsum(shotVar_sorted);

cumF0 = cumF0(M_list);
cumShotVar = cumShotVar(M_list);

% Compute alpha_min heatmap
% alphaMin will be K x M:
% rows = number of frames averaged
% columns = number of pixels averaged
alphaMin = NaN(numel(K_list), numel(M_list));

for kk = 1:numel(K_list)
    K = K_list(kk);

    % Equal temporal weights h(t)=1:
    % signal = alpha * K * sum(F0_i)
    % noise std = sqrt(K * sum(shotVar_i))
    alphaMin(kk,:) = dcrit .* sqrt(K .* cumShotVar) ./ (K .* cumF0);
end

%% Plot low bound ∆F/F heatmap

% Single-pixel, single-frame noise variance (M = 1, K = 1).
% shotVar_vec is in camera counts^2 per pixel per frame.
% Dividing by F0^2 converts it to variance of dF/F for one pixel, one frame.
noiseVar_counts2 = shotVar_vec;                       % counts^2 / pixel / frame
noiseVar_dFF = shotVar_vec ./ F0_vec.^2;              % fractional (dF/F)^2

fprintf('\n--- Single-pixel, single-frame noise (n = %d ROI pixels) ---\n', nPixAvailable);
fprintf('Shot-noise variance [counts^2]: median %.4g | mean %.4g | IQR %.4g - %.4g\n', ...
    median(noiseVar_counts2), mean(noiseVar_counts2), ...
    prctile(noiseVar_counts2, 25), prctile(noiseVar_counts2, 75));
fprintf('Shot-noise std      [counts]  : median %.4g\n', median(sqrt(noiseVar_counts2)));
fprintf('F0                  [counts]  : median %.4g\n', median(F0_vec));
fprintf('Noise variance of dF/F        : median %.4g  (std = %.3g %%dF/F)\n', ...
    median(noiseVar_dFF), 100*median(sqrt(noiseVar_dFF)));
fprintf('Best-ranked pixel: var %.4g counts^2, F0 %.4g counts, %.3g %%dF/F std\n', ...
    shotVar_sorted(1), F0_sorted(1), 100*sqrt(shotVar_sorted(1))/F0_sorted(1));
fprintf('Minimum detectable dF/F at M=1, K=1 (d'' = %g): %.3g %%\n\n', ...
    dcrit, 100*alphaMin(1,1));

figure(336); clf;
% alphaMin is fractional dF/F.
% For display, convert to percent.
alphaPct = 100 * alphaMin;

pcolor(M_list, K_list, log10(alphaPct));
clim([-2 1])
shading interp
axis xy;
xlabel('Spatial averaging area (µm^2)');
ylabel('Temporal averaging window (ms)');

colormap(turbo);
cb = colorbar;
cb.Ticks=[-2:1:1];
cb.TickLabels=10.^[-2:1:1];
ylabel(cb, 'Minimum detectable \DeltaF/F (%)');
hold on; axis tight;

contourLevels = [0.05 0.2 0.5];

[C, h] = contour(M_list, K_list, alphaPct, contourLevels, ...
    'LineColor', 'k', ...
    'LineWidth', 1.5);

clabel(C, h, ...
    'Color', 'k', ...
    'FontSize', 11, ...
    'LabelSpacing', 500);
set(gca, 'FontSize', 16,'xscale','log','xtick',[1 10 100 1000]);
%scatter(squeeze(sum(Result.ftprnt>100000,[1 2])),20,'rd','filled')
hold off;

L_um = 110;       % length constant, um
tau_ms = 15;      % membrane time constant, ms
dt_ms = 1;        % frame interval, ms. Change this to your actual imaging dt.
pix_um = 1.17;       % pixel size, um/pixel. Change this to your actual pixel size.
dendWidth_um=5;

D_um2_per_ms = L_um^2 / tau_ms;

K_overlay = K_list;
t_ms = K_overlay * dt_ms;

% Axial spread along dendrite
d_um = sqrt(2 * D_um2_per_ms .* t_ms);

% Two-sided spread along a thin dendrite
area_dend_um2 = 2 .* d_um .* dendWidth_um;

% Convert area to approximate number of pixels
M_dend_spread = area_dend_um2 ./ pix_um^2;

hold on;
plot(M_dend_spread, K_overlay, 'r--', 'LineWidth', 2);

% Label the overlay with the diffusive-spread formula and the value of D
% used here (D = \lambda^2/\tau_m from the cable equation).
text(1500, K_overlay(round(end*0.6)), ...
    sprintf([' dendritic cable spread\n' ...
             ' d = (2Dt)^{1/2}\n' ...
             ' D = \\lambda^2/\\tau_m = %.0f \\mum^2/ms'], D_um2_per_ms), ...
    'Color', 'r', ...
    'FontSize', 14, ...
    'Interpreter', 'tex', ...
    'HorizontalAlignment','right', ...
    'VerticalAlignment','middle');

%xlim([10 max(M_list)])

set_font('Arial');
set_fontsize(18);
set_figsize(200,115);
hold off;

%% Histogram of per-pixel shot-noise variance inside roiMask
% Sanity check for the uniform-variance approximation below: how wide is the
% spread of shotVar across ROI pixels, and is the median representative?
maxShotVar = prctile(shotVar_vec,99);       
medShotVar = median(shotVar_vec, 'omitnan');          % counts^2 / pixel / frame
medNoiseVar_dFF = median(noiseVar_dFF, 'omitnan');    % (dF/F)^2 / pixel / frame

figure(337); clf;

subplot(1,3,1);
histogram(shotVar_vec, 50);
hold on; xline(medShotVar, 'r-', 'LineWidth', 2);
xlabel('Shot-noise variance (counts^2)');
ylabel('# pixels');
title(sprintf('median = %.4g', medShotVar));
axis square; hold off;

subplot(1,3,2);
histogram(log10(shotVar_vec), 50);
hold on; xline(log10(medShotVar), 'r-', 'LineWidth', 2);
xlabel('log_{10} variance (counts^2)');
ylabel('# pixels');
title(sprintf('IQR = %.3g - %.3g', ...
    prctile(shotVar_vec,25), prctile(shotVar_vec,75)));
axis square; hold off;

subplot(1,3,3);
histogram(100*sqrt(noiseVar_dFF), 50);
hold on; xline(100*sqrt(medNoiseVar_dFF), 'r-', 'LineWidth', 2);
xlabel('Single-pixel noise (% \DeltaF/F)');
ylabel('# pixels');
title(sprintf('median = %.3g %%', 100*sqrt(medNoiseVar_dFF)));
axis square; hold off;

set_font('Arial');
set_fontsize(14);

fprintf('\n--- Uniform-variance approximation (n = %d ROI pixels) ---\n', nPixAvailable);
fprintf('Median shot variance [counts^2]: %.4g (mean %.4g, ratio mean/median %.3g)\n', ...
    medShotVar, mean(shotVar_vec), mean(shotVar_vec)/medShotVar);
fprintf('Spread: 10th %%ile %.4g | 90th %%ile %.4g (90/10 = %.3g)\n', ...
    prctile(shotVar_vec,10), prctile(shotVar_vec,90), ...
    prctile(shotVar_vec,90)/prctile(shotVar_vec,10));

%% Calculate low bound threshold ∆F/F  --  single (median) variance version
% Same d' calculation as above, but every ROI pixel is assigned the SAME
% noise variance, medShotVar, instead of its own empirical shotVar.
% With uniform variance the optimal pixel ordering is simply brightest-F0-first.

[F0_sorted_u, sortIdx_u] = sort(F0_vec, 'descend');

cumF0_u = cumsum(F0_sorted_u);
cumF0_u = cumF0_u(M_list);
cumF0_u = cumF0_u(:);

% sum of variances over M pixels = M * medShotVar
cumShotVar_u = maxShotVar .* M_list(:);   % column, same shape as cumF0_u

alphaMin_u = NaN(numel(K_list), numel(M_list));

for kk = 1:numel(K_list)
    K = K_list(kk);
    % signal = alpha * K * sum(F0_i); noise std = sqrt(K * M * medShotVar)
    alphaMin_u(kk,:) = dcrit .* sqrt(K .* cumShotVar_u) ./ (K .* cumF0_u);
end

fprintf('Minimum detectable dF/F at M=1, K=1, uniform variance (d'' = %g): %.3g %%\n', ...
    dcrit, 100*alphaMin_u(1,1));
fprintf('Per-pixel-variance version gave %.3g %% (ratio uniform/per-pixel = %.3g)\n\n', ...
    100*alphaMin(1,1), alphaMin_u(1,1)/alphaMin(1,1));

%% Plot low bound ∆F/F heatmap -- single (median) variance version

figure(338); clf;
alphaPct_u = 100 * alphaMin_u;

pcolor(M_list, K_list, log10(alphaPct_u));
clim([-2 1])
shading interp
axis xy;
xlabel('Spatial averaging area (µm^2)');
ylabel('Temporal averaging window (ms)');
title(sprintf('Uniform variance = %.4g counts^2/pixel/frame', medShotVar));

colormap(turbo);
cb = colorbar;
cb.Ticks=[-2:1:1];
cb.TickLabels=10.^[-2:1:1];
ylabel(cb, 'Minimum detectable \DeltaF/F (%)');
hold on; axis tight;

[C_u, h_u] = contour(M_list, K_list, alphaPct_u, contourLevels, ...
    'LineColor', 'k', ...
    'LineWidth', 1.5);

clabel(C_u, h_u, ...
    'Color', 'k', ...
    'FontSize', 11, ...
    'LabelSpacing', 500);
set(gca, 'FontSize', 16,'xscale','log','xtick',[1 10 100 1000]);

% Same dendritic cable-spread overlay as figure 336
plot(M_dend_spread, K_overlay, 'r--', 'LineWidth', 2);
text(1500, K_overlay(round(end*0.6)), ...
    sprintf([' dendritic cable spread\n' ...
             ' d = (2Dt)^{1/2}\n' ...
             ' D = \\lambda^2/\\tau_m = %.0f \\mum^2/ms'], D_um2_per_ms), ...
    'Color', 'r', ...
    'FontSize', 14, ...
    'Interpreter', 'tex', ...
    'HorizontalAlignment','right', ...
    'VerticalAlignment','middle');

set_font('Arial');
set_fontsize(18);
set_figsize(200,115);
hold off;

%% Impose artificial signal.
if ~isfile(fullfile(fpath{f},'injectionResult.mat'))
    winLen = 1001;                                  % background window length [frames], same as old [4850:5850]
    nReps  = 40;                                     % # random background windows per condition, used for detection probability
    validStarts = 1 : (numel(frm2read) - winLen);     % valid start indices (relative to frm2read) for a full window
    rng(1);                                           % reproducible random frame-window draws
    repStarts = validStarts(randi(numel(validStarts), 1, nReps));

    [voltSTA]=get_STA(Result.normTraces(1,:)./RobustdFFfit.ScaleFactor(1),max(Result.SpClass([1 2],:),[],1)>0,30,30);
    voltSTA=voltSTA-median(voltSTA(1:20));
    SpHeight=max(voltSTA);
    InterDenddist=Result.interDendDist;
    InterDenddist=max(cat(3,InterDenddist,InterDenddist'),[],3);

    % alpha_injection_test.m
    % Estimate detection limit by injecting synthetic alpha-function events
    % into mov_res and running the existing detection algorithm.
    %
    % Requires in workspace:
    %   mov_res      : base movie [Y x X x T]
    %   Result.swc   : dendrite skeleton [N x 2], columns = [col(x), row(y)] in pixels

    % Parameters
    lengths   = [1:20:385];         % dendrite segment lengths [µm]
    %lengths   = [100 200];         % dendrite segment lengths [µm]
    alphas    = [0.005 0.01 0.02 0.04 0.08 0.14];  % amplitude fractions (dF)
    %alphas    = [0.005 0.04];  % amplitude fractions (dF)
    taus      = [1:5:70];          % alpha time constants [frames]

    % Options
    [nY, nX, ~] = size(mov_res);
    nT         = winLen;
    t0         = round(nT / 2);    % injection frame (relative to each random window)
    tol        = 7;                % detection tolerance [frames]
    s          = -1;               % sign: -1 for Voltron2 depolarization
    px_per_um  = 1/1.17;                % pixels per µm
    gauss_sigma = 2.0;             % Gaussian width [px] for skeleton footprint

    % N nominated injection points [row, col] — sampled randomly each rep
    nominated_centers = [78 217;   % anchor 1
                         90 301;   % anchor 2
                         95 357];  % anchor 3  <- adjust to your dendrite locations
    N_anchors = size(nominated_centers, 1);

    % Step 2: Baseline image F0
    F0 = mean(mov_mc, 3);   % [Y x X]

    % Step 3: Skeleton and arc-length
    skel_col = Result.SWC(:, 1);
    skel_row = Result.SWC(:, 2);
    [xx, yy] = meshgrid(1:nX, 1:nY);
    d_arc = [0; cumsum(sqrt(diff(skel_col).^2 + diff(skel_row).^2))];

    % Precompute per-anchor: anchor skeleton index and ROI_injected
    neighbor_dist_um = 35;
    anchor_idx_all      = zeros(1, N_anchors);
    ROI_injected_all    = cell(1, N_anchors);
    ROI_check_all       = cell(1, N_anchors);

    for iN = 1:N_anchors
        sc = nominated_centers(iN, :);
        dist_to_center = sqrt((skel_row - sc(1)).^2 + (skel_col - sc(2)).^2);
        [~, anchor_idx_all(iN)] = min(dist_to_center);
        ROI_injected_n = find(Result.ftprnt(sc(1), sc(2), :) > 0);
        ROI_neighbors_n = find(InterDenddist(ROI_injected_n, :) * Pixelsize(f) <= neighbor_dist_um);
        ROI_injected_all{iN} = ROI_injected_n;
        ROI_check_all{iN}    = union(ROI_injected_n, ROI_neighbors_n);
    end

    % Step 4: Precompute spatial footprints S_all [nY x nX x nL x N_anchors]
    nL = length(lengths);
    S_all = zeros(nY, nX, nL, N_anchors);

    for iN = 1:N_anchors
        arc_anchor_n = d_arc(anchor_idx_all(iN));
        for iL = 1:nL
            ell_px   = lengths(iL) * px_per_um;
            arc_sigma = ell_px / 2;
            arc_dist_n = abs(d_arc - arc_anchor_n);

            S = zeros(nY, nX);
            for k = 1:numel(skel_row)
                w  = exp(-arc_dist_n(k)^2 / (2 * arc_sigma^2));
                dr = yy - skel_row(k);
                dc = xx - skel_col(k);
                S  = S + w * exp(-(dr.^2 + dc.^2) / (2 * gauss_sigma^2));
            end
            if max(S(:)) > 0, S = S / max(S(:)); end
            S_all(:,:,iL,iN) = S;
        end
        fprintf('Precomputed footprints for anchor %d/%d\n', iN, N_anchors);
    end
    %
    % Footprint preview — show all length masks overlaid on F0
    % figure('Name', 'Injection footprints');
    % tiledlayout(1, nL, 'TileSpacing', 'compact', 'Padding', 'compact');
    %
    % for iL = 1:nL
    %     nexttile;
    %     imagesc(F0); colormap(gray); axis image off; hold on;
    %     Overlay footprint in red
    %     S_show = S_all(:,:,iL);
    %     rgb = cat(3, ones(nY,nX), zeros(nY,nX), zeros(nY,nX));
    %     h_ov = imagesc(rgb);
    %     set(h_ov, 'AlphaData', S_show * 0.6);
    %     Draw full skeleton
    %     scatter(skel_col, skel_row, 10,'c', 'filled');
    %     Mark anchor
    %     plot(skel_col(anchor_idx), skel_row(anchor_idx), 'y*', 'MarkerSize', 8);
    %     title(sprintf('\\ell = %d µm', lengths(iL)));
    %     hold off;
    % end
    % sgtitle('Injection footprints along dendrite skeleton (red = S, cyan = full skeleton)');

    % Pre-allocate results
    nL        = numel(lengths);
    nA        = numel(alphas);
    nTau      = numel(taus);
    P_detect  = nan(nL, nA, nTau);

    % -------------------------------------------------------------------------
    % Main loop
    % -------------------------------------------------------------------------
    nCond    = nL * nA * nTau;
    cond_idx = 0;
    injectionResult = {};

    warning('off', 'signal:findpeaks:largeMinPeakHeight');

    for iL = 1:nL

        for iA = 1:nA
            alpha = alphas(iA);

            for iTau = 1:nTau
                tau = taus(iTau);
                cond_idx = cond_idx + 1;
                fprintf('Condition %d/%d  ell=%.1f um  alpha=%.3f  tau=%d frames  (%d reps x %d anchors)\n', ...
                    cond_idx, nCond, lengths(iL), alpha, tau, nReps, N_anchors);

                % Step 5: Alpha waveform h(t; tau) = (t/tau)*exp(1 - t/tau)
                % (same for every repeat — depends only on t0, tau, nT)
                t_rel = (0:nT-1)' - t0;
                h = zeros(nT, 1);
                mask = t_rel > 0;
                tr = t_rel(mask);
                h(mask) = (tr / tau) .* exp(1 - tr / tau);

                if max(h) > 0, h = h / max(h); end
                h3 = reshape(h, 1, 1, nT);   % [1 x 1 x nT] for broadcast
                [~, t_shift]=max(h);
                t_shift=t_shift-t0;

                % Repeat over random background windows AND random anchor points
                peak_hit      = false(1, nReps);
                spike_hit     = false(1, nReps);   % spike at any checked ROI
                spike_hit_inj = false(1, nReps);   % spike only at ROI_injected

                for iRep = 1:nReps
                    % Sample random anchor for this rep
                    iN_rep       = randi(N_anchors);
                    S            = S_all(:,:,iL,iN_rep);
                    ROI_injected = ROI_injected_all{iN_rep};
                    ROI_check    = ROI_check_all{iN_rep};

                    % Step 6: build injected movie for this rep's anchor + background
                    % I_inj(x,y,t) = I_base(x,y,t) + s*alpha*F0(x,y)*S(x,y)*h(t)
                    spatial = s * alpha * F0 .* S;
                    inj     = spatial .* h3;

                    injectfrm_r      = repStarts(iRep) + (0:winLen-1);
                    mov_res2inject_r = mov_res(:,:,injectfrm_r);
                    normFactor_r     = mean(Result.SpikeHeight_fit(frm2read(injectfrm_r)));
                    SpikeInj_r       = Result.spike(:,frm2read(injectfrm_r));

                    mov_inj = mov_res2inject_r + inj;

                    OriginalVtrace=tovec(mov_res2inject_r)'*tovec(Result.ftprnt);
                    OriginalVtrace=-OriginalVtrace'./normFactor_r;
                    Origin_tr=OriginalVtrace./get_threshold(OriginalVtrace,1);
                    Originspike_tr=find_spike_bh(Origin_tr-movmedian(Origin_tr,20,2),2.5,1.5);

                    injectedVtrace=tovec(mov_inj)'*tovec(Result.ftprnt);
                    injectedVtrace=-injectedVtrace'./normFactor_r;
                    inject_tr=injectedVtrace./get_threshold(injectedVtrace,1);

                    spike_tr=find_spike_bh(inject_tr-movmedian(inject_tr,20,2),2.5,1.5);
                    spike_tr=max(cat(3,SpikeInj_r,spike_tr),[],3);
                    spike_tr2=(spike_tr-SpikeInj_r-Originspike_tr)>0;   % injection-evoked spikes only
                    dFFtrace=injectedVtrace./RobustdFFfit.ScaleFactor;
                    subthreshold_tr=get_subthreshold(dFFtrace,spike_tr,7,20);

                    % Peak/spike detection across the checked ROIs
                    peak_found      = false;
                    spike_found     = false;
                    spike_found_inj = false;
                    h_active = find(h > 0.1);

                    for n = ROI_check(:)'
                        FilterTr = subthreshold_tr(n,:) / SpHeight;
                        FilterTr = FilterTr - median(FilterTr, 2,'omitnan');

                        [~, locs] = findpeaks(FilterTr, 'MinPeakHeight', 0.2, ...
                            'MinPeakDistance', 50, 'MinPeakProminence', 0.1);
                        if any(ismember(locs, h_active)) && ismember(n, ROI_injected)
                            peak_found = true;
                        end

                        sp_n = spike_tr2(n, :);
                        tWin = max(1, t0-tol+t_shift) : min(nT, t0+tol+t_shift);
                        if any(sp_n(tWin))
                            spike_found = true;
                            if ismember(n, ROI_injected)
                                spike_found_inj = true;
                            end
                        end
                    end
                    peak_hit(iRep)      = peak_found;
                    spike_hit(iRep)     = spike_found;
                    spike_hit_inj(iRep) = spike_found_inj;
                end

                either_hit = peak_hit | spike_hit_inj;   % P(peak OR spike at injected ROI)

                injectionResult(end+1, :) = {lengths(iL), alpha, tau, t0, ...
                    mean(peak_hit), mean(spike_hit), mean(spike_hit_inj), mean(either_hit), nReps};
            end
        end
    end

    fieldNames = {'ell_um','alpha','tau','t0', ...
        'peak_detect_prob','spike_detect_prob','spike_detect_prob_injROI','either_detect_prob','nReps'};
    injectionResult = cell2table(injectionResult, 'VariableNames', fieldNames);
    save(fullfile(fpath{f},'injectionResult.mat'),'injectionResult','lengths','taus','alphas','S_all','-v7.3');
    disp('Injection result is saved');
end
%% Step 10: 3D scatter plot
load(fullfile(fpath{f},'injectionResult.mat'));
ell_vec   = (injectionResult.ell_um);
alpha_vec = (injectionResult.alpha);
tau_vec   = (injectionResult.tau);
peak_vec  = vec2cmap(injectionResult.peak_detect_prob,'jet');
spike_vec = (injectionResult.spike_detect_prob);

% Plot 2: 2D heatmap (Length x Tau) per alpha, peak+spike overlaid, with cable spread
% Cable spread line: ell_spread(tau) = sqrt(2 * D * tau_ms)
% taus in frames -> ms using dt_ms

D_um2_per_ms  = L_um^2 / tau_ms;   % reuse L_um, tau_ms from earlier section
tau_ms_vec    = taus * dt_ms;       % injection taus -> ms
ell_cable     = sqrt(2 * D_um2_per_ms .* tau_ms_vec);   % [um], one per tau

% Plot 2: Detection probability heatmap (Length x Tau), one column-pair per alpha
% Each alpha gets 2 tiles: left = peak prob, right = spike prob
nAlpha = numel(alphas);
alphas2show = 1:nAlpha-1;
det_types  = {'peak_detect_prob','spike_detect_prob_injROI','either_detect_prob'};
det_labels = {'Peak (ROI injected)','Spike (ROI injected)','Peak or Spike (ROI injected)'};
det_types_str={'peaks','spikes','either'};

figure(17); clf;
tiledlayout(2, length(alphas2show), 'TileSpacing', 'compact', 'Padding', 'compact');

for iDet = 1:2
    for iA = alphas2show
        nexttile; hold on;
        prob_map = nan(numel(taus), numel(lengths));
        for iTau = 1:numel(taus)
            for iL = 1:numel(lengths)
                idx = injectionResult.ell_um == lengths(iL) & ...
                    injectionResult.tau    == taus(iTau)  & ...
                    injectionResult.alpha  == alphas(iA);
                if ~any(idx), continue; end
                prob_map(iTau, iL) = injectionResult.(det_types{iDet})(idx);
            end
        end
        imagesc(lengths/2, taus, prob_map, [0 1]);
        colormap(turbo); axis xy; axis tight;
        plot(ell_cable, taus, 'w--', 'LineWidth', 2);
        if iA==alphas2show(1)
        ylabel('\tau (ms)');
        else
            set(gca,'Ytick',[]);
        end
        if iDet==1
            title(sprintf('%2.1f%% ∆F/F', alphas(iA)*100));
            if iA == length(alphas2show)
                cb = colorbar; cb.Label.String = sprintf('Peak detection \n probability');
            end
        end
        if iDet==2
            xlabel('\sigma_L (µm)'); 
            if iA == length(alphas2show)
                cb = colorbar; cb.Label.String = sprintf('Spike detection \n probability');
            end
        end
        xlim([0 max(lengths)/2]);
        hold off;
    end
end
set_font('Arial'); set_fontsize(16);
set_figsize(230,130);

%% Plot 2b: Detection probability heatmap (Length x Alpha), one column-pair per tau
% Each tau gets 2 tiles: top row = peak prob, bottom row = spike prob
showTau = 1:3;%numel(taus);

figure(18); clf;
tiledlayout(2, length(showTau), 'TileSpacing', 'compact', 'Padding', 'compact');

for iDet = 1:2
    for iTau = showTau
        nexttile; hold on;
        prob_map2 = nan(numel(alphas), numel(lengths));
        for iA = 1:numel(alphas)
            for iL = 1:numel(lengths)
                idx = injectionResult.ell_um == lengths(iL) & ...
                    injectionResult.tau    == taus(iTau)  & ...
                    injectionResult.alpha  == alphas(iA);
                if ~any(idx), continue; end
                prob_map2(iA, iL) = injectionResult.(det_types{iDet})(idx);
            end
        end
        imagesc(lengths/2, alphas * 100, prob_map2, [0 1]);
        colormap(turbo); axis xy; axis tight;
        xline(ell_cable(iTau), 'w--', 'LineWidth', 2);
       if iDet==2
        xlabel('\sigma_L (\mum)');
       end
       if iTau == 1
           ylabel('\alpha (% \DeltaF/F)');
       end
        set(gca, 'YScale', 'log', 'YTick', alphas*100, ...
            'YTickLabel', arrayfun(@(a) sprintf('%.1f',a), alphas*100, 'UniformOutput', false));
        title(sprintf('\\tau = %d ms', taus(iTau)));
        hold off;
    end
end
set_font('Arial'); set_fontsize(14);
set_figsize(148,130);

%% Plot 2c: Detection of ANY event (peak OR spike), one heatmap per alpha
% Same (sigma_L x tau) axes as figure 17, but collapsing the two detectors
% into either_detect_prob = P(peak OR spike at the injected ROI) per rep.

hasEither = ismember('either_detect_prob', injectionResult.Properties.VariableNames);
if ~hasEither
    warning(['either_detect_prob not found — plotting min(peak+spike,1) as an ' ...
        'upper bound. Re-run the injection loop for exact values.']);
end

figure(22); clf;
tiledlayout(1, numel(alphas2show), 'TileSpacing', 'compact', 'Padding', 'compact');

for iA = alphas2show
    nexttile; hold on;
    prob_map_either = nan(numel(taus), numel(lengths));
    for iTau = 1:numel(taus)
        for iL = 1:numel(lengths)
            idx = injectionResult.ell_um == lengths(iL) & ...
                injectionResult.tau    == taus(iTau)  & ...
                injectionResult.alpha  == alphas(iA);
            if ~any(idx), continue; end
            if hasEither
                prob_map_either(iTau, iL) = mean(injectionResult.either_detect_prob(idx), 'omitnan');
            else
                prob_map_either(iTau, iL) = min(mean(injectionResult.peak_detect_prob(idx), 'omitnan') + ...
                    mean(injectionResult.spike_detect_prob_injROI(idx), 'omitnan'), 1);
            end
        end
    end
    imagesc(lengths/2, taus, prob_map_either, [0 1]);
    colormap(turbo); axis xy; axis tight;
    plot(ell_cable, taus, 'w--', 'LineWidth', 2);
    title(sprintf('%2.1f%% ∆F/F', alphas(iA)*100));
    xlabel('\sigma_L (µm)');
    if iA == alphas2show(1)
        ylabel('\tau (ms)');
    else
        set(gca,'Ytick',[]);
    end
    if iA == alphas2show(end)
        cb = colorbar; cb.Label.String = sprintf('Peak or spike \n detection probability');
    end
    xlim([0 max(lengths)/2]);
    hold off;
end
set_font('Arial'); set_fontsize(16);
set_figsize(230,75);

%% Plot 2e: Minimum detectable amplitude (mV) at a fixed detection probability
% Figure 22 fixes the amplitude and shows P(detect) over (sigma_L, tau). This
% inverts it: fix the detection probability and show, for every (sigma_L, tau),
% the SMALLEST injected amplitude that reaches it, in mV.
%
% For each (sigma_L, tau) the tested amplitudes give a curve P(alpha). P must
% rise with alpha, but finite reps make it jitter, so it is forced
% non-decreasing (cummax) before the crossing is located. The threshold is then
% interpolated in log(alpha) -- the alpha grid is roughly geometric (0.5 1 2 4
% 8 14 %) -- which gives a continuous map instead of one quantised to the six
% simulated amplitudes.
%
% Two edge cases are kept distinct rather than both collapsing to NaN:
%   * already detected at the SMALLEST alpha: the true threshold lies below the
%     simulated range. The pixel is censored at alpha(1) and dotted white, so a
%     floor of 4.2 mV is not mistaken for a measurement.
%   * still undetected at the LARGEST alpha: left NaN (grey) = not reachable
%     with the amplitudes simulated.
pTargets     = [0.75];      % detection probabilities to map
contourMV    = [40 56 67];          % contour levels drawn on every tile, mV
pSmooth      = 1.5;                 % Gaussian sigma, in GRID CELLS, applied to P over
                                    % (sigma_L, tau) BEFORE the threshold is read off. 0 = off.
                                    % Each P is a binomial estimate from nReps, so one lucky
                                    % rep at the smallest amplitude flips a pixel to "already
                                    % detected" and it is drawn at the floor -- those are the
                                    % isolated dark squares. The true detection surface is
                                    % smooth in both sigma_L and tau, so averaging neighbours
                                    % removes that sampling noise at its source; smoothing the
                                    % THRESHOLD map instead would blur real structure and
                                    % leave the censored outliers behind. sigma = 1 cell pools
                                    % ~9 neighbours, i.e. ~3x lower binomial SE (~9x the reps).
mvper1perdFF = 100/12;              % 12% dF/F per 100 mV, Pojeong's paper
detVar       = 'either_detect_prob';% peak OR spike at the injected ROI
assert(ismember(detVar, injectionResult.Properties.VariableNames), ...
       '%s not found in injectionResult -- re-run the injection loop', detVar);

tolA = 1e-12; tolT = 1e-12; tolL = 1e-9;
nL = numel(lengths); nT = numel(taus); nA = numel(alphas);
alpha_mV = alphas(:) * 100 * mvper1perdFF;      % tested amplitudes, mV
logA     = log(alpha_mV);

Pgrid = nan(nT, nL, nA);                        % P(tau, sigma_L, alpha)
for iTau = 1:nT
    for iL = 1:nL
        for iA = 1:nA
            idx = abs(injectionResult.ell_um - lengths(iL)) < tolL & ...
                  abs(injectionResult.tau    - taus(iTau))  < tolT & ...
                  abs(injectionResult.alpha  - alphas(iA))  < tolA;
            if any(idx)
                Pgrid(iTau,iL,iA) = mean(injectionResult.(detVar)(idx), 'omitnan');
            end
        end
    end
end

% Which (tau, sigma_L) cells were simulated at all. Recorded BEFORE smoothing,
% because the NaN-aware kernel below would otherwise fill an unsimulated cell
% from its neighbours and invent a threshold where nothing was ever run.
hadData = squeeze(any(isfinite(Pgrid), 3));
if ~all(hadData(:))
    fprintf(['[detect] %d of %d (tau, sigma_L) cells have no simulated reps at any ' ...
             'amplitude -> NaN (grey)\n'], nnz(~hadData), numel(hadData));
end

% ---- smooth P over the (sigma_L, tau) plane, one amplitude at a time -------
% NaN-aware: the kernel is renormalised by the weight that actually landed on
% simulated pixels, so edges and gaps are not pulled toward zero.
if pSmooth > 0
    hw = max(1, ceil(2*pSmooth));
    [gx, gy] = meshgrid(-hw:hw, -hw:hw);
    Kp = exp(-(gx.^2 + gy.^2) / (2*pSmooth^2));  Kp = Kp / sum(Kp(:));
    for iA = 1:nA
        A = Pgrid(:,:,iA);  v = isfinite(A);  A(~v) = 0;
        den = conv2(double(v), Kp, 'same');
        S   = conv2(A, Kp, 'same') ./ den;
        S(den <= 0 | ~hadData) = NaN;     % never fill a cell that was not simulated
        Pgrid(:,:,iA) = S;
    end
    fprintf('[detect] P smoothed over (sigma_L, tau) with a %.2g-cell Gaussian before inversion\n', pSmooth);
end

thrMV  = nan(nT, nL, numel(pTargets));
censLo = false(nT, nL, numel(pTargets));
for ip = 1:numel(pTargets)
    pT = pTargets(ip);
    for iTau = 1:nT
        for iL = 1:nL
            p = squeeze(Pgrid(iTau,iL,:));
            if all(isnan(p)), continue; end
            p = cummax(fillmissing(p, 'linear', 'EndValues', 'nearest'));
            if p(1) >= pT                                  % threshold below range
                thrMV(iTau,iL,ip) = alpha_mV(1);  censLo(iTau,iL,ip) = true;  continue
            end
            k = find(p >= pT, 1, 'first');
            if isempty(k), continue; end                   % never reached -> NaN
            if p(k) > p(k-1)
                w = (pT - p(k-1)) / (p(k) - p(k-1));
                thrMV(iTau,iL,ip) = exp(logA(k-1) + w*(logA(k) - logA(k-1)));
            else
                thrMV(iTau,iL,ip) = alpha_mV(k);
            end
        end
    end
end

figure(23); clf;
tiledlayout(1, numel(pTargets), 'TileSpacing','compact','Padding','compact');
climMV = [min(thrMV(:),[],'omitnan') max(thrMV(:),[],'omitnan')];
if ~(climMV(2) > climMV(1)), climMV = climMV(1)*[1 1+eps]; end
for ip = 1:numel(pTargets)
    nexttile; hold on;
    M = thrMV(:,:,ip);
    imagesc(lengths/2, taus, M, 'AlphaData', isfinite(M));
    set(gca,'Color',[.75 .75 .75]);                        % NaN = not reachable
    axis xy; axis tight; colormap(turbo); clim(climMV); set(gca,'ColorScale','log');
    [Cc,hc] = contour(lengths/2, taus, M, sort(contourMV), 'k-', 'LineWidth', 1);
    if ~isempty(Cc), clabel(Cc, hc, 'FontSize', 9, 'Color', 'k'); end
    plot(ell_cable, taus, 'k--', 'LineWidth', 2);
    [yc, xc] = find(censLo(:,:,ip));                       % censored at the floor
    if ~isempty(yc), plot(lengths(xc)/2, taus(yc), 'w.', 'MarkerSize', 4); end
    title(sprintf('P(detect) \\geq %g%%', 100*pTargets(ip)));
    xlabel('\sigma_L (µm)');
    if ip == 1, ylabel('\tau (ms)'); else, set(gca,'YTick',[]); end
    if ip == numel(pTargets)
        cb = colorbar; cb.Label.String = sprintf('Minimum detectable\namplitude (mV)');
    end
    xlim([0 max(lengths)/2]);
    hold off;
end
set_font('Arial'); set_fontsize(16);
set_figsize(230,75);

fprintf('\n[detect] amplitude grid: %s %%dF/F  =  %s mV\n', ...
        num2str(alphas*100, '%.2g '), num2str(alpha_mV', '%.1f '));
for ip = 1:numel(pTargets)
    M = thrMV(:,:,ip);
    nCens    = nnz(censLo(:,:,ip));
    nNoData  = nnz(~hadData);                            % nothing simulated there
    nUnreach = nnz(~isfinite(M) & hadData);              % simulated, target never reached
    fprintf('[detect] P >= %3g%%: threshold %.1f-%.1f mV over the grid\n', ...
            100*pTargets(ip), min(M(:),[],'omitnan'), max(M(:),[],'omitnan'));
    fprintf('         NaN (grey): %d px never reach %g%% even at %.1f mV, %d px not simulated\n', ...
            nUnreach, 100*pTargets(ip), alpha_mV(end), nNoData);
    fprintf('         censored:   %d px already detected at the %.1f mV floor (white dots)\n', ...
            nCens, alpha_mV(1));
    % same number read along the cable-spread line, one value per tau
    onCable = arrayfun(@(iTau) interp1(lengths(:)/2, M(iTau,:)', ell_cable(iTau), ...
                                       'linear', NaN), (1:nT)');
    fprintf('         on the cable-spread line: %.1f-%.1f mV (median %.1f)\n', ...
            min(onCable,[],'omitnan'), max(onCable,[],'omitnan'), median(onCable,'omitnan'));
end

%% Plot 2d: Detection probability at cable-spread length
% Same box overlay
% x = tau
% y = detection probability
mvper1perdFF=100/12; %12% per 100 mV, Pojeong's paper
sigma_L_fixed_um = 56;                  % fixed cable-spread length sigma_L [um]
ell_fixed_um     = 2 * sigma_L_fixed_um; % full injection length [um]
lengths_vec = lengths(:);                % nLength x 1
[~, nearest_iL] = min(abs(lengths_vec - ell_fixed_um));         % closest tested length (scalar)
nearest_ell_um   = repmat(lengths_vec(nearest_iL), 1, numel(taus));  % same length for all tau
nearest_sigma_um = nearest_ell_um ./ 2;          % = sigma_L (fixed)

% Build detection probability matrices
prob_cable = struct();

tol_alpha = 1e-12; tol_tau   = 1e-12; tol_ell   = 1e-9;
det_types_2c     = {'peak_detect_prob','spike_detect_prob_injROI'};
det_types_str_2c = {'peaks','spikes'};
for iDet = 1:numel(det_types_2c)
    pm = nan(numel(alphas), numel(taus));
    for iTau = 1:numel(taus)
        ell_use = nearest_ell_um(iTau);
        for iA = 1:numel(alphas)
            idx = abs(injectionResult.ell_um - ell_use)      < tol_ell & ...
                abs(injectionResult.tau    - taus(iTau))   < tol_tau & ...
                abs(injectionResult.alpha  - alphas(iA))   < tol_alpha;
            if ~any(idx)
                continue;
            end
            pm(iA, iTau) = mean(injectionResult.(det_types_2c{iDet})(idx), 'omitnan');
        end
    end
    prob_cable.(det_types_str_2c{iDet}) = pm;
end

% P(peak OR spike): use saved per-rep boolean OR if available, otherwise
% fall back to inclusion-exclusion upper bound (peaks + spikes, clamped).
if ismember('either_detect_prob', injectionResult.Properties.VariableNames)
    pm = nan(numel(alphas), numel(taus));
    for iTau = 1:numel(taus)
        ell_use = nearest_ell_um(iTau);
        for iA = 1:numel(alphas)
            idx = abs(injectionResult.ell_um - ell_use)      < tol_ell & ...
                abs(injectionResult.tau    - taus(iTau))   < tol_tau & ...
                abs(injectionResult.alpha  - alphas(iA))   < tol_alpha;
            if ~any(idx), continue; end
            pm(iA, iTau) = mean(injectionResult.either_detect_prob(idx), 'omitnan');
        end
    end
    prob_cable.either = pm;
else
    warning('either_detect_prob not found — using min(peaks+spikes,1) as upper bound. Re-run injection loop to get exact values.');
    prob_cable.either = min(prob_cable.peaks + prob_cable.spikes, 1);
end

figure(21); clf;
showInd    = tau_ms_vec > 0 & tau_ms_vec < 120;
showtau    = tau_ms_vec(showInd);
alpha_pct  = alphas(:) * 100;
cmap_alpha = flipud(copper(numel(alphas)));

tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% --- Left tile: peak (circle) and spike (triangle) separately ---
nexttile; hold on;
for iA = 1:numel(alphas)
    c = cmap_alpha(iA,:);
    plot(showtau, prob_cable.peaks(iA, showInd), ...
        'Color', c, 'MarkerFaceColor', c, 'MarkerEdgeColor', c, ...
        'LineWidth', 1.8, 'MarkerSize', 5, 'Marker', 'o');
    plot(showtau, prob_cable.spikes(iA, showInd), ...
        'Color', c, 'MarkerFaceColor', c, 'MarkerEdgeColor', c, ...
        'LineWidth', 1.8, 'MarkerSize', 5, 'Marker', '^');
end
xlabel('\tau (ms)');
ylabel('P_{detect}');
%title('Peak  \circ  /  Spike  \Delta');
ylim([0 1.15]); xlim([0 23]); box off;

% --- Right tile: P(peak OR spike) ---
nexttile; hold on;
for iA = 1:numel(alphas)
    c = cmap_alpha(iA,:);
    plot(showtau, prob_cable.either(iA, showInd), ...
        'Color', c, 'MarkerFaceColor', c, 'MarkerEdgeColor', c, ...
        'LineWidth', 1.8, 'MarkerSize', 5, 'Marker', 's');
end
xlabel('\tau (ms)');
%ylabel(sprintf('P_{detect} (\\tau, \\alpha, \\sigma_L = %3.0f µm)', nearest_sigma_um(1)));
ylabel(sprintf('Probability of \n peak or spike detection'));
ylim([0 1.15]); xlim([0 70]); box off;
legend(num2str([alpha_pct alpha_pct*mvper1perdFF],'%2.1f%% ∆F/F\\approx%3.0f mV'), 'Location', 'northoutside', 'Box', 'off','NumColumns',1);

set_font('Arial');
set_fontsize(16);
set_figsize(160, 120);



%% Plot 3: Gallery of injection waveforms (1 row: 3 spatial + 1 temporal panel)
% Left 3 tiles: spatial footprint overlaid on F0 for 3 representative lengths
% Right tile (spanning 2): all 5 taus overlaid, normalized waveform

% Pick 3 representative length indices: small / mid / large
L_gallery_idx = [2 4];   % -> lengths [1 61 141] um

cmap_tau_gal = nebula(numel(taus));   % one color per tau

figure(19); clf;
tiledlayout(length(L_gallery_idx),1, 'TileSpacing', 'loose', 'Padding', 'loose');

% --- 3 spatial footprints ---
for k = 1:length(L_gallery_idx)
    iL = L_gallery_idx(k);
    nexttile([1 1]); hold on;
    imagesc(F0); colormap(gray); axis image off;
    S_show = S_all(:,:,iL);
    rgb_ov = cat(3, ones(nY,nX), zeros(nY,nX), zeros(nY,nX));
    h_ov = imagesc(rgb_ov);
    set(h_ov, 'AlphaData', S_show * 0.9);
    scatter(skel_col(anchor_idx), skel_row(anchor_idx)+20,40 ,'y^', 'filled');
    title(sprintf('2 * \\sigma_L = %2.0f µm', lengths(iL)));
    set(gca,'ydir','reverse')
    hold off;
end
drawScaleBar(100/Pixelsize(f),'horizontal','color','w','linewidth',2)
set_font('Arial'); set_fontsize(16);
set_figsize(130, 109);

figure(20); clf;
% --- Temporal waveforms: all 5 taus overlaid ---
nT_show   = 300;
t_rel_show = (0:nT_show-1)' - 50;   % start 50 frames before onset
hold all;
for iTau = 1:3:numel(taus)
    tau_k  = taus(iTau);
    h_wave = zeros(nT_show, 1);
    mask_w = t_rel_show > 0;
    tr_w   = t_rel_show(mask_w);
    h_wave(mask_w) = (tr_w / tau_k) .* exp(1 - tr_w / tau_k);
    if max(h_wave) > 0, h_wave = h_wave / max(h_wave); end
    plot(t_rel_show, h_wave, 'Color', cmap_tau_gal(iTau,:), 'LineWidth', 2, ...
        'DisplayName', sprintf('\\tau = %d ms', tau_k));
end
xline(0, 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
xlabel('Time relative to onset (ms)');
ylabel('h(\it{t})');
legend('Location', 'northeast');
xlim([t_rel_show(1) t_rel_show(end)]);
ylim([0 1.05]);
grid on; box off;
hold off;

set_font('Arial'); set_fontsize(16);
set_figsize(127, 100);