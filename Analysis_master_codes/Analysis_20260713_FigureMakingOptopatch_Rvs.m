% Analysis_20260713_FigureMakingOptopatch_Rvs
% =====================================================================
% Optopatch figure-making master script (reorganized).
%
% Data model
%   OpResult : {WvfCode, RepeatID, NeuronID} cell array. Each cell is the
%              per-recording Result struct plus derived fields (dendaxis,
%              Rheobase, ScaleFactor, State, ...). WvfCode encodes the
%              stimulation region x waveform (see the dictionaries below).
%   OpTable  : one row per loaded recording, a human-readable index into
%              OpResult (NeuronID / StimRegion / Waveform / RepeatID /
%              State / FileInd / path / hasRobustdFF). Built in Section 1.
%
% Normalization
%   Voltage traces are converted to dF/F by dividing by the robust dF/F
%   scale factor, RobustdFFfit.ScaleFactor (loaded per recording from
%   <fpath>/RobustdFFfit.mat), replacing the older ./F0_PCA convention.
%   See Analysis_20260523_RobustdFF_optopatch.m for how it is produced.
%
% Table of contents
%   Section 0  Paths, metadata, constants
%   Section 1  Build OpResult + OpTable (ScaleFactor-normalized; recordings
%              selected by the spreadsheet isUse flag, gated on RobustdFFfit)
%   Section 2  Kymograph interactive labeling
%   Section 3  Load / save pre-built data
%   ---- Figures ----
%   Fig 4  Short-pulse Soma vs Dendrite  (LabelMat_SP)
%   Fig 5  Pulse & Ramp                  (LabelMat)
%   Fig S  Soma vs Dendrite stim (complex-spike fraction)
%   Misc   Representative neurons / traces / interactive browser
%   Archive  Retired exploratory blocks (kept commented for reference)
% =====================================================================

%% Section 0. Paths, metadata, constants
clear
clc;
cd '/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/Statistics_Optopatch_Prism';
[~, ~, raw] = xlsread(['/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/' ...
    'Prism_OptopatchData_Arrangement.xlsx'], 'Sheet1', 'C5:Z192');

save_dir='/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Figures/invivoPrism/FigureOptopatch';
data_dir='/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/Statistics_Optopatch_Prism';

% --- Spreadsheet columns (raw column 1 == worksheet column C) ---
isUse=cell2mat(raw(:,19));
fpath=raw(:,1);
Mouse=cell2mat(raw(:,2));
CamType=raw(:,3);
NeuronInd=cell2mat(raw(:,5));
StimROI=raw(:,6);
StimWfn=raw(:,7);
State=raw(:,9);                 % Awake / Iso brain state (col K)
StructureData=raw(:,10);
isGoodCell=cell2mat(raw(:,11));
PixelSize=cell2mat(cellfun(@(x) (str2num(num2str(x))),raw(:,12),'UniformOutput',false));
refROI=cellfun(@(x) (str2num(num2str(x))),raw(:,14),'UniformOutput',false);
maintrunkROI=cellfun(@(x) (str2num(num2str(x))),raw(:,15),'UniformOutput',false);
TimeSegFrame=cell2mat(cellfun(@(x) (str2num(num2str(x))),raw(:,11),'UniformOutput',false));

place_bin=150; time_segment=15000; overlap=200;
alignedMovFN = {'STA_Mat_SS','STA_Mat_CS','STA_Mat_dSP'};
bound=6;
title_str={'Basal','Apical','Peri-Soma'};
[~, unqInd] = unique([Mouse NeuronInd] ,'row');
set(0,'DefaultFigureWindowStyle','docked')

% --- Stimulation label dictionaries ---
% StimROI / StimWfn spreadsheet strings are matched against these keys.
StimROI_Ind = {'Soma','Distal Dend','WF'};     % region keys
StimWfn_Ind = {'Ramp Stim','Short Pulse'};     % waveform keys
RegionNames   = {'Soma','DistalDend','WF'};    % tidy region names
WaveformNames = {'Ramp','ShortPulse'};         % tidy waveform names
% WvfCode: 1=Soma/Ramp 2=Soma/Pulse 3=Dend/Ramp 4=Dend/Pulse 5=WF/Ramp 6=WF/Pulse
%   region  = ceil(WvfCode/2)      -> RegionNames
%   waveform= 2-mod(WvfCode,2)     -> WaveformNames (odd=Ramp, even=Pulse)
wvf_patterns = [1 0 0 1 0; 1 0 0 0 1; 0 1 0 1 0; 0 1 0 0 1; 0 0 1 1 0; 0 0 1 0 1];
wvf_values   = [1; 2; 3; 4; 5; 6];
%% Section 1. Build OpResult + OpTable (RobustdFFfit.ScaleFactor-normalized)
% Recording selection is driven by the spreadsheet "Use?" flag (isUse, col U):
% every recording flagged isUse==1 is a candidate. Neurons are the unique
% (Mouse, NeuronInd) cells among those recordings, ordered by first appearance
% in the sheet; that order defines NeuronID (the 3rd OpResult index). A
% recording is actually loaded only if <fpath>/RobustdFFfit.mat is ready;
% others are skipped and listed in OpTable_skipped.
useNeurons = unique([Mouse(isUse>0) NeuronInd(isUse>0)], 'rows', 'stable');
nNeuron = size(useNeurons,1);

OpResult = cell(6, 1, nNeuron);
OpEntry = struct([]);        % accumulates one struct per loaded recording
skipEntry = struct([]);      % recordings skipped for missing RobustdFFfit
kLoaded = 0; kSkipped = 0;

for g2 = 1:nNeuron                            % neuron (3rd OpResult index)
    SameCellInd = find(Mouse==useNeurons(g2,1) & NeuronInd==useNeurons(g2,2));
    isSoma = ~cellfun(@isempty,strfind(StimROI(SameCellInd), StimROI_Ind{1}));
    isDend = ~cellfun(@isempty,strfind(StimROI(SameCellInd), StimROI_Ind{2}));
    isWF   = ~cellfun(@isempty,strfind(StimROI(SameCellInd), StimROI_Ind{3}));
    isRamp = ~cellfun(@isempty,strfind(StimWfn(SameCellInd), StimWfn_Ind{1}));
    isSP   = ~cellfun(@isempty,strfind(StimWfn(SameCellInd), StimWfn_Ind{2}));

    ROIwvf_ind = [isSoma isDend isWF isRamp isSP];
    % Load only recordings the sheet marks isUse==1 that also match a pattern.
    validind = find(sum(ROIwvf_ind,2)>=2 & isUse(SameCellInd)>0);

    g = ones(1,6);                            % per-WvfCode repeat counter
    for j = 1:length(validind)
        f2read = SameCellInd(validind(j));
        wfn = wvf_values(ismember(wvf_patterns, ROIwvf_ind(validind(j),:), 'rows'));
        if isempty(wfn); continue; end

        % --- Gate: require the robust dF/F fit to be ready ---
        robustFile = fullfile(fpath{f2read}, 'RobustdFFfit.mat');
        if ~isfile(robustFile)
            kSkipped = kSkipped + 1;
            skipEntry(kSkipped).NeuronID = g2;
            skipEntry(kSkipped).MouseID  = Mouse(f2read);
            skipEntry(kSkipped).FileInd  = f2read;
            skipEntry(kSkipped).WvfCode  = wfn;
            skipEntry(kSkipped).fpath    = string(fpath{f2read});
            fprintf('[skip] no RobustdFFfit: neuron %d, file %d (%s)\n', g2, f2read, fpath{f2read});
            continue;
        end

        load(fullfile(fpath{f2read},'OP_Result.mat'),'Result');
        RobustdFFfit = importdata(robustFile);

        r = g(wfn);                           % repeat index for this WvfCode
        Result.fpath        = fpath{f2read};
        Result.fileInd      = f2read;
        Result.pixelsize    = PixelSize(f2read);
        Result.maintrunkROI = maintrunkROI{f2read};
        Result.NeuronIndex  = g2;
        Result.MouseNumber  = Mouse(f2read);
        Result.State        = State{f2read};          % Awake / Iso
        Result.ScaleFactor  = RobustdFFfit.ScaleFactor;   % robust dF/F normalizer

        % --- Spike-triggered-average sign check / auto-correction ---
        % Some recordings have their normTraces polarity inverted, which makes
        % the bAP read as a negative deflection (and yields negative transmission
        % ratios downstream). Build the spike-triggered average of the NORMALIZED
        % trace (normTraces./ScaleFactor, the quantity every figure uses), pick
        % the strongest-responding ROI, and require its peri-spike peak to be
        % positive; if it is negative, flip the sign of normTraces.
        Result.SignFlipped = false;
        staWin = 10;
        if isfield(Result,'spike') && sum(Result.spike(1,:)>0) >= 3
            nrm = Result.normTraces ./ Result.ScaleFactor;
            STAsc = get_STA(nrm, Result.spike(1,:), staWin, staWin);
            STAsc = STAsc - median(STAsc(:,1:staWin), 2);      % pre-spike baseline
            pk = max(abs(STAsc), [], 2); pk(~isfinite(pk)) = -Inf;
            [~, refROIsc] = max(pk);                            % dominant spiking ROI
            [~, tpk] = max(abs(STAsc(refROIsc,:)));
            if isfinite(STAsc(refROIsc,tpk)) && STAsc(refROIsc,tpk) < 0
                Result.normTraces = -Result.normTraces;
                Result.SignFlipped = true;
                fprintf('[sign-flip] neuron %d, file %d (%s): normTraces polarity corrected.\n', ...
                    g2, f2read, fpath{f2read});
            end
        end

        % Signed distance-from-soma axis (um).
        Dsign = ones(1,size(Result.interDendDist,2));
        Dsign(Result.dist_order(1:find(Result.dist_order==1)-1)) = -1;
        Result.dendaxis = Result.interDendDist(1,:).*Dsign * Result.pixelsize;

        % Optical rheobase (ramp waveforms only: Soma/Dend/WF ramps).
        if ismember(wfn,[1 3 5])
            bwBlue = bwlabel(Result.Blue>0);
            bluePeriod = regionprops(bwBlue,'Area'); bluePeriod = cat(1, bluePeriod.Area);
            blueRampN = find(bluePeriod>2500);
            RampSpike = find(Result.spike(1,:).*(bwBlue==blueRampN));
            if ~isempty(RampSpike)
                Rheobase = mean(Result.Blue(RampSpike(1:3)));
            else
                PulseSpike = find(Result.spike(1,:).*(bwBlue>0));
                Rheobase = mean(Result.Blue(PulseSpike(1:3)));
            end
            Result.Rheobase = Rheobase;
            Result.RheobaseBlue = Result.Blue/Rheobase;
        end

        OpResult{wfn, r, g2} = Result;

        % --- OpTable row (index into OpResult) ---
        kLoaded = kLoaded + 1;
        OpEntry(kLoaded).NeuronID     = g2;
        OpEntry(kLoaded).MouseID      = Mouse(f2read);
        OpEntry(kLoaded).OrigNeuronInd= NeuronInd(f2read);
        OpEntry(kLoaded).FileInd      = f2read;
        OpEntry(kLoaded).WvfCode      = wfn;
        OpEntry(kLoaded).StimRegion   = string(RegionNames{ceil(wfn/2)});
        OpEntry(kLoaded).Waveform     = string(WaveformNames{2-mod(wfn,2)});
        OpEntry(kLoaded).RepeatID     = r;
        OpEntry(kLoaded).State        = string(State{f2read});
        OpEntry(kLoaded).PixelSize    = PixelSize(f2read);
        OpEntry(kLoaded).hasRobustdFF = true;
        OpEntry(kLoaded).SignFlipped  = Result.SignFlipped;
        OpEntry(kLoaded).fpath        = string(fpath{f2read});

        g(wfn) = g(wfn) + 1;
    end
end

% Assemble the readable index tables.
if isempty(OpEntry)
    warning('No recordings loaded (no isUse-flagged recording has RobustdFFfit.mat).');
    OpTable = table();
else
    OpTable = struct2table(OpEntry);
    OpTable.StimRegion = categorical(OpTable.StimRegion);
    OpTable.Waveform   = categorical(OpTable.Waveform);
    OpTable.State      = categorical(OpTable.State);
    OpTable = sortrows(OpTable, {'NeuronID','WvfCode','RepeatID'});
end
if kSkipped>0
    OpTable_skipped = struct2table(skipEntry);
else
    OpTable_skipped = table();
end
nFlipped = 0; if ~isempty(OpEntry), nFlipped = sum([OpEntry.SignFlipped]); end
fprintf('OpResult built: %d recordings loaded, %d skipped (missing RobustdFFfit), %d sign-corrected.\n', ...
    kLoaded, kSkipped, nFlipped);

%% Section 2. Kymograph labels (currently USE ALL; interactive labeling TODO)
% Builds OpResultLabel (one row per blue-stim kymograph) to match the current
% OpResult, then marks every kymograph as used (Kymo2use all true). When you
% are ready to hand-pick kymographs again, run InteractivelabelImages (see
% below) instead of the all-ones line. Must be run after Section 1 so indices
% into OpResult stay consistent.
KymoCell=[]; OpResultLabel=[];
[Protocol Repeat Neuron] = ind2sub(size(OpResult), find(cellfun(@(x) ~isempty(x),OpResult)));
g=1;
for img=1:length(Protocol)
    RstTemp=OpResult{Protocol(img),Repeat(img),Neuron(img)};
    bwBlue=bwlabel(RstTemp.Blue>0);       % label ALL blue stim periods; filter (e.g. last-6 / >100ms) downstream
    for b=1:max(bwBlue)
        KymoFrame=find(bwBlue==b);
        KymoCell{g} = RstTemp.normTraces./RstTemp.ScaleFactor;
        KymoCell{g} = KymoCell{g}(RstTemp.dist_order,KymoFrame);
        OpResultLabel(g,:)=[Protocol(img),Repeat(img),Neuron(img),b];
        g=g+1;
    end
end
%Kymo2use = InteractivelabelImages(KymoCell);   % TODO: redo interactive labeling here
Kymo2use = ones(1,numel(KymoCell));              % for now: use ALL kymographs
fieldName={'Stim_Region_Waveform','Repeat','Neuron','BlueStim_Order'};
if ~istable(OpResultLabel)                       % avoid nested table-of-tables on re-run
    OpResultLabel=array2table(OpResultLabel,'VariableNames',fieldName);
end

%% Section 3. Save / load pre-built data (LEGACY LOAD - off by default)
% The default workflow is Sections 1-2 (fresh build + use-all kymo labels).
% The saved bundle below uses the OLD neuron ordering and predates
% ScaleFactor/isUse, so loading it would overwrite the Section 1 OpResult with
% incompatible data. Leave the load commented unless you deliberately want the
% legacy bundle (in which case skip Sections 1-2 and also convert OpResultLabel
% to a table with the array2table line).
%--- Save the freshly built data (uncomment to overwrite) ---
save(fullfile(data_dir,'Prism_optogenetic_Results_20260713.mat'),'OpResult','Kymo2use','OpResultLabel','OpTable','-v7.3')
%--- Legacy load (leave commented for the Sections 1-2 workflow) ---
load(fullfile(data_dir,'Prism_optogenetic_Results_20260713.mat'))
fieldName={'Stim_Region_Waveform','Repeat','Neuron','BlueStim_Order'};
% Only convert if the loaded OpResultLabel is a plain numeric array; the saved
% bundle already stores it as a table, and re-running array2table on a table
% produces a nested table-of-tables (breaks OpResultLabel.Stim_Region_Waveform).
if ~istable(OpResultLabel)
    OpResultLabel=array2table(OpResultLabel,'VariableNames',fieldName);
end
%% Show representative neuron
figure(20); clf; tiledlayout(2,1,'padding','compact');
nexttile([1 1])
i=72;
load(fullfile(fpath{i},'OP_Result.mat'))
imshow2(Result.Structure',[0 4]); hold all
plot(Result.bluePatt{1, 1}(:,1),Result.bluePatt{1, 1}(:,2),'color',[0 0.6 1])
view([180 90])

nexttile([1 1])
i=74;
load(fullfile(fpath{i},'OP_Result.mat'))
imshow2(Result.Structure',[0 4]); hold all
plot(Result.bluePatt{1, 1}(:,1),Result.bluePatt{1, 1}(:,2),'color',[0 0.6 1])
set(gca,'ydir','reverse')
view([180 90])

%% Show STA trace
i=[73]; nTau=[-15:15];
load(fullfile(fpath{i},'OP_Result.mat'))
RobustdFFfit=importdata(fullfile(fpath{i},'RobustdFFfit.mat'));
taxis=[1:size(Result.normTraces,2)]/1000;

noi=[1:size(Result.ftprnt,3)];
noi_dist=ismember(Result.dist_order,noi);

Dsign=ones(1,size(Result.interDendDist,2));
Dsign(Result.dist_order(1:find(Result.dist_order==1)-1))=-1;
dendaxis=Result.interDendDist(1,Result.dist_order(noi_dist)).*Dsign(Result.dist_order(noi_dist));

normResidue=SeeResiduals(permute(Result.normTraces,[3 1 2]),Result.mc);
normResidue=SeeResiduals(permute(Result.normTraces,[3 1 2]),Result.mc.^2);
normResidue=SeeResiduals(permute(Result.normTraces,[3 1 2]),Result.mc(:,1).*Result.mc(:,end));
normResidue=squeeze(normResidue);
%F0PCA=get_F0PCA(normResidue);

normTr=normResidue./RobustdFFfit.ScaleFactor;
normTr=pcafilterTrace(normTr,[1:15]);
normTr=normTr(Result.dist_order(noi_dist),:);

STAtr=get_STA(normTr,Result.spike(1,:),-nTau(1),nTau(end));
STAtr=STAtr-median(STAtr(:,1:5),2);
figure(31); clf;
l=plot(nTau,STAtr');
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(turbo(33),2));

% %% Stimulation subthreshold figure;
% bound=6;
% STAmovie=[]; g=1; BlueonSpike=[]; BlueBoundary=[];
% for i=[73 75]
%     load(fullfile(fpath{i},'OP_Result.mat'))
%     mov=readBinMov_BHL(fpath{i});
%     mov_res= mov-mean(mov,3);
%     bkg = zeros(2, size(mov,3));
%     bkg(1,:) = linspace(-1, 1, size(mov,3));  % linear term
%     bkg(2,:) = linspace(-1, 1, size(mov,3)).^2;  % quadratic term
%     mov_res = SeeResiduals(mov_res,Result.mc);
%     mov_res = SeeResiduals(mov_res,Result.mc.^2);
%     mov_res = SeeResiduals(mov_res,Result.mc(:,1).*Result.mc(:,end));
%     mov_res= SeeResiduals(mov_res,bkg,1);
%
%     bwBlue=bwlabel(Result.Blue);
%     spvec=Result.spike;
%     [~, normTrMat spind]=get_STA(Result.normTraces(1,:),Result.spike(1,:),1,1);
%     [~, shift]=max(normTrMat,[],3);
%     spTime=spind+shift-2;
%     spvec_shifted=ind2vec(size(Result.normTraces,2),spTime,1);
%
%     BlueonSpike=[];
%     for b=1:max(bwBlue)
%         sptime_tmp=find(spvec_shifted(1,find(bwBlue==b,1)+[0:50])>0,1);
%         BlueonSpike=[BlueonSpike sptime_tmp+find(bwBlue==b,1)-1];
%     end
%
%     statmp = get_STA(tovec(mov_res),ind2vec(size(mov_res,3),BlueonSpike,1),50,50);
%     STAmovie(:,:,:,g)=toimg(statmp,size(mov_res,1),size(mov_res,2));
%
%     BlueBoundary{g}=(Result.blueDMDimg);
%     g=g+1;
% end
%
% mov_filt=imresize(pcafilt(imresize(mov_res,1/5),3),5);
% lag1mean=abs(sqrt(mean(mov_filt(:,:,2:end).*mov_filt(:,:,1:end-1),3)));
% mov_res_shrink=imresize(movmedian(mov_res,10,3),1/5);
% [V D]=get_eigvector(tovec(mov_res_shrink),10);
% F0img=sqrt(sum((V.^2).*D(1:10)',2));
% F0img=toimg(F0img,size(mov_res_shrink,1),size(mov_res_shrink,2));
% F0img=imresize(F0img,5);
% %F0img=get_F0img(mov_res);
%
% STAmov_norm=-STAmovie./F0img;
% STAmov_norm=STAmov_norm-mean(STAmov_norm(:,:,[1:15],:),3);
%
% Rfixed = imref2d([size(Result.ref_im,1) size(Result.ref_im,2)]);
% inverseTform = invert(Result.tform);
% %revertedStruct = imwarp(Result.Structure, inverseTform,'OutputView',Rfixed);
% revertedStruct= Result.Structure;
% revertedStruct(revertedStruct==0)=prctile(revertedStruct(:),30);
% revertedStruct=mat2gray(revertedStruct);
%
% STAmovie_normStr=[];
% crop_roi=[94 6  300  159];
% STAmov_norm_crop=STAmov_norm(crop_roi(2):crop_roi(2)+crop_roi(4),crop_roi(1):crop_roi(1)+crop_roi(3),:,:);
% revertedStruct_crop=revertedStruct(crop_roi(2):crop_roi(2)+crop_roi(4),crop_roi(1):crop_roi(1)+crop_roi(3));
% cax_sub=[0.004 0.01];
% cax_sp=[0.004 0.025];
% STAmov_norm_crop_filt=[];
% for spclass_ind=1:2
%     %STAmovie_norm{spclass_ind}=imgaussfilt3(STAnorm_sub./F_refImg,[1.5 1.5 0.1]);%.*SkelDend(bound:end-bound,bound:end-bound);
%     STAmov_norm_crop_filt{spclass_ind}=imgaussfilt(pcafilt(STAmov_norm_crop(:,:,:,spclass_ind),15),4);
%     colorSTA=grs2rgb(tovec(STAmov_norm_crop_filt{spclass_ind}),colormap('turbo'),cax_sub(1),cax_sub(2));
%     colorSTA=reshape(colorSTA,size(STAmov_norm_crop,1),size(STAmov_norm_crop,2),size(STAmov_norm_crop,3),[]);
%     colorSTA=permute(colorSTA,[1 2 4 3]);
%
%     colorSTA2=grs2rgb(tovec(STAmov_norm_crop_filt{spclass_ind}),colormap('turbo'),cax_sp(1),cax_sp(2));
%     colorSTA2=reshape(colorSTA2,size(STAmov_norm_crop,1),size(STAmov_norm_crop,2),size(STAmov_norm_crop,3),[]);
%     colorSTA2=permute(colorSTA2,[1 2 4 3]);
%
%     STAmovie_normStr{spclass_ind}=colorSTA.*revertedStruct_crop*3;%.*SkelDend(bound:end-bound,bound:end-bound);
%     STAmovie_normStr2{spclass_ind}=colorSTA2.*revertedStruct_crop*3;%.*SkelDend(bound:end-bound,bound:end-bound);
% end
%
% % sptype={'SomStim','ddStim'};
% % cax=[-0.005,0.02];
% % writeMov4d(fullfile(save_dir,['STA_dFFStructgrsrgb_ShortPulse']),[imrotate(STAmovie_normStr{1},90) imrotate(STAmovie_normStr{2},90)],[-50:50],10,1,cax)
%
% figure(21); clf; tiledlayout(2,3)
% t_show=[37:41]; t_show_spike=[51:53]; ax1=[];
% for stimROI=1:2
%     subimage=grs2rgb(mean(STAmov_norm_crop_filt{stimROI}(:,:,t_show),3),colormap(turbo),cax_sub(1),cax_sub(2));
%     subimage=subimage.*revertedStruct_crop*3;
%     spimage=grs2rgb(max(STAmov_norm_crop_filt{stimROI}(:,:,t_show_spike),[],3),colormap(turbo),cax_sp(1),cax_sp(2));
%     spimage=spimage.*revertedStruct_crop*3;
%
%     ax1=[ax1 nexttile([1 1])];
%     imshow2(imrotate([Result.Structure(crop_roi(2):crop_roi(2)+crop_roi(4),crop_roi(1):crop_roi(1)+crop_roi(3))],90),[]); hold all
%     Blbd=bwboundaries(imrotate(BlueBoundary{stimROI}(crop_roi(2):crop_roi(2)+crop_roi(4),crop_roi(1):crop_roi(1)+crop_roi(3)),90));
%     plot(Blbd{1}(:,2),Blbd{1}(:,1),'color',[0 0.6 1])
%     nexttile([1 1]);
%     imshow2(imrotate([subimage],90),[]); hold all
%     cb=colorbar;
%     cb.Ticks=[0 0.5 1];
%     cb.TickLabels=[cax_sub(1) mean(cax_sub) cax_sub(2)];
%     cb.Label.String = '\DeltaF/F';
%
%     nexttile([1 1]);
%     imshow2(imrotate([spimage],90),[]); hold all
%     cb=colorbar;
%     cb.Ticks=[0 0.5 1];
%     cb.TickLabels=[cax_sp(1) mean(cax_sp) cax_sp(2)];
%     cb.Label.String = '\DeltaF/F';
%     %plot(BlueBoundary{1}(:,2),BlueBoundary{1}(:,1),'color',[0 0.5 1])
% end
% colormap(turbo);
% colormap(ax1(1),gray);
% colormap(ax1(2),gray);
% %%
% for n=1:size(OpResult,3)
%     [nonemptyr nonemptyc]=find(cellfun(@(x) ~isempty(x), OpResult(:,:,n)));
%     catTr=[];
%     if ~isempty(nonemptyr)
%     for e=1:length(nonemptyr)
%         normsubTrtmp=get_subthreshold(OpResult{nonemptyr(e),nonemptyc(e),n}.normTraces,OpResult{nonemptyr(e),nonemptyc(e),n}.spike(1,:),7,17);
%         catTr=[catTr normsubTrtmp];
%     end
%     [F0PCA nPC]=get_F0PCA(catTr);
%     if nPC==1
%     [F0PCA nPC]=get_F0PCA(catTr,5);
%     end
%     for e=1:length(nonemptyr)
%        OpResult{nonemptyr(e),nonemptyc(e),n}.F0_PCA=F0PCA;
%     end
%     end
% end
%% Short pulses Soma vs Dend
nTau=[60 100];
SP_STA=[]; distFromSoma=[]; SP_STA_std=[];
dendriteaxis_bin=[-280:50:600];
dendriteaxis_bin2=[-300 -50 50 140 460];
perisomadist=[-70 70]; cmap=[0 0 0; 1 0 0];  cmap_light=[0.7 0.7 0.7; 1 0.7 0.7];
prespikesubtime=[-10:-3]; dax=[];
LabelMat_SP=[]; AlignedbAPall_SP=[]; dax_all=[]; g=1;
for n=1:size(OpResult,3) %neuron
    for wvf=[2 4] %stim region
        for rep=1:size(OpResult,2) %repeat 
            if ~isempty(OpResult{wvf,rep,n})
                normTr=OpResult{wvf,rep,n}.normTraces./OpResult{wvf,rep,n}.ScaleFactor;
                %noi=OpResult{wvf,rep,n}.maintrunkROI; % only main trunks
                noi=[1:size(OpResult{wvf,rep,n}.normTraces,1)]; % all ROIs
                dax=OpResult{wvf,rep,n}.dendaxis(noi);
                MouseInd=Mouse(OpResult{wvf,rep,n}.fileInd); fileInd=OpResult{wvf,rep,n}.fileInd;
                spikeheight=get_STA(mean(normTr(dax>dendriteaxis_bin2(2) & dax <= dendriteaxis_bin2 (3),:),1,'omitnan'),OpResult{wvf,rep,n}.spike(1,:),0,0);
                %normTr=normTr./spikeheight;

                nROI=length(noi);
                nTime=size(OpResult{wvf,rep,n}.normTraces,2);
                noi_dist=find(ismember(OpResult{wvf,rep,n}.dist_order,noi));
                dOrder=OpResult{wvf,rep,n}.dist_order(noi_dist);

                [~, blueOff2, blueOff3] = get_blueoffTrace(zeros(1,nTime), OpResult{wvf,rep,n}.Blue>0, 25, 15);
                blueDialated=(1-blueOff2);  bwBlue=bwlabel(blueDialated>0);
                [bb blueonset]=unique(bwBlue);  blueonset(bb==0)=[];
                som_spike=max(OpResult{wvf,rep,n}.spike([1],:),[],1).*(blueDialated>0);
                som_spike_blueIndexed=max(som_spike,[],1).*(bwBlue);

                SubLarge=get_subthreshold(normTr, max([som_spike; bwBlue>0],[],1),50,500); %subthreshold
                SubTr=get_subthreshold(normTr, max([som_spike],[],1),7,17); %subthreshold
                normTrSub=normTr-SubLarge;

                [STAbAP, AlignedbAP_SP, Sptime]=get_STA(normTrSub,som_spike>0,nTau(1),nTau(2));
                AlignedbAP_SP=permute(AlignedbAP_SP,[1 3 2]);
                AlignedbAP_SP_preSubSubtracted=AlignedbAP_SP-mean(AlignedbAP_SP(:,nTau(1)+prespikesubtime,:),2,'omitnan');
                AlignedbAPall_SP=[AlignedbAPall_SP; squeeze(num2cell(AlignedbAP_SP(dOrder,:,:),[1 2]))];
                dax_all{g,1}=dax; dax_all{g,2}=dax(OpResult{wvf,rep,n}.dist_order);
                sptimefromblueonset=[]; sporderinblue=[]; BluePulseN=[];
                for s=1:length(Sptime)
                    sptimefromblueonset(s)=Sptime(s)-blueonset(som_spike_blueIndexed(Sptime(s)));
                    som_spike_blueIndexed_sub=find(som_spike_blueIndexed==bwBlue(Sptime(s)));
                    sporderinblue(s)=find(som_spike_blueIndexed_sub==Sptime(s));
                    BluePulseN(s)=bwBlue(Sptime(s));
                end

                [~, tmax]=max(STAbAP(:,nTau(1)+[-2:4]),[],2);
                tmax=tmax+nTau(1)-3;
                [AUCbAP, AUCrawbAP_SP]=get_AUC(AlignedbAP_SP,repmat(tmax,1,1,size(AlignedbAP_SP,3)),2,3);
                [~, AUCrawbAP_SP_short, ~, AmprawbAP_SP_short]=get_AUC(AlignedbAP_SP,repmat(tmax,1,1,size(AlignedbAP_SP,3)),1,1);
                [~, AUCrawbAP_SP_subtracted]=get_AUC(AlignedbAP_SP_preSubSubtracted,repmat(tmax,1,1,size(AlignedbAP_SP,3)),2,3);

                AUCbAP_SP_cell = num2cell(AUCrawbAP_SP, 1); % Change this line to use AUC vs AUC raw for Tx rate
                AUCbAP_SP_cell = cellfun(@(x) [dax' x(:,1)],AUCbAP_SP_cell,'UniformOutput',false);
                AUCshortbAP_SP_cell = num2cell(AUCrawbAP_SP_short, 1);
                AUCshortbAP_SP_cell = cellfun(@(x) [dax' x(:,1)],AUCshortbAP_SP_cell,'UniformOutput',false);
                AmpbAP_SP_cell = num2cell(AmprawbAP_SP_short, 1);
                AmpbAP_SP_cell = cellfun(@(x) [dax' x(:,1)],AmpbAP_SP_cell,'UniformOutput',false);
                KinkbAP_SP_cell = num2cell(AUCrawbAP_SP_subtracted, 1);
                KinkbAP_SP_cell = cellfun(@(x) [dax' x(:,1)],KinkbAP_SP_cell,'UniformOutput',false);

                [~, AUCbAPbin_SP_cell, dcenter, dbin_index] = zscore_binning(AUCbAP_SP_cell, dendriteaxis_bin2);
                [~, AUCshortbAPbin_SP_cell] = zscore_binning(AUCshortbAP_SP_cell, dendriteaxis_bin2);
                [~, AmpbAPbin_SP_cell] = zscore_binning(AmpbAP_SP_cell, dendriteaxis_bin2);
                [~, KinkbAPbin_SP_cell] = zscore_binning(KinkbAP_SP_cell, dendriteaxis_bin2);
                emptycell=cellfun(@isempty,AUCbAPbin_SP_cell);

                preSub_SP_Soma=squeeze(mean(AlignedbAP_SP(find(dbin_index{2,1}),nTau(1)+prespikesubtime,:),[1 2],'omitnan'));
                preSub_SP_Apical=squeeze(mean(AlignedbAP_SP(find(dbin_index{4,1}),nTau(1)+prespikesubtime,:),[1 2],'omitnan'));
                trnsmit_SP_bAPAUC=NaN(length(Sptime),1);
                trnsmit_SP_bAPAUC(find(~emptycell))=cellfun(@(x,y) x(4,2)/y(2,2),AUCbAPbin_SP_cell(~emptycell),AUCshortbAPbin_SP_cell(~emptycell)); %Divided by short AUC

                bAPAUC_SP_Apical=cellfun(@(x) x(4,2),AUCbAPbin_SP_cell(~emptycell));
                bAPAUC_SP_Soma=cellfun(@(x) x(2,2),AUCshortbAPbin_SP_cell(~emptycell));
                bAPKink_SP_Soma=cellfun(@(x) x(2,2),KinkbAPbin_SP_cell(~emptycell));
                bAPKink_SP_Apical=cellfun(@(x) x(4,2),KinkbAPbin_SP_cell(~emptycell));
                bAPAmp_SP_Soma=cellfun(@(x) x(2,2),AmpbAPbin_SP_cell(~emptycell));
                bAPAmp_SP_Apical=cellfun(@(x) x(4,2),AmpbAPbin_SP_cell(~emptycell));
                ISI=[NaN diff(Sptime)]';

                catlab=[repmat(n,length(Sptime),1) Sptime' sptimefromblueonset' ...
                    BluePulseN' repmat(MouseInd,length(Sptime),1) repmat(fileInd,length(Sptime),1) repmat(wvf,length(Sptime),1) repmat(rep,length(Sptime),1)...
                    sporderinblue' trnsmit_SP_bAPAUC bAPAUC_SP_Apical bAPAUC_SP_Soma ISI preSub_SP_Apical preSub_SP_Soma bAPKink_SP_Soma bAPKink_SP_Apical, bAPAmp_SP_Soma, bAPAmp_SP_Apical];
                LabelMat_SP=[LabelMat_SP; catlab];  %Session number, Neuron #, Spike Frame, Stimulation Pattern, Rheobase, Spike Order, Mouse#
                g=g+1;
            end
        end
    end
end
fieldName={'NeuronID','SpikeFrame','SpikeFrameBlueonset','StimOrder','MouseID','FileID','StimRegion','RepeatID',...
    'SpikeOrder','TransmissionRatio','AUC_apical','AUC_Soma','ISI','PreSub_apical','PreSub_soma','Kink_soma','Kink_apical','Amp_Soma','Amp_apical'};
LabelMat_SP=array2table(LabelMat_SP,'VariableNames',fieldName);

for n=unique(LabelMat_SP.NeuronID)'
    nindex=find(LabelMat_SP.NeuronID==n);
    LabelMat_SP.AUC_apical_norm(nindex)=zscore(LabelMat_SP.AUC_apical(nindex));
end

%% Monitor kymographs of negative transmission-ratio spikes (QC before figure 4)
% Each row of LabelMat_SP is 1:1 with AlignedbAPall_SP (both appended per spike
% in the Short-pulse section), so a spike with TransmissionRatio < 0 has bAP
% kymograph AlignedbAPall_SP{idx}. Rows are dist-ordered ROIs (y = distance from
% soma), columns are peri-spike time. Shown as paged montages so each flagged
% event can be eyeballed; title lists TX ratio and the apical/soma AUCs that
% produce it (TX < 0 usually means a negative apical AUC).
nTauSP=[60 100]; tSP=-nTauSP(1):nTauSP(2);     % must match the Short-pulse section
regionName=containers.Map([2 4],{'Soma','Dend'});
perPage=12;                                     % kymographs per figure window

negIdx=find(LabelMat_SP.TransmissionRatio<0);
fprintf('%d / %d short-pulse spikes have TransmissionRatio < 0.\n', numel(negIdx), height(LabelMat_SP));
if numel(AlignedbAPall_SP)~=height(LabelMat_SP)
    warning('AlignedbAPall_SP (%d) and LabelMat_SP (%d) are out of sync; re-run the Short-pulse section.', ...
        numel(AlignedbAPall_SP), height(LabelMat_SP));
end

if ~isempty(negIdx)
    allv=cell2mat(cellfun(@(c) c(:), AlignedbAPall_SP(negIdx),'UniformOutput',false));
    caxSP=prctile(allv,[2 98]);                 % shared robust color scale
    for pg=1:ceil(numel(negIdx)/perPage)
        figure(140+pg); clf;
        tiledlayout('flow','Padding','compact','TileSpacing','compact');
        pageIdx=negIdx((pg-1)*perPage+1:min(pg*perPage,numel(negIdx)));
        for ii=pageIdx'
            sr=LabelMat_SP.StimRegion(ii); rp=LabelMat_SP.RepeatID(ii); nid=LabelMat_SP.NeuronID(ii);
            kymo=AlignedbAPall_SP{ii};
            nexttile;
            imagesc(tSP,1:size(kymo,1),kymo,caxSP); colormap(turbo); hold on;
            xline(0,'w:');
            if ~isempty(OpResult{sr,rp,nid})
                set_kymoYtick(OpResult{sr,rp,nid}.dendaxis(OpResult{sr,rp,nid}.dist_order));
            end
            title(sprintf('#%d N%d %s\nTX=%.2f Aap=%.2f As=%.2f', ii, nid, regionName(sr), ...
                LabelMat_SP.TransmissionRatio(ii), LabelMat_SP.AUC_apical(ii), LabelMat_SP.AUC_Soma(ii)), 'FontSize',8);
            xlabel('peri-spike (ms)');
        end
    end
end

% --- Interactive one-by-one browser (uncomment to step through with Enter) ---
% for ii=negIdx'
%     sr=LabelMat_SP.StimRegion(ii); rp=LabelMat_SP.RepeatID(ii); nid=LabelMat_SP.NeuronID(ii);
%     figure(139); clf; imagesc(tSP,1:size(AlignedbAPall_SP{ii},1),AlignedbAPall_SP{ii},caxSP); colormap(turbo); hold on; xline(0,'w:');
%     set_kymoYtick(OpResult{sr,rp,nid}.dendaxis(OpResult{sr,rp,nid}.dist_order));
%     title(sprintf('#%d N%d %s TX=%.2f', ii, nid, regionName(sr), LabelMat_SP.TransmissionRatio(ii)));
%     xlabel('peri-spike (ms)'); input('next?\n');
% end

%% Show spike order and TX ratio (figure 4)

timebin=[0 5 20 40 60 100];
maxTXratio=6; stimstr={'Soma','Dendrite'}; g=1; cmap_sporder=gray(10); cmap_sporder=cmap_sporder(1:2:10,:);
minTXratio=-1;
% figure(14): paired Soma (black) vs Apical-dendrite (gray) transmission
% ratio, grouped by evoked spike order (1st/2nd/3rd). One point per neuron
% (mean over its spikes); gray lines connect the same neuron across the two
% stimulation sites. Significance: signrank within each order (soma vs dend)
% and across orders for soma stimulation.
orders2show = 1:3;
pairGap = 1; groupGap = 1.6; bw = 0.6; jitter = 0.10;
faceCol = {[1 1 1],[0.6 0.6 0.6]};   edgeCol = {[0 0 0],[0.5 0.5 0.5]};
faceA   = [1 0.35];                   ptCol   = {[0.2 0.2 0.2],[0.45 0.45 0.45]};
neurons = unique(LabelMat_SP.NeuronID)';
Vsoma = nan(numel(neurons),numel(orders2show));   % per-neuron mean, soma (wvf 2)
Vdend = nan(numel(neurons),numel(orders2show));   % per-neuron mean, dend (wvf 4)
xS = zeros(1,numel(orders2show)); xD = zeros(1,numel(orders2show));

figure(14); clf; hold on;
for si = orders2show
    xc = (si-1)*(pairGap+groupGap);
    xS(si) = xc+1; xD(si) = xc+1+pairGap;
    for k = 1:numel(neurons)
        n = neurons(k);
        iS = LabelMat_SP.NeuronID==n & LabelMat_SP.StimRegion==2 & abs(LabelMat_SP.TransmissionRatio)<maxTXratio & (LabelMat_SP.TransmissionRatio)>minTXratio & LabelMat_SP.SpikeOrder==si;
        iD = LabelMat_SP.NeuronID==n & LabelMat_SP.StimRegion==4 & abs(LabelMat_SP.TransmissionRatio)<maxTXratio & (LabelMat_SP.TransmissionRatio)>minTXratio & LabelMat_SP.SpikeOrder==si;
        if any(iS), Vsoma(k,si) = mean(LabelMat_SP.TransmissionRatio(iS),'omitnan'); end
        if any(iD), Vdend(k,si) = mean(LabelMat_SP.TransmissionRatio(iD),'omitnan'); end
    end
    % paired connecting lines (neurons measured at both sites)
    for k = 1:numel(neurons)
        if ~isnan(Vsoma(k,si)) && ~isnan(Vdend(k,si))
            plot([xS(si) xD(si)],[Vsoma(k,si) Vdend(k,si)],'-','Color',[0.82 0.82 0.82],'LineWidth',1);
        end
    end
    % boxes + jittered points for each stimulation site
    dataCond = {Vsoma(:,si), Vdend(:,si)}; xpos = [xS(si) xD(si)];
    for cc = 1:2
        v = dataCond{cc}; v = v(~isnan(v)); xx = xpos(cc);
        if isempty(v), continue; end
        q1 = prctile(v,25); med = prctile(v,50); q3 = prctile(v,75); iqrv = q3-q1;
        lo = max(min(v),q1-1.5*iqrv); hi = min(max(v),q3+1.5*iqrv);
        patch([xx-bw/2 xx+bw/2 xx+bw/2 xx-bw/2],[q1 q1 q3 q3],faceCol{cc}, ...
            'FaceAlpha',faceA(cc),'EdgeColor',edgeCol{cc},'LineWidth',1);
        plot([xx-bw/2 xx+bw/2],[med med],'-','Color',edgeCol{cc},'LineWidth',2);
        plot([xx-bw/10 xx+bw/20],[hi hi],'-','Color',edgeCol{cc},'LineWidth',1);
        plot([xx-bw/10 xx+bw/20],[lo lo],'-','Color',edgeCol{cc},'LineWidth',1);
        plot([xx xx],[q3 hi],'-','Color',edgeCol{cc},'LineWidth',1);
        plot([xx xx],[q1 lo],'-','Color',edgeCol{cc},'LineWidth',1);
        plot(xx+jitter*(2*rand(size(v))-1),v,'o','MarkerSize',6, ...
            'MarkerFaceColor',ptCol{cc},'MarkerEdgeColor','none');
    end
end

% --- axis cosmetics ---
ymin = min([Vsoma(:);Vdend(:)],[],'omitnan');
ymax = max([Vsoma(:);Vdend(:)],[],'omitnan');
set(gca,'XTick',(xS+xD)/2,'XTickLabel',counting_string(orders2show));
xlabel('Evoked spike order'); ylabel('Transmission ratio'); box off;
xlim([xS(1)-0.9 xD(end)+0.9]);
% text(0.98,0.99,'Soma stimulation','Units','normalized','HorizontalAlignment','right', ...
%     'VerticalAlignment','top','Color',edgeCol{1},'FontWeight','bold');
% text(0.98,0.92,'Apical dendrite stimulation','Units','normalized','HorizontalAlignment','right', ...
%     'VerticalAlignment','top','Color',edgeCol{2});

% --- significance: soma vs dend within each order (local star) ---
starOf = @(p) repmat('*',1,(p<0.05)+(p<0.01)+(p<0.001));
for si = orders2show
    cs = ~isnan(Vsoma(:,si)) & ~isnan(Vdend(:,si));
    if sum(cs) >= 2
        p = signrank(Vsoma(cs,si),Vdend(cs,si));
        if p < 0.05
            yl = max([Vsoma(cs,si);Vdend(cs,si)]) + 0.15;
            plot([xS(si) xS(si) xD(si) xD(si)],[yl yl+0.05 yl+0.05 yl],'k-','LineWidth',1);
            text((xS(si)+xD(si))/2,yl+0.06,starOf(p),'HorizontalAlignment','center', ...
                'VerticalAlignment','bottom','FontSize',12);
        end
    end
end
% --- significance: soma across spike orders (stacked brackets on top) ---
yb = ymax + 0.5; step = 0.32; lvl = 0; somaPairs = [1 2; 1 3; 2 3];
for r = 1:size(somaPairs,1)
    a = somaPairs(r,1); b = somaPairs(r,2);
    cs = ~isnan(Vsoma(:,a)) & ~isnan(Vsoma(:,b));
    if sum(cs) >= 2
        p = signrank(Vsoma(cs,a),Vsoma(cs,b));
        if p < 0.05
            yy = yb + lvl*step;
            plot([xS(a) xS(a) xS(b) xS(b)],[yy yy+0.05 yy+0.05 yy],'k-','LineWidth',1);
            text((xS(a)+xS(b))/2,yy+0.06,starOf(p),'HorizontalAlignment','center', ...
                'VerticalAlignment','bottom','FontSize',12);
            lvl = lvl+1;
        end
    end
end
ylim([0.9 4.5]);

% --- report N (points / neurons / mice) ---
mouseOf = arrayfun(@(n) LabelMat_SP.MouseID(find(LabelMat_SP.NeuronID==n,1)), neurons);
fprintf('\nfigure(14) transmission-ratio summary (|TX|<%.1f):\n', maxTXratio);
regname = {'Soma','Dend'}; Vcond = {Vsoma, Vdend};
for cc = 1:2
    for si = orders2show
        has = ~isnan(Vcond{cc}(:,si)); wvfc = 2*cc;    % 2->Soma, 4->Dend
        nEv = sum(LabelMat_SP.StimRegion==wvfc & LabelMat_SP.SpikeOrder==si & abs(LabelMat_SP.TransmissionRatio)<maxTXratio);
        fprintf('  %-4s order %d: %2d points, %2d mice, %3d spikes\n', ...
            regname{cc}, si, sum(has), numel(unique(mouseOf(has))), nEv);
    end
end
anyHas = any(~isnan([Vsoma Vdend]),2);
nPointTot = sum(~isnan(Vsoma(:))) + sum(~isnan(Vdend(:)));   % total plotted dots
nNeuronTot = sum(anyHas);
nMouseTot  = numel(unique(mouseOf(anyHas)));
nEvTot = sum(ismember(LabelMat_SP.StimRegion,[2 4]) & ismember(LabelMat_SP.SpikeOrder,orders2show) & abs(LabelMat_SP.TransmissionRatio)<maxTXratio);
fprintf('  TOTAL: %d points, %d neurons, %d mice, %d spikes\n', nPointTot, nNeuronTot, nMouseTot, nEvTot);
text(0.02,0.98,sprintf('%d neurons, %d mice', nNeuronTot, nMouseTot), 'Units','normalized', ...
    'HorizontalAlignment','left','VerticalAlignment','top','Color',[0.2 0.2 0.2]);

set_fontsize(16); set_font('Arial');
set_figsize(100,120);
% 
% figure(17); clf; g=1; cmap_stimregion=[0 0 0; 0.5 0.5 0.5];
% tiledlayout(1,3,'Padding','tight');
% sp2show_sporder_meanN=NaN(max(LabelMat_SP.NeuronID),3,4);
% usedN=[];
% for wvf=[2 4]
%     for sporder=1:3
%         for n=unique(LabelMat_SP.NeuronID)'
%             show_ind=find(LabelMat_SP.NeuronID==n & LabelMat_SP.StimRegion==wvf & abs(LabelMat_SP.TransmissionRatio)<maxTXratio & LabelMat_SP.SpikeOrder==sporder); %Pulse & Soma
%             if isempty(show_ind)
%                 sp2show_sporder_meanN(n,sporder,wvf)=NaN;
%             else
%                 sp2show_sporder_meanN(n,sporder,wvf)=mean(LabelMat_SP.TransmissionRatio(show_ind),'omitnan');
%                 usedN=[usedN; [LabelMat_SP.NeuronID(show_ind) LabelMat_SP.MouseID(show_ind)]];
%             end
%         end
%         length(unique(LabelMat_SP.MouseID(find(LabelMat_SP.StimRegion==wvf & abs(LabelMat_SP.TransmissionRatio)<maxTXratio & LabelMat_SP.SpikeOrder==sporder))))
%     end
% end
% for sporder=1:3
%     nexttile([1 1]);
%     p=Boxplot_wPoints2(permute(sp2show_sporder_meanN(:,sporder,[2 4]),[1 3 2]),cmap_stimregion);
%     p123_soma=get_pValue(sp2show_sporder_meanN(:,[1 2 3],[2]),1);
%     drawPValueLines(p,0,'TextYOffset',0.1,'XCoord',[1:2]); box off;
%     set(gca,'XTick',[]);
%     xlabel(counting_string(sporder)); ylabel('Transmission ratio'); box off;
% end

%% Show spike order and bAP apical AUC (figure 4)

timebin=[0 5 20 40 60 100];
maxTXratio=6; stimstr={'Soma','Dendrite'}; g=1; cmap_sporder=gray(10); cmap_sporder=cmap_sporder(1:2:10,:);
figure(14); clf; tiledlayout(2,1);
for wvf=[2 4]
    nexttile([1 1]);
    % [M S t_center N SubSet]=binning_data({[LabelMat_SP.SpikeFrameBlueonset(show_ind) LabelMat_SP.TransmissionRatio(show_ind)]},timebin);
    % %p=get_pValue(SubSet,0);
    % scatter_heatmap(LabelMat_SP.SpikeFrameBlueonset(show_ind),LabelMat_SP.TransmissionRatio(show_ind),100,200);
    sp2show_sporder=[];
    for sporder=1:4
        % show_ind=find(LabelMat_SP.StimRegion==wvf & abs(LabelMat_SP.TransmissionRatio)<maxTXratio & LabelMat_SP.SpikeOrder==sporder); %Pulse & Soma
        show_ind=find(LabelMat_SP.NeuronID==5 & LabelMat_SP.StimRegion==wvf & abs(LabelMat_SP.TransmissionRatio)<maxTXratio & LabelMat_SP.SpikeOrder==sporder); %Pulse & Soma
        sp2show_sporder{sporder}=LabelMat_SP.AUC_apical_norm(show_ind);
    end
    p=Violin_wPoints(sp2show_sporder,cmap_sporder(1:4,:));
    % hold all;
    % errorbar(t_center, M, S ./ sqrt(cellfun(@sum,N)'), 'Color', [1 0 0], 'LineWidth', 1.5);
    %xlim([0 100]); ylim([0 4]);
    %drawPValueLines(p,0,'TextYOffset',0.1,'XCoord',t_center);
    xlabel('Time after blue onset (ms)'); ylabel('Transmission ratio'); box off;
    title([stimstr{g} ' Stim. (Pulse)']); caxis([0 0.0112]);
    colormap([gen_colormap([0 0 0; parula(10)])])
    g=g+1;
end

figure(17); clf; g=1; cmap_stimregion=[0 0 0; 0.5 0.5 0.5];
tiledlayout(1,3,'Padding','tight');
sp2show_sporder_meanN=NaN(max(LabelMat_SP.NeuronID),3,4);
usedN=[];
for wvf=[2 4]
    for sporder=1:3
        for n=unique(LabelMat_SP.NeuronID)'
            show_ind=find(LabelMat_SP.NeuronID==n & LabelMat_SP.StimRegion==wvf & abs(LabelMat_SP.TransmissionRatio)<maxTXratio & LabelMat_SP.SpikeOrder==sporder); %Pulse & Soma
            if isempty(show_ind)
                sp2show_sporder_meanN(n,sporder,wvf)=NaN;
            else
                sp2show_sporder_meanN(n,sporder,wvf)=mean(LabelMat_SP.AUC_apical_norm(show_ind),'omitnan');
                usedN=[usedN; [LabelMat_SP.NeuronID(show_ind) LabelMat_SP.MouseID(show_ind)]];
            end
        end
        length(unique(LabelMat_SP.MouseID(find(LabelMat_SP.StimRegion==wvf & abs(LabelMat_SP.TransmissionRatio)<maxTXratio & LabelMat_SP.SpikeOrder==sporder))))
    end
end
for sporder=1:3
    nexttile([1 1]);
    p=Boxplot_wPoints2(permute(sp2show_sporder_meanN(:,sporder,[2 4]),[1 3 2]),cmap_stimregion);
    p123_soma=get_pValue(sp2show_sporder_meanN(:,[1 2 3],[2]),1);
    drawPValueLines(p,0,'TextYOffset',0.1,'XCoord',[1:2]); box off;
    set(gca,'XTick',[]);
    xlabel(counting_string(sporder)); ylabel('bAP apical AUC (z-score)'); box off;
    ylim([-1.5 2])
end

%% Kink vs Amp vs AUC
figure(15); clf;
show_ind_somaStim=find(LabelMat_SP.StimRegion==2 & abs(LabelMat_SP.TransmissionRatio)<maxTXratio & LabelMat_SP.SpikeOrder==1 & LabelMat_SP.NeuronID==1);
show_ind_dendStim=find(LabelMat_SP.StimRegion==4 & abs(LabelMat_SP.TransmissionRatio)<maxTXratio & LabelMat_SP.SpikeOrder==1 & LabelMat_SP.NeuronID==1);
nexttile([1 1])
scatter(LabelMat_SP.PreSub_soma(show_ind_somaStim),LabelMat_SP.Kink_soma(show_ind_somaStim),'filled','MarkerFaceAlpha',0.5); hold all
scatter(LabelMat_SP.PreSub_soma(show_ind_dendStim),LabelMat_SP.Kink_soma(show_ind_dendStim),'filled','MarkerFaceAlpha',0.5); hold all
xlabel('Prespike sub. (soma)'); ylabel('Kink amplitude (soma)'); legend({'Soma stim.','Apical stim.'});
nexttile([1 1])
scatter(LabelMat_SP.PreSub_apical(show_ind_somaStim),LabelMat_SP.Kink_apical(show_ind_somaStim),'filled','MarkerFaceAlpha',0.5); hold all
scatter(LabelMat_SP.PreSub_apical(show_ind_dendStim),LabelMat_SP.Kink_apical(show_ind_dendStim),'filled','MarkerFaceAlpha',0.5); hold all
xlabel('Prespike sub. (apical)'); ylabel('Kink amplitude (apical)'); legend({'Soma stim.','Apical stim.'});

nexttile([1 1])
scatter(LabelMat_SP.PreSub_apical(show_ind_somaStim),LabelMat_SP.AUC_apical(show_ind_somaStim),'filled','MarkerFaceAlpha',0.5); hold all
scatter(LabelMat_SP.PreSub_apical(show_ind_dendStim),LabelMat_SP.AUC_apical(show_ind_dendStim),'filled','MarkerFaceAlpha',0.5); hold all
xlabel('Prespike sub. (apical)'); ylabel('AUC (apical)'); legend({'Soma stim.','Apical stim.'});

nexttile([1 1])
scatter(LabelMat_SP.PreSub_apical(show_ind_somaStim),LabelMat_SP.Amp_apical(show_ind_somaStim),'filled','MarkerFaceAlpha',0.5); hold all
scatter(LabelMat_SP.PreSub_apical(show_ind_dendStim),LabelMat_SP.Amp_apical(show_ind_dendStim),'filled','MarkerFaceAlpha',0.5); hold all
xlabel('Prespike sub. (apical)'); ylabel('Amplitude (apical)'); legend({'Soma stim.','Apical stim.'});

nexttile([1 1])
scatter(LabelMat_SP.Kink_soma(show_ind_somaStim),LabelMat_SP.Kink_apical(show_ind_somaStim),'filled','MarkerFaceAlpha',0.5); hold all
scatter(LabelMat_SP.Kink_soma(show_ind_dendStim),LabelMat_SP.Kink_apical(show_ind_dendStim),'filled','MarkerFaceAlpha',0.5); hold all
xlabel('Kink amplitude (soma)'); ylabel('Kink amplitude (apical)'); legend({'Soma stim.','Apical stim.'});

nexttile([1 1])
scatter(LabelMat_SP.Amp_Soma(show_ind_somaStim),LabelMat_SP.AUC_apical(show_ind_somaStim),'filled','MarkerFaceAlpha',0.5); hold all
scatter(LabelMat_SP.Amp_Soma(show_ind_dendStim),LabelMat_SP.AUC_apical(show_ind_dendStim),'filled','MarkerFaceAlpha',0.5); hold all
xlabel('Amplitude (soma)'); ylabel('AUC (apical)'); legend({'Soma stim.','Apical stim.'});

nexttile([1 1])
scatter(LabelMat_SP.PreSub_soma(show_ind_somaStim),LabelMat_SP.PreSub_apical(show_ind_somaStim),'filled','MarkerFaceAlpha',0.5); hold all
scatter(LabelMat_SP.PreSub_soma(show_ind_dendStim),LabelMat_SP.PreSub_apical(show_ind_dendStim),'filled','MarkerFaceAlpha',0.5); hold all
xlabel('Prespike sub. (soma)'); ylabel('Prespike sub. (apical)'); legend({'Soma stim.','Apical stim.'});

nexttile([1 1])
scatter(LabelMat_SP.AUC_apical(show_ind_somaStim),LabelMat_SP.Amp_apical(show_ind_somaStim),'filled','MarkerFaceAlpha',0.5); hold all
scatter(LabelMat_SP.AUC_apical(show_ind_dendStim),LabelMat_SP.Amp_apical(show_ind_dendStim),'filled','MarkerFaceAlpha',0.5); hold all
xlabel('AUC (apical)'); ylabel('Amplitude (apical)'); legend({'Soma stim.','Apical stim.'});

nexttile([1 1])
scatter(LabelMat_SP.Kink_apical(show_ind_somaStim),LabelMat_SP.AUC_apical(show_ind_somaStim),'filled','MarkerFaceAlpha',0.5); hold all
scatter(LabelMat_SP.Kink_apical(show_ind_dendStim),LabelMat_SP.AUC_apical(show_ind_dendStim),'filled','MarkerFaceAlpha',0.5); hold all
xlabel('Kink amplitude (apical)'); ylabel('AUC (apical)'); legend({'Soma stim.','Apical stim.'});
%% Show Representative STA (figure 4)
colormap(turbo); cmap_ExTr=gen_colormap(Plasma,10);
figure(14); clf;
figure(15); clf; tiledlayout(2,1,'padding','compact');
figure(16); clf; tiledlayout(2,1,'padding','tight','TileSpacing','compact');
show_neuron=11; ax2=[]; cax=[-0.01 0.12]; t2show=[-22 -10 1 2 3];
SWCpoints=OpResult{1,1,show_neuron}.SWC; SWCpoints(1,3)=70;

dax=OpResult{1,1,show_neuron}.dendaxis;
ROIvec=OpResult{1,1,show_neuron}.maintrunkROI; ftprint2show=OpResult{1,1,show_neuron}.ftprnt(:,:,ROIvec); ftprintall=OpResult{1,1,show_neuron}.ftprnt;
nROI=size(OpResult{1,1,show_neuron}.ftprnt,3); ftprint2show_trace=get_coord(ftprint2show);
somROI=1; dendROI=7;

% figure(16): precompute each SWC point's ROI (footprint) membership ONCE, so
% every peri-spike time point can be drawn in a single axis per stim region
% (2 subplots total) with an x offset, instead of one tile per time point.
% Mirrors showScaleScatter's point-in-footprint masking.
ftO16 = ftprintall(:,:,OpResult{1,1,show_neuron}.dist_order);
sz16 = [size(ftO16,1) size(ftO16,2)];
swc16 = SWCpoints; if size(swc16,2)<4, swc16(:,4)=1; end; swc16(swc16(:,4)==0,4)=1;
ptX16 = swc16(:,2); ptY16 = swc16(:,1); ptSz16 = swc16(:,4)+2; ptROI16 = zeros(size(swc16,1),1); ptSz16(1)=40;
valid16 = ptY16>0.5 & ptY16<sz16(2)-0.5 & ptX16>0.5 & ptX16<sz16(1)-0.5;
lin16 = nan(size(swc16,1),1); lin16(valid16) = sub2ind(sz16, round(ptX16(valid16)), round(ptY16(valid16)));
for r = 1:size(ftO16,3)
    m = ftO16(:,:,r)>0; hit = false(size(swc16,1),1); hit(valid16) = m(lin16(valid16)); ptROI16(hit) = r;
end
xStep16 = 0.97*(max(ptX16)-min(ptX16));   % horizontal spacing between time points

figure(14);
nexttile([1 1])
boundary_s=bwboundaries(ftprint2show(:,:,somROI)); boundary_ap=bwboundaries(ftprint2show(:,:,dendROI));
boundary_blue=OpResult{2,1,show_neuron}.bluePatt;
boundary_blue2=OpResult{4,1,show_neuron}.bluePatt;
showScaleScatter([1:nROI], SWCpoints, ftprintall,repmat([0 0 0],256,1)); hold all
plot(boundary_s{1}(:,1),boundary_s{1}(:,2),'color',cmap_ExTr(1,:),'linewidth',2);
plot(boundary_ap{1}(:,1),boundary_ap{1}(:,2),'color',cmap_ExTr(end,:),'linewidth',2);
plot(boundary_blue{1}(:,1),boundary_blue{1}(:,2),'color',[0 0.5 1],'linewidth',2);
plot(boundary_blue2{1}(:,1),boundary_blue2{1}(:,2),'color',[0 0.5 1],'linewidth',2);
plot(ftprint2show_trace(:,2),ftprint2show_trace(:,1),'r','linewidth',2);
drawScaleBar(100/OpResult{1,1,show_neuron}.pixelsize,'vertical')
view([180 -90])

for wvf=[2 4]

    figure(15);
    ax2=[ax2 nexttile([1 1])];
    sp2mean=LabelMat_SP.SpikeOrder==1 & LabelMat_SP.StimRegion==4 & LabelMat_SP.NeuronID==5;
    staSel=LabelMat_SP.SpikeOrder==1 & LabelMat_SP.StimRegion==wvf & LabelMat_SP.NeuronID==show_neuron;
    STA2showall=mean(cat(3,AlignedbAPall_SP{staSel}),3,'omitnan');
    stimNm={'Soma','Dend'};
    fprintf('figure(15) STA kymo: neuron %d, %s stim (wvf %d), 1st-order spikes averaged = %d\n', ...
        show_neuron, stimNm{wvf/2}, wvf, sum(staSel));
    STA2show=STA2showall(ismember(OpResult{1,1,show_neuron}.dist_order,ROIvec)',:);
    STA2show_interp=interp1(dax(ROIvec),STA2show,linspace(dax(ROIvec(1)),dax(ROIvec(end)),10));
    imagesc([-nTau(1):nTau(2)],ROIvec,STA2show_interp,cax);
    %imagesc([-nTau(1):nTau(2)],[1:size(STA2showall,1)],STA2show_interp,cax);
    set_kymoYtick(dax(ROIvec)); ylabel('Distance (µm)');

    % ax2=[ax2 nexttile([1 1])];
    % plot([-nTau(1):nTau(2)],STA2show(1,:),'color',cmap_ExTr(1,:)); hold all;
    % plot([-nTau(1):nTau(2)],STA2show(end,:),'color',cmap_ExTr(end,:));
    % ylim(cax); xlim([-50 50]);
    % axis off;
    % drawScaleBar(1,'vertical');
    % drawScaleBar(50,'horizontal');
    if wvf==4
        xlabel('Time (ms)');
    end
     cb=colorbar;
     cb.Label.String = '∆F/F';

    figure(16);
    nexttile([1 1]); hold on;
    for t=1:length(t2show)
        V2show=STA2showall(:,t2show(t)+nTau(1));
        cols16=vec2cmap(V2show,'turbo',cax); xoff16=(t-1)*xStep16;
        g0=ptROI16==0;
        scatter(ptX16(g0)+xoff16, ptY16(g0), ptSz16(g0), [0.5 0.5 0.5], 'filled');
        for r=1:max(ptROI16)
            sel=ptROI16==r;
            if any(sel), scatter(ptX16(sel)+xoff16, ptY16(sel), ptSz16(sel), cols16(r,:), 'filled'); end
        end
    end
    set(gca,'xdir','reverse');
    colormap(turbo); caxis(cax); axis equal tight off; view([180 -90]);
end
figure(15);
linkaxes(ax2,'x'); xlim([-50 50]);
set_font('Arial'); set_fontsize(24);
set_figsize(200,200);

figure(16);
set_figsize(200,162);
%% Show Short pulse stimulation voltage trace (Raw trace)
figure(21); clf; tiledlayout(12,6,'padding','tight'); ax1=[]; ax3=[];
dendriteaxis_bin=[-100:25:300]; cax=[-0.3 2]*0.01; t_show=[3 13]; cmap=[118 85 157; 194 102 56]/256;
SomaROI=[1]; ApicalROI=[6 7]; g=1;
ShowT_zoom=[315]+[-30:60];
for i=[73 75];
    load(fullfile(fpath{i},'OP_Result.mat'))
    RobustdFFfit=importdata(fullfile(fpath{i},'RobustdFFfit.mat'));
    taxis=[1:size(Result.normTraces,2)]/1000;

    noi=maintrunkROI{i};
    noi_dist=ismember(Result.dist_order,noi);

    Dsign=ones(1,size(Result.interDendDist,2));
    Dsign(Result.dist_order(1:find(Result.dist_order==1)-1))=-1;
    dendaxis=Result.interDendDist(1,Result.dist_order(noi_dist)).*Dsign(Result.dist_order(noi_dist));

    normResidue=permute(Result.normTraces,[3 1 2]);
    normResidue=SeeResiduals_tiling(normResidue,Result.mc,1,5000);
    normResidue=SeeResiduals_tiling(normResidue,Result.mc.^2,1,5000);
    normResidue=SeeResiduals_tiling(normResidue,Result.mc(:,1).*Result.mc(:,end),1,5000);
    normResidue=squeeze(normResidue);

    %F0PCA=get_F0PCA(normResidue);

    %normTr=normResidue./F0PCA;
    normTr=normResidue./RobustdFFfit.ScaleFactor;
    %normTr=normTr-movprc(normTr,200,30,2);
    %subTr=get_blueoffTrace(normTr(1,:),Result.Blue,50,70);
    %normTr=normTr-subTr;
    %normTr=normTr-movmedian(normTr,100,2);
    normTr=pcafilterTrace(normTr,[1:15]);
    normTr=normTr(Result.dist_order(noi_dist),:);
    %normTr=get_bandstop(normTr,1000,[49 80]);

    ax3=[ax3 nexttile((g-1)*30+1,[5 4])];
    plot(taxis,mean(normTr(SomaROI,:),1,'omitnan'),'color',cmap(1,:),'linewidth',1.5); hold all
    plot(taxis,mean(normTr(ApicalROI,:),1,'omitnan'),'color',cmap(2,:),'linewidth',1.5);
    %ylim([-0.015 0.03])
    axis off

    ax1=[ax1 nexttile((g-1)*30+5,[2 2])];
    [~, sortind]=sort(dendaxis,'ascend');
    [Xq, Yq] = meshgrid([taxis], dendriteaxis_bin);
    normTr_interp = interp2([taxis], dendaxis(sortind), normTr(sortind,:), Xq, Yq, 'linear');
    normTr_interpShow=normTr_interp(:,ShowT_zoom);
    normTr_interpShow=normTr_interpShow-mean(normTr_interpShow(:,1:10),2);
    pcolor(taxis(ShowT_zoom),dendriteaxis_bin,normTr_interpShow);
    caxis(cax); colormap(turbo(256)); shading flat;
    cb=colorbar;
    cb.Label.String='\DeltaF/F';
    set(gca,'XTick',t_show,'XTickLabel',t_show-t_show(1),'YDir','reverse')
    ylim([0 250])

    ax1=[ax1 nexttile((g-1)*30+17,[3 2])];
    normTr_show=[mean(normTr(SomaROI,ShowT_zoom),1,'omitnan'); mean(normTr(ApicalROI,ShowT_zoom),1,'omitnan')];
    normTr_show=normTr_show-mean(normTr_show(:,1:10),2);
    l=plot(taxis(ShowT_zoom),normTr_show','linewidth',1.5);
    arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(cmap,2)); axis tight off;
    ylim(cax)
    g=g+1;
end

ax2=nexttile([2 4]);
plot(taxis, Result.Blue);
linkaxes([ax2 ax3],'x');  xlim(t_show); axis off;
ax1=[ax1 nexttile(65,[2 2])];
plot(taxis(ShowT_zoom), Result.Blue(ShowT_zoom)); axis tight off;
linkaxes([ax1],'x')
%% Pulse & Ramp
nTau=[-20 20]; nTauPeriSp=[-3:7]; nTauPeriSp2=[-3:3];
g=ones(1,4); SP_STA=[]; distFromSoma=[]; LabelMat=[]; AlignedbAPall=[];
dendriteaxis_bin=[-200 -100 -50 50 100:50:500];
dendriteaxis_bin2=[-300 -50 50 140 460];
timebin=[0 30 35 70 100 200:100:500];
Rheobin=[1 1.3 1.6 2 3 4 5];
perisomadist=[-30 30]; cmap=[0 0 0; 1 0 0];  cmap_light=[0.7 0.7 0.7; 1 0.7 0.7];
Kymos=[]; AUCbAPall=[];
% Guard: if OpResultLabel was accidentally passed through array2table twice it
% becomes a nested table-of-tables, and OpResultLabel.Stim_Region_Waveform==wvf
% then returns a table (breaking find). Un-nest to a flat numeric-column table.
% No-op on a healthy OpResultLabel.
while istable(OpResultLabel) && any(varfun(@istable,OpResultLabel,'OutputFormat','uniform'))
    OpResultLabel = OpResultLabel.Variables;
end

SpikeAmp=[]; SpikeAUC=[]; dax=[]; g=1;
for wvf=[1]; % wvf=1 : soma stimulation, wvf=3 : Dendrite stimulation, wvf =5 : WF stimulation
    BlueStimN=OpResultLabel(find(Kymo2use'>0 & OpResultLabel.Stim_Region_Waveform==wvf),[3 4 2]); %Neuron#, bwBlue, Repeat
    fprintf('[Pulse&Ramp] wvf %d: %d blue-stim sessions\n', wvf, size(BlueStimN,1));
    for n=1:size(BlueStimN,1) % session
        cax=[];
        Neuron=BlueStimN.Neuron(n);
        rep=BlueStimN.Repeat(n);
        BluePulseN=BlueStimN.BlueStim_Order(n);
        fprintf('  wvf %d  session %2d/%2d  (neuron %d, rep %d, pulse %d)\n', ...
            wvf, n, size(BlueStimN,1), Neuron, rep, BluePulseN);
        noi=[1:size(OpResult{wvf,rep,Neuron}.normTraces,1)];
        dax{n,wvf}=OpResult{wvf,rep,Neuron}.dendaxis(noi);
        perisomaInd=find(dax{n,wvf}'>perisomadist(1) & dax{n,wvf}'<perisomadist(2));

        if ~isempty(OpResult{wvf,rep,Neuron})
            MouseInd=OpResult{wvf,rep,Neuron}.MouseNumber;
            fileInd=OpResult{wvf,rep,Neuron}.fileInd;
            normTr=OpResult{wvf,rep,Neuron}.normTraces./OpResult{wvf,rep,Neuron}.ScaleFactor; %Read traces
            %noi=OpResult{wvf,rep,Neuron}.maintrunkROI;
            normTr=normTr(noi,:);
            nROI=length(noi); nTime=size(OpResult{wvf,rep,Neuron}.normTraces,2); %set ROI and time
            dOrder=OpResult{wvf,rep,Neuron}.dist_order(noi);
            Classvec = get_Class2index(OpResult{wvf,rep,Neuron}.SpClass([1 2],:));
            som_spike = Classvec .* OpResult{wvf,rep,Neuron}.spike(1,:);
            bwBlue=bwlabel(OpResult{wvf,rep,Neuron}.Blue);   % all-periods labels (match OpResultLabel BlueStim_Order)
            % Ramp (wvf 1/3): drop the front pilot stims, keep only the last 6
            % stimulation periods. BluePulseN indexes all periods, so a pilot is
            % any period at or before max-6.
            if ismember(wvf,[1 3]) && BluePulseN <= max(bwBlue)-6
                continue;
            end

            [STAbAP]=get_STA(normTr,max(som_spike>0,[],1),-nTau(1),nTau(2));
            [~, tmax]=max(STAbAP(:,-nTau(1)+[-2:4]),[],2);
            tmax=tmax-nTau(1)-3;

            orderVector = get_BurstOrder(max(som_spike,[],1), 20);
            SpikeinBlue=max(som_spike,[],1).*(bwBlue==BluePulseN);
            SpikeinBlue_ind=find(SpikeinBlue);
            if isempty(SpikeinBlue_ind)   % no spike evoked in this blue period -> nothing to average
                fprintf('    (no spike in blue period %d; skipped)\n', BluePulseN);
                continue;
            end
            BlueOnset=find(bwBlue==BluePulseN,1);
            SubLarge=get_subthreshold(normTr, max([som_spike; bwBlue>0],[],1),200,1000); %subthreshold
            %SubLarge=get_subthreshold(normTr, max([som_spike],[],1), 7,50);
            normTrSub=normTr-SubLarge;
            kymoFrames=unique(find(bwBlue==BluePulseN)'+[-50:50]);
            kymoFrames(kymoFrames<1 | kymoFrames>nTime)=[];   % clip ±50 window to trace bounds
            Kymos{n,wvf}=normTr(OpResult{wvf,rep,Neuron}.dist_order,kymoFrames);

            [~, AlignedbAP, Sptime]=get_STA(normTrSub,SpikeinBlue>0,-nTau(1),nTau(2));
            % RheobaseBlue only exists for ramp sessions (Rheobase computed from
            % the ramp). For pulse sessions, express this session's blue in
            % rheobase units using the neuron's Rheobase from its soma-ramp
            % (wvf==1) session; NaN if no soma-ramp Rheobase is available.
            if isfield(OpResult{wvf,rep,Neuron},'RheobaseBlue')
                RheoBaseAmp=OpResult{wvf,rep,Neuron}.RheobaseBlue(Sptime);
            else
                somaRheobase=NaN;
                for rr=1:size(OpResult,2)
                    if ~isempty(OpResult{1,rr,Neuron}) && isfield(OpResult{1,rr,Neuron},'Rheobase')
                        somaRheobase=OpResult{1,rr,Neuron}.Rheobase; break;
                    end
                end
                RheoBaseAmp=OpResult{wvf,rep,Neuron}.Blue(Sptime)/somaRheobase;
            end
            IsCS=OpResult{wvf,rep,Neuron}.CStrace(Sptime);
            BurstOrder=orderVector(Sptime);
            Sptime_inBlue=Sptime-BlueOnset+1;
            AlignedbAP=permute(AlignedbAP,[1 3 2]);
            mat2cat=squeeze(num2cell(AlignedbAP(OpResult{wvf,rep,Neuron}.dist_order,:,:), [1 2]))';
            AlignedbAPall=[AlignedbAPall mat2cat];
            [AUCbAP, AUCrawbAP kinkAmp kink_raw]=get_AUC(AlignedbAP,repmat(tmax,1,1,size(AlignedbAP,3)),2,3);
            [AUCbAP_short, AUCrawbAP_short]=get_AUC(AlignedbAP,repmat(tmax,1,1,size(AlignedbAP,3)),1,1);

            AUCbAP_cell = num2cell(AUCrawbAP, 1); % Change this line to use AUC vs AUC raw for Tx rate
            AUCbAP_cell = cellfun(@(x) [dax{n,wvf}' x(:,1)],AUCbAP_cell,'UniformOutput',false);
            AUCbAPall=[AUCbAPall AUCbAP_cell];
            AmpbAP_cell = num2cell(AUCrawbAP_short, 1);
            AmpbAP_cell = cellfun(@(x) [dax{n,wvf}' x(:,1)],AmpbAP_cell,'UniformOutput',false);

            [~, AUCbAPbin_cell, dcenter] = zscore_binning(AUCbAP_cell, dendriteaxis_bin2);
            [~, AmpbAPbin_cell, dcenter] = zscore_binning(AmpbAP_cell, dendriteaxis_bin2);
            emptycell=cellfun(@isempty,AUCbAPbin_cell);

            centroid_bAPAUC=sum(dax{n,wvf}'.*AUCrawbAP,1,'omitnan')./sum(AUCrawbAP,1,'omitnan');
            trnsmit_bAPAUC=NaN(length(Sptime),1);
            trnsmit_bAPAUC(find(~emptycell))=cellfun(@(x,y) x(4,2)/y(2,2),AUCbAPbin_cell(~emptycell),AmpbAPbin_cell(~emptycell)); %Divided by short AUC
            trnsmit_bAPAUC2=NaN(length(Sptime),1);
            trnsmit_bAPAUC2(find(~emptycell))=cellfun(@(x) x(4,2)/x(2,2),AUCbAPbin_cell(~emptycell)); %Divided by same length AUC
            bAPAUC_Apical=cellfun(@(x) x(4,2),AUCbAPbin_cell(~emptycell));
            bAPAUC_Soma=cellfun(@(x) x(2,2),AmpbAPbin_cell(~emptycell));
            ISI=[NaN diff(Sptime)]';

            [SpType, ~]=find(som_spike(:,Sptime));

            %AUCbAP=zscore(AUCbAP,[],2);
            %AUCrawbAP=AUCrawbAP./mean(AUCrawbAP(perisomaInd,:),1,'omitnan');
            ztemp = num2cell(AUCrawbAP, 1);
            ztemp = cellfun(@(x) [dax{n,wvf}' x(:,1)],ztemp,'UniformOutput',false);
            SpikeAUC=[SpikeAUC ztemp];
            dax{n,wvf}=dax{n,wvf}(OpResult{wvf,rep,Neuron}.dist_order);

            catlab=[repmat(n,length(Sptime),1) repmat(Neuron,length(Sptime),1) Sptime_inBlue' ...
                repmat(BluePulseN,length(Sptime),1) RheoBaseAmp' [1:length(Sptime)]'  repmat(MouseInd,length(Sptime),1) repmat(fileInd,length(Sptime),1) repmat(wvf,length(Sptime),1) ...
                IsCS' BurstOrder' repmat(rep,length(Sptime),1) SpType trnsmit_bAPAUC bAPAUC_Apical bAPAUC_Soma ISI];
            LabelMat=[LabelMat; catlab];  %Session number, Neuron #, Spike Frame, Stimulation Pattern, Rheobase, Spike Order, Mouse#
        end
    end
end
fieldName={'SessionID','NeuronID','SpikeFrame','StimOrder','Rheobase','SpikeOrder','MouseID','FileID','StimRegion',...
    'InCS','BurstOrder','RepeatID','SpikeType','TransmissionRatio','AUC_apical','Amp_Soma','ISI'};
LabelMat=array2table(LabelMat,'VariableNames',fieldName);

NeuronIDs_unq=unique(LabelMat.NeuronID)';
NormConstant=NaN(1,max(NeuronIDs_unq));
NormConstant_vec=NaN(size(LabelMat,1),1);
NormConstantAA=NaN(1,max(NeuronIDs_unq));
NormConstant_vecAA=NaN(size(LabelMat,1),1);
for n=NeuronIDs_unq
    show_ind=find(LabelMat.StimOrder~=6 & LabelMat.SpikeFrame>=200 & LabelMat.StimRegion==1 ...
                & abs(LabelMat.TransmissionRatio)<7 & LabelMat.NeuronID==n); %Pulse & Soma
    NormConstant(n)=median(LabelMat.TransmissionRatio(show_ind));
    NormConstant_vec(LabelMat.NeuronID==n,1)=NormConstant(n);

    NormConstantAA(n)=mean(LabelMat.AUC_apical(show_ind));
    NormConstant_vecAA(LabelMat.NeuronID==n,1)=NormConstantAA(n);
end
LabelMat.TransmissionRatio_norm=LabelMat.TransmissionRatio./NormConstant_vec;
LabelMat.AUC_apical_norm=LabelMat.AUC_apical./NormConstant_vecAA;
%% Soma normalized transmission ratio over time (Figure 5D)
% The single best panel: soma stimulation (wvf==1), normalized transmission
% ratio vs time after blue onset. One dot per pulse spike (density-coloured)
% with the binned mean +/- SEM, plus pairwise statistical tests between time
% bins (ranksum; full pairwise matrix printed, adjacent-bin significance drawn).
timebin=[0 30 35 85 110 200:100:500]; maxTXratio=10;
badNeurons=[1 13 22 32];                 % neurons flagged noisy by the figure(27) diagnostic
figure(23); clf; hold on;
show_ind=find(LabelMat.StimOrder~=6 & LabelMat.StimRegion==1 & abs(LabelMat.TransmissionRatio)<maxTXratio ...
    & ~ismember(LabelMat.NeuronID, badNeurons) & LabelMat.TransmissionRatio_norm>0.3);
Xt=LabelMat.SpikeFrame(show_ind); Y=LabelMat.TransmissionRatio_norm(show_ind);
scatter_density(Xt, Y, 50, parula(256));
cb=colorbar;
cb.Label.String='Relative density';
br=binning_data({[Xt Y]}, timebin);
errorbar(br.centers, br.mean, br.std./sqrt(cellfun(@numel,br.values)), 'Color',[1 0 0],'LineWidth',2.5);
xlabel('Time after blue onset (ms)'); ylabel('Normalized transmission ratio'); box off;
xlim([timebin(1) timebin(end)]);

% --- statistical test between time bins (pairwise ranksum on pooled spikes) ---
pBins = get_pValue(br.values, 0);                 % prints + returns full pairwise matrix
pAdj = nan(size(pBins));                           % draw only adjacent-bin comparisons
for b=1:numel(br.centers)-1
    pAdj(b,b+1)=pBins(b,b+1); pAdj(b+1,b)=pBins(b,b+1);
end
drawPValueLines(pBins, 0.1, 'XCoord', br.centers);
ylim([0.3 2]);
set_fontsize(14); set_font('Arial'); set_figsize(110,100);

%% Per-neuron transmission ratio vs time -- soma stim, normalized & un-normalized (figure 27)
% Transmission ratio as a function of time after blue onset, traced separately
% for each neuron (soma wvf==1, pulse spikes; un-binned). Top = normalized,
% bottom = un-normalized; each neuron is one coloured trace (spikes sorted by
% time). Noisy neurons (robust isoutlier on the normalized spread OR mean) are
% drawn thicker and returned in badNeurons for filtering.
wvf=1; maxTXratio=7;
show_ind=find(LabelMat.StimOrder~=6 & LabelMat.StimRegion==wvf & abs(LabelMat.TransmissionRatio)<maxTXratio);
Xt=LabelMat.SpikeFrame(show_ind);
Ynorm=LabelMat.TransmissionRatio_norm(show_ind); Yraw=LabelMat.TransmissionRatio(show_ind);
Nid=LabelMat.NeuronID(show_ind); Mid=LabelMat.MouseID(show_ind);
neurons=unique(Nid)'; nN=numel(neurons);
muNorm=nan(nN,1); sdNorm=nan(nN,1); muRaw=nan(nN,1); sdRaw=nan(nN,1); nSpk=zeros(nN,1); mouseN=zeros(nN,1);
for k=1:nN
    sel=Nid==neurons(k);
    muNorm(k)=mean(Ynorm(sel),'omitnan'); sdNorm(k)=std(Ynorm(sel),'omitnan');
    muRaw(k) =mean(Yraw(sel),'omitnan');  sdRaw(k) =std(Yraw(sel),'omitnan');
    nSpk(k)=nnz(sel); mouseN(k)=Mid(find(sel,1));
end
isBad = isoutlier(sdNorm) | isoutlier(muNorm);    % robust: median +/- 3*scaled MAD
badNeurons=neurons(isBad);

figure(27); clf;
cmap=distinguishable_colors(nN);
Ysets={Ynorm,Yraw}; ylabels={'Normalized transmission ratio','Transmission ratio'};
for pnl=1:2
    for k=1:nN
        nexttile([1 1]);
        sel=Nid==neurons(k);
        [xs,si]=sort(Xt(sel)); ys=Ysets{pnl}(sel); ys=ys(si);
        h=plot(xs, ys, '.-', 'Color',cmap(k,:), 'MarkerSize',7, 'LineWidth',0.5+isBad(k));
        title(neurons(k));
    end
end
set_fontsize(11);
%% Show TX ratio over time (figure 5)
timebin=[0 30 35 85 110 200:100:500];
Rheobin=[1 1.3 1.6 2 3 4 5];
maxTXratio=7; stimstr={'Soma','Dendrite'}; g=1;
figure(23); clf; tiledlayout(1,2);
for wvf=[1]
    nexttile([1 1]);
    show_ind=find(LabelMat.StimOrder~=6 & LabelMat.StimRegion==wvf & abs(LabelMat.TransmissionRatio)<maxTXratio); %Pulse & Soma
    binRes=binning_data({[LabelMat.SpikeFrame(show_ind) LabelMat.TransmissionRatio_norm(show_ind)]},timebin);
    M=binRes.mean; S=binRes.std; t_center=binRes.centers; SubSet=binRes.values; N=cellfun(@numel,binRes.values);
    p=get_pValue(SubSet,0);
    scatter_density(LabelMat.SpikeFrame(show_ind),LabelMat.TransmissionRatio_norm(show_ind),40);
    %scatter_heatmap(LabelMat.SpikeFrame(show_ind),LabelMat.TransmissionRatio(show_ind),100,200);
    hold all;
    errorbar(t_center, M, S ./ sqrt(N), 'Color', [1 0 0], 'LineWidth', 1.5);
    xlim([0 500]); ylim([0 2.5]);
    %set(gca,'xscale','log')
    %drawPValueLines(p,0,'TextYOffset',0.1,'XCoord',t_center);
    xlabel('Time after blue onset (ms)'); ylabel('Normalized transmission ratio'); box off;
    title([stimstr{g} ' Stim. (Pulse)']); caxis([0 0.0012]);
    colormap([gen_colormap([0 0 0; parula(10)])])

    nexttile([1 1]);
    show_ind=find(LabelMat.StimOrder==6 & LabelMat.StimRegion==wvf & abs(LabelMat.TransmissionRatio)<maxTXratio); %Pulse & Soma
    binRes=binning_data({[LabelMat.Rheobase(show_ind) LabelMat.TransmissionRatio_norm(show_ind)]},Rheobin);
    M=binRes.mean; S=binRes.std; Rheo_center=binRes.centers; SubSet=binRes.values; N=cellfun(@numel,binRes.values);
    pR=get_pValue(SubSet,0);
    scatter_density(LabelMat.Rheobase(show_ind),LabelMat.TransmissionRatio_norm(show_ind),40);
    %scatter_heatmap(LabelMat.Rheobase(show_ind),LabelMat.TransmissionRatio(show_ind),500,200);
    hold all;
    errorbar(Rheo_center, M, S ./ sqrt(N), 'Color', [1 0 0], 'LineWidth', 1.5);
    xlim([0.5 5]); ylim([0 4]); box off;
    %drawPValueLines(pR,0,'TextYOffset',0.1,'XCoord',Rheo_center);
    xlabel('Rheobase'); ylabel('Transmission ratio')
    title([stimstr{g} ' Stim. (Ramp)']); caxis([0 0.195]);
    colormap([gen_colormap([parula(10)])])
    g=g+1;
end

% Show TX ratio, spike order
TXratio_sporder=[]; g2=1;
for stimregion2show=[1 3];
    for sporder=1:20
        g=1;
        [unqLab, refind]=unique([LabelMat.SessionID LabelMat.StimRegion],'row');
        for n=refind(unqLab(:,2)==stimregion2show,1)' %soma targeted
            sessioninterested=unique(LabelMat.SessionID(n));
            stimregioninterested=unique(LabelMat.StimRegion(n));
            ind2norm=find(LabelMat.SessionID==sessioninterested & LabelMat.StimRegion==stimregioninterested);
            show_ind=find(LabelMat.StimOrder~=6 & LabelMat.StimRegion==stimregioninterested & abs(LabelMat.TransmissionRatio_norm)<maxTXratio & LabelMat.SpikeOrder==sporder & LabelMat.SessionID==sessioninterested); %Pulse & Soma
            % TXratio_sporder(sporder,g,1)=mean(LabelMat.TransmissionRatio(show_ind)./median(LabelMat.TransmissionRatio(ind2norm),'omitnan'),'omitnan');
            % TXratio_sporder(sporder,g,2)=std(LabelMat.TransmissionRatio(show_ind)./median(LabelMat.TransmissionRatio(ind2norm),'omitnan'),'omitnan');
            TXratio_sporder(sporder,g,1)=mean(LabelMat.TransmissionRatio_norm(show_ind),'omitnan');
            TXratio_sporder(sporder,g,2)=std(LabelMat.TransmissionRatio_norm(show_ind),'omitnan');
            g=g+1;
        end
    end
    g2=g2+1;
end
figure(24); clf; cmap_order=gray(20); show_sporder=[1:10];
nexttile([1 1]);
TXratio_sporder(TXratio_sporder>5 | TXratio_sporder<0)=NaN;
TXratio_sporder_mean=mean(TXratio_sporder(show_sporder,:,1)',1,'omitnan');
TXratio_sporder_sem=std(TXratio_sporder(show_sporder,:,1)',0,1,'omitnan')./sqrt(sum(~isnan(TXratio_sporder(show_sporder,:,1)),2)');
%plot(show_sporder,TXratio_sporder(show_sporder,:,1)','color',[0.7 0.7 0.7]); hold all;
errorbar_shade(show_sporder,TXratio_sporder_mean,TXratio_sporder_sem,[1 0 0]);
%Violin_wPoints(TXratio_sporder(show_sporder,:,1)',cmap_order(show_sporder,:));
p_values = get_pValue(TXratio_sporder(show_sporder,:,1)', 1);
p_values_show=p_values;
p_values_show(setdiff([1:size(p_values,1)],[3 4]),:)=1;
p_values_show=min(cat(3,p_values_show,p_values_show'),[],3);
drawPValueLines(p_values_show,0,'StepHeight',0.03,'TextYOffset',0.02);
xlim([show_sporder(1)-0.5 show_sporder(end)+0.5]); ylim([0.85 1.55]);
set(gca,'xtick',show_sporder,'XTickLabel',counting_string(show_sporder))
ylabel(sprintf('Normalized\n transmission ratio'));  xlabel('Evoked spike order');
box off; set_font('Arial'); set_fontsize(16); set_figsize(120,120);
%% Show bAP apical AUC over time (figure 5)
timebin=[0 30 35 85 110 200:100:500];
Rheobin=[1 1.3 1.6 2 3 4 5];
maxTXratio=7; stimstr={'Soma','Dendrite'}; g=1;
figure(23); clf; tiledlayout(1,2);
for wvf=[1]
    nexttile([1 1]);
    show_ind=find(LabelMat.StimOrder~=6 & LabelMat.StimRegion==wvf & abs(LabelMat.TransmissionRatio)<maxTXratio); %Pulse & Soma
    binRes=binning_data({[LabelMat.SpikeFrame(show_ind) LabelMat.AUC_apical_norm(show_ind)]},timebin);
    M=binRes.mean; S=binRes.std; t_center=binRes.centers; SubSet=binRes.values; N=cellfun(@numel,binRes.values);
    p=get_pValue(SubSet,0);
    scatter_density(LabelMat.SpikeFrame(show_ind),LabelMat.AUC_apical_norm(show_ind),40);
    %scatter_heatmap(LabelMat.SpikeFrame(show_ind),LabelMat.TransmissionRatio(show_ind),100,200);
    hold all;
    errorbar(t_center, M, S ./ sqrt(N), 'Color', [1 0 0], 'LineWidth', 1.5);
    xlim([0 500]); ylim([0 2.5]);
    %set(gca,'xscale','log')
    %drawPValueLines(p,0,'TextYOffset',0.1,'XCoord',t_center);
    xlabel('Time after blue onset (ms)'); ylabel('Normalized bAP apical AUC'); box off;
    title([stimstr{g} ' Stim. (Pulse)']); caxis([0 0.0012]);
    colormap([gen_colormap([0 0 0; parula(10)])])

    nexttile([1 1]);
    show_ind=find(LabelMat.StimOrder==6 & LabelMat.StimRegion==wvf & abs(LabelMat.TransmissionRatio)<maxTXratio); %Pulse & Soma
    binRes=binning_data({[LabelMat.Rheobase(show_ind) LabelMat.AUC_apical_norm(show_ind)]},Rheobin);
    M=binRes.mean; S=binRes.std; Rheo_center=binRes.centers; SubSet=binRes.values; N=cellfun(@numel,binRes.values);
    pR=get_pValue(SubSet,0);
    scatter_density(LabelMat.Rheobase(show_ind),LabelMat.AUC_apical_norm(show_ind),40);
    %scatter_heatmap(LabelMat.Rheobase(show_ind),LabelMat.TransmissionRatio(show_ind),500,200);
    hold all;
    errorbar(Rheo_center, M, S ./ sqrt(N), 'Color', [1 0 0], 'LineWidth', 1.5);
    xlim([0.5 5]); ylim([0 4]); box off;
    %drawPValueLines(pR,0,'TextYOffset',0.1,'XCoord',Rheo_center);
    xlabel('Rheobase'); ylabel('Transmission ratio')
    title([stimstr{g} ' Stim. (Ramp)']); caxis([0 0.195]);
    colormap([gen_colormap([parula(10)])])
    g=g+1;
end

% Show TX ratio, spike order
TXratio_sporder=[]; g2=1;
for stimregion2show=[1 3];
    for sporder=1:13
        g=1;
        [unqLab, refind]=unique([LabelMat.SessionID LabelMat.StimRegion],'row');
        for n=refind(unqLab(:,2)==stimregion2show,1)' %soma targeted
            sessioninterested=unique(LabelMat.SessionID(n));
            stimregioninterested=unique(LabelMat.StimRegion(n));
            ind2norm=find(LabelMat.SessionID==sessioninterested & LabelMat.StimRegion==stimregioninterested);
            show_ind=find(LabelMat.StimOrder~=6 & LabelMat.StimRegion==stimregioninterested & abs(LabelMat.TransmissionRatio)<maxTXratio & LabelMat.SpikeOrder==sporder & LabelMat.SessionID==sessioninterested); %Pulse & Soma
            % TXratio_sporder(sporder,g,1)=mean(LabelMat.TransmissionRatio(show_ind)./median(LabelMat.TransmissionRatio(ind2norm),'omitnan'),'omitnan');
            % TXratio_sporder(sporder,g,2)=std(LabelMat.TransmissionRatio(show_ind)./median(LabelMat.TransmissionRatio(ind2norm),'omitnan'),'omitnan');
            TXratio_sporder(sporder,g,1)=mean(LabelMat.AUC_apical_norm(show_ind),'omitnan');
            TXratio_sporder(sporder,g,2)=std(LabelMat.AUC_apical_norm(show_ind),'omitnan');
            g=g+1;
        end
    end
    g2=g2+1;
end
figure(24); clf; cmap_order=gray(10); show_sporder=[1:7];
nexttile([1 1]);
TXratio_sporder(TXratio_sporder>5 | TXratio_sporder<0)=NaN;
TXratio_sporder_mean=mean(TXratio_sporder(show_sporder,:,1)',1,'omitnan');
TXratio_sporder_sem=std(TXratio_sporder(show_sporder,:,1)',0,1,'omitnan')./sqrt(sum(~isnan(TXratio_sporder(show_sporder,:,1)),2)');
%plot([1:4],TXratio_sporder(show_sporder,:,1)','color',[0.7 0.7 0.7]); hold all;
%errorbar([1:4],TXratio_sporder_mean,TXratio_sporder_sem,'r');
Violin_wPoints(TXratio_sporder(show_sporder,:,1)',cmap_order(show_sporder,:));
p_values = get_pValue(TXratio_sporder(show_sporder,:,1)', 1);
drawPValueLines(p_values,0.1);
xlim([show_sporder(1)-0.5 show_sporder(end)+0.5]); ylim([0 2.5]);
set(gca,'xtick',[1:5],'XTickLabel',counting_string(1:5))
ylabel('Normalized bAP apical AUC');  xlabel('Evoked spike order');

%% Representative Pulse stim (Figure 5a)
% show example bAP amplitude/AUC image
ShowSession=[105 1]; ExampleRep=1;
ExampleSession=ShowSession(1,1); wvf=ShowSession(1,2);
ExampleNeuron=unique(LabelMat.NeuronID(find(LabelMat.SessionID==ExampleSession & LabelMat.StimRegion==wvf)));
[nROI nTime]=size(OpResult{1,1,ExampleNeuron}.normTraces);
%ROIvec=setdiff([1:nROI],[]); 
ROIvec=OpResult{1,1,ExampleNeuron}.maintrunkROI;
ftprint2show=OpResult{1,1,ExampleNeuron}.ftprnt(:,:,ROIvec); ftprintall=OpResult{1,1,ExampleNeuron}.ftprnt;
nROI=size(OpResult{1,1,ExampleNeuron}.ftprnt,3); ftprint2show_trace=get_coord(ftprint2show);
somROI=1; dendROI=[7 8]; caxshow=[-0.02 0.17];
SWCpoints=OpResult{1,1,ExampleNeuron}.SWC; SWCpoints(1,3)=70; cmap_ExTr=gen_colormap(Plasma,10);
time2show=[2980:3550];
% 
% figure(25); clf;
% nexttile([1 1])
% boundary_s=bwboundaries(max(ftprint2show(:,:,somROI),[],3)); boundary_ap=bwboundaries(max(ftprint2show(:,:,dendROI),[],3));
% boundary_blue=bwboundaries(OpResult{1,1,ExampleNeuron}.BlueDMDimg);
% showScaleScatter([1:nROI], SWCpoints, ftprintall,repmat([0 0 0],256,1)); hold all
% plot(boundary_s{1}(:,1),boundary_s{1}(:,2),'color',cmap_ExTr(1,:),'linewidth',2);
% plot(boundary_ap{1}(:,1),boundary_ap{1}(:,2),'color',cmap_ExTr(end,:),'linewidth',2);
% plot(boundary_blue{1}(:,1),boundary_blue{1}(:,2),'color',[0 0.5 1],'linewidth',2);
% plot(ftprint2show_trace(:,2),ftprint2show_trace(:,1),'r','linewidth',2);
% drawScaleBar(100/OpResult{1,1,ExampleNeuron}.pixelsize,'vertical')
% view([210 90])

nROI=size(OpResult{1,1,ExampleNeuron}.normTraces,1);
ftprint=OpResult{1,1,ExampleNeuron}.ftprnt(:,:,OpResult{1,1,ExampleNeuron}.dist_order)>0;
dax=OpResult{1,1,ExampleNeuron}.dendaxis;

normTr=OpResult{wvf,1,ExampleNeuron}.normTraces./OpResult{wvf,1,ExampleNeuron}.ScaleFactor;
normTr_res=normTr(ROIvec,time2show);
normTr_res=pcafilterTrace(normTr_res,setdiff([1:9],[3]));
normTr_res=normTr_res-median(normTr_res(:,1:35),2);
Tr2show=[mean(normTr_res(somROI,:),1); mean(normTr_res(dendROI,:),1)];
normTr_interp=interp1(dax(ROIvec),normTr_res,linspace(dax(ROIvec(1)),dax(ROIvec(end)),16));
dax2show=linspace(dax(ROIvec(1)),dax(ROIvec(end)),10);

figure(27); clf; tiledlayout(7,1,'TileSpacing','tight'); % show example bAP amplitude/AUC image
ax1=nexttile([3 1]);
normTr_interp=normTr_interp-median(normTr_interp(:,1:35),2);
t_ax=[1:size(normTr_interp,2)]*0.001;
imagesc(t_ax, dax2show, normTr_interp,caxshow); colormap(ax1,'turbo'); cb=colorbar; cb.Label.String = '∆F/F';
set(gca,'xtick',[]); ylabel('Distance (µm)');

ax2=nexttile([4 1]);
l=plot(t_ax,Tr2show([1 2],:)','linewidth',1.5);
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(cmap_ExTr([1 10],:),2));
ylabel('∆F/F'); xlabel('Time (ms)');
axis tight off;
ylim([-0.05 0.2]);
drawScaleBar(0.1, 'horizontal','color','k','position',[0.570 -0.02]);
drawScaleBar(0.1, 'vertical','color','k','position',[0.570 0.08]);
linkaxes([ax1 ax2],'x'); box off;
set_font('Arial'); set_fontsize(16);
set_figsize(200,100);
%%
ss=1; t_start=[3000 3290]-time2show(1)+1; t_showlength=90;
t_ax=[1:t_showlength+1];
figure(30); clf;
tiledlayout(5,length(t_start),'padding','compact','TileSpacing','tight'); ax2=[];
for t=1:length(t_start)
    ax1=nexttile(t,[2 1]);
    kymos2show=normTr_interp(:,t_start(t):t_start(t)+t_showlength);
    %kymos2show=pcafilterTrace(kymos2show,[1:5]);
    %kymos2show=kymos2show(ShowROI,:);
    imagesc(t_ax, dax2show, kymos2show,caxshow); colormap(ax1,'turbo');
    Trcat=Tr2show(:,t_start(t):t_start(t)+t_showlength);
    set(gca,'xtick',[]);
if t==1
    ylabel('Distance (µm)');
else
    set(gca,'ytick',[]);
    cb=colorbar; cb.Label.String = '∆F/F';
end
    ax2=[ax2 nexttile(t+length(t_start)*2,[3 1])];
    l=plot(t_ax,Trcat([1 2],:)','linewidth',2);
    arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(cmap_ExTr([1 10],:),2));
    ylabel('∆F/F'); xlabel('Time (ms)');
    axis tight off;
    ylim([-0.05 0.2]);
end

drawScaleBar(10, 'horizontal','color','k','position',[93 0]);
drawScaleBar(0.1, 'vertical','color','k','position',[93 0]);

set_font('Arial'); set_fontsize(16.5);
set_figsize(330,120);

%% Pad the spike-group AUC maps into ONE axis with an x-offset (instead of
% separate figures 32-34). Precompute each SWC point's ROI membership once
% (mirrors showScaleScatter), then scatter a shifted copy per spike group.


Splist=find(LabelMat.SessionID==ExampleSession & LabelMat.StimRegion==1 & LabelMat.NeuronID==ExampleNeuron);
showSP={[1],[2 3 4],[13 14 15]}; excludeROI=[15 35]; cax_AUC=[0.02 0.4]/6;
AUC2show=cell2mat(cellfun(@(x) x(:,2),AUCbAPall(Splist),'UniformOutput',false));
ftprintall2show=ftprintall;
ftprintall2show(:,:,excludeROI)=[];
AUC2show=AUC2show;%./AUC2show(1,:);
AUC2show(excludeROI,:)=[];
%tiledlayout(1,9,'Padding','tight'); 

sz32=[size(ftprintall2show,1) size(ftprintall2show,2)];
swc32=SWCpoints; if size(swc32,2)<4, swc32(:,4)=1; end; swc32(swc32(:,4)==0,4)=1;
ptX32=swc32(:,2); ptY32=swc32(:,1); ptSz32=swc32(:,4); ptROI32=zeros(size(swc32,1),1);
ptSz32(1)=40; ptSz32=ptSz32+4;
valid32=ptY32>0.5 & ptY32<sz32(2)-0.5 & ptX32>0.5 & ptX32<sz32(1)-0.5;
lin32=nan(size(swc32,1),1); lin32(valid32)=sub2ind(sz32, round(ptX32(valid32)), round(ptY32(valid32)));
for r=1:size(ftprintall2show,3)
    m=ftprintall2show(:,:,r)>0; hit=false(size(swc32,1),1); hit(valid32)=m(lin32(valid32)); ptROI32(hit)=r;
end
% Rotate the SWC display coordinates in-plane (membership is unchanged -- it is
% computed on the original coords above). Tune rotAngle32 to orient the neuron.
rotAngle32=-20;   % degrees
ct32=mean([ptX32 ptY32],1);
Rm32=[cosd(rotAngle32) -sind(rotAngle32); sind(rotAngle32) cosd(rotAngle32)];
XYr32=(Rm32*([ptX32 ptY32]-ct32)')' + ct32;
ptXr32=XYr32(:,1); ptYr32=XYr32(:,2);
xStep32=1.02*(max(ptXr32)-min(ptXr32));   % spacing between spike groups (rotated frame)

figure(32); clf; hold on;
for s=1:length(showSP)
    showvec=mean(AUC2show(:,showSP{s}),2,'omitnan')/6;
    cols=vec2cmap(showvec,'turbo',cax_AUC); xoff=(s-1)*xStep32;
    g0=ptROI32==0;
    scatter(ptXr32(g0)-xoff, ptYr32(g0), ptSz32(g0), [0.5 0.5 0.5], 'filled');
    for r=1:max(ptROI32)
        sel=ptROI32==r;
        if any(sel), scatter(ptXr32(sel)-xoff, ptYr32(sel), ptSz32(sel), cols(r,:), 'filled'); end
    end
end
colormap(turbo); caxis(cax_AUC); axis equal tight off;
set(gca,'ydir','reverse','xdir','reverse');
cb=colorbar; cb.Label.String='mean ∆F/F';
set_font('Arial'); set_fontsize(16);
set_figsize(150,90);
% 
% figure(32); clf;
% set(gcf,'Color','w')
% imshow2(Fcat,[]);
% colormap(turbo); 
% cb=colorbar; 
% cb.Ticks=[0 1]; caxis([0 1]);
% cb.TickLabels=num2cell(cax_AUC)';
% cb.Label.String = 'AUC (A.U.)';


%% Soma stim VS Dendrite stim (Figure S)
SomRP_path=[15 54 72 86 92 4 76]; %22
DenRP_path=[13 55 74 88 91 3 79]; %23
validresultind=find(cell2mat(cellfun(@(x) ~isempty(x),OpResult,'UniformOutput',false)));
finterest=cell2mat(cellfun(@(x) x.fileInd,OpResult(validresultind),'UniformOutput',false));
finterest_SomRP=validresultind(ismember(finterest,SomRP_path));
finterest_DenRP=validresultind(ismember(finterest,DenRP_path));
[finterest_SomRPw, finterest_SomRPr, finterest_SomRPn]=ind2sub(size(OpResult),finterest_SomRP);
[finterest_DenRPw, finterest_DenRPr, finterest_DenRPn]=ind2sub(size(OpResult),finterest_DenRP);

%path_cat=[SomRP_path; DenRP_path];
rheo_bin=[0:1.3:10];
CSfraction=NaN(size(finterest_SomRP,1),2);
for n=1:size(finterest_SomRP,1) %neuron
    resultdat=OpResult{finterest_SomRPw(n),finterest_SomRPr(n),finterest_SomRPn(n)};
    spikeInfoInd=(LabelMat.NeuronID==finterest_SomRPn(n) & LabelMat.RepeatID==finterest_SomRPr(n) & LabelMat.StimRegion==finterest_SomRPw(n));
    CSfraction(n,1)=sum(LabelMat.InCS(spikeInfoInd))/sum(spikeInfoInd);

    spikeInfoInd=(LabelMat.NeuronID==finterest_DenRPn(n) & LabelMat.RepeatID==finterest_DenRPr(n) & LabelMat.StimRegion==finterest_DenRPw(n));
    CSfraction(n,2)=sum(LabelMat.InCS(spikeInfoInd))/sum(spikeInfoInd);
end

showneuron=4;
result2show=OpResult{finterest_SomRPw(showneuron),finterest_SomRPr(showneuron),finterest_SomRPn(showneuron)};
ROIvec=result2show.maintrunkROI;
normTR=result2show.normTraces./result2show.ScaleFactor;
[nROI nTime]=size(result2show.normTraces);
%[S]=get_STA(normTR,result2show.spike(1,:),0,0);
%normTR=normTR./S(1);
somROI=1; dendROI=10;
normTR2show_S=mean(normTR(somROI,:),1,'omitnan'); normTR2show_D=mean(normTR(dendROI,:),1,'omitnan');
%normTR2show_S=get_bandstop(normTR2show_S,1000,[70 140]);
normTR2show_S_CS=normTR2show_S; normTR2show_S_CS(result2show.CStrace==0)=NaN;
normTR2show_D_CS=normTR2show_D; normTR2show_D_CS(result2show.CStrace==0)=NaN;
load(fullfile(fpath{result2show.fileInd},"output_data.mat"))
[~,~,~,DMDimg]=get_blueDMDPatt(Device_Data);
ROIbigger = double(Device_Data{1, 3}.ROI([1 3 2 4])); imgsize=ROIbigger(:,3:4);
ROIbigger(:,1:2) = ROIbigger(:,1:2)-imgsize/2; ROIbigger(:,3:4) = imgsize*2;
DMDimg=imcrop(DMDimg,ROIbigger);
blueboundary=cell2mat(bwboundaries(imgaussfilt(DMDimg,4)>0.5));
blueboundary=blueboundary-imgsize([2 1])/2;

figure(33); clf; scaleoffset=6; tiledlayout(2,7);
ax1=nexttile(1,[1 2]);
SWC2show=result2show.SWC; SWC2show(:,4)=3; SWC2show(1,4)=20;
Ftprnts=result2show.ftprnt;
showScaleScatter(ones(nROI),SWC2show,Ftprnts,repmat([0 0 0],256,1),[0 1]); hold all;
plot(blueboundary(:,1),blueboundary(:,2),'color',[0 0.6 1],'linewidth',1.5);
drawScaleBar(100/resultdat.pixelsize,'vertical')

nexttile(3,[2 3]);
plot(normTR2show_S,'color',[0.4 0.4 0.4]); hold all
plot(normTR2show_S_CS,'color',[1 0 0]); hold all
axis off;
drawScaleBar(1000,'horizontal')

result2show=OpResult{finterest_DenRPw(showneuron),finterest_DenRPr(showneuron),finterest_DenRPn(showneuron)};
ROIvec=result2show.maintrunkROI;
normTR=result2show.normTraces./result2show.ScaleFactor;
load(fullfile(fpath{result2show.fileInd},"output_data.mat"))
[~,~,~,DMDimg]=get_blueDMDPatt(Device_Data);
ROIbigger = double(Device_Data{1, 3}.ROI([1 3 2 4])); imgsize=ROIbigger(:,3:4);
ROIbigger(:,1:2) = ROIbigger(:,1:2)-imgsize/2; ROIbigger(:,3:4) = imgsize*2;
DMDimg=imcrop(DMDimg,ROIbigger);
blueboundary=cell2mat(bwboundaries(imgaussfilt(DMDimg,4)>0.5));
blueboundary=blueboundary-imgsize([2 1])/2;

%[S]=get_STA(normTR,result2show.spike(1,:),0,0);
%normTR=normTR./S(1);
somROI=1; dendROI=10;
normTR2show_S=mean(normTR(somROI,:),1,'omitnan'); normTR2show_D=mean(normTR(dendROI,:),1,'omitnan');
normTR2show_S=get_bandstop(normTR2show_S,1000,[20 50]);
normTR2show_S_CS=normTR2show_S; normTR2show_S_CS(result2show.CStrace==0)=NaN;
normTR2show_D_CS=normTR2show_D; normTR2show_D_CS(result2show.CStrace==0)=NaN;

ax2=nexttile(8,[1 2]);
showScaleScatter(ones(nROI),SWC2show,Ftprnts,repmat([0 0 0],256,1),[0 1]); hold all;
plot(blueboundary(:,1),blueboundary(:,2),'color',[0 0.6 1],'linewidth',1.5);

nexttile(3,[2 3]);
plot(normTR2show_S+scaleoffset,'color',[0.4 0.4 0.4]); hold all
plot(normTR2show_S_CS+scaleoffset,'color',[1 0 0]); hold all
axis off;
nexttile(6,[2 2])
p=Boxplot_wPoints2(CSfraction,[0 0 0]); box off;
drawPValueLines(p,0,'Format','text');
set(gca,'xtick',[1 2],'xticklabel',{'Soma stim.','Dendrite stim.'});
ylabel('Fraction of complex spikes');

p = linkprop([ax1 ax2], {'CameraPosition', 'CameraTarget', ...
                       'CameraUpVector', 'CameraViewAngle'});
set_fontsize(14);
%% [ARCHIVE] Earlier dendrite/soma-stim complex-spike exploration (commented)


    % Ramp_period=find(bwblue==max(bwblue));
    % Rheobase_frame=find(max(Result.spike(:,Ramp_period),[],1))+Ramp_period(1)-1;
    % Rheobase_blue=Result.Blue(round(mean(Rheobase_frame(1:3))));
    % %Rheobase_blue=1;
    % 
    % firstpulst_bw=max(bwblue)-5;
    % for b=1:5
    % pulse_period=find(bwblue==firstpulst_bw+b-1);    
    % pulseBlue_rb(b,f,r)=median(Result.Blue(pulse_period))/Rheobase_blue;
    % pulseFR(b,:,f,r)=sum(Result.SpClass(1:2,pulse_period),2)';
    % end
    % 
    % RampFR(:,:,f,r)=NaN(length(rheo_bin),2);
    % Ramp_period=find(bwblue==max(bwblue));    
    % binnedRheo_ramp=ceil(Result.Blue(Ramp_period)/Rheobase_blue/(rheo_bin(2)-rheo_bin(1)));
    % for b=binnedRheo_ramp
    % bin_frm=find(binnedRheo_ramp==b);
    % RampFR(b,:,f,r)=sum(Result.SpClass(1:2,bin_frm+Ramp_period(1)-1),2)';
    % end
    % end
% 
% %% Update of dendStim/SomaStim (2025.03.20)
% 
% 
% SomRP_path=[15 22 54 72 86 92 4];
% DenRP_path=[13 23 55 74 88 91 3];
% path_cat=[SomRP_path; DenRP_path];
% rheo_bin=[0:1.3:10];
% pulseFR=[]; RampFR=[]; pulseBlue_rb=[];
% SpikeMatInfo=[];
% for f=1:size(path_cat,2)
%     for r=1:2 %soma and dendrite
%     SpikeMatInfo{f,r}=[];   
%     i=path_cat(r,f);
%     load(fullfile(fpath{i},'OP_Result'))
%     bwblue=bwlabel(Result.Blue);
%     Ramp_period=(bwblue==max(bwblue));
%     Rheobase_frame=find(max(Result.spike(:,Ramp_period),[],1))+find(Ramp_period,1)-1;
%     Rheobase_blue=Result.Blue(round(mean(Rheobase_frame(1:3))));
%     %Rheobase_blue=1;
% 
%     SP_trace=[Result.spike(1,:).*(1-Result.CStrace);Result.spike(1,:).*Result.CStrace];
%     SP_trace=SP_trace.*(1-Ramp_period).*double(Result.Blue>0);
%     BlueOnset=find((double(Result.Blue(2:end)>0)-double(Result.Blue(1:end-1)>0))>0);
% 
%     for stype=1:2 %SS and CS
%     Sptime=find(SP_trace(stype,:));
%     SptimeMod=mod(Sptime-BlueOnset(1),BlueOnset(2)-BlueOnset(1));
%     SpikeMatInfo{f,r}=[SpikeMatInfo{f,r}; [Result.Blue(Sptime)'/Rheobase_blue SptimeMod' ones(length(Sptime),1)*stype]];
%     end
%     end
% end
% 
% 
% 
% 
% figure(15); clf; cmap=[0.4 0.4 0.4; 1 0 0]; tiledlayout(1,2);
% RheoEdge=[0 1 15];
% TimeEdge=[0:100:500];
% % for r=1:2
% %     nexttile([1 1])
% % ShowDat=cell2mat(SpikeMatInfo(:,r));
% % for stype=1:2
% %     scatter(ShowDat(ShowDat(:,3)==stype,2),ShowDat(ShowDat(:,3)==stype,1),10,cmap(stype,:),'o','filled'); hold all
% % end
% % end
% 
% BinRheoRatio=cellfun(@(x) histcounts(x(x(:,3)==2,1),RheoEdge)./histcounts(x(:,1),RheoEdge),SpikeMatInfo,'UniformOutput',false);
% BinRheoRatio_weak=cellfun(@(x) x(1),BinRheoRatio,'UniformOutput',false);
% BinRheoRatio_strong=cellfun(@(x) x(2),BinRheoRatio,'UniformOutput',false);
% nexttile([1 1]);
% Boxplot_wPoints(cell2mat(BinRheoRatio_weak(:,1)),cell2mat(BinRheoRatio_weak(:,2)),[1 0 0],'pairwise')
% set(gca,'XTick',[1 2],'XTickLabel',{'Soma Stim.','Dend. Stim.'})
% ylabel('Fraction of complex spike')
% nexttile([1 1]);
% p=Boxplot_wPoints(cell2mat(BinRheoRatio_strong(:,1)),cell2mat(BinRheoRatio_strong(:,2)),[1 0 0],'pairwise')
% set(gca,'XTick',[1 2],'XTickLabel',{'Soma Stim.','Dend. Stim.'})
% ylabel('Fraction of complex spike')
% 
% 
% 
%     firstpulst_bw=max(bwblue)-5;
%     for b=1:5
%     pulse_period=find(bwblue==firstpulst_bw+b-1);    
%     pulseBlue_rb(b,f,r)=median(Result.Blue(pulse_period))/Rheobase_blue;
%     pulseFR(b,:,f,r)=sum(Result.SpClass(1:2,pulse_period),2)';
%     end
% 
%     RampFR(:,:,f,r)=NaN(length(rheo_bin),2);
%     Ramp_period=find(bwblue==max(bwblue));    
%     binnedRheo_ramp=ceil(Result.Blue(Ramp_period)/Rheobase_blue/(rheo_bin(2)-rheo_bin(1)));
%     for b=binnedRheo_ramp
%     bin_frm=find(binnedRheo_ramp==b);
%     RampFR(b,:,f,r)=sum(Result.SpClass(1:2,bin_frm+Ramp_period(1)-1),2)';
%     end
%     end
% end
% 
% 
% figure(27); clf; cmap=distinguishable_colors(6); bin_width=1.3;
% tiledlayout(6,2);
% load(fullfile(fpath{SomRP_path(5)},'OP_Result'));
% tr_somStim=rescale(Result.normTraces(1,:));
% tr_cssomStim=Result.CStrace;
% load(fullfile(fpath{DenRP_path(6)},'OP_Result'));
% tr_ddStim=rescale(Result.normTraces(1,:));
% tr_csddStim=Result.CStrace;
% ax1=nexttile([3 2]);
% plot(tr_somStim,'color','k'); hold all
% show_cs=tr_somStim.*tr_cssomStim; show_cs(tr_cssomStim==0)=NaN;
% plot(show_cs,'color',[0.7 0 0.2]); hold all
% plot(tr_ddStim+1,'color','k'); hold all
% show_cs=tr_ddStim.*tr_csddStim; show_cs(tr_csddStim==0)=NaN;
% axis off
% text(-500,0.7,'Soma Stimulation','FontSize',13)
% text(-500,1.7,'Distal dendrite Stimulation','FontSize',13)
% text(9500,2,'Complex spikes','FontSize',13,'color',[0.7 0 0.2])
% plot(show_cs+1,'color',[0.7 0 0.2]); hold all
% ax2=nexttile([1 2]);
% plot(Result.Blue)
% linkaxes([ax1 ax2],'x')
% 
% nexttile([2 1])
% CS_fracPulseMatSom=squeeze(pulseFR(:,2,:,1)./sum(pulseFR(:,:,:,1),2));
% CS_fracPulseMatDD=squeeze(pulseFR(:,2,:,2)./sum(pulseFR(:,:,:,2),2));
% markers = {'o', 's', 'd', '^', 'v', '>', '<', 'p', 'h', '+', '*'};
% for f=1:size(pulseFR,3)
% plot(squeeze(pulseBlue_rb(:,f,1))',CS_fracPulseMatSom(:,f),'marker',markers{f},'color',cmap(1,:),'linestyle','none'); hold all
% plot(squeeze(pulseBlue_rb(:,f,2))',CS_fracPulseMatDD(:,f),'marker',markers{f},'color',cmap(2,:),'linestyle','none'); hold all
% end
%     clear Ms Ss
% for b=1:max(max(ceil(pulseBlue_rb(:,:,1)/bin_width)))
%     Ms(b)=mean(CS_fracPulseMatSom(b==ceil(pulseBlue_rb(:,:,1)/bin_width)),'omitnan');
%     Ss(b)=std(CS_fracPulseMatSom(b==ceil(pulseBlue_rb(:,:,1)/bin_width)),'omitnan');
% end
% clear Md Sd
% for b=1:max(max(ceil(pulseBlue_rb(:,:,2)/bin_width)))
%     Md(b)=mean(CS_fracPulseMatDD(b==ceil(pulseBlue_rb(:,:,2)/bin_width)),'omitnan');
%     Sd(b)=std(CS_fracPulseMatDD(b==ceil(pulseBlue_rb(:,:,2)/bin_width)),'omitnan');
% end
% errorbar([1:max(max(ceil(pulseBlue_rb(:,:,1)/bin_width)))],Ms,Ss,'color',cmap(1,:),'LineWidth',1.5); hold all
% errorbar([1:max(max(ceil(pulseBlue_rb(:,:,2)/bin_width)))],Md,Sd,'color',cmap(2,:),'LineWidth',1.5);
% 
% xlabel('Optical Rheobase')
% ylabel('Fraction of complex spike')
% title('Pulse Stimulation')
% 
% nexttile([2 1])
% CS_fracRampMatSom=squeeze(RampFR(:,2,:,1)./sum(RampFR(:,:,:,1),2));
% CS_fracRampMatDD=squeeze(RampFR(:,2,:,2)./sum(RampFR(:,:,:,2),2));
% % for f=1:size(RampFR,3)
% % plot(rheo_bin,CS_fracRampMatSom,'color',cmap(1,:)); hold all
% % plot(rheo_bin,CS_fracRampMatDD,'color',cmap(2,:))
% % end
% errorbar(rheo_bin,mean(CS_fracRampMatSom,2,'omitnan'),std(CS_fracRampMatSom,0,2,'omitnan'),'color',cmap(1,:),'LineWidth',1.5); hold all
% errorbar(rheo_bin,mean(CS_fracRampMatDD,2,'omitnan'),std(CS_fracRampMatDD,0,2,'omitnan'),'color',cmap(2,:),'LineWidth',1.5)
% xlabel('Optical Rheobase')
% ylabel('Fraction of complex spike')
% title('Ramp Stimulation')
% 
















%% Representative Pulse ramp stim, soma & apical dend stim (Figure 4)
% show example bAP amplitude/AUC image

ExampleNeuron=6; 
ROIvec=setdiff(OpResult{1,1,ExampleNeuron}.maintrunkROI, []); 
ftprint2show=OpResult{1,1,ExampleNeuron}.ftprnt(:,:,ROIvec); ftprintall=OpResult{1,1,ExampleNeuron}.ftprnt;
nROI=size(OpResult{1,1,ExampleNeuron}.ftprnt,3); ftprint2show_trace=get_coord(ftprint2show);
caxshow=[-2.5 6.5];
%SWCpoints=OpResult{1,1,ExampleNeuron}.SWC; SWCpoints(1,3)=70; 
cmap_ExTr=gen_colormap(Plasma,10);
nROI=size(OpResult{1,1,ExampleNeuron}.normTraces,1);
ftprint=OpResult{1,1,ExampleNeuron}.ftprnt(:,:,OpResult{1,1,ExampleNeuron}.dist_order)>0;
dax=OpResult{1,1,ExampleNeuron}.dendaxis;
% 
% figure(25); clf;
% nexttile([1 1])
% boundary_s=bwboundaries(max(ftprint2show(:,:,somROI),[],3)); boundary_ap=bwboundaries(max(ftprint2show(:,:,dendROI),[],3));
% boundary_blue=bwboundaries(OpResult{1,1,ExampleNeuron}.BlueDMDimg);
% showScaleScatter([1:nROI], SWCpoints, ftprintall,repmat([0 0 0],256,1)); hold all
% plot(boundary_s{1}(:,1),boundary_s{1}(:,2),'color',cmap_ExTr(1,:),'linewidth',2);
% plot(boundary_ap{1}(:,1),boundary_ap{1}(:,2),'color',cmap_ExTr(end,:),'linewidth',2);
% plot(boundary_blue{1}(:,1),boundary_blue{1}(:,2),'color',[0 0.5 1],'linewidth',2);
% plot(ftprint2show_trace(:,2),ftprint2show_trace(:,1),'r','linewidth',2);
% drawScaleBar(100/OpResult{1,1,ExampleNeuron}.pixelsize,'vertical')
% view([210 90])

figure(29); clf;
tiledlayout(5,1); % show example bAP amplitude/AUC image
for wvf=[1 3]
normTr=OpResult{wvf,1,ExampleNeuron}.normTraces./OpResult{wvf,1,ExampleNeuron}.ScaleFactor;
normTr=normTr(ROIvec,:);
%normTr=pcafilterTrace(normTr,[1 2 4 6 7 8 9]);
%sublarge=get_subthreshold(normTr,OpResult{wvf,ExampleRep,ExampleNeuron}.Blue>0,300,1000);
Tr2show=[mean(normTr(1,:),1); mean(normTr(end,:),1)];
Tr2showCS=Tr2show;
Tr2showCS(:,OpResult{wvf,1,ExampleNeuron}.CStrace==0)=NaN;
normTr_interp=interp1(dax(ROIvec),normTr,linspace(dax(ROIvec(1)),dax(ROIvec(end)),10));

ax1=nexttile([2 1]);
plot(Tr2show'+[0 5],'color',[0.4 0.4 0.4]); hold all
plot(Tr2showCS'+[0 5],'color',[1 0 0]); hold all
end


%% [ARCHIVE] Orphan duplicate STA / representative blocks (stale variables)
% Two near-identical leftover blocks that reference variables not defined in
% this scope (show_neuron, STA2show_interp, ShowSession, ShowROI, ...). The
% live versions are "Show Representative STA (figure 4)" and "Representative
% Pulse stim (Figure 5a)". Kept commented to preserve history.
%{
for wvf=[2 4]

    figure(15);
    ax2=[ax2 nexttile([1 1])];
    sp2mean=LabelMat_SP.SpikeOrder==1 & LabelMat_SP.StimRegion==4 & LabelMat_SP.NeuronID==5;
    STA2showall=mean(cat(3,AlignedbAPall_SP{LabelMat_SP.SpikeOrder==1 & LabelMat_SP.StimRegion==wvf & LabelMat_SP.NeuronID==5}),3,'omitnan');
    STA2show=STA2showall(ismember(OpResult{1,1,show_neuron}.dist_order,ROIvec)',:);

    imagesc([-nTau(1):nTau(2)],ROIvec,STA2show_interp,cax);
    set_kymoYtick(dax(ROIvec));

    ax2=[ax2 nexttile([1 1])];
    plot([-nTau(1):nTau(2)],STA2show(1,:),'color',cmap_ExTr(1,:)); hold all;
    plot([-nTau(1):nTau(2)],STA2show(end,:),'color',cmap_ExTr(end,:));
    ylim(cax); xlim([-50 50]);
    axis off;
    drawScaleBar(1,'vertical');
    drawScaleBar(50,'horizontal');
    xlabel('Peri-spike time (ms)')

    figure(16);
    for t=1:length(t2show)
        nexttile([1 1]);
        V2show=STA2showall(:,t2show(t)+nTau(1));
        %V2showcolor=vec2cmap(V2show,turbo,cax);
        showScaleScatter(V2show, SWCpoints, ftprintall(:,:,OpResult{1,1,show_neuron}.dist_order),'turbo',cax);
        axis tight;
        view([180 -90])
    end
end
figure(15);
linkaxes(ax2,'x'); xlim([-50 50]);


figure(25); clf; cmap_ExTr=gen_colormap(Plasma,10);
for ss=1:size(ShowSession,1)
    nexttile([1 1])
    ExampleSession=ShowSession(ss,1); wvf=ShowSession(ss,2);
    ExampleNeuron=unique(LabelMat.NeuronID(find(LabelMat.SessionID==ExampleSession & LabelMat.StimRegion==wvf)));
    ExampleRep=unique(LabelMat.RepeatID(find(LabelMat.SessionID==ExampleSession & LabelMat.StimRegion==wvf)));
    BlueBin=OpResult{wvf,1,ExampleNeuron}.BlueDMDimg;
    SWCpoints=OpResult{1,1,ExampleNeuron}.SWC;
    SWCpoints(1,3)=50;
    ROIvec=[1:nROI]; %ROIvec(ShowROI([35 37 39]))=2; ROIvec(ShowROI([3]))=1;
    ftprint_s=max(ftprint(:,:,ShowROI(SomaROI)),[],3); ftprint_ap=max(ftprint(:,:,ShowROI(ApicalROI)),[],3);
    boundary_s=bwboundaries(ftprint_s); boundary_ap=bwboundaries(ftprint_ap); boundary_blue=bwboundaries(BlueBin);
    showScaleScatter(ROIvec, SWCpoints, ftprint,repmat([0 0 0],256,1)); hold all
    plot(boundary_s{1}(:,1),boundary_s{1}(:,2),'color',cmap_ExTr(1,:),'linewidth',2);
    plot(boundary_ap{1}(:,1),boundary_ap{1}(:,2),'color',cmap_ExTr(end,:),'linewidth',2);
    plot(boundary_blue{1}(:,1),boundary_blue{1}(:,2),'color',[0 0.5 1],'linewidth',2);
    drawScaleBar(100/OpResult{1,1,ExampleNeuron}.pixelsize,'vertical')
    view([180 90])
end


%%

for wvf=[2 4]

    figure(15);
    ax2=[ax2 nexttile([1 1])];
    sp2mean=LabelMat_SP.SpikeOrder==1 & LabelMat_SP.StimRegion==4 & LabelMat_SP.NeuronID==5;
    STA2showall=mean(cat(3,AlignedbAPall_SP{LabelMat_SP.SpikeOrder==1 & LabelMat_SP.StimRegion==wvf & LabelMat_SP.NeuronID==5}),3,'omitnan');
    STA2show=STA2showall(ismember(OpResult{1,1,show_neuron}.dist_order,ROIvec)',:);
    
    imagesc([-nTau(1):nTau(2)],ROIvec,STA2show_interp,cax);
    set_kymoYtick(dax(ROIvec));

    ax2=[ax2 nexttile([1 1])];
    plot([-nTau(1):nTau(2)],STA2show(1,:),'color',cmap_ExTr(1,:)); hold all;
    plot([-nTau(1):nTau(2)],STA2show(end,:),'color',cmap_ExTr(end,:));
    ylim(cax); xlim([-50 50]);
    axis off;
    drawScaleBar(1,'vertical');
    drawScaleBar(50,'horizontal');
    xlabel('Peri-spike time (ms)')

    figure(16);
    for t=1:length(t2show)
        nexttile([1 1]);
        V2show=STA2showall(:,t2show(t)+nTau(1));
        %V2showcolor=vec2cmap(V2show,turbo,cax);
        showScaleScatter(V2show, SWCpoints, ftprintall(:,:,OpResult{1,1,show_neuron}.dist_order),'turbo',cax);
        axis tight;
        view([180 -90])
    end
end
figure(15);
linkaxes(ax2,'x'); xlim([-50 50]);


figure(25); clf; cmap_ExTr=gen_colormap(Plasma,10);
for ss=1:size(ShowSession,1)
    nexttile([1 1])
    ExampleSession=ShowSession(ss,1); wvf=ShowSession(ss,2);
    ExampleNeuron=unique(LabelMat.NeuronID(find(LabelMat.SessionID==ExampleSession & LabelMat.StimRegion==wvf)));
    ExampleRep=unique(LabelMat.RepeatID(find(LabelMat.SessionID==ExampleSession & LabelMat.StimRegion==wvf)));
    BlueBin=OpResult{wvf,1,ExampleNeuron}.BlueDMDimg;
    SWCpoints=OpResult{1,1,ExampleNeuron}.SWC;
    SWCpoints(1,3)=50;
    ROIvec=[1:nROI]; %ROIvec(ShowROI([35 37 39]))=2; ROIvec(ShowROI([3]))=1;
    ftprint_s=max(ftprint(:,:,ShowROI(SomaROI)),[],3); ftprint_ap=max(ftprint(:,:,ShowROI(ApicalROI)),[],3);
    boundary_s=bwboundaries(ftprint_s); boundary_ap=bwboundaries(ftprint_ap); boundary_blue=bwboundaries(BlueBin);
    showScaleScatter(ROIvec, SWCpoints, ftprint,repmat([0 0 0],256,1)); hold all
    plot(boundary_s{1}(:,1),boundary_s{1}(:,2),'color',cmap_ExTr(1,:),'linewidth',2);
    plot(boundary_ap{1}(:,1),boundary_ap{1}(:,2),'color',cmap_ExTr(end,:),'linewidth',2);
    plot(boundary_blue{1}(:,1),boundary_blue{1}(:,2),'color',[0 0.5 1],'linewidth',2);
    drawScaleBar(100/OpResult{1,1,ExampleNeuron}.pixelsize,'vertical')
    view([180 90])
end
%}

%% Transmission rate by spike order (soma stim)
figure(41); clf;
ShowInd=LabelMat.StimRegion==1;
TransmissionRate_order=[];
for ord=1:5
    TransmissionRate_order{ord}=transmissionRate(ShowInd & LabelMat.SpikeOrder==ord);
end
p_values = Violin_wPoints(TransmissionRate_order, autumn(5));
%% Interactive per-neuron trace browser (dendrite stim; press Enter to advance)
figure(222);
Neurons=unique(LabelMat.NeuronID);
for n=Neurons'
    clf;
    tiledlayout('vertical');
    nexttile([1 1])
    imshow2(OpResult{1,1,n}.ref_im,[])
    title([num2str(n)]);
    for wvf=[3]
        maxr=sum(1-cellfun(@isempty,OpResult(wvf,:,n)));
        for r=1:maxr
            ax2=[ax2 nexttile([1 1])];
            ntr=OpResult{wvf,r,n}.normTraces./OpResult{wvf,r,n}.ScaleFactor;
            %ntr=ntr(OpResult{wvf,r,n}.dist_order,:);
            %imagesc(ntr,[-0.01 0.02]);
            plot(ntr(1,:))
            title(['stim region: ' num2str(wvf)]);
        end
    end
    colormap(turbo);
    %    linkaxes(ax2,'x');
    ss=input('yes?\n');
end
%% Ramp (Soma Stim): normalized spike amplitude vs optical rheobase (all ROIs)
nTau=[-20 20]; nTauPeriSp=[-3:7]; nTauPeriSp2=[-3:3]; conductionspeed=170; %um/ms
g=ones(1,4); SP_STA=[]; distFromSoma=[]; LabelMat=[];
dendriteaxis_bin=[-200 -100 -50 50 100:50:500];
dendriteaxis_bin2=[-200 -50 50 150 600];
timebin=[0 20 40 60 100 200:100:500];
perisomadist=[-50 50]; cmap=[0 0 0; 1 0 0];  cmap_light=[0.7 0.7 0.7; 1 0.7 0.7];
%BlueStimN=[2 5 1;3 3 1;3 4 1;4 2 1;4 3 1;4 4 1;5 3 1;5 5 1;6 4 1;6 5 1;6 3 2;6 4 2;6 5 2;7 1 1;7 2 1;7 3 1;];
BlueStimN=[2 6 1;3 6 1;4 6 1;5 6 1;6 6 1;6 6 2;7 6 1];

SpikeAmp=[]; SpikeAUC=[]; g=1;
for n=1:size(BlueStimN,1) % neuron
    cax=[];
    Neuron=BlueStimN(n,1);
    rep=BlueStimN(n,3);
    BluePulseN=BlueStimN(n,2);
    wvf=1;
    if ~isempty(OpResult{wvf,rep,Neuron})

        normTr=OpResult{wvf,rep,Neuron}.normTraces./OpResult{wvf,rep,Neuron}.ScaleFactor;
        noi=[1:size(OpResult{wvf,rep,Neuron}.normTraces,1)];
        %noi=OpResult{wvf,rep,Neuron}.maintrunkROI;
        normTr=normTr(noi,:);

        nROI=length(noi);
        nTime=size(OpResult{wvf,rep,Neuron}.normTraces,2);
        dOrder=OpResult{wvf,rep,Neuron}.dist_order(noi);
        %[~, dOrder]=sort(OpResult{wvf,rep,Neuron}.dist_order(noi),'ascend');
        Classvec = get_Class2index(OpResult{wvf,rep,Neuron}.SpClass([1 2],:));
        som_spike = Classvec .* OpResult{wvf,rep,Neuron}.spike(1,:);

        bwBlue=bwlabel(OpResult{wvf,rep,Neuron}.Blue);
        bwBlue((bwBlue-(max(bwBlue)-5))<0)=0;
        bwBlue=bwlabel(bwBlue>0);

        [STAbAP]=get_STA(normTr,max(som_spike>0,[],1),-nTau(1),nTau(2));
        [~, tmax]=max(STAbAP(:,-nTau(1)+[-2:4]),[],2);
        tmax=tmax-nTau(1)-3;

        SpikeinBlue=max(som_spike,[],1).*(bwBlue==BluePulseN);
        SpikeinBlue_ind=find(SpikeinBlue);
        BlueOnset=find(bwBlue==BluePulseN,1);
        SubLarge=get_subthreshold(normTr, max([som_spike; bwBlue>0],[],1), 400,3000);
        normTrSub=normTr-SubLarge;

        Rheobase=mean(OpResult{wvf,rep,Neuron}.Blue(SpikeinBlue_ind(1:2)));


        Dsign=ones(1,size(OpResult{wvf,rep,Neuron}.interDendDist,2));
        Dsign(OpResult{wvf,rep,Neuron}.dist_order(1:find(OpResult{wvf,rep,Neuron}.dist_order==1)-1))=-1;
        dendaxis=OpResult{wvf,rep,Neuron}.interDendDist(1,:).*Dsign;
        dendaxis=dendaxis(noi);
        distFromSoma=dendaxis*OpResult{wvf,rep,Neuron}.pixelsize;
        perisomaInd=find(distFromSoma'>perisomadist(1) & distFromSoma'<perisomadist(2));

        SpikeMat=permute(reshape(normTr(:,SpikeinBlue_ind'+nTauPeriSp),nROI,[],length(nTauPeriSp)),[1 3 2]); %1: ROI, 2:time, 3:event
        %SpikeAmp{g}=[[NaN SpikeinBlue_ind-BlueOnset]; [distFromSoma' squeeze(max(SpikeMat,[],2))]];
        RheoBase_trace=OpResult{wvf,rep,Neuron}.Blue(SpikeinBlue_ind)/Rheobase;
        SpikeAmp{g}=[[NaN RheoBase_trace]; [distFromSoma' squeeze(max(SpikeMat,[],2))]]; %RheoBase
        SpikeAmp{g}(2:end,2:end)=SpikeAmp{g}(2:end,2:end)./mean(SpikeAmp{g}(perisomaInd+1,2:end),[1]);
        SpikeAmp{g}(2:end,:)=SpikeAmp{g}(dOrder+1,:);
        [~, reorder_dist]=sort(SpikeAmp{g}(2:end,1),'ascend');
        SpikeAmp{g}(2:end,:)=SpikeAmp{g}(reorder_dist,:);

        SPdelay=round(abs(distFromSoma)/conductionspeed);
        sub=[];
        for r=1:size(normTr,1)
            sp_vec=ind2vec(size(normTr,2),find(som_spike)+SPdelay(r),1,0);
            [~, sub(r,:)]=get_subthreshold(normTr(r,:),sp_vec,2*(2+SPdelay(r))+1,15);
        end
        normTr_sub=normTr-sub;

        SpikeMat_sub=permute(reshape(normTr_sub(:,SpikeinBlue_ind'+nTauPeriSp),nROI,[],length(nTauPeriSp)),[1 3 2]); %1: ROI, 2:time, 3:event
        AUC=squeeze(sum(SpikeMat_sub,2,'omitnan'));
        %
        % for r=1:size(SpikeMat,1)
        %     for s=1:size(SpikeMat,3)
        %     AUC(r,s)=get_AUC(squeeze(SpikeMat(r,:,s)),-nTauPeriSp(1)+1+SPdelay(r),3,3);
        %     end
        % end

        %SpikeAUC{g}=[[NaN SpikeinBlue_ind-BlueOnset]; [distFromSoma' AUC]];
        SpikeAUC{g}=[[NaN RheoBase_trace]; [distFromSoma' AUC]]; %RheoBase
        SpikeAUC{g}(2:end,2:end)=SpikeAUC{g}(2:end,2:end)./mean(SpikeAUC{g}(perisomaInd+1,2:end),[1]);
        SpikeAUC{g}(2:end,:)=SpikeAUC{g}(dOrder+1,:);

        g=g+1;
    end
end

figure(14); clf; %tiledlayout(1,2);
nexttile([1 1]);
[binnedZ binX binY]=show3Dbinning(SpikeAmp, 11, 8, 'image'); hold all
shading flat
allX=[]; allY=[]; allZ=[];
for n = 1:numel(SpikeAmp)
    timeAxis = SpikeAmp{n}(1, 2:end);
    xAxis = SpikeAmp{n}(2:end, 1);
    zData = SpikeAmp{n}(2:end, 2:end);
    [X, Y] = meshgrid(timeAxis, xAxis);
    allX = [allX; X(:)];
    allY = [allY; Y(:)];
    allZ = [allZ; zData(:)];
end
%scatter3(allX,allY,allZ,20,[0.7 0.7 0.7],'filled')
set(gca, 'YDir', 'reverse');
colormap("turbo")
%zlim([0.3 1.6])
xlabel('Optical Rheobase')
ylabel('Distance from soma (\mum)')
title('Normalized Spike Amplitude')


%% Representative pulse/Ramp

normTr=OpResult{1,1,7}.normTraces./OpResult{1,1,7}.ScaleFactor;
normTr=normTr(OpResult{1,1,7}.dist_order,:);
normTr=pcafilterTrace(normTr,[1 2 3 4 9 10]);
cax=[-0.005 0.02];

Dsign=ones(1,size(OpResult{1,1,7}.interDendDist,2));
Dsign(OpResult{1,1,7}.dist_order(1:find(OpResult{1,1,7}.dist_order==1)-1))=-1;
dendaxis=OpResult{1,1,7}.interDendDist(1,:).*Dsign;
dendaxis=dendaxis(OpResult{1,1,7}.dist_order);

SomaROI=[6 10 16];
distalROI=[29 31 33];
dendriteaxis_bin=[-100:10:220];
taxis=[-99:600];
[~, sortind]=sort(dendaxis,'ascend');
[Xq, Yq] = meshgrid([taxis], dendriteaxis_bin);
normTr_interp = interp2([taxis], dendaxis(sortind), normTr(sortind,2901:3600), Xq, Yq, 'linear');
Sptime=find(OpResult{1,1,7}.spike(1,2901:3600)>0)-100;

figure(14); clf; tiledlayout(4,1);
ax1=nexttile([1 1]);
% pcolor(Xq,Yq,normTr_interp); hold all
plot(taxis,mean(normTr(SomaROI,2901:3600),1)); hold all
plot(taxis,mean(normTr(distalROI,2901:3600),1)); hold all
%plot(Sptime,normTr(SomaROI,Sptime+3000),'kv','markersize',3); axis tight off;
%set(gca,'YDir','reverse')
%caxis([cax])
%shading flat
%colormap(turbo);
ax2=nexttile([1 1]);
plot(taxis,OpResult{1,1,7}.Blue(2901:3600),'color',[0 0.6 1]); axis off;
linkaxes([ax1 ax2],'x')

N=7; wvf=3; rep=2;
SomaROI=[6 10 16];
distalROI=[29 31 33];
normTr=OpResult{wvf,rep,N}.normTraces./OpResult{wvf,rep,N}.ScaleFactor;
normTr=normTr(OpResult{wvf,rep,N}.dist_order,:);
%normTr=pcafilterTrace(normTr,[1 2 3 4 9 10]);
cax=[-0.005 0.02];

Dsign=ones(1,size(OpResult{wvf,rep,N}.interDendDist,2));
Dsign(OpResult{wvf,rep,N}.dist_order(1:find(OpResult{wvf,rep,N}.dist_order==1)-1))=-1;
dendaxis=OpResult{wvf,rep,N}.interDendDist(1,:).*Dsign;
dendaxis=dendaxis(OpResult{wvf,rep,N}.dist_order);

dendriteaxis_bin=[-100:10:220];
taxis=[1:9999]; time_show=[1:9999];
[~, sortind]=sort(dendaxis,'ascend');
[Xq, Yq] = meshgrid([taxis], dendriteaxis_bin);
normTr_interp = interp2([taxis], dendaxis(sortind), normTr(sortind,time_show), Xq, Yq, 'linear');
Sptime=find(OpResult{wvf,rep,N}.spike(1,time_show)>0)-100;
%figure(15); clf; tiledlayout(2,1);
ax3=nexttile([1 1]);
plot(taxis,mean(normTr(SomaROI,time_show),1)); hold all
plot(taxis,mean(normTr(distalROI,time_show),1)); hold all
%pcolor(Xq,Yq,normTr_interp); hold all
%plot(Sptime,-110,'k','marker','v','markersize',3); axis tight off;
%set(gca,'YDir','reverse')
%caxis([cax]);
%shading flat
%colormap(turbo);
ax4=nexttile([1 1]);
plot(taxis,OpResult{wvf,rep,N}.Blue(time_show),'color',[0 0.6 1]); axis off;
linkaxes([ax4 ax3],'x')

%% Dendrite-stim complex-spike example traces
figure(15); clf;
% show dStim, CS examples
BlueStimN=[6 1;6 2;7 1]; scale=0.03;
wvf=3;
for n=1:size(BlueStimN,1) % neuron
    cax=[];
    Neuron=BlueStimN(n,1);
    rep=BlueStimN(n,2);

    normTr=OpResult{wvf,rep,Neuron}.normTraces./OpResult{wvf,rep,Neuron}.ScaleFactor;
    tax=[1:size(normTr,2)];
    CStr=normTr(1,:);
    CStr(OpResult{wvf,rep,Neuron}.CStrace==0)=NaN;
    plot(tax,normTr(1,:)+scale*(n-1),'k');
    hold all
    plot(tax,CStr+scale*(n-1),'color',[0.7 0 0.1]);
end
%% Ramp (Soma Stim): normalized spike amplitude vs optical rheobase (main-trunk ROIs)
nTau=[-70:70]; nTauPeriSp=[-3:7]; nTauPeriSp2=[-3:3]; conductionspeed=170; %um/ms
g=ones(1,4); SP_STA=[]; distFromSoma=[];
dendriteaxis_bin=[-80:80:600];
perisomadist=[-50 50]; cmap=[0 0 0; 1 0 0];  cmap_light=[0.7 0.7 0.7; 1 0.7 0.7];
%BlueStimN=[2 5 1;3 3 1;3 4 1;4 2 1;4 3 1;4 4 1;5 3 1;5 5 1;6 4 1;6 5 1;6 3 2;6 4 2;6 5 2;7 1 1;7 2 1;7 3 1;];
BlueStimN=[2 6 1;3 6 1;4 6 1;5 6 1;6 6 1;6 6 2;7 6 1];

SpikeAmp=[]; SpikeAUC=[]; g=1;
for n=1:size(BlueStimN,1) % neuron
    cax=[];
    Neuron=BlueStimN(n,1);
    rep=BlueStimN(n,3);
    BluePulseN=BlueStimN(n,2);
    wvf=1;
    if ~isempty(OpResult{wvf,rep,Neuron})

        normTr=OpResult{wvf,rep,Neuron}.normTraces./OpResult{wvf,rep,Neuron}.ScaleFactor;
        %noi=[1:size(OpResult{wvf,rep,Neuron}.normTraces,1)];
        noi=OpResult{wvf,rep,Neuron}.maintrunkROI;
        normTr=normTr(noi,:);

        nROI=length(noi);
        nTime=size(OpResult{wvf,rep,Neuron}.normTraces,2);
        [~, dOrder]=sort(OpResult{wvf,rep,Neuron}.dist_order(noi),'ascend');
        som_spike=max(OpResult{wvf,rep,Neuron}.SpClass([1 2],:),[],1);
        bwBlue=bwlabel(OpResult{wvf,rep,Neuron}.Blue);
        bwBlue((bwBlue-(max(bwBlue)-5))<0)=0;
        bwBlue=bwlabel(bwBlue>0);
        SpikeinBlue=som_spike.*(bwBlue==BluePulseN);
        SpikeinBlue_ind=find(SpikeinBlue);
        BlueOnset=find(bwBlue==BluePulseN,1);
        Rheobase=mean(OpResult{wvf,rep,Neuron}.Blue(SpikeinBlue_ind(1:2)));


        normTr=normTr-prctile(normTr(:,BlueOnset+[-500:-200]),30,2);

        Dsign=ones(1,size(OpResult{wvf,rep,Neuron}.interDendDist,2));
        Dsign(OpResult{wvf,rep,Neuron}.dist_order(1:find(OpResult{wvf,rep,Neuron}.dist_order==1)-1))=-1;
        dendaxis=OpResult{wvf,rep,Neuron}.interDendDist(1,:).*Dsign;
        dendaxis=dendaxis(noi);
        distFromSoma=dendaxis*OpResult{wvf,rep,Neuron}.pixelsize;
        perisomaInd=find(distFromSoma'>perisomadist(1) & distFromSoma'<perisomadist(2));

        SpikeMat=permute(reshape(normTr(:,SpikeinBlue_ind'+nTauPeriSp),nROI,[],length(nTauPeriSp)),[1 3 2]); %1: ROI, 2:time, 3:event
        %SpikeAmp{g}=[[NaN SpikeinBlue_ind-BlueOnset]; [distFromSoma' squeeze(max(SpikeMat,[],2))]];
        RheoBase_trace=OpResult{wvf,rep,Neuron}.Blue(SpikeinBlue_ind)/Rheobase;
        SpikeAmp{g}=[[NaN RheoBase_trace]; [distFromSoma' squeeze(max(SpikeMat,[],2))]]; %RheoBase
        SpikeAmp{g}(2:end,2:end)=SpikeAmp{g}(2:end,2:end)./mean(SpikeAmp{g}(perisomaInd+1,2:end),[1]);
        SpikeAmp{g}(2:end,:)=SpikeAmp{g}(dOrder+1,:);
        [~, reorder_dist]=sort(SpikeAmp{g}(2:end,1),'ascend');
        SpikeAmp{g}(2:end,:)=SpikeAmp{g}(reorder_dist,:);

        SPdelay=round(abs(distFromSoma)/conductionspeed);
        sub=[];
        for r=1:size(normTr,1)
            sp_vec=ind2vec(size(normTr,2),find(som_spike)+SPdelay(r),1,0);
            [~, sub(r,:)]=get_subthreshold(normTr(r,:),sp_vec,2*(2+SPdelay(r))+1,15);
        end
        normTr_sub=normTr-sub;

        SpikeMat_sub=permute(reshape(normTr_sub(:,SpikeinBlue_ind'+nTauPeriSp),nROI,[],length(nTauPeriSp)),[1 3 2]); %1: ROI, 2:time, 3:event
        AUC=squeeze(sum(SpikeMat_sub,2,'omitnan'));
        %
        % for r=1:size(SpikeMat,1)
        %     for s=1:size(SpikeMat,3)
        %     AUC(r,s)=get_AUC(squeeze(SpikeMat(r,:,s)),-nTauPeriSp(1)+1+SPdelay(r),3,3);
        %     end
        % end

        %SpikeAUC{g}=[[NaN SpikeinBlue_ind-BlueOnset]; [distFromSoma' AUC]];
        SpikeAUC{g}=[[NaN RheoBase_trace]; [distFromSoma' AUC]]; %RheoBase
        SpikeAUC{g}(2:end,2:end)=SpikeAUC{g}(2:end,2:end)./mean(SpikeAUC{g}(perisomaInd+1,2:end),[1]);
        SpikeAUC{g}(2:end,:)=SpikeAUC{g}(dOrder+1,:);

        g=g+1;
    end
end

figure(14); clf; %tiledlayout(1,2);
nexttile([1 1]);
[binnedZ binX binY]=show3Dbinning(SpikeAmp, 11, 8, 'image'); hold all
shading flat
allX=[]; allY=[]; allZ=[];
for n = 1:numel(SpikeAmp)
    timeAxis = SpikeAmp{n}(1, 2:end);
    xAxis = SpikeAmp{n}(2:end, 1);
    zData = SpikeAmp{n}(2:end, 2:end);
    [X, Y] = meshgrid(timeAxis, xAxis);
    allX = [allX; X(:)];
    allY = [allY; Y(:)];
    allZ = [allZ; zData(:)];
end
%scatter3(allX,allY,allZ,20,[0.7 0.7 0.7],'filled')
set(gca, 'YDir', 'reverse');
colormap("turbo")
%zlim([0.3 1.6])
xlabel('Optical Rheobase')
ylabel('Distance from soma (\mum)')
title('Normalized Spike Amplitude')
%view(0,90)
%
% nexttile([1 1]);
% [binnedZ binX binY]=show3Dbinning(SpikeAUC, 11, 8, 'image'); hold all
% shading flat
% allX=[]; allY=[]; allZ=[];
%   for n = 1:numel(SpikeAUC)
%         timeAxis = SpikeAUC{n}(1, 2:end);
%         xAxis = SpikeAUC{n}(2:end, 1);
%         zData = SpikeAUC{n}(2:end, 2:end);
%         [X, Y] = meshgrid(timeAxis, xAxis);
%         allX = [allX; X(:)];
%         allY = [allY; Y(:)];
%         allZ = [allZ; zData(:)];
%     end
% %scatter3(allX,allY,allZ,20,[0.8 0.8 0.8],'filled')
% set(gca, 'YDir', 'reverse');
% colormap("turbo")
% %zlim([0.3 1.4])
% xlabel('Optical Rheobase')
% ylabel('Distance from soma (\mum)')
% title('Normalized Spike AUC')
% %view(0,90)

% %% Pulse (WF Stim)
% nTau=[-70:70]; nTauPeriSp=[-3:7]; nTauPeriSp2=[-3:3]; conductionspeed=170; %um/ms
% g=ones(1,4); SP_STA=[]; distFromSoma=[];
% dendriteaxis_bin=[-80:80:600];
% perisomadist=[-50 50]; cmap=[0 0 0; 1 0 0];  cmap_light=[0.7 0.7 0.7; 1 0.7 0.7];
% BlueStimN=[10 3 1;10 4 1;10 5 1;11 4 1;11 5 1;12 3 1;12 4 1;12 5 1;13 4 1;13 5 1; 14 5 1;15 4 1;15 5 1]; % Neuron number, Blue Pulse # , N th session
% %BlueStimN=[2 6 1;3 6 1;4 6 1;5 6 1;6 6 1;6 6 2;7 6 1];
%
% SpikeAmp=[]; SpikeAUC=[]; g=1;
% for n=1:size(BlueStimN,1) % neuron
%     cax=[];
%     Neuron=BlueStimN(n,1);
%     rep=BlueStimN(n,3);
%     BluePulseN=BlueStimN(n,2);
%     wvf=5;
%     if ~isempty(OpResult{wvf,rep,Neuron})
%
%         normTr=OpResult{wvf,rep,Neuron}.normTraces./OpResult{wvf,rep,Neuron}.ScaleFactor;
%         %noi=[1:size(OpResult{wvf,rep,Neuron}.normTraces,1)];
%         noi=OpResult{wvf,rep,Neuron}.maintrunkROI;
%         normTr=normTr(noi,:);
%
%         nROI=length(noi);
%         nTime=size(OpResult{wvf,rep,Neuron}.normTraces,2);
%         [~, dOrder]=sort(OpResult{wvf,rep,Neuron}.dist_order(noi),'ascend');
%         som_spike=max(OpResult{wvf,rep,Neuron}.SpClass([1 2],:),[],1);
%         bwBlue=bwlabel(OpResult{wvf,rep,Neuron}.Blue);
%         bwBlue((bwBlue-(max(bwBlue)-5))<0)=0;
%         bwBlue=bwlabel(bwBlue>0);
%         SpikeinBlue=som_spike.*(bwBlue==BluePulseN);
%         SpikeinBlue_ind=find(SpikeinBlue);
%         BlueOnset=find(bwBlue==BluePulseN,1);
%         normTr=normTr-prctile(normTr(:,BlueOnset+[-500:-200]),30,2);
%
%         Dsign=ones(1,size(OpResult{wvf,rep,Neuron}.interDendDist,2));
%         Dsign(OpResult{wvf,rep,Neuron}.dist_order(1:find(OpResult{wvf,rep,Neuron}.dist_order==1)-1))=-1;
%         dendaxis=OpResult{wvf,rep,Neuron}.interDendDist(1,:).*Dsign;
%         dendaxis=dendaxis(noi);
%         distFromSoma=dendaxis*OpResult{wvf,rep,Neuron}.pixelsize;
%         perisomaInd=find(distFromSoma'>perisomadist(1) & distFromSoma'<perisomadist(2));
%
%         SpikeMat=permute(reshape(normTr(:,SpikeinBlue_ind'+nTauPeriSp),nROI,[],length(nTauPeriSp)),[1 3 2]); %1: ROI, 2:time, 3:event
%         SpikeAmp{g}=[[NaN SpikeinBlue_ind-BlueOnset]; [distFromSoma' squeeze(max(SpikeMat,[],2))]];
%         SpikeAmp{g}(2:end,2:end)=SpikeAmp{g}(2:end,2:end)./mean(SpikeAmp{g}(perisomaInd+1,2:end),[1]);
%         SpikeAmp{g}(2:end,:)=SpikeAmp{g}(dOrder+1,:);
%
%         SPdelay=round(abs(distFromSoma)/conductionspeed);
%         sub=[];
%         for r=1:size(normTr,1)
%             sp_vec=ind2vec(size(normTr,2),find(som_spike)+SPdelay(r),1,0);
%         [~, sub(r,:)]=get_subthreshold(normTr(r,:),sp_vec,2*(2+SPdelay(r))+1,15);
%         end
%         normTr_sub=normTr-sub;
%
%         SpikeMat_sub=permute(reshape(normTr_sub(:,SpikeinBlue_ind'+nTauPeriSp),nROI,[],length(nTauPeriSp)),[1 3 2]); %1: ROI, 2:time, 3:event
%         AUC=squeeze(sum(SpikeMat_sub,2,'omitnan'));
%         %
%         % for r=1:size(SpikeMat,1)
%         %     for s=1:size(SpikeMat,3)
%         %     AUC(r,s)=get_AUC(squeeze(SpikeMat(r,:,s)),-nTauPeriSp(1)+1+SPdelay(r),3,3);
%         %     end
%         % end
%
%         SpikeAUC{g}=[[NaN SpikeinBlue_ind-BlueOnset]; [distFromSoma' AUC]];
%         SpikeAUC{g}(2:end,2:end)=SpikeAUC{g}(2:end,2:end)./mean(SpikeAUC{g}(perisomaInd+1,2:end),[1]);
%         SpikeAUC{g}(2:end,:)=SpikeAUC{g}(dOrder+1,:);
%
%         g=g+1;
%     end
% end
%
% figure(13); clf; tiledlayout(1,2);
% nexttile([1 1]);
% [binnedZ binX binY]=show3Dbinning(SpikeAmp, 9, 6, 'image'); hold all
% allX=[]; allY=[]; allZ=[];
%   for n = 1:numel(SpikeAmp)
%         timeAxis = SpikeAmp{n}(1, 2:end);
%         xAxis = SpikeAmp{n}(2:end, 1);
%         zData = SpikeAmp{n}(2:end, 2:end);
%         [X, Y] = meshgrid(timeAxis, xAxis);
%         allX = [allX; X(:)];
%         allY = [allY; Y(:)];
%         allZ = [allZ; zData(:)];
%     end
% scatter3(allX,allY,allZ,20,[0.7 0.7 0.7],'filled')
% set(gca, 'YDir', 'reverse');
% colormap("turbo")
% %zlim([0.3 1.6])
% xlabel('Time after blue onset (ms)')
% ylabel('Distance from soma (\mum)')
% title('Normalized Spike Amplitude')
% %view(0,90)
%
% nexttile([1 1]);
% [binnedZ binX binY]=show3Dbinning(SpikeAUC, 9, 7, 'image'); hold all
% allX=[]; allY=[]; allZ=[];
%   for n = 1:numel(SpikeAUC)
%         timeAxis = SpikeAUC{n}(1, 2:end);
%         xAxis = SpikeAUC{n}(2:end, 1);
%         zData = SpikeAUC{n}(2:end, 2:end);
%         [X, Y] = meshgrid(timeAxis, xAxis);
%         allX = [allX; X(:)];
%         allY = [allY; Y(:)];
%         allZ = [allZ; zData(:)];
%     end
% %scatter3(allX,allY,allZ,20,[0.8 0.8 0.8],'filled')
% set(gca, 'YDir', 'reverse');
% colormap("turbo")
% %zlim([0.3 1.4])
% xlabel('Time after blue onset (ms)')
% ylabel('Distance from soma (\mum)')
% title('Normalized Spike AUC')
% %view(0,90)

% %% Ramp (WF Stim)
% nTau=[-70:70]; nTauPeriSp=[-2:3]; nTauPeriSp2=[-3:3]; conductionspeed=170; %um/ms
% g=ones(1,4); SP_STA=[]; distFromSoma=[];
% dendriteaxis_bin=[-80:80:600];
% perisomadist=[-10 10]; cmap=[0 0 0; 1 0 0];  cmap_light=[0.7 0.7 0.7; 1 0.7 0.7];
% BlueStimN=[10 6 1;12 6 1;13 6 1;14 6 1];
%
% SpikeAmp=[]; SpikeAUC=[]; g=1;
% for n=1:size(BlueStimN,1) % neuron
%     cax=[];
%     Neuron=BlueStimN(n,1);
%     rep=BlueStimN(n,3);
%     BluePulseN=BlueStimN(n,2);
%     wvf=5;
%     if ~isempty(OpResult{wvf,rep,Neuron})
%
%         normTr=OpResult{wvf,rep,Neuron}.normTraces./OpResult{wvf,rep,Neuron}.ScaleFactor;
%         %noi=[1:size(OpResult{wvf,rep,Neuron}.normTraces,1)];
%         noi=OpResult{wvf,rep,Neuron}.maintrunkROI;
%         normTr=normTr(noi,:);
%
%         nROI=length(noi);
%         nTime=size(OpResult{wvf,rep,Neuron}.normTraces,2);
%         [~, dOrder]=sort(OpResult{wvf,rep,Neuron}.dist_order(noi),'ascend');
%         som_spike=max(OpResult{wvf,rep,Neuron}.SpClass([1 2],:),[],1);
%         bwBlue=bwlabel(OpResult{wvf,rep,Neuron}.Blue);
%         bwBlue((bwBlue-(max(bwBlue)-5))<0)=0;
%         bwBlue=bwlabel(bwBlue>0);
%         SpikeinBlue=som_spike.*(bwBlue==BluePulseN);
%         SpikeinBlue_ind=find(SpikeinBlue);
%         BlueOnset=find(bwBlue==BluePulseN,1);
%         Rheobase=mean(OpResult{wvf,rep,Neuron}.Blue(SpikeinBlue_ind(1:2)));
%
%         normTr=normTr-prctile(normTr(:,BlueOnset+[-500:-200]),30,2);
%
%         Dsign=ones(1,size(OpResult{wvf,rep,Neuron}.interDendDist,2));
%         Dsign(OpResult{wvf,rep,Neuron}.dist_order(1:find(OpResult{wvf,rep,Neuron}.dist_order==1)-1))=-1;
%         dendaxis=OpResult{wvf,rep,Neuron}.interDendDist(1,:).*Dsign;
%         dendaxis=dendaxis(noi);
%         distFromSoma=dendaxis*OpResult{wvf,rep,Neuron}.pixelsize;
%         perisomaInd=find(distFromSoma'>perisomadist(1) & distFromSoma'<perisomadist(2));
%
%         SpikeMat=permute(reshape(normTr(:,SpikeinBlue_ind'+nTauPeriSp),nROI,[],length(nTauPeriSp)),[1 3 2]); %1: ROI, 2:time, 3:event
%         %SpikeAmp{g}=[[NaN SpikeinBlue_ind-BlueOnset]; [distFromSoma' squeeze(max(SpikeMat,[],2))]];
%         RheoBase_trace=OpResult{wvf,rep,Neuron}.Blue(SpikeinBlue_ind)/Rheobase;
%         SpikeAmp{g}=[[NaN RheoBase_trace]; [distFromSoma' squeeze(max(SpikeMat,[],2))]]; %RheoBase
%         SpikeAmp{g}(2:end,2:end)=SpikeAmp{g}(2:end,2:end)./mean(SpikeAmp{g}(perisomaInd+1,2:end),[1]);
%         %SpikeAmp{g}(2:end,:)=SpikeAmp{g}(dOrder+1,:);
%
%         SPdelay=round(abs(distFromSoma)/conductionspeed);
%         sub=[];
%         for r=1:size(normTr,1)
%             sp_vec=ind2vec(size(normTr,2),find(som_spike)+SPdelay(r),1,0);
%         [~, sub(r,:)]=get_subthreshold(normTr(r,:),sp_vec,2*(2+SPdelay(r))+1,15);
%         end
%         normTr_sub=normTr-sub;
%
%         SpikeMat_sub=permute(reshape(normTr_sub(:,SpikeinBlue_ind'+nTauPeriSp),nROI,[],length(nTauPeriSp)),[1 3 2]); %1: ROI, 2:time, 3:event
%         AUC=squeeze(sum(SpikeMat_sub,2,'omitnan'));
%         %
%         % for r=1:size(SpikeMat,1)
%         %     for s=1:size(SpikeMat,3)
%         %     AUC(r,s)=get_AUC(squeeze(SpikeMat(r,:,s)),-nTauPeriSp(1)+1+SPdelay(r),3,3);
%         %     end
%         % end
%
%         %SpikeAUC{g}=[[NaN SpikeinBlue_ind-BlueOnset]; [distFromSoma' AUC]];
%         SpikeAUC{g}=[[NaN RheoBase_trace]; [distFromSoma' AUC]]; %RheoBase
%         SpikeAUC{g}(2:end,2:end)=SpikeAUC{g}(2:end,2:end)./mean(SpikeAUC{g}(perisomaInd+1,2:end),[1]);
%         %SpikeAUC{g}(2:end,:)=SpikeAUC{g}(dOrder+1,:);
%
%         g=g+1;
%     end
% end
%
% figure(14); clf; tiledlayout(1,2);
% nexttile([1 1]);
% [binnedZ binX binY]=show3Dbinning(SpikeAmp, 11, 8, 'image'); hold all
% allX=[]; allY=[]; allZ=[];
%   for n = 1:numel(SpikeAmp)
%         timeAxis = SpikeAmp{n}(1, 2:end);
%         xAxis = SpikeAmp{n}(2:end, 1);
%         zData = SpikeAmp{n}(2:end, 2:end);
%         [X, Y] = meshgrid(timeAxis, xAxis);
%         allX = [allX; X(:)];
%         allY = [allY; Y(:)];
%         allZ = [allZ; zData(:)];
%     end
% %scatter3(allX,allY,allZ,20,[0.7 0.7 0.7],'filled')
% set(gca, 'YDir', 'reverse');
% colormap("turbo")
% %zlim([0.3 1.6])
% xlabel('Optical Rheobase')
% ylabel('Distance from soma (\mum)')
% title('Normalized Spike Amplitude')
% %view(0,90)
%
% nexttile([1 1]);
% [binnedZ binX binY]=show3Dbinning(SpikeAUC, 11, 8, 'image'); hold all
% allX=[]; allY=[]; allZ=[];
%   for n = 1:numel(SpikeAUC)
%         timeAxis = SpikeAUC{n}(1, 2:end);
%         xAxis = SpikeAUC{n}(2:end, 1);
%         zData = SpikeAUC{n}(2:end, 2:end);
%         [X, Y] = meshgrid(timeAxis, xAxis);
%         allX = [allX; X(:)];
%         allY = [allY; Y(:)];
%         allZ = [allZ; zData(:)];
%     end
% %scatter3(allX,allY,allZ,20,[0.8 0.8 0.8],'filled')
% set(gca, 'YDir', 'reverse');
% colormap("turbo")
% %zlim([0.3 1.4])
% xlabel('Optical Rheobase')
% ylabel('Distance from soma (\mum)')
% title('Normalized Spike AUC')
% %view(0,90)

%% Dendrite targeted stimulation



