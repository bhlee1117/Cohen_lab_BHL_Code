% Run this after running 2nd session of Analysis_20250331_FigureMakingSeeSaw

nBin=100; nShuffles=2500;
%foi2=[15 16 17 18 19 20 21 22 23 24 25 27];
for f=foi(1:end)
    f
    load(fullfile(fpath{f},'PC_Result.mat'),'Result');

    nTime=size(Subthreshold{f},2);
    nROI=size(Subthreshold{f},1);

    [SubV SubD subTrace]=get_eigvector(Subthreshold{f}(:,sum(isnan(Subthreshold{f}),1)==0)',nROI);
    subTrace_onFrame=NaN(nROI,nTime);
    subTrace_onFrame(:,sum(isnan(Subthreshold{f}),1)==0)=subTrace';
    subTrace_onFrame(:,BlueStim{f}>0)=NaN;

    BlueOffFrame=BlueStim{f}==0;
    RunTime=VRtrack{f}(end,:)>0.001;
    StimOnLap=unique(VRtrack{f}(8,:).*(BlueStim{f}>0));
    TrainingEnvLap=unique(VRtrack{f}(8,:).*(VRtrack{f}(2,:)==1));
    StimOfflap=setdiff(unique(VRtrack{f}(8,:)),[StimOnLap TrainingEnvLap]);
    FrameOfInterest= BlueOffFrame & RunTime & ismember(VRtrack{f}(8,:),StimOfflap);

    % [Result.ShuffleFR_percentile Result.ShuffleFR_zscore] = position_selectivity(allSpikeMat{f}(1,FrameOfInterest),VRtrack{f}(5,FrameOfInterest),VRtrack{f}(8,FrameOfInterest) ...
    %     ,'nBins',nBin,'nShuffles',nShuffles,'Show','off');

    % [Result.SI_FRreal Result.SI_FRnull] = SpatialInformation_shuffling(allSpikeMat{f}(1,FrameOfInterest), VRtrack{f}(5,FrameOfInterest), VRtrack{f}(8,FrameOfInterest), nShuffles, 150);
    % for dclass=1:4
    % [Result.SI_FRClassReal(dclass,1) Result.SI_FRClassnull(dclass,:)] = SpatialInformation_shuffling(allSpikeClassVecMat{f}(dclass,FrameOfInterest), VRtrack{f}(5,FrameOfInterest), VRtrack{f}(8,FrameOfInterest), nShuffles, 150);
    % end
    % 
    % for pc=1:10
    %     pc
    %     % [Result.ShuffleEigTr_percentile(pc,:) Result.ShuffleEigTr_zscore(pc,:)] = position_selectivity(subTrace_onFrame(pc,FrameOfInterest),VRtrack{f}(5,FrameOfInterest),VRtrack{f}(8,FrameOfInterest) ...
    %     %     ,'nBins',nBin,'nShuffles',nShuffles,'Show','off');
    % 
    %     Result.MI_EigTrReal(pc,1) = mi_cont_cont(subTrace_onFrame(pc,FrameOfInterest), VRtrack{f}(5,FrameOfInterest));
    %     Result.MI_EigTrFR(pc,:) = mi_discrete_cont(subTrace_onFrame(pc,FrameOfInterest), allSpikeMat{f}(1,FrameOfInterest));
    % 
    %     for dclass=1:4
    %         [Result.MI_EigTrFRdClass(dclass,pc)] = mi_discrete_cont(subTrace_onFrame(pc,FrameOfInterest), allSpikeClassVecMat{f}(dclass,FrameOfInterest));
    %     end
    % end

    for v=1:nROI
        v
        % [Result.ShuffleV_percentile(v,:) Result.ShuffleV_zscore(v,:)] = position_selectivity(Subthreshold{f}(v,FrameOfInterest),VRtrack{f}(5,FrameOfInterest),VRtrack{f}(8,FrameOfInterest) ...
        %     ,'nBins',nBin,'nShuffles',nShuffles,'Show','off');

        % Mutual information between subthreshold and VR track
        Result.MI_SubVReal(v,:) = mi_cont_cont(Subthreshold{f}(v,FrameOfInterest), VRtrack{f}(5,FrameOfInterest));
        % Mutual information between spike and subthreshold
        Result.MI_SubVFR(v,:) = mi_discrete_cont(Subthreshold{f}(v,FrameOfInterest), allSpikeMat{f}(1,FrameOfInterest));
        for dclass=1:4
            [Result.MI_SubVFRdClass(dclass,v)] = mi_discrete_cont(Subthreshold{f}(v,FrameOfInterest), allSpikeClassVecMat{f}(dclass,FrameOfInterest));
        end
    end

    save(fullfile(fpath{f},'PC_Result.mat'),'Result','-v7.3');
    try
        load(fullfile(fpath{f},'PC_Result.mat'),'Result');
    catch
        save(fullfile(fpath{f},'PC_Result.mat'),'Result','-v7.3');
    end
    load(fullfile(fpath{f},'PC_Result.mat'),'Result');
end

