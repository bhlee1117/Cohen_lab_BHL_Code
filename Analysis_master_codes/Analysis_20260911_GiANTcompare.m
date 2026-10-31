%% GIAnT (SILo) vs. extractGluSNFR3 — glutamate extraction benchmark
%
%  Implements GluSNFR_GIAnT_Benchmark_Runbook.md against the actual code:
%    master pipeline : Analysis_master_codes/Analysis_20260713_SD_VoltronST_plus_Glu.m
%    extractor A     : Cohen_lab_BHL_Code/extractGluSNFR3.m
%    extractor B     : GIAnT-MATLAB/source_extraction/SILo.m
%
%  Part A  real-movie head-to-head on one artifact-corrected 1000-frame movie
%  Part B  paired ground-truth simulation (200x200x1000 real background + injected sources)
%  Part C  end-to-end check: inject BEFORE artifact regression, then run everything
%
%  ------------------------------------------------------------------------
%  FACTS ESTABLISHED BY READING THE THREE SOURCES (runbook section 14)
%  ------------------------------------------------------------------------
%  Movie loading / orientation
%    readBinMov_times(fn, nRow, nCol, frameList) -> nRow x nCol x T uint16 (Y x X x T).
%    Glutamate (blue) camera: Device_Data{4}.ROI([2 4]) = [nCol nRow];
%    exposuretime = Device_Data{4}.exposuretime (seconds).  Chunks are
%    'mc2%02d.bin' with motion traces 'mc2Trace%02d.mat' (var mcTrace2, field .xymean).
%
%  Artifact regression (master Section 2, verbatim operations)
%    meanF   = mean over space
%    y_fit   = expfitDM_2((1:T)', meanF, (1:T)', [10000 100])       % bleach
%    filt    = mov ./ y_fit                                          % bleach-normalised
%    filt_lw = movmedian(filt, round(5/exposuretime), 3)             % 5 s slow baseline
%    sub     = filt - filt_lw                                        % high-pass residual
%    [V,~,~] = get_eigvector(tovec(imgaussfilt(sub,1))', 40)
%    sub     = SeeResiduals(sub, V(:,nPC))                           % PC regression
%    sub     = SeeResiduals(sub, mcTrace2.xymean)                    % motion regression
%    filt2   = sub + filt_lw                                         % re-add slow baseline
%    The master script picks nPC interactively; here it is a frozen cfg value
%    (cfg.regressPCs) so the benchmark is reproducible and identical for both
%    methods.  PCA is refit inside every movie (runbook leakage warning).
%
%  extractGluSNFR3 contract
%    R = extractGluSNFR3(mov2_sub, mov2_mc_filt, exposuretime2, opts [, S_glu_in])
%    mov2_sub     : high-pass residual movie (Y x X x T)     -> detection + T_glu
%    mov2_mc_filt : bleach-normalised fluorescence movie     -> F0_glu, dFF_glu
%    Returns .S_glu (Y x X x K, L2-normalised, non-negative z-score footprints),
%            .T_glu (K x T), .F0_glu (K x T), .dFF_glu = T_glu ./ mean(F0_glu,2),
%            .synCentroid (K x 2 [x y]), .synEventCount, .eventList, .opts.
%    Internal preprocessing: imgaussfilt(sigma=1) + robust per-pixel MAD z-score.
%    It does NOT detrend or normalise on its own — the caller supplies both movies.
%
%  GIAnT contract
%    Entry point is SILo(pathToTrialTable, params); it is a FILE-BASED pipeline,
%    there is no in-memory movie function.  The adapter below writes a minimal
%    GIAnT workspace and reads experiment_summary.h5 back:
%      <run>/motion_correction/<stem>_REGISTERED.h5      movie, /mov = Y x X x T
%      <run>/motion_correction/<stem>_ALIGNMENTDATA.h5   numChannels, frametime,
%                                                        alignHz, motionC/R,
%                                                        motionDSc/r, recNegErr, meanIM
%      <run>/trial_table.h5                              datadr, savedr, filename,
%                                                        motion_correction.{fn_reg_ds,
%                                                        fn_raw, fn_adata,
%                                                        registration_failed}
%    loadAndProcessTrialAsync accepts '.h5' registered movies (single root
%    dataset) and derives the alignment filename by truncating at '_REGISTERED',
%    hence the naming above.  Outputs read back:
%      sources/spatial/profiles   K x 1 x Y x X
%      sources/spatial/coords     K x 3, 0-indexed [z y x]
%      sources/temporal/{dF_ls,dF_denoised,events,F0}  K x nChan x T
%      sources/temporal/SNR       K x 1
%    Internal preprocessing: 5th-percentile baseline subtraction, movmedian F0,
%    variance-stabilising photon scaling, exponential matched filter (tau_s),
%    difference-of-Gaussians, then constrained-NMF source extraction.  It needs
%    POSITIVE FLUORESCENCE COUNTS — it cannot take a mean-subtracted residual
%    or a dF/F movie.
%
%  Common representation (runbook section 4)
%    Because extractGluSNFR3 wants a high-pass residual and GIAnT wants photon
%    counts, no single array satisfies both.  The least-processed common
%    representation is the artifact-corrected movie IN CAMERA COUNTS:
%        movCleanF = (sub + filt_lw) .* y_fit
%    Both adapters receive a byte-identical copy of movCleanF and each applies
%    only its own documented internal transform, recorded in
%    result.extra.preproc.  This deviation from "byte-identical input to the
%    algorithm core" is unavoidable and is reported (runbook section 11).
%
%  Frame rate  : 1/Device_Data{4}.exposuretime, read from the recording.
%  Pixel size  : NOT stored in Prism_V2+Glu_Data_Arrangement.xlsx (its columns
%                are Data Location / Quality / Mouse / Ind / EndFrame / Morphology
%                / Misc / VR data / VR_nth) and not in Device_Data.
%                cfg.pixelSize_um MUST be set by hand in makeBenchmarkConfig.
%
%  Byung Hun Lee, 2026-09-11

%% ========================================================================
%  Section 0 : Setup, paths, configuration, provenance
% ========================================================================
clear; clc;
set(0,'DefaultFigureWindowStyle','docked')

giantRoot = '/Users/bhlee1117/Documents/GitHub/GIAnT-MATLAB';
labRoot   = '/Users/bhlee1117/Documents/GitHub/Cohen_lab_BHL_Code';
addpath(genpath(labRoot));
addpath(genpath(giantRoot));

% Remove any shadow of the built-in 'graph' (old cvx/sdpt3 ships graph.m and
% breaks the centroid clustering inside extractGluSNFR3) — same guard as the
% master script.
gshadow = which('graph','-all');
for gi = 1:numel(gshadow)
    if ~contains(gshadow{gi}, fullfile('toolbox','matlab'))
        rmpath(fileparts(gshadow{gi}));
        fprintf('Removed shadowing graph.m folder from path: %s\n', fileparts(gshadow{gi}));
    end
end

cfg = makeBenchmarkConfig();
rng(cfg.seed,'twister');

outRoot = fullfile(cfg.outParent, ['GIAnT_vs_extractGluSNFR3_' datestr(now,'yyyymmdd_HHMMSS')]);
mkdir(outRoot);
mkdir(fullfile(outRoot,'real_movie'));
mkdir(fullfile(outRoot,'simulation'));
mkdir(fullfile(outRoot,'figures'));
mkdir(fullfile(outRoot,'results'));
mkdir(fullfile(outRoot,'giant_runs'));
cfg.outRoot = outRoot;
cfg.pixelSize_um=0.7312; %µm
fprintf('Output root: %s\n', outRoot);

prov = collectProvenance(cfg, labRoot, giantRoot);
save(fullfile(outRoot,'results','config_and_provenance.mat'),'cfg','prov','-v7.3');

%% ========================================================================
%  Section 1 : Load and standardize the real glutamate movie
% ========================================================================
%  The benchmark segment is chosen by a FROZEN rule written down before any
%  extraction was run, and recorded in makeBenchmarkConfig:
%    recording : the first sheet entry graded Quality == 'Good' that has
%                glutamate chunks on disk (readGluMeta asserts the grade)
%    tile      : the first chunk, mc201.bin
%    frames    : cfg.realStartFrame onwards, far enough in that the 488 nm
%                modulation onset and the initial bleach transient are excluded
%  It was NOT chosen after looking at either method's output.
%
%  Only the frames actually needed are read: the 1000 benchmark frames plus
%  cfg.contextWindows x (5 s movmedian window) of context on each side, so the
%  bleach fit and the slow-baseline filter are not edge-limited at the window
%  boundaries.  Reading a whole 644x480 x 5000-frame tile in double is 12 GB;
%  this plus cfg.realCrop keeps the working set manageable.

meta = readGluMeta(cfg);
cfg.frameRate_Hz = meta.frameRate_Hz;
cfg.exposuretime = meta.exposuretime;
assert(isfinite(cfg.pixelSize_um) && cfg.pixelSize_um > 0, ...
    ['cfg.pixelSize_um is not set.  It is not recorded in the data-arrangement ' ...
     'sheet or in Device_Data; set it by hand in makeBenchmarkConfig() from the ' ...
     'objective/camera calibration for this recording (%s).'], meta.fpath);
cfg = derivePixelUnits(cfg);

win = planReadWindow(meta, cfg);
fprintf('Reading frames %d-%d of tile %d (%d frames, benchmark window %d-%d)\n', ...
    win.first, win.last, cfg.realTile, win.nFrames, win.benchFirst, win.benchLast);
fprintf('  estimated working set: %.1f GB\n', win.estGB);

[rawWin, meta] = readGluFrames(meta, win, cfg);

assert(ndims(rawWin) == 3, 'Movie must be Y x X x T.');
assert(all(isfinite(rawWin(:))), 'Movie contains NaN or Inf values.');

% QC on the raw benchmark window, before any regression
qc0 = win.benchIdx;
rawQC.meanImg     = mean(rawWin(:,:,qc0), 3);
rawQC.stdImg      = std(double(rawWin(:,:,qc0)), 0, 3);
rawQC.globalTrace = squeeze(mean(rawWin(:,:,qc0), [1 2], 'omitnan'));
rawQC.window      = win;
save(fullfile(outRoot,'real_movie','raw_QC.mat'),'rawQC','meta','win','-v7.3');

%% ========================================================================
%  Section 2 : Artifact regression (master-pipeline logic, non-interactive)
% ========================================================================
%  With cfg.regressPCsInteractive = true this is the one place that asks which
%  PCs to regress out (the master script's "PCs to regress out" prompt).
[cleanWin, artifactModel, artifactQC] = regressArtifactWrapper(rawWin, meta, cfg);
assert(isequal(size(cleanWin.F), size(rawWin)), 'Regression changed the movie size.');

%  Freeze the choice.  Part C regresses the movie again once per repeat, and
%  every one of those must use the same artifact model as Part A — otherwise the
%  attenuation measurement and the two extractors would each be looking at a
%  differently corrected movie.  Clearing the flag also stops the prompt from
%  reappearing mid-run.
cfg.regressPCs            = artifactModel.regressPCs;
cfg.regressPCsInteractive = true;
fprintf('Artifact model frozen: regressing PCs %s + motion (x,y) for the rest of the run.\n', ...
    mat2str(cfg.regressPCs));
prov.cfg = cfg;
prov.regressPCsChosenInteractively = artifactModel.pcsChosenInteractively;
save(fullfile(outRoot,'results','config_and_provenance.mat'),'cfg','prov','-v7.3');

% Crop the frozen benchmark window out of each preprocessed representation
movCleanF = single(cleanWin.F      (:,:,qc0));  % camera counts, artifact-corrected
movSub    = single(cleanWin.sub    (:,:,qc0));  % high-pass residual
movLw     = single(cleanWin.filt_lw(:,:,qc0));  % slow baseline (bleach-normalised)
yFit      = cleanWin.y_fit(qc0);
clear cleanWin rawWin

artifactQC = augmentArtifactQC(artifactQC, movCleanF, movSub, cfg);
makeArtifactQCFigures(artifactQC, cfg, fullfile(outRoot,'figures'));
save(fullfile(outRoot,'real_movie','artifact_model.mat'), ...
    'artifactModel','artifactQC','yFit','-v7.3');

%% ========================================================================
%  Section 3 : Normalize the input contract
% ========================================================================
%  movCleanF (camera counts) is the common input.  F0 is defined exactly as the
%  master pipeline defines it — the temporal mean of the bleach-normalised
%  slow-baseline movie, rescaled to counts — NOT the runbook's placeholder 20th
%  percentile, because the pipeline already has an F0 definition.

F0img = single(mean(double(movLw) .* reshape(yFit(:),1,1,[]), 3));
F0img = max(F0img, eps('single'));
cfg.F0definition = 'mean_t( movmedian(mov./expfit, 5 s) .* expfit )   [master-pipeline baseline]';
save(fullfile(outRoot,'real_movie','F0.mat'),'F0img','cfg','-v7.3');
save(fullfile(outRoot,'real_movie','movCleanF.mat'),'movCleanF','-v7.3');
clear movLw movSub

inputForBoth = movCleanF;      % byte-identical copy handed to both adapters

%% ========================================================================
%  Section 4 : Run both extractors on the real movie
% ========================================================================
resExtract = runExtractGluSNFR3(inputForBoth, cfg.extractParams, cfg, 'realA_extract');
resGIAnT   = runGIAnT          (inputForBoth, cfg.giantParams,   cfg, 'realA_giant');

save(fullfile(outRoot,'real_movie','resExtract.mat'),'resExtract','-v7.3');
save(fullfile(outRoot,'real_movie','resGIAnT.mat'),  'resGIAnT',  '-v7.3');
fprintf('\nReal movie: extractGluSNFR3 %d components in %.1f s | GIAnT %d in %.1f s\n', ...
    size(resExtract.footprints,3), resExtract.runtime_s, ...
    size(resGIAnT.footprints,3),   resGIAnT.runtime_s);

%% ========================================================================
%  Section 5 : Compare components in the real movie
% ========================================================================
stdExtract = standardizeFootprints(resExtract, cfg);
stdGIAnT   = standardizeFootprints(resGIAnT,   cfg);

realMatches = matchComponents(stdExtract, stdGIAnT, cfg);
realMetrics = compareRealComponents(movCleanF, F0img, stdExtract, stdGIAnT, realMatches, cfg);

makeRealMovieFigures(movCleanF, stdExtract, stdGIAnT, realMatches, realMetrics, ...
                     cfg, fullfile(outRoot,'figures'));
save(fullfile(outRoot,'results','realMovieComparison.mat'), ...
    'realMatches','realMetrics','stdExtract','stdGIAnT','-v7.3');

%% ========================================================================
%  Section 6 : Unmatched components — did either method miss a real signal?
% ========================================================================
%  A component found by one method and not the other is either a real source
%  the other method MISSED, or a FALSE POSITIVE.  A footprint picture cannot
%  tell those apart, so for every unmatched component this pulls the trace the
%  MOVIE ITSELF supports over that footprint's support — the same model-free
%  common-mask dF/F used for matched pairs in Section 5.  Events present in
%  that trace mean the other method missed a real source; a flat trace means
%  the component is a false positive.
%
%  This is evidence, not proof: the common-mask trace is derived from the same
%  movie, so a signal too weak to show up there is invisible to this test too.
%  It bounds the misses from below.  Part B measures misses against genuine
%  ground truth; this is the only handle available on the real movie, where no
%  ground truth exists.
%
%  Section 6 reports counts and shows examples; it never labels a component a
%  miss or a false positive -- that judgement is yours, from the figure and
%  results/unmatched_components.csv.
%
%  It also records WHY each component went unmatched, which the counts hide:
%  no counterpart nearby at all, a counterpart close enough but below the IoU
%  rule, or a qualifying counterpart that the one-to-one assignment gave to a
%  different source.  Those three mean very different things — only the first
%  is a detection miss.

unmatchedReport = analyzeUnmatched(movCleanF, F0img, stdExtract, stdGIAnT, realMatches, cfg);
printMissSummary(unmatchedReport, cfg);

makeUnmatchedSummaryFigure(movCleanF, stdExtract, stdGIAnT, unmatchedReport, cfg, ...
                           fullfile(outRoot,'figures'));
writeIf(unmatchedReport.table, fullfile(outRoot,'results','unmatched_components.csv'));
save(fullfile(outRoot,'results','unmatchedComponents.mat'),'unmatchedReport','-v7.3');

%% ========================================================================
%  Section 7 : Part B — ground-truth simulation
% ========================================================================
[bgF, bgF0, bgInfo] = selectFixedBackgroundCrop(movCleanF, F0img, cfg);
bgStrata            = stratifyBackground(bgF, bgF0, cfg);
save(fullfile(outRoot,'simulation','background.mat'),'bgF','bgF0','bgInfo','bgStrata','-v7.3');

simManifest = makeSimulationManifest(cfg);
fprintf('Simulation manifest: %d movies, %d injected sources total.\n', ...
    height(simManifest), sum(simManifest.nSources));
writetable(simManifest, fullfile(outRoot,'simulation','simulation_manifest.csv'));

simRows = table();
simFail = struct('movie_id',{},'method',{},'message',{});
for i = 1:height(simManifest)
    row = simManifest(i,:);
    fprintf('\n=== Simulated movie %d/%d  [%s]  amp=%.3f  sigma=%.2f um  n=%d  seed=%d ===\n', ...
        i, height(simManifest), row.dataset, row.amplitude_dFF, row.sigma_um, ...
        row.nSources, row.seed);

    [movSim, truth] = injectSyntheticSignals(bgF, bgF0, bgStrata, row, cfg);

    [outE, failE] = tryRun(@() runExtractGluSNFR3(movSim, cfg.extractParams, cfg, ...
                        sprintf('sim%04d_extract', row.movie_id)), size(movSim));
    [outG, failG] = tryRun(@() runGIAnT(movSim, cfg.giantParams, cfg, ...
                        sprintf('sim%04d_giant', row.movie_id)), size(movSim));

    if ~isempty(failE)
        simFail(end+1) = struct('movie_id',row.movie_id,'method','extractGluSNFR3','message',failE); %#ok<SAGROW>
    end
    if ~isempty(failG)
        simFail(end+1) = struct('movie_id',row.movie_id,'method','GIAnT','message',failG); %#ok<SAGROW>
    end

    simRows = [simRows; scoreAgainstTruth(outE, truth, row, cfg, 'extractGluSNFR3'); ...
                        scoreAgainstTruth(outG, truth, row, cfg, 'GIAnT')]; %#ok<AGROW>

    % Raw outputs are preserved, not just the metrics (runbook section 11)
    save(fullfile(outRoot,'simulation',sprintf('sim%04d_raw.mat', row.movie_id)), ...
        'row','truth','outE','outG','failE','failG','-v7.3');
    writetable(simRows, fullfile(outRoot,'results','simulation_rows.csv'));
end
save(fullfile(outRoot,'results','simulationRows.mat'),'simRows','simManifest','simFail','-v7.3');

%% ========================================================================
%  Section 8 : Part C — end-to-end (inject before artifact regression)
% ========================================================================
%  A prespecified subset of Part B conditions is injected into the RAW
%  fluorescence window, then pushed through regressArtifactWrapper and both
%  extractors.  Reported separately from extractor performance: it measures how
%  much true glutamate signal the artifact regression itself removes.
endToEnd = runEndToEndSubset(meta, win, bgInfo, bgStrata, cfg);
save(fullfile(outRoot,'results','endToEnd.mat'),'endToEnd','-v7.3');

%% ========================================================================
%  Section 9 : Summarize, figure, save
%% ========================================================================
results = summarizeBenchmark(realMetrics, simRows, endToEnd, simFail, cfg);
makeSimulationFigures(simRows, results, cfg, fullfile(outRoot,'figures'));
saveBenchmark(results, resExtract, resGIAnT, artifactModel, artifactQC, cfg);
printDecisionSummary(results);


%% ========================================================================
%  LOCAL FUNCTIONS
%% ========================================================================

function cfg = makeBenchmarkConfig()
% Every constant lives here.  Anything that could bias the comparison is frozen
% before the first run.

cfg = struct;
cfg.seed      = 1;
cfg.outParent = '/Volumes/BHL18TB_D2/GIAnT_benchmark';
if ~exist(cfg.outParent,'dir'), cfg.outParent = pwd; end

% ---- Real-movie selection (FROZEN before any extraction) ----------------
%  Selection rule, written down before any extraction was run: the FIRST
%  recording in the sheet whose Quality column reads 'Good' and that has
%  glutamate chunks on disk.  Today that is index 14 (20260712/160418BHLm501_VR1);
%  the index is pinned so the benchmark is reproducible if the sheet changes,
%  and readGluMeta asserts the quality still reads 'Good'.
cfg.xlsPath        = '/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/Prism_V2+Glu_Data_Arrangement.xlsx';
cfg.xlsRange       = 'C5:AA45';   % col 1 = Data Location, col 2 = Quality
cfg.realFileIdx    = 15;
cfg.requireQuality = 'Good';
cfg.realTile       = 1;      % mc2<tile>.bin
cfg.realStartFrame = 1501;   % skip the 488 nm modulation onset + initial bleach
cfg.nRealFrames    = 1000;
cfg.contextWindows = 1.0;    % frames of context each side = this x the 5 s window

% Spatial crop applied BEFORE regression, so the regression, both extractors
% and every metric see exactly the same field.  [] = the full 644x480 field,
% which the runbook prefers but which costs roughly 5 full copies of the read
% window: at ~330 Hz that is >30 GB, and extractGluSNFR3's per-pixel findpeaks
% loop over 309k pixels runs for tens of minutes.  The default below keeps the
% field comfortably larger than the 200x200 simulation crop while fitting in a
% normal workstation.  Format: [y0 y1 x0 x1].
cfg.realCrop    = [100 256 550 750];
cfg.memWarn_GB  = 24;        % warn if the planned working set exceeds this

% ---- Acquisition metadata ----------------------------------------------
cfg.frameRate_Hz = NaN;      % filled from Device_Data{4}.exposuretime
cfg.exposuretime = NaN;
cfg.pixelSize_um = NaN;      % <<<< MUST BE SET BY HAND (see header note)

% ---- Artifact regression (master Section 2) -----------------------------
cfg.nEig         = 40;       % get_eigvector components computed
cfg.regressPCs   = [1 2 3];      % PCs regressed out
cfg.bleachTau    = [10000 100];
cfg.slowWindow_s = 5;        % movmedian window, seconds
cfg.smoothForPCA = 1;        % imgaussfilt sigma before PCA

% Pick the PCs by eye, the way the master script does (it prompts with
% "PCs to regress out" after showing imshow2_patch(u_sub)).
%
% The prompt fires ONCE, on the Part A movie.  Section 2 then writes the answer
% back into cfg.regressPCs and clears this flag, so every later regression —
% the six Part C re-runs in particular — reuses exactly the same component set.
% Letting it re-prompt per movie would put a different artifact model in front
% of each extractor run and make the comparison meaningless.
%
% Leave this false for a reproducible batch run; the chosen set is recorded in
% artifactModel.regressPCs and in results/config_and_provenance.mat either way.
cfg.regressPCsInteractive = true;

% ---- Simulation grid ----------------------------------------------------
cfg.simSize         = [200 200 1000];
cfg.amplitudes_dFF  = [0.005 0.01 0.02 0.05 0.10];
cfg.sigma_um        = [0.8 1.5 2.5];
cfg.tauRise_s       = [0.01 0.03 0.06];
cfg.tauDecay_s      = [0.05 0.15 0.40];
cfg.sourcesPerMovie = [1 5 20];
cfg.nRepeats        = 20;
cfg.nEmptyControls  = 10;    % no-injection movies -> false-positive rate
cfg.eventRate_Hz    = 0.5;   % mean event rate per synthetic source
cfg.pilotOnly       = true;  % reduced repeats for the pilot (runbook section 8);
                             % set false once the pilot grid is frozen
cfg.pilotRepeats    = 3;

% ---- Matching rules (FROZEN before examining any comparison) ------------
cfg.maxCentroidDistance_um = 3;
cfg.minFootprintIoU        = 0.10;
cfg.footprintThreshFrac    = 0.20;   % binary support = fp >= 0.2*max(fp)
cfg.minComponentArea_um2   = 0.30;   % biologically equivalent limits, both methods
cfg.maxComponentArea_um2   = 20;
cfg.boundaryMargin_um      = 1.5;

% ---- Temporal scoring (used ONLY for event metrics) ---------------------
cfg.eventZ             = 4;      % threshold in robust SD of the trace
cfg.eventTol_s         = 0.03;   % timing tolerance against ground truth
cfg.eventRefrac_s      = 0.02;
cfg.maxEventsPerSource = 10;
% Window used ONLY to measure an event's amplitude once its frame is found on
% the raw trace (a centred mean, which is far less biased than a max).  Event
% DETECTION never smooths — see detectEvents.
cfg.eventAmpWindow_s   = 0.02;

% ---- Unmatched-component adjudication (Section 6) -----------------------
% How many events the model-free common-mask trace must show before an
% unmatched component counts as a real source the other method missed.
% Frozen with the other thresholds: raising it after seeing the result would
% let the miss count be tuned.
cfg.unmatchedMinEvents      = 2;
cfg.unmatchedExamplesPerSide = 5;     % examples shown per direction in fig06

% ---- extractGluSNFR3 parameters (defaults; UI disabled) -----------------
%  Everything except the doPlot/ui* flags is the function's own default, so
%  this is the published/default operating point (runbook section 5).
ep = struct();
ep.smoothSigma   = 1;
ep.zSeedMin      = 4;
ep.zThrFoot      = 2.5;
ep.minPeakDist_s = 0.2;
ep.minPeakProm   = 0.75;
ep.minArea       = 2;
ep.eventMaxArea  = 30;
ep.nmfReps       = 5;
ep.r_min         = 3;
ep.r_max         = 3.5;
ep.bound         = 10;
ep.minEventsSyn  = 1;
ep.footSigma     = 0.75;
ep.usePeakWeight = true;
ep.doPlot        = false;    % the benchmark must be non-interactive
ep.uiSeedThresh  = false;
ep.uiClusterSize = false;
ep.uiSynFilter   = false;
ep.plotCentroids = false;
ep.showNMFDiag   = false;
cfg.extractParams = ep;

% ---- GIAnT / SILo parameters (setParams('SILo') defaults; GUI disabled) --
gp = struct();
gp.isSLAP2                = false;
gp.sigma_px               = 1;
gp.nmfIter                = 2;
gp.dXY                    = 3;
gp.photonScale            = [];    % estimated independently per movie
gp.lambda                 = 0.5;
gp.phi                    = 0.1;
gp.denoiseWindow_s        = 0.2;
gp.baselineWindow_Glu_s   = 4;
gp.baselineWindow_Ca_s    = 4;
gp.activityChannel        = 1;
gp.tau_s                  = 0.03;
gp.tau2_s                 = 0.15;
gp.VIF                    = 1;
gp.peakth                 = 3;
gp.minPeakDistance        = 1;
gp.nWorkers               = 1;     % single trial: a process pool is pure overhead
gp.drawUserRois           = false; % no GUI
gp.motionThresh           = 2;
gp.nanThresh              = 0.33;
gp.discardInitial_s       = 0;
gp.includeIntegrationROIs = false;
cfg.giantParams = gp;

cfg.giantKeepWorkspace = false;    % delete the movie copy after reading results
end

% -------------------------------------------------------------------------
function cfg = derivePixelUnits(cfg)
p = cfg.pixelSize_um;
cfg.maxCentroidDistance_px = cfg.maxCentroidDistance_um / p;
cfg.minComponentArea_px    = cfg.minComponentArea_um2  / p^2;
cfg.maxComponentArea_px    = cfg.maxComponentArea_um2  / p^2;
cfg.sigma_px_grid          = cfg.sigma_um / p;

% Boundary margin.  extractGluSNFR3 ALREADY drops every synapse whose
% footprint falls outside a border of opts.bound pixels (its Step 5), while
% GIAnT applies no edge exclusion at all and is prone to activity-image edge
% artifacts.  If the shared margin were looser than opts.bound the comparison
% would silently penalise GIAnT with edge false positives that extractGluSNFR3
% never had the chance to report.  Take the stricter of the two so both methods
% are judged on the same interior field.
marginFromUm = cfg.boundaryMargin_um / p;
cfg.boundaryMargin_px = max(marginFromUm, cfg.extractParams.bound);
cfg.boundaryMargin_um_effective = cfg.boundaryMargin_px * p;
if cfg.boundaryMargin_px > marginFromUm
    fprintf(['Boundary margin raised from %.1f px (%.1f um) to %.1f px to match ' ...
             'extractGluSNFR3''s internal opts.bound = %d px.\n'], ...
        marginFromUm, cfg.boundaryMargin_um, cfg.boundaryMargin_px, cfg.extractParams.bound);
end
end

% -------------------------------------------------------------------------
function prov = collectProvenance(cfg, labRoot, giantRoot)
prov = struct();
prov.matlabVersion = version;
prov.toolboxes     = ver;
prov.timestamp     = datestr(now,'yyyy-mm-dd HH:MM:SS');
prov.cfg           = cfg;
prov.gitLab        = gitHash(labRoot);
prov.gitGIAnT      = gitHash(giantRoot);
prov.rngState      = rng;
try
    prov.hostname = char(java.net.InetAddress.getLocalHost.getHostName);
catch
    prov.hostname = 'unknown';
end
end

function h = gitHash(dr)
h = 'unavailable';
try
    [st, out] = system(sprintf('git -C "%s" rev-parse HEAD', dr));
    if st == 0
        h = strtrim(out);
        [~, dirtyOut] = system(sprintf('git -C "%s" status --porcelain', dr));
        if ~isempty(strtrim(dirtyOut)), h = [h ' (dirty)']; end
    end
catch
end
end

% =========================================================================
%  MOVIE LOADING
% =========================================================================

function meta = readGluMeta(cfg)
%READGLUMETA  Geometry, frame rate and file paths of the glutamate camera.
[~,~,raw] = xlsread(cfg.xlsPath, 'Sheet1', cfg.xlsRange);
fpath   = raw(:,1)';
quality = raw(:,2)';
assert(cfg.realFileIdx <= numel(fpath), 'realFileIdx %d exceeds the sheet.', cfg.realFileIdx);
meta.fpath   = fpath{cfg.realFileIdx};
meta.quality = quality{cfg.realFileIdx};
if ~isempty(cfg.requireQuality)
    assert(ischar(meta.quality) && strcmpi(strtrim(meta.quality), cfg.requireQuality), ...
        ['Recording %d is graded "%s", not "%s".  The frozen selection rule is the ' ...
         'first Quality=="Good" recording with glutamate chunks; re-check the sheet ' ...
         'before changing cfg.realFileIdx.'], cfg.realFileIdx, ...
        string(meta.quality), cfg.requireQuality);
end

S = load(fullfile(meta.fpath,'output_data.mat'));    % Device_Data
Device_Data = S.Device_Data;

meta.nCol         = double(Device_Data{4}.ROI(2));   % blue (glutamate) camera
meta.nRow         = double(Device_Data{4}.ROI(4));
meta.exposuretime = Device_Data{4}.exposuretime;
meta.frameRate_Hz = 1 / meta.exposuretime;
meta.DAQ_rate     = Device_Data{1,2}.Counter_Inputs(1,1).rate;

meta.binFile = fullfile(meta.fpath, sprintf('mc2%02d.bin', cfg.realTile));
meta.mcFile  = fullfile(meta.fpath, sprintf('mc2Trace%02d.mat', cfg.realTile));
assert(exist(meta.binFile,'file')==2, 'Missing glutamate chunk: %s', meta.binFile);
assert(exist(meta.mcFile, 'file')==2, 'Missing motion trace: %s',    meta.mcFile);

d = dir(meta.binFile);
meta.nFramesTile = floor(d.bytes / (meta.nRow*meta.nCol*2));
end

function win = planReadWindow(meta, cfg)
%PLANREADWINDOW  Benchmark frames plus enough context for the bleach fit and
%the 5 s movmedian, clipped to the tile.
lwWin = max(3, round(cfg.slowWindow_s / meta.exposuretime));
ctx   = round(cfg.contextWindows * lwWin);

benchFirst = cfg.realStartFrame;
benchLast  = cfg.realStartFrame + cfg.nRealFrames - 1;
assert(benchLast <= meta.nFramesTile, ...
    'Tile has %d frames; benchmark window ends at %d.', meta.nFramesTile, benchLast);

win.first = max(1, benchFirst - ctx);
win.last  = min(meta.nFramesTile, benchLast + ctx);
win.frames  = win.first:win.last;
win.nFrames = numel(win.frames);
win.benchFirst = benchFirst;
win.benchLast  = benchLast;
win.benchIdx   = (benchFirst - win.first + 1) : (benchLast - win.first + 1);
win.lwWin      = lwWin;
win.contextFrames = ctx;

if isempty(cfg.realCrop)
    H = meta.nRow;  W = meta.nCol;
else
    H = cfg.realCrop(2)-cfg.realCrop(1)+1;
    W = cfg.realCrop(4)-cfg.realCrop(3)+1;
end
win.H = H;  win.W = W;
% Peak working set: the read window, filt, filt_lw, sub, and the vectorised
% copy get_eigvector builds internally — roughly five double copies.
win.estGB = H*W*win.nFrames*8*5 / 1e9;
if win.estGB > cfg.memWarn_GB
    warning('Benchmark:Memory', ...
        ['The planned artifact regression needs roughly %.0f GB (%dx%d x %d frames). ' ...
         'Shrink cfg.realCrop or cfg.contextWindows if this machine has less.'], ...
        win.estGB, H, W, win.nFrames);
end
end

function [mov, meta] = readGluFrames(meta, win, cfg)
%READGLUFRAMES  Read the planned frames (and optional spatial crop) as double.
mov = double(readBinMov_times(meta.binFile, meta.nRow, meta.nCol, win.frames));
if ~isempty(cfg.realCrop)
    r = cfg.realCrop;
    mov = mov(r(1):r(2), r(3):r(4), :);
end
M = load(meta.mcFile);
meta.mcTrace = M.mcTrace2.xymean(win.frames, :);
meta.readWindow = win;
end

% =========================================================================
%  ARTIFACT REGRESSION
% =========================================================================

function [out, model, qc] = regressArtifactWrapper(movIn, meta, cfg)
%REGRESSARTIFACTWRAPPER  The master pipeline's artifact regression, frozen.
%
%  Operations and their order are taken verbatim from
%  Analysis_20260713_SD_VoltronST_plus_Glu.m Section 2; the only change is that
%  the number of PCs to regress is cfg.regressPCs instead of a keyboard prompt,
%  so the operation is deterministic and identical for both methods.
%
%  out.F       artifact-corrected movie in CAMERA COUNTS (mean preserved)
%  out.sub     high-pass residual after PC + motion regression
%  out.filt    bleach-normalised movie  (mov ./ y_fit)
%  out.filt_lw slow baseline of filt    (5 s movmedian)
%  out.y_fit   bleach curve
%
%  Mean fluorescence is PRESERVED: out.F = (sub + filt_lw) .* y_fit and sub is a
%  residual about zero, so the temporal mean of out.F tracks that of movIn up to
%  the regressed components.  Nothing downstream re-centres the movie.
%
%  Leakage: the PCA is refit inside this call, on this movie only.  No ground
%  truth and no extraction result is used to choose the component count.

T  = size(movIn,3);
tS = tic;

meanF = squeeze(mean(movIn, [1 2]));
y_fit = expfitDM_2((1:T)', meanF(:), (1:T)', cfg.bleachTau);
y_fit = y_fit(:);

filt    = movIn ./ reshape(y_fit,1,1,[]);
lwWin   = max(3, round(cfg.slowWindow_s / meta.exposuretime));
filt_lw = movmedian(filt, lwWin, 3);
sub     = filt - filt_lw;

[V, D, u_sub] = get_eigvector(tovec(imgaussfilt(sub, cfg.smoothForPCA))', cfg.nEig);
u_img = reshape(u_sub, size(movIn,1), size(movIn,2), []);

pcIdx = cfg.regressPCs(cfg.regressPCs <= size(V,2));
answered = false;
if isfield(cfg,'regressPCsInteractive') && cfg.regressPCsInteractive
    [pcIdx, answered] = promptForPCs(u_img, V, D, pcIdx, meta);
end
sub   = SeeResiduals(sub, V(:,pcIdx));

mcReg = meta.mcTrace(1:T,:);
sub   = SeeResiduals(sub, mcReg);

out.F       = (sub + filt_lw) .* reshape(y_fit,1,1,[]);
out.sub     = sub;
out.filt    = filt;
out.filt_lw = filt_lw;
out.y_fit   = y_fit;

model = struct();
model.method            = 'master Analysis_20260713 Section 2 (bleach / 5 s movmedian / PCA / motion)';
model.bleachTau         = cfg.bleachTau;
model.y_fit             = y_fit;
model.slowWindow_frames = lwWin;
model.smoothForPCA      = cfg.smoothForPCA;
model.nEig              = cfg.nEig;
model.regressPCs        = pcIdx;
model.pcTraces          = V(:,1:min(cfg.nEig,size(V,2)));
model.pcEigenvalues     = D;
model.pcImages          = u_img;
model.motionRegressors  = mcReg;
model.meanPreserved     = true;
model.pcsChosenInteractively = answered;   % true only if a human actually answered
model.runtime_s         = toc(tS);

qc = struct();
qc.varBefore      = var(filt, 0, 3);
qc.varAfter       = var(sub + filt_lw, 0, 3);
qc.globalRaw      = meanF(:);
qc.globalFit      = y_fit;
qc.regressors     = [V(:,pcIdx) mcReg];
qc.regressorNames = [arrayfun(@(k) sprintf('PC%d',k), pcIdx, 'UniformOutput', false), {'mcX','mcY'}];
qc.frameRate_Hz   = meta.frameRate_Hz;
end

% -------------------------------------------------------------------------
function [pcIdx, answered] = promptForPCs(u_img, V, D, pcIdxDefault, meta)
%PROMPTFORPCS  Show the PC spatial maps and ask which to regress out.
%  Same judgement call the master script's `input("PCs to regress out")` asks
%  for, with the PC temporal traces and the eigenvalue spectrum added because
%  an artifact component is usually obvious in its trace, not just its map.
%
%  Returns pcIdxDefault unchanged if the session cannot prompt (-batch, -nodisplay,
%  a parallel worker), so a headless run degrades to the frozen value instead of
%  hanging on a blocked input().

pcIdx    = pcIdxDefault;
answered = false;
if ~usejava('desktop') || ~feature('ShowFigureWindows')
    warning('Benchmark:NoInteractivePrompt', ...
        ['cfg.regressPCsInteractive is true but this session cannot prompt ' ...
         '(batch/headless).  Using the frozen cfg.regressPCs = %s.'], mat2str(pcIdx));
    return
end

nShow = min(size(V,2), 20);

fh = figure('Name','Artifact PCs — choose components to regress out', ...
            'Position',[40 40 1500 850]);
tiledlayout(fh, 3, 1, 'TileSpacing','compact','Padding','compact');

% (1) spatial maps, as the master script shows them
nexttile;
try
    imshow2_patch(u_img(:,:,1:nShow));
catch
    montage(rescale(u_img(:,:,1:nShow)), 'Size', [ceil(nShow/5) 5]);
end
title(sprintf('PC spatial maps 1-%d   (%s)', nShow, meta.fpath), 'Interpreter','none');

% (2) temporal traces, offset
nexttile;
tv = (0:size(V,1)-1) / meta.frameRate_Hz;
plot(tv, rescale2(V(:,1:nShow)', 2)' + (1:nShow));
set(gca,'YTick',1:nShow,'YTickLabel',arrayfun(@(k) sprintf('PC%d',k), 1:nShow, 'UniformOutput',false));
xlabel('Time (s)'); ylabel('component'); axis tight;
title('PC temporal traces (these are the regressors)');

% (3) eigenvalue spectrum
nexttile;
semilogy(D(1:nShow) / sum(D), 'o-'); hold on;
plot(pcIdxDefault, D(pcIdxDefault)/sum(D), 'ro', 'MarkerSize',10, 'LineWidth',1.5);
xlabel('component'); ylabel('variance fraction'); grid on;
title('Eigenvalue spectrum (red = current default)');
drawnow;

fprintf('\n------------------------------------------------------------\n');
fprintf('Choose the PCs to regress out of the glutamate movie.\n');
fprintf('  default : %s\n', mat2str(pcIdxDefault));
fprintf('  examples: 1:3   [1 2 5]   1   []  (regress none)\n');
fprintf('  Enter alone keeps the default.  %d components available.\n', size(V,2));
fprintf('------------------------------------------------------------\n');

while true
    raw = input('PCs to regress out: ', 's');
    if isempty(strtrim(raw))
        pcIdx = pcIdxDefault;
        break
    end
    try
        % eval so that range syntax like 1:3 works, as in the master script
        val = eval(['[' raw ']']);
    catch
        fprintf(2,'  Could not parse "%s".  Try again.\n', raw);
        continue
    end
    if isempty(val)
        pcIdx = [];
        break
    end
    val = unique(round(val(:)'));
    if any(val < 1) || any(val > size(V,2))
        fprintf(2,'  Components must be between 1 and %d.  Try again.\n', size(V,2));
        continue
    end
    pcIdx = val;
    break
end

answered = true;
if isempty(pcIdx)
    fprintf('Regressing out NO principal components.\n');
else
    fprintf('Regressing out PCs %s (%.1f%% of the residual variance).\n', ...
        mat2str(pcIdx), 100*sum(D(pcIdx))/sum(D));
end
close(fh);
end

% -------------------------------------------------------------------------
function qc = augmentArtifactQC(qc, movCleanF, movSub, cfg)
%  Residual-vs-regressor correlations and power spectra in vessel-rich and
%  vessel-poor ROIs (runbook section 3).
[H,W,T] = size(movSub);
sub2 = reshape(double(movSub), H*W, T);
gTrace = mean(sub2, 1)';

nReg = size(qc.regressors,2);
qc.residualRegressorCorr = nan(1,nReg);
for k = 1:nReg
    r = qc.regressors(1:T,k);
    c = corrcoef(gTrace, r);
    qc.residualRegressorCorr(k) = c(1,2);
end

% Vessel-rich / vessel-poor proxy: the darkest and brightest deciles of the
% mean image (vessels are the dark, high-variance structures in this prep).
mImg = mean(double(movCleanF),3);
maskVessel = mImg <= prctile(mImg(:), 10);
maskTissue = mImg >= prctile(mImg(:), 90);

trV = mean(sub2(maskVessel(:),:), 1)';
trT = mean(sub2(maskTissue(:),:), 1)';
[qc.psdVessel, qc.psdF] = simplePSD(trV, qc.frameRate_Hz);
qc.psdTissue            = simplePSD(trT, qc.frameRate_Hz);
qc.maskVessel = maskVessel;
qc.maskTissue = maskTissue;
qc.cfgSnapshot = cfg;
end

function [p, f] = simplePSD(x, fs)
x = double(x(:));  x = x - mean(x);
n = numel(x);
X = fft(x .* hann(n));
p = abs(X(1:floor(n/2)+1)).^2 / (fs*n);
f = (0:floor(n/2))' * fs / n;
end

% =========================================================================
%  ADAPTERS — common output contract
%    result.name, .footprints (Y x X x K), .traces (K x T), .tracesDFF (K x T),
%    .params, .runtime_s, .extra, .failed, .message
% =========================================================================

function result = runExtractGluSNFR3(movCleanF, params, cfg, tag)
%RUNEXTRACTGLUSNFR3  Adapter.
%  extractGluSNFR3 takes two movies: the high-pass residual it detects on, and
%  the bleach-normalised fluorescence it computes F0 from.  Both are derived
%  here from the common input with the SAME operations the master pipeline uses,
%  and the derivation is recorded in result.extra.preproc.

tStart = tic;
mov = double(movCleanF);
T   = size(mov,3);

% --- documented internal preprocessing -----------------------------------
meanF   = squeeze(mean(mov,[1 2]));
y_fit   = expfitDM_2((1:T)', meanF(:), (1:T)', cfg.bleachTau);
filt    = mov ./ reshape(y_fit(:),1,1,[]);
lwWin   = max(3, round(cfg.slowWindow_s * cfg.frameRate_Hz));
filt_lw = movmedian(filt, lwWin, 3);
sub     = filt - filt_lw;

opts = params;
opts.diagSaveDir = fullfile(cfg.outRoot,'figures',['extract_' tag]);
opts.diagPrefix  = tag;

R = extractGluSNFR3(sub, filt, cfg.exposuretime, opts);

result = struct();
result.name       = 'extractGluSNFR3';
result.footprints = double(R.S_glu);
result.traces     = double(R.T_glu);
if isfield(R,'dFF_glu') && ~isempty(R.dFF_glu)
    result.tracesDFF = double(R.dFF_glu);
else
    result.tracesDFF = [];
end
result.params    = R.opts;
result.runtime_s = toc(tStart);
result.extra     = R;
result.extra.preproc = struct( ...
    'description', ['bleach fit expfitDM_2([10000 100]) -> divide; 5 s movmedian ' ...
                    'slow baseline -> subtract; sub = detection movie, filt = F0 movie'], ...
    'y_fit', y_fit(:), 'slowWindow_frames', lwWin);
result.failed  = false;
result.message = '';
end

% -------------------------------------------------------------------------
function result = runGIAnT(movCleanF, params, cfg, tag)
%RUNGIANT  Adapter around SILo.
%  SILo is file-based, so the movie is written into a minimal GIAnT workspace
%  (registered .h5 + _ALIGNMENTDATA.h5 + trial_table.h5), SILo is run on it, and
%  experiment_summary.h5 is read back into the common output contract.

tStart  = tic;
runDir  = fullfile(cfg.outRoot,'giant_runs',tag);
mocoDir = fullfile(runDir,'motion_correction');
if exist(runDir,'dir'), rmdir(runDir,'s'); end
mkdir(mocoDir);

stem      = 'bench_trial01';
regFn     = [stem '_REGISTERED.h5'];
adataFn   = [stem '_ALIGNMENTDATA.h5'];
regPath   = fullfile(mocoDir, regFn);
adataPath = fullfile(mocoDir, adataFn);

mov = single(movCleanF);
[H,W,T] = size(mov);

% GIAnT subtracts a 5th-percentile baseline and estimates a photon scale from
% the variance-mean relationship, so it needs non-negative counts.  The
% artifact-corrected movie is already in counts; clip the few pixels the
% regression can push slightly negative and record how many.
negFrac = mean(mov(:) < 0);
mov = max(mov, 0);

% ---- registered movie (single root dataset, as loadAndProcessTrialAsync expects)
h5create(regPath, '/mov', [H W T], 'Datatype','single', ...
         'ChunkSize',[H W min(T,64)], 'Deflate',1);
h5write(regPath, '/mov', mov);

% ---- alignment data ------------------------------------------------------
aData = struct();
aData.numChannels = 1;
aData.frametime   = 1 / cfg.frameRate_Hz;    % SILo: analyzeHz = 1/frametime
aData.alignHz     = cfg.frameRate_Hz;
aData.motionDSc   = zeros(1,T);
aData.motionDSr   = zeros(1,T);
aData.motionC     = zeros(1,T);              % movie is already motion-corrected
aData.motionR     = zeros(1,T);
aData.recNegErr   = zeros(1,T);              % -> no motion censoring
aData.meanIM      = reshape(single(mean(mov,3)), 1, H, W);   % channels x rows x cols
saveStructToH5(aData, adataPath);

% ---- trial table ---------------------------------------------------------
% datadr points at motion_correction and filename at the registered movie, so
% verifyFiles' "source data file exists" check passes without a dummy file; for
% non-SLAP2 data that source file is never read (fn_raw is used instead).
tt = struct();
tt.datadr        = mocoDir;
tt.savedr        = runDir;
tt.filename      = {regFn};
tt.true_trial_ix = 1;
tt.epoch         = 1;
tt.motion_correction = struct( ...
    'fn_reg_ds', {{regFn}}, ...
    'fn_adata',  {{adataFn}}, ...
    'fn_raw',    {{regFn}}, ...
    'registration_failed', uint8(0));
ttPath = fullfile(runDir,'trial_table.h5');
saveStructToH5(tt, ttPath);

% ---- run SILo ------------------------------------------------------------
result = struct();
result.name    = 'GIAnT';
result.failed  = false;
result.message = '';
usedParams     = params;
try
    usedParams = SILo(ttPath, params);
    G = readGiantSummary(fullfile(runDir,'source_extraction','experiment_summary.h5'), usedParams);
    result.footprints = G.footprints;
    result.traces     = G.dF_ls;
    result.tracesDFF  = G.dFF;
    result.extra      = G;
catch ME
    result.failed     = true;
    result.message    = ME.message;
    result.footprints = zeros(H,W,0);
    result.traces     = zeros(0,T);
    result.tracesDFF  = zeros(0,T);
    result.extra      = struct('error', getReport(ME,'extended','hyperlinks','off'));
    fprintf(2,'GIAnT run "%s" FAILED: %s\n', tag, ME.message);
end
result.params    = usedParams;
result.runtime_s = toc(tStart);
result.extra.preproc = struct( ...
    'description', ['movie written as photon-count registered H5; SILo applies its own ' ...
                    '5th-pct baseline subtraction, movmedian F0, photon scaling, ' ...
                    'matched filter and DoG internally'], ...
    'negativePixelFractionClipped', negFrac, ...
    'runDir', runDir);

if ~cfg.giantKeepWorkspace && ~result.failed
    try, delete(regPath); catch, end   % the movie copy dominates the run dir
end
end

% -------------------------------------------------------------------------
function G = readGiantSummary(h5file, params)
%READGIANTSUMMARY  experiment_summary.h5 -> common contract.
%  Dataset shapes are MATLAB size() (GIAnT writes row_major = 0):
%    profiles   K x 1 x Y x X
%    coords     K x 3, 0-indexed [z y x]
%    dF_ls etc. K x nChan x T

assert(exist(h5file,'file')==2, 'GIAnT wrote no experiment_summary.h5 (%s)', h5file);
p = '/Path1';

prof = h5read(h5file, [p '/sources/spatial/profiles']);      % K x 1 x Y x X
K = size(prof,1);  Y = size(prof,3);  X = size(prof,4);
G.footprints = zeros(Y,X,K);
for k = 1:K
    G.footprints(:,:,k) = squeeze(prof(k,1,:,:));
end
G.coords = h5read(h5file, [p '/sources/spatial/coords']) + 1;   % -> 1-indexed [z y x]

ch = params.activityChannel;
G.dF_ls       = pickChannel(h5read(h5file, [p '/sources/temporal/dF_ls']),       ch, K);
G.dF_denoised = pickChannel(h5read(h5file, [p '/sources/temporal/dF_denoised']), ch, K);
G.events      = pickChannel(h5read(h5file, [p '/sources/temporal/events']),      ch, K);
G.F0          = pickChannel(h5read(h5file, [p '/sources/temporal/F0']),          ch, K);
try
    G.SNR = h5read(h5file, [p '/sources/temporal/SNR']);
catch
    G.SNR = nan(K,1);
end

F0safe = G.F0;
F0safe(~isfinite(F0safe) | F0safe <= 0) = NaN;
G.dFF = G.dF_ls ./ F0safe;

G.actIM  = h5read(h5file, [p '/visualizations/act_im']);
G.meanIM = h5read(h5file, [p '/visualizations/mean_im']);
G.h5file = h5file;
end

function M = pickChannel(A, ch, K)
%  A is K x nChan x T -> K x T for one channel, without chaining an index onto
%  a function call (illegal in MATLAB).
M = reshape(A(:, ch, :), K, []);
end

% =========================================================================
%  FOOTPRINT STANDARDIZATION AND MATCHING
% =========================================================================

function S = standardizeFootprints(result, cfg)
%STANDARDIZEFOOTPRINTS  Runbook 6.1.  Applied identically to both methods.
%   1. negative weights -> 0 (neither method defines a negative weight)
%   2. peak-normalised footprint for display
%   3. binary support at cfg.footprintThreshFrac * max, largest connected part
%   4. area limits in um^2 and a boundary margin, equal for both methods

fp = result.footprints;
fp(~isfinite(fp)) = 0;
fp = max(fp, 0);
[H,W,K] = size(fp);

bw     = false(H,W,K);
normFp = zeros(H,W,K);
cent   = nan(K,2);      % [x y], 1-indexed pixels
area   = zeros(K,1);
[rr,cc] = ndgrid(1:H, 1:W);
for k = 1:K
    f = fp(:,:,k);
    m = max(f(:));
    if m <= 0, continue; end
    normFp(:,:,k) = f / m;
    b = f >= cfg.footprintThreshFrac * m;
    CC = bwconncomp(b, 8);
    if CC.NumObjects > 1                       % keep the part holding the maximum
        [~, imax] = max(f(:));
        for c = 1:CC.NumObjects
            if any(CC.PixelIdxList{c} == imax)
                b = false(H,W);  b(CC.PixelIdxList{c}) = true;
                break
            end
        end
    end
    bw(:,:,k) = b;
    area(k)   = nnz(b);
    w = f .* b;
    if sum(w(:)) > 0
        cent(k,:) = [sum(w(:).*cc(:)) sum(w(:).*rr(:))] / sum(w(:));
    end
end

bm      = cfg.boundaryMargin_px;
tooSmall = area < cfg.minComponentArea_px;
tooLarge = area > cfg.maxComponentArea_px;
inBounds = all(isfinite(cent),2) & ...
           cent(:,1) > bm & cent(:,1) < W - bm & ...
           cent(:,2) > bm & cent(:,2) < H - bm;
keep = ~tooSmall & ~tooLarge & inBounds;

S = struct();
S.name         = result.name;
S.keepIdx      = find(keep);
S.footprints   = fp(:,:,keep);
S.normFp       = normFp(:,:,keep);
S.bw           = bw(:,:,keep);
S.centroid_xy  = cent(keep,:);
S.area_px      = area(keep);
S.area_um2     = area(keep) * cfg.pixelSize_um^2;
S.traces       = pickRows(result.traces,    keep);
S.tracesDFF    = pickRows(result.tracesDFF, keep);
S.nBefore      = K;
S.nAfter       = nnz(keep);
S.runtime_s    = result.runtime_s;
S.failed       = result.failed;
S.message      = result.message;
S.rejectReason = struct('tooSmall', nnz(tooSmall), 'tooLarge', nnz(tooLarge), ...
                        'atBoundary', nnz(~inBounds));
S.imSize = [H W];
end

function out = pickRows(M, keep)
if isempty(M), out = []; else, out = M(keep,:); end
end

% -------------------------------------------------------------------------
function M = matchComponents(A, B, cfg)
%MATCHCOMPONENTS  One-to-one Hungarian assignment between two component sets.
%  cost = 0.5*(d/dmax) + 0.5*(1 - IoU).  Pairs failing the preregistered rule
%  (d <= dmax AND IoU >= minIoU) are made unselectable, so a match is never
%  forced onto a pair that does not qualify.

nA = size(A.centroid_xy,1);
nB = size(B.centroid_xy,1);
M = struct('pairs',zeros(0,2),'cost',[],'dist_um',[],'iou',[],'spatialCorr',[], ...
           'unmatchedA',(1:nA)','unmatchedB',(1:nB)','nA',nA,'nB',nB, ...
           'distMatrix',[],'iouMatrix',[],'corrMatrix',[]);
if nA == 0 || nB == 0, return; end

[D, IOU, SC] = pairwiseSpatial(A, B, cfg);

cost    = 0.5*(D/cfg.maxCentroidDistance_um) + 0.5*(1-IOU);
allowed = D <= cfg.maxCentroidDistance_um & IOU >= cfg.minFootprintIoU;
cost(~allowed) = 1e6;

pairs = matchpairs(cost, 1.0);      % a pair must beat cost 1 per unmatched side
if ~isempty(pairs)
    ok = false(size(pairs,1),1);
    for k = 1:size(pairs,1)
        ok(k) = allowed(pairs(k,1), pairs(k,2));
    end
    pairs = pairs(ok,:);
end

M.pairs      = pairs;
M.distMatrix = D;  M.iouMatrix = IOU;  M.corrMatrix = SC;
if ~isempty(pairs)
    idx = sub2ind([nA nB], pairs(:,1), pairs(:,2));
    M.cost = cost(idx);  M.dist_um = D(idx);  M.iou = IOU(idx);  M.spatialCorr = SC(idx);
end
M.unmatchedA = setdiff((1:nA)', pairs(:,1));
M.unmatchedB = setdiff((1:nB)', pairs(:,2));
end

function [D, IOU, SC] = pairwiseSpatial(A, B, cfg)
nA = size(A.centroid_xy,1);
nB = size(B.centroid_xy,1);
D   = inf(nA,nB);
IOU = zeros(nA,nB);
SC  = nan(nA,nB);
for i = 1:nA
    for j = 1:nB
        d = hypot(A.centroid_xy(i,1)-B.centroid_xy(j,1), ...
                  A.centroid_xy(i,2)-B.centroid_xy(j,2)) * cfg.pixelSize_um;
        D(i,j) = d;
        if d > 3*cfg.maxCentroidDistance_um, continue; end   % skip the pixel work
        IOU(i,j) = footprintIoU(A.bw(:,:,i), B.bw(:,:,j));
        SC(i,j)  = localSpatialCorr(A.normFp(:,:,i), B.normFp(:,:,j), ...
                                    A.bw(:,:,i) | B.bw(:,:,j));
    end
end
end

function v = footprintIoU(a, b)
u = nnz(a | b);
if u == 0, v = 0; else, v = nnz(a & b) / u; end
end

function c = localSpatialCorr(a, b, mask)
%  Correlation inside the local bounding box of the union support.
c = NaN;
if ~any(mask(:)), return; end
[rr,cc] = find(mask);
r = max(1,min(rr)-2) : min(size(a,1), max(rr)+2);
k = max(1,min(cc)-2) : min(size(a,2), max(cc)+2);
x = double(reshape(a(r,k),[],1));
y = double(reshape(b(r,k),[],1));
if std(x)==0 || std(y)==0, return; end
cm = corrcoef(x,y);  c = cm(1,2);
end

% =========================================================================
%  REAL-MOVIE COMPARISON
% =========================================================================

function R = compareRealComponents(movCleanF, F0img, A, B, M, cfg)
%COMPAREREALCOMPONENTS  Runbook 6.3: native traces AND common-mask traces.
%  The common-mask trace is extracted from the union support of the two
%  footprints with the same background subtraction and normalisation for both,
%  so it separates footprint differences from trace-estimation differences.

R = struct();
R.nA = M.nA;  R.nB = M.nB;  R.nMatched = size(M.pairs,1);
R.runtime  = struct('extractGluSNFR3', A.runtime_s, 'GIAnT', B.runtime_s);
R.area_um2 = struct('extractGluSNFR3', A.area_um2,  'GIAnT', B.area_um2);

np = R.nMatched;
rows = struct('idxA',{},'idxB',{},'dist_um',{},'iou',{},'spatialCorr',{}, ...
              'areaRatio',{},'pearsonNative',{},'spearmanNative',{},'lagNative_s',{}, ...
              'corrToCommon_extract',{},'corrToCommon_giant',{}, ...
              'ampRatio',{},'snrRatio',{},'noise_extract',{},'noise_giant',{});

R.commonTraces = cell(np,1);
R.nativeTraces = cell(np,2);
for k = 1:np
    i = M.pairs(k,1);  j = M.pairs(k,2);
    mask = A.bw(:,:,i) | B.bw(:,:,j);

    common  = commonMaskDFF(movCleanF, F0img, mask, cfg);
    nativeA = nativeTrace(A, i);
    nativeB = nativeTrace(B, j);

    [pN, sN, lagN] = traceAgreement(nativeA, nativeB, cfg);
    pCA = traceAgreement(nativeA, common, cfg);
    pCB = traceAgreement(nativeB, common, cfg);

    ampA = robustPeakAmp(nativeA, cfg);  ampB = robustPeakAmp(nativeB, cfg);
    snrA = traceSNR(nativeA, cfg);       snrB = traceSNR(nativeB, cfg);

    rows(k).idxA = i;                      rows(k).idxB = j;
    rows(k).dist_um = M.dist_um(k);        rows(k).iou = M.iou(k);
    rows(k).spatialCorr = M.spatialCorr(k);
    rows(k).areaRatio = A.area_px(i) / max(B.area_px(j), eps);
    rows(k).pearsonNative = pN;            rows(k).spearmanNative = sN;
    rows(k).lagNative_s = lagN;
    rows(k).corrToCommon_extract = pCA;    rows(k).corrToCommon_giant = pCB;
    rows(k).ampRatio = ampA / max(ampB, eps);
    rows(k).snrRatio = snrA / max(snrB, eps);
    rows(k).noise_extract = robustSD(nativeA);
    rows(k).noise_giant   = robustSD(nativeB);

    R.commonTraces{k}   = common;
    R.nativeTraces{k,1} = nativeA;
    R.nativeTraces{k,2} = nativeB;
end
if np > 0
    R.per = struct2table(rows);
else
    R.per = table();
end

R.summary = struct( ...
    'medianDist_um',        medianOr(R.per,'dist_um'), ...
    'medianIoU',            medianOr(R.per,'iou'), ...
    'medianSpatialCorr',    medianOr(R.per,'spatialCorr'), ...
    'medianPearsonNative',  medianOr(R.per,'pearsonNative'), ...
    'medianAreaRatio',      medianOr(R.per,'areaRatio'), ...
    'medianCorrToCommon_extract', medianOr(R.per,'corrToCommon_extract'), ...
    'medianCorrToCommon_giant',   medianOr(R.per,'corrToCommon_giant'), ...
    'matchRateA', R.nMatched / max(R.nA,1), ...
    'matchRateB', R.nMatched / max(R.nB,1));

% Blinded random sample of unmatched components for visual inspection
rng(cfg.seed + 7,'twister');
R.blindedUnmatched = blindSample(M);
end

function v = medianOr(T, name)
if isempty(T) || ~ismember(name, T.Properties.VariableNames)
    v = NaN;
else
    v = median(T.(name), 'omitnan');
end
end

function tr = nativeTrace(S, k)
if ~isempty(S.tracesDFF)
    tr = S.tracesDFF(k,:);
else
    tr = S.traces(k,:);
end
tr = double(tr(:))';
end

function d = commonMaskDFF(movCleanF, F0img, mask, cfg)
%COMMONMASKDFF  Mean over the mask, slow baseline removed, expressed as dF/F.
[H,W,T] = size(movCleanF);
flat = reshape(double(movCleanF), H*W, T);
x    = mean(flat(mask(:), :), 1)';
f0   = mean(double(F0img(mask)));
win  = max(3, round(cfg.slowWindow_s * cfg.frameRate_Hz));
d    = ((x - movmedian(x, win)) / max(f0, eps))';
end

function [p, s, lag_s] = traceAgreement(a, b, cfg)
a = double(a(:));  b = double(b(:));
n = min(numel(a), numel(b));
a = a(1:n);  b = b(1:n);
good = isfinite(a) & isfinite(b);
p = NaN;  s = NaN;  lag_s = NaN;
if nnz(good) < 10 || std(a(good)) == 0 || std(b(good)) == 0, return; end
p = corr(a(good), b(good), 'type','Pearson');
s = corr(a(good), b(good), 'type','Spearman');
aa = a;  bb = b;  aa(~good) = 0;  bb(~good) = 0;
maxLag = max(1, round(0.1 * cfg.frameRate_Hz));
[xc, lags] = xcorr(aa - mean(aa), bb - mean(bb), maxLag, 'coeff');
[~, im] = max(xc);
lag_s = lags(im) / cfg.frameRate_Hz;
end

function a = robustPeakAmp(tr, cfg)
[~, pk] = detectEvents(tr, cfg);
if isempty(pk), a = NaN; else, a = mean(pk); end
end

function s = traceSNR(tr, cfg)
n = robustSD(tr);
[~, pk] = detectEvents(tr, cfg);
if isempty(pk) || ~isfinite(n) || n == 0, s = NaN; else, s = mean(pk) / n; end
end

function n = robustSD(tr)
%ROBUSTSD  Baseline noise of a trace whose events are strictly positive-going.
%  A plain MAD is inflated by the events themselves: with a 150 ms decay at a
%  few hundred Hz, a handful of events already covers a third of the samples
%  and the MAD grows with them.  Estimating from the BELOW-median half only is
%  immune to that (for symmetric noise the half-MAD equals 0.6745*sigma, the
%  same constant), as long as fewer than half the samples are event samples.
tr = double(tr(:));  tr = tr(isfinite(tr));
if isempty(tr), n = NaN; return; end
m  = median(tr);
lo = m - tr(tr < m);
if isempty(lo)
    n = 1.4826 * median(abs(tr - m));   % degenerate (e.g. constant trace)
else
    n = 1.4826 * median(lo);
end
end

function [frames, peaks] = detectEvents(tr, cfg)
%DETECTEVENTS  Threshold at cfg.eventZ robust SD with a refractory period.
%  Used ONLY for temporal scoring, never for detection or spatial scoring.
%
%  Detection runs on the RAW trace, deliberately.  Smoothing first would buy
%  SNR, but it is a matched-filter-like operation that NEITHER method performs
%  in its own output, so a smoothing scorer does work on both methods' behalf
%  and can mask the very kinetic differences the benchmark exists to measure
%  (GIAnT filters with an assumed tau internally; extractGluSNFR3 does not).
%  Keeping the scorer neutral costs sensitivity on noisy traces — expect more
%  spurious peaks where per-frame SNR is near 1 — but that cost is identical
%  for both methods, so the paired comparison stays fair.
%
%  Height alone is not enough: the decay tail of a real event sits above the
%  threshold for many frames and every noise ripple on it is a local maximum.
%  Requiring the same z in PROMINENCE keeps one detection per event.
%
%  Amplitudes are measured as a short centred MEAN, not a max.  Taking max over
%  a window is badly biased upwards — measured at +34% to +93% depending on
%  amplitude and decay, worst at low amplitude — because it selects the largest
%  of several noise samples.  A centred mean measures within about +/-3%.  This
%  is amplitude MEASUREMENT, not detection: the peak locations come from the
%  raw trace above.
tr = double(tr(:));
tr(~isfinite(tr)) = 0;
frames = [];  peaks = [];

n = robustSD(tr);
if ~isfinite(n) || n == 0, return; end
thr  = median(tr) + cfg.eventZ * n;
minD = max(1, round(cfg.eventRefrac_s * cfg.frameRate_Hz));

% A trace with no suprathreshold sample is a normal, informative outcome here,
% not a problem worth warning about.
ws = warning('off','signal:findpeaks:largeMinPeakHeight');
restoreWarn = onCleanup(@() warning(ws)); %#ok<NASGU>
[~, frames] = findpeaks(tr, ...
    'MinPeakHeight',     thr, ...
    'MinPeakProminence', cfg.eventZ * n, ...
    'MinPeakDistance',   minD);

frames = frames(:)';
peaks  = zeros(1, numel(frames));
hw     = max(1, round(cfg.eventAmpWindow_s * cfg.frameRate_Hz / 2));
for k = 1:numel(frames)
    lo = max(1, frames(k)-hw);
    hi = min(numel(tr), frames(k)+hw);
    peaks(k) = mean(tr(lo:hi));
end
end

function B = blindSample(M)
nShow = 12;
ua = M.unmatchedA(:);  ub = M.unmatchedB(:);
pick = [ones(numel(ua),1) ua; 2*ones(numel(ub),1) ub];
B = struct('order',[],'key',zeros(0,2),'note', ...
    'inspect in the given order without reading the key, then unblind');
if isempty(pick), return; end
ord  = randperm(size(pick,1), min(nShow, size(pick,1)));
B.order = ord(:);
B.key   = pick(ord,:);   % [method (1 = extractGluSNFR3, 2 = GIAnT), component index]
end

% =========================================================================
%  UNMATCHED-COMPONENT ADJUDICATION
% =========================================================================

function U = analyzeUnmatched(movCleanF, F0img, A, B, M, cfg)
%ANALYZEUNMATCHED  For every component one method found and the other did not,
%extract the trace the DATA supports at that location and classify why it went
%unmatched.
%
%  Answers "did the other method miss a real signal here?" on the real movie,
%  where there is no ground truth.  For each unmatched component:
%    nativeTrace  the discovering method's own trace
%    dataTrace    common-mask dF/F straight from the movie over that footprint's
%                 support — model-free, identical recipe for both methods
%  and a reason code:
%    noNearbyComponent  nothing from the other method within maxCentroidDistance
%                       -> a genuine detection difference
%    nearbyButLowIoU    a counterpart is close but the footprints disagree
%                       -> both found it; they disagree about its extent
%    assignedElsewhere  a qualifying counterpart existed but the one-to-one
%                       assignment gave it to a different source
%                       -> a matching artefact, not a detection difference
%  Only noNearbyComponent is evidence of a miss.

rows = struct('candidate',{},'method',{},'idx',{},'x_px',{},'y_px',{}, ...
              'area_um2',{},'reason',{}, ...
              'nearestOther_idx',{},'nearestOther_dist_um',{},'nearestOther_iou',{}, ...
              'nativeSNR',{},'nativeEventCount',{}, ...
              'commonSNR',{},'commonEventCount',{},'commonPeakDFF',{}, ...
              'dataSupportsSignal',{});
nativeTraces = {};
dataTraces   = {};

sets = {struct('S',A,'other',B,'name',"extractGluSNFR3",'idxList',M.unmatchedA(:)','isA',true), ...
        struct('S',B,'other',A,'name',"GIAnT",          'idxList',M.unmatchedB(:)','isA',false)};

n = 0;
for s = 1:2
    cur = sets{s};
    for i = cur.idxList
        n = n + 1;
        mask = cur.S.bw(:,:,i);

        native = nativeTrace(cur.S, i);
        dataTr = commonMaskDFF(movCleanF, F0img, mask, cfg);

        % nearest counterpart from the other method
        if isempty(M.distMatrix)
            dmin = Inf;  jbest = NaN;  ioubest = 0;
        else
            if cur.isA
                dvec   = M.distMatrix(i,:);
                iouvec = M.iouMatrix(i,:);
            else
                dvec   = M.distMatrix(:,i)';
                iouvec = M.iouMatrix(:,i)';
            end
            if isempty(dvec)
                dmin = Inf;  jbest = NaN;  ioubest = 0;
            else
                [dmin, jbest] = min(dvec);
                ioubest = iouvec(jbest);
            end
        end

        if ~isfinite(dmin) || dmin > cfg.maxCentroidDistance_um
            reason = "noNearbyComponent";
        elseif ioubest < cfg.minFootprintIoU
            reason = "nearbyButLowIoU";
        else
            reason = "assignedElsewhere";
        end

        [~, nativePk] = detectEvents(native, cfg);
        [~, dataPk]   = detectEvents(dataTr, cfg);

        rows(n).candidate = n;
        rows(n).method    = cur.name;
        rows(n).idx       = i;
        rows(n).x_px      = cur.S.centroid_xy(i,1);
        rows(n).y_px      = cur.S.centroid_xy(i,2);
        rows(n).area_um2  = cur.S.area_um2(i);
        rows(n).reason    = reason;
        rows(n).nearestOther_idx     = jbest;
        rows(n).nearestOther_dist_um = dmin;
        rows(n).nearestOther_iou     = ioubest;
        rows(n).nativeSNR        = traceSNR(native, cfg);
        rows(n).nativeEventCount = numel(nativePk);
        rows(n).commonSNR        = traceSNR(dataTr, cfg);
        rows(n).commonEventCount = numel(dataPk);
        if isempty(dataPk), rows(n).commonPeakDFF = NaN; else, rows(n).commonPeakDFF = max(dataPk); end

        % The adjudication: does the movie itself show events here?
        rows(n).dataSupportsSignal = double(numel(dataPk) >= cfg.unmatchedMinEvents);

        nativeTraces{n} = native; %#ok<AGROW>
        dataTraces{n}   = dataTr; %#ok<AGROW>
    end
end

U = struct();
if n == 0
    U.table = table();
else
    U.table = struct2table(rows, 'AsArray', true);
end
U.nativeTraces = nativeTraces;
U.dataTraces   = dataTraces;

U.note = ['dataSupportsSignal = the common-mask dF/F over this footprint shows ' ...
          'at least cfg.unmatchedMinEvents events; reason = why it went unmatched'];
end

% -------------------------------------------------------------------------
function printMissSummary(U, cfg)
%PRINTMISSSUMMARY  What each method found alone, and whether the movie shows
%events over those pixels.  Reports counts only — calling a component a miss or
%a false positive is left to the person inspecting the blinded gallery.
fprintf('\n---------- Unmatched components on the real movie ----------\n');
if isempty(U.table)
    fprintf('Every component matched; neither method reported a unique source.\n');
    return
end
T = U.table;
methods = ["extractGluSNFR3","GIAnT"];
for m = 1:2
    sel = T.method == methods(m);
    other = methods(3-m);
    fprintf('%s found %d component(s) with no match in %s:\n', methods(m), nnz(sel), other);
    for r = ["noNearbyComponent","nearbyButLowIoU","assignedElsewhere"]
        k = sel & T.reason == r;
        if ~any(k), continue; end
        fprintf('   %-18s %3d   (data shows events in %d of them)\n', ...
            r, nnz(k), nnz(k & T.dataSupportsSignal == 1));
    end
    % Only noNearbyComponent is a unilateral result.  A nearbyButLowIoU
    % component was found by BOTH methods — they disagree about its extent, so
    % it is neither a miss by the other method nor a false positive by this one.
    unilateral = sel & T.reason == "noNearbyComponent";
    withSig = unilateral & T.dataSupportsSignal == 1;
    fprintf('   => %d unique to %s, movie shows >= %d event(s) over those pixels', ...
        nnz(withSig), methods(m), cfg.unmatchedMinEvents);
    if any(withSig)
        fprintf(' (median peak %.3f dF/F, median SNR %.1f)', ...
            median(T.commonPeakDFF(withSig),'omitnan'), median(T.commonSNR(withSig),'omitnan'));
    end
    fprintf('\n');
    noSig = unilateral & T.dataSupportsSignal == 0;
    fprintf('   => %d unique to %s, movie shows fewer events over those pixels\n', ...
        nnz(noSig), methods(m));
    amb = sel & T.reason ~= "noNearbyComponent";
    if any(amb)
        fprintf(['   => %d component(s) both methods found but the matching rule ' ...
                 'rejected (footprint disagreement, not a detection difference)\n'], nnz(amb));
    end
end
fprintf(['Counts use >= %d event(s) at %g robust SD on the raw trace (no smoothing).\n' ...
         'These are descriptions, not verdicts; every component is listed in\n' ...
         'results/unmatched_components.csv for you to judge.\n'], ...
    cfg.unmatchedMinEvents, cfg.eventZ);
fprintf('------------------------------------------------------------\n');
end

% -------------------------------------------------------------------------
function makeUnmatchedSummaryFigure(movCleanF, A, B, U, cfg, figDir)
%MAKEUNMATCHEDSUMMARYFIGURE  The readable one: a few examples per direction.
%
%  Top block  : found by extractGluSNFR3, not by GIAnT
%  Bottom block: found by GIAnT, not by extractGluSNFR3
%
%  LABELLED, unlike the blinded gallery, because this figure exists to show
%  what the two methods disagree about — not to score it.  The quantitative
%  record stays the blinded gallery plus unmatched_components.csv.
%
%  Only reason == "noNearbyComponent" is shown: those are the real detection
%  differences.  Examples are the strongest available evidence (highest
%  common-mask SNR), which is stated on the figure so it is not mistaken for a
%  random sample.
if isempty(U.table), return; end
T = U.table;
n = cfg.unmatchedExamplesPerSide;

methods = ["extractGluSNFR3","GIAnT"];
picks = cell(1,2);
for m = 1:2
    sel = find(T.method == methods(m) & T.reason == "noNearbyComponent");
    if isempty(sel), picks{m} = []; continue; end
    [~, ord] = sort(T.commonSNR(sel), 'descend', 'MissingPlacement','last');
    picks{m} = sel(ord(1:min(n, numel(ord))));
end
if isempty(picks{1}) && isempty(picks{2}), return; end

rowsAll = [picks{1}; picks{2}];
nRows   = numel(rowsAll);
mImg    = mean(double(movCleanF),3);
tv      = (0:size(movCleanF,3)-1) / cfg.frameRate_Hz;

pad  = 15;
tile = 2*pad + 1;
sep  = 3;
colE = [0.10 0.35 0.85];    % extractGluSNFR3
colG = [0.95 0.45 0.05];    % GIAnT

% Outlines of EVERY component of each method, once, over the whole field; each
% crop below is just a window onto these.  Showing both methods inside the same
% crop is the point: an unmatched component reads as an outline of one colour
% with no outline of the other beside it.
perimA = false(size(mImg));
for k = 1:size(A.bw,3), perimA = perimA | bwperim(A.bw(:,:,k)); end
perimB = false(size(mImg));
for k = 1:size(B.bw,3), perimB = perimB | bwperim(B.bw(:,:,k)); end

strip  = ones(nRows*(tile+sep) - sep, tile, 3);
labels = strings(nRows,1);
tileTop = zeros(nRows,1);
for r = 1:nRows
    k  = rowsAll(r);
    xr = max(1,round(T.x_px(k))-pad) : min(size(mImg,2), round(T.x_px(k))+pad);
    yr = max(1,round(T.y_px(k))-pad) : min(size(mImg,1), round(T.y_px(k))+pad);
    nx = numel(xr);  ny = numel(yr);

    g = zeros(tile);  g(1:ny,1:nx) = mat2gray(mImg(yr,xr));
    pa = false(tile); pa(1:ny,1:nx) = perimA(yr,xr);
    pb = false(tile); pb(1:ny,1:nx) = perimB(yr,xr);

    % this component's own outline, thickened so it stands out from neighbours
    own = false(tile);
    ownFull = imdilate(bwperim(pickBw(A,B,T.method(k),T.idx(k))), ones(2));
    own(1:ny,1:nx) = ownFull(yr,xr);

    R = g; G = g; Bl = g;
    R(pa) = colE(1)*0.75; G(pa) = colE(2)*0.75; Bl(pa) = colE(3)*0.75;
    R(pb) = colG(1)*0.75; G(pb) = colG(2)*0.75; Bl(pb) = colG(3)*0.75;
    if T.method(k) == "extractGluSNFR3", oc = colE; else, oc = colG; end
    R(own) = oc(1); G(own) = oc(2); Bl(own) = oc(3);

    r0 = (r-1)*(tile+sep) + 1;
    strip(r0:r0+tile-1,:,:) = cat(3, R, G, Bl);
    tileTop(r) = r0;
    labels(r)  = sprintf('%d.  %s #%d', r, shortName(T.method(k)), T.idx(k));
end

% ---- traces as z-scores, one common offset step --------------------------
%  Each trace is (x - median) / robustSD of that trace, so components with
%  different baseline noise are directly comparable by eye.  robustSD uses the
%  below-median half, so positive-going events do not inflate their own
%  denominator.  Detection still runs on the trace in its own units; only the
%  display is normalised.
zTr = cell(nRows,1);
pp  = zeros(nRows,1);
for r = 1:nRows
    d  = double(U.nativeTraces{rowsAll(r)});
    sd = robustSD(d);
    if ~isfinite(sd) || sd == 0, sd = 1; end
    zTr{r} = (d - median(d,'omitnan')) / sd;
    pp(r)  = prctile(zTr{r},99.5) - prctile(zTr{r},0.5);
end
step = 1.15 * max(pp);
if ~isfinite(step) || step <= 0, step = 1; end

fh = figure('Name','Unmatched components — examples','Color','w', ...
            'Position',[40 40 1250 max(380, 95*nRows + 110)]);

% ---- left: one crop per example, BOTH methods' footprints drawn in ----------
axS = axes(fh, 'Position',[0.012 0.06 0.080 0.84]);
image(axS, strip); axis(axS,'image'); axis(axS,'off'); hold(axS,'on'); cleanAxes(axS);
for r = 1:nRows
    text(axS, 2, tileTop(r) + 5, sprintf('%d', r), ...
        'Color','w', 'FontSize', 8, 'FontWeight','bold');
end
title(axS, {'ROI', '\color[rgb]{0.10,0.35,0.85}extract  \color[rgb]{0.95,0.45,0.05}GIAnT'}, ...
    'FontSize', 8);

axT = axes(fh, 'Position',[0.175 0.06 0.810 0.84]); cleanAxes(axT);
hold(axT,'on');
for r = 1:nRows
    k    = rowsAll(r);
    dRaw = double(U.nativeTraces{k});
    z    = zTr{r};
    y0   = (nRows - r) * step;          % r = 1 at the top, matching the ROI strip
    if T.method(k) == "extractGluSNFR3", col = colE; else, col = colG; end
    plot(axT, tv, z + y0, 'Color', col, 'LineWidth', 0.5);

    % Markers sit ON the plotted trace: detection returns the frame plus an
    % amplitude ESTIMATE (a windowed mean), and plotting the marker at that
    % estimate puts it off the curve.  Use the displayed trace's own value at
    % the detected frame instead.
    ef = detectEvents(dRaw, cfg);
    if ~isempty(ef)
        plot(axT, tv(ef), z(ef) + y0, 'v', ...
             'Color',[0.75 0 0], 'MarkerFaceColor',[0.75 0 0], 'MarkerSize',4);
    end
end
axis(axT,'tight');
% YTick must ascend, so flip the labels to keep row 1 at the top
set(axT, 'YTick', (0:nRows-1)*step, 'YTickLabel', flipud(labels), 'FontSize', 8, ...
         'TickLabelInterpreter','none', 'Box','off');
xlabel(axT, 'time (s)', 'FontSize', 9);

% dF/F scale bar
xl = xlim(axT);
xb   = xl(2) - 0.02*diff(xl);
yb   = -0.95*step;                     % clear of the bottom trace's excursions
zbar = round(0.5*step);
if zbar < 1, zbar = 0.5*step; end
plot(axT, [xb xb], [yb yb+zbar], 'k', 'LineWidth', 2);
text(axT, xb - 0.008*diff(xl), yb + zbar/2, sprintf('%g\\sigma', zbar), ...
     'HorizontalAlignment','right', 'FontSize', 8);
ylim(axT, [yb - 0.10*step, (nRows-1)*step + 0.7*step]);

sgtitle(fh, sprintf(['Components found by one method and not the other ' ...
                     '(top %d per side by common-mask SNR)   ' ...
                     'blue = extractGluSNFR3, orange = GIAnT   ' ...
                     'traces z-scored, peaks marked from the raw trace'], n), ...
        'FontSize', 10, 'Interpreter','none');
safeExport(fh, fullfile(figDir,'fig06_unmatched_examples.png'));
end

function b = pickBw(A, B, method, idx)
if method == "extractGluSNFR3", b = A.bw(:,:,idx); else, b = B.bw(:,:,idx); end
end

function s = shortName(m)
if m == "extractGluSNFR3", s = 'extract'; else, s = 'GIAnT'; end
end

function cleanAxes(ax)
%  Hide the interaction toolbar; exportgraphics otherwise renders it as marks
%  in the corner of every panel.
try, ax.Toolbar.Visible = 'off'; catch, end
try, disableDefaultInteractivity(ax); catch, end
end

% =========================================================================
%  SYNTHETIC SIGNAL GENERATION
% =========================================================================

function [bgF, bgF0, info] = selectFixedBackgroundCrop(movCleanF, F0img, cfg)
%SELECTFIXEDBACKGROUNDCROP  Frozen rule, defined before any simulation ran:
%  slide a simSize window over the field and take the one whose mean image
%  maximises median(brightness) * IQR(brightness) — bright tissue that also
%  contains dark vascular structure.  Never chosen from extraction performance.

sz = cfg.simSize(1:2);
[H,W,T] = size(movCleanF);
assert(H >= sz(1) && W >= sz(2), ...
    'Field %dx%d is smaller than the %dx%d simulation crop.', H, W, sz(1), sz(2));
assert(T >= cfg.simSize(3), 'Movie has %d frames, simulation needs %d.', T, cfg.simSize(3));

mImg = double(mean(movCleanF,3));
step = 25;
best = -Inf;  bestYX = [1 1];
for y0 = 1:step:(H - sz(1) + 1)
    for x0 = 1:step:(W - sz(2) + 1)
        p = mImg(y0:y0+sz(1)-1, x0:x0+sz(2)-1);
        score = median(p(:)) * iqr(p(:));
        if score > best, best = score;  bestYX = [y0 x0]; end
    end
end
y0 = bestYX(1);  x0 = bestYX(2);
bgF  = single(movCleanF(y0:y0+sz(1)-1, x0:x0+sz(2)-1, 1:cfg.simSize(3)));
bgF0 = single(F0img    (y0:y0+sz(1)-1, x0:x0+sz(2)-1));
info = struct('y0',y0,'x0',x0,'size',cfg.simSize,'score',best, ...
              'rule','max median(meanImg)*iqr(meanImg) on a 25 px grid');
end

function S = stratifyBackground(bgF, bgF0, cfg)
%STRATIFYBACKGROUND  Location classes for stratified injection (runbook 7.6).
mImg = double(mean(bgF,3));
vImg = var(double(bgF), 0, 3);
[H,W] = size(mImg);

S = struct();
S.meanImg = mImg;  S.varImg = vImg;  S.F0 = double(bgF0);

% Vessels: dark, locally structured regions of the mean image.
vesselMask = imopen(mImg < prctile(mImg(:),15), strel('disk',2));
S.vesselMask      = vesselMask;
S.distToVessel_px = bwdist(vesselMask);

S.highVar = vImg >= prctile(vImg(:), 70);
S.lowVar  = vImg <= prctile(vImg(:), 30);
S.bright  = mImg >= prctile(mImg(:), 70);
S.dim     = mImg <= prctile(mImg(:), 30);

edge = round(cfg.boundaryMargin_px) + 12;      % footprints must fit inside
S.valid = false(H,W);
S.valid(edge+1:H-edge, edge+1:W-edge) = true;
S.valid = S.valid & ~vesselMask & isfinite(mImg) & mImg > 0;

S.classes = struct( ...
    'mixed',      S.valid, ...
    'lowVar',     S.valid & S.lowVar, ...
    'highVar',    S.valid & S.highVar, ...
    'nearVessel', S.valid & S.distToVessel_px <= 5, ...
    'farVessel',  S.valid & S.distToVessel_px >= 15, ...
    'bright',     S.valid & S.bright, ...
    'dim',        S.valid & S.dim, ...
    'center',     S.valid & centerMask(H,W,true), ...
    'edge',       S.valid & centerMask(H,W,false));
S.classNames = fieldnames(S.classes);
end

function m = centerMask(H, W, wantCenter)
[rr,cc] = ndgrid(1:H, 1:W);
d = hypot(rr-(H+1)/2, cc-(W+1)/2) / (min(H,W)/2);
if wantCenter, m = d <= 0.5; else, m = d >= 0.7; end
end

% -------------------------------------------------------------------------
function T = makeSimulationManifest(cfg)
%MAKESIMULATIONMANIFEST  One row per simulated movie; datasets 1-6 of runbook
%section 8.  Under cfg.pilotOnly the repeat count is reduced but every
%amplitude — including the lowest — is kept, so the detection threshold stays
%bracketed from below and above.

nRep = cfg.nRepeats;
if cfg.pilotOnly, nRep = cfg.pilotRepeats; end

ds = strings(0,1);  amp = []; sg = []; tr = []; td = []; ns = [];
lc = strings(0,1);  pd = []; to = [];

    function push(dsName, a, s_, r_, d_, n_, class_, pdist, tover, reps)
        for q = 1:reps
            ds(end+1,1) = dsName;   amp(end+1,1) = a;    sg(end+1,1) = s_; %#ok<AGROW>
            tr(end+1,1) = r_;       td(end+1,1) = d_;    ns(end+1,1) = n_; %#ok<AGROW>
            lc(end+1,1) = class_;   pd(end+1,1) = pdist; to(end+1,1) = tover; %#ok<AGROW>
        end
    end

% 1. single-source sensitivity: amplitude x size x kinetics
for a = cfg.amplitudes_dFF
    for s_ = cfg.sigma_um
        for ki = 1:numel(cfg.tauRise_s)
            push("single_source", a, s_, cfg.tauRise_s(ki), cfg.tauDecay_s(ki), ...
                 1, "mixed", NaN, NaN, nRep);
        end
    end
end
% 2. multi-source specificity
for a = cfg.amplitudes_dFF(2:4)
    for n_ = cfg.sourcesPerMovie(2:end)
        push("multi_source", a, cfg.sigma_um(2), cfg.tauRise_s(2), cfg.tauDecay_s(2), ...
             n_, "mixed", NaN, NaN, nRep);
    end
end
% 2b. empty controls — their detections are counted as false positives
push("empty_control", 0, cfg.sigma_um(2), cfg.tauRise_s(2), cfg.tauDecay_s(2), ...
     0, "none", NaN, NaN, cfg.nEmptyControls);
% 3. source separation: distance x temporal overlap
for dist = [0.5 1 2 3 5]
    for ov = [0 1]
        push("separation", cfg.amplitudes_dFF(4), cfg.sigma_um(2), cfg.tauRise_s(2), ...
             cfg.tauDecay_s(2), 2, "mixed", dist, ov, nRep);
    end
end
% 4. background difficulty (stratified single sources)
for class_ = ["lowVar","highVar","nearVessel","farVessel","bright","dim","center","edge"]
    for a = cfg.amplitudes_dFF([2 4])
        push("background", a, cfg.sigma_um(2), cfg.tauRise_s(2), cfg.tauDecay_s(2), ...
             1, class_, NaN, NaN, nRep);
    end
end
% 5. kernel mismatch — kinetics swept while each method keeps its fixed kernel
for ri = 1:numel(cfg.tauRise_s)
    for di = 1:numel(cfg.tauDecay_s)
        push("kernel_mismatch", cfg.amplitudes_dFF(4), cfg.sigma_um(2), ...
             cfg.tauRise_s(ri), cfg.tauDecay_s(di), 1, "mixed", NaN, NaN, ...
             max(2, round(nRep/4)));
    end
end

keep = (td > tr) | (ns == 0);      % tauDecay must exceed tauRise
T = table(ds(keep), amp(keep), sg(keep), tr(keep), td(keep), ns(keep), ...
          lc(keep), pd(keep), to(keep), ...
    'VariableNames', {'dataset','amplitude_dFF','sigma_um','tauRise_s','tauDecay_s', ...
                      'nSources','location_class','pairDistance_um','temporalOverlap'});
T.movie_id = (1:height(T))';
T.seed     = cfg.seed*1000 + T.movie_id;
T = movevars(T, 'movie_id', 'Before', 1);
end

% -------------------------------------------------------------------------
function [movSim, truth] = injectSyntheticSignals(bgF, bgF0, bgStrata, row, cfg)
%INJECTSYNTHETICSIGNALS  Add ground-truth sources to a real residual background.
%      F_sim = F_bg + F0 .* A .* S(x,y) .* T(t)
%  with S a peak-normalised subpixel-centred Gaussian and T a peak-normalised
%  difference-of-exponentials train.  The SAME array is handed to both methods.

rng(row.seed, 'twister');
[H,W,T] = size(bgF);
movSim = single(bgF);

truth = struct();
truth.frameRate_Hz  = cfg.frameRate_Hz;
truth.x = [];  truth.y = [];  truth.sigma_px = [];  truth.amplitude_dFF = [];
truth.eventFrames = {};
truth.locationClass = strings(0,1);
truth.footprint = zeros(H,W,0,'single');
truth.bw        = false(H,W,0);
truth.trace     = zeros(0,T,'single');
truth.kernel    = [];
truth.peakFrames = {};
truth.kernelPeakOffset = 0;
truth.verifiedPeakDFF = [];
truth.minSeparation_px = NaN;

if row.nSources == 0
    return
end

sigma_px = row.sigma_um / cfg.pixelSize_um;
k = glutamateKernel(cfg.frameRate_Hz, row.tauRise_s, row.tauDecay_s, ...
                    min(8*row.tauDecay_s, T/cfg.frameRate_Hz/4));
truth.kernel = k;
kLen = numel(k);
% A difference-of-exponentials peaks some frames AFTER the release event, by an
% amount that depends on tauRise/tauDecay.  Event scoring must compare detected
% peak times against the expected PEAK frames, not the onsets, or the
% kernel-mismatch condition would be scored as a timing failure that is really
% just the kernel's own delay.
[~, kPeakOffset] = max(k);
kPeakOffset = kPeakOffset - 1;
truth.kernelPeakOffset = kPeakOffset;

% ---- choose locations ----------------------------------------------------
className = char(row.location_class);
if ~isfield(bgStrata.classes, className)
    className = 'mixed';
end
cands = find(bgStrata.classes.(className));
assert(~isempty(cands), 'No valid injection sites for class %s', className);

if row.dataset == "separation"
    c0 = cands(randi(numel(cands)));
    [y0,x0] = ind2sub([H W], c0);
    dpx = row.pairDistance_um / cfg.pixelSize_um;
    th  = 2*pi*rand;
    xs = [x0 - dpx/2*cos(th), x0 + dpx/2*cos(th)] + (rand(1,2)-0.5);
    ys = [y0 - dpx/2*sin(th), y0 + dpx/2*sin(th)] + (rand(1,2)-0.5);
    xs = min(max(xs, 6), W-5);  ys = min(max(ys, 6), H-5);
else
    % Individually resolvable sources.  Rejection sampling can fail when the
    % requested count, footprint size and valid area are in tension, so the
    % minimum separation is relaxed in steps rather than aborting the sweep.
    % The separation actually achieved is recorded in truth.minSeparation_px.
    minSep0 = max(4, 3*sigma_px);
    minSepFloor = max(2, 2*sigma_px);
    minSep = minSep0;
    xs = [];  ys = [];
    while true
        xs = [];  ys = [];
        for tries = 1:20000
            c = cands(randi(numel(cands)));
            [yy,xx] = ind2sub([H W], c);
            xx = xx + rand - 0.5;  yy = yy + rand - 0.5;
            if isempty(xs) || min(hypot(xs-xx, ys-yy)) >= minSep
                xs(end+1) = xx;  ys(end+1) = yy; %#ok<AGROW>
            end
            if numel(xs) == row.nSources, break; end
        end
        if numel(xs) == row.nSources || minSep <= minSepFloor, break; end
        minSep = max(minSepFloor, minSep * 0.8);
    end
    assert(numel(xs) == row.nSources, ...
        ['Could not place %d sources of sigma %.1f px in movie %d even at the ' ...
         'minimum separation of %.1f px (%d valid pixels in class %s).  Reduce ' ...
         'cfg.sourcesPerMovie or enlarge cfg.simSize.'], ...
        row.nSources, sigma_px, row.movie_id, minSepFloor, numel(cands), className);
    truth.minSeparation_px = minSep;
end

nS = numel(xs);
truth.x = xs(:);  truth.y = ys(:);
truth.sigma_px      = repmat(sigma_px, nS, 1);
truth.amplitude_dFF = repmat(row.amplitude_dFF, nS, 1);
truth.locationClass = repmat(string(className), nS, 1);
truth.footprint = zeros(H,W,nS,'single');
truth.bw        = false(H,W,nS);
truth.trace     = zeros(nS,T,'single');
truth.eventFrames = cell(nS,1);
truth.peakFrames  = cell(nS,1);

for s = 1:nS
    Sxy = gaussianFootprint(H, W, xs(s), ys(s), sigma_px);

    if row.dataset == "separation" && row.temporalOverlap == 1 && s == 2
        ev = truth.eventFrames{1};                       % perfectly overlapping pair
    else
        nEv = max(1, min(cfg.maxEventsPerSource, ...
                    poissrnd(cfg.eventRate_Hz * T / cfg.frameRate_Hz)));
        lo = kLen + 1;  hi = T - kLen;                   % margin = kernel duration
        ev = sort(randperm(hi-lo+1, min(nEv, hi-lo+1)) + lo - 1);
        if numel(ev) > 1                                 % keep isolated events isolated
            ev = ev([true, diff(ev) > kLen]);
        end
    end
    truth.eventFrames{s} = ev;
    truth.peakFrames{s}  = min(ev + kPeakOffset, T);

    train = zeros(T,1,'single');
    train(ev) = 1;
    tr = conv(train, k, 'full');
    tr = tr(1:T);
    truth.trace(s,:) = tr;

    signalDFF = row.amplitude_dFF .* reshape(Sxy,H,W,1) .* reshape(tr,1,1,T);
    movSim = movSim + single(double(bgF0) .* double(signalDFF));

    truth.footprint(:,:,s) = Sxy;
    truth.bw(:,:,s) = Sxy >= cfg.footprintThreshFrac;    % Sxy is peak-normalised
end

% ---- verify numerically that the injected peak equals the requested dF/F ----
truth.verifiedPeakDFF = nan(nS,1);
for s = 1:nS
    yi = min(max(round(truth.y(s)),1),H);
    xi = min(max(round(truth.x(s)),1),W);
    pk = max(truth.trace(s,:));
    if pk > 0
        truth.verifiedPeakDFF(s) = truth.amplitude_dFF(s) * ...
            double(truth.footprint(yi,xi,s)) * double(pk);
    end
end
end

function S = gaussianFootprint(Y, X, x0, y0, sigma_px)
[xx,yy] = meshgrid(1:X, 1:Y);
S = exp(-((xx-x0).^2 + (yy-y0).^2) / (2*sigma_px^2));
S = single(S ./ max(S(:)));
end

function k = glutamateKernel(frameRate_Hz, tauRise_s, tauDecay_s, duration_s)
assert(tauDecay_s > tauRise_s, 'tauDecay_s must exceed tauRise_s.');
t = (0:ceil(duration_s*frameRate_Hz)-1)' / frameRate_Hz;
k = exp(-t/tauDecay_s) - exp(-t/tauRise_s);
k(k < 0) = 0;
k = single(k ./ max(k));
end

% =========================================================================
%  GROUND-TRUTH SCORING
% =========================================================================

function Tout = scoreAgainstTruth(result, truth, row, cfg, methodName)
%SCOREAGAINSTTRUTH  Detection, spatial and temporal metrics, one row per true
%source (or one summary row for an empty-control movie, with nSourcesTrue = 0
%so its detections count as false positives).
%
%  A detection is a true positive only if it passes BOTH preregistered spatial
%  criteria: centroid error <= cfg.maxCentroidDistance_um AND IoU >= cfg.minFootprintIoU.
%  Assignment is one-to-one (Hungarian), never nearest-neighbour per component.

S = standardizeFootprints(result, cfg);
nTrue = size(truth.footprint,3);
nDet  = size(S.footprints,3);

[pairs, dpair, ipair, cpair] = matchToTruth(S, truth, cfg);
TP = size(pairs,1);
FN = nTrue - TP;
FP = nDet  - TP;

% Sensitivity check: the same matching with the shared area/boundary filter
% switched off.  The filter is biologically motivated and applied identically
% to both methods, but the two produce very differently shaped footprints
% (extractGluSNFR3 accumulates thresholded z-score blobs, typically ~90 px
% above 0.2*max; GIAnT returns Gaussian-smoothed NMF weights, typically ~15 px)
% so it can still bite them unequally.  Reporting both keeps that visible
% instead of burying it in the primary number.
%
% NOTE when reading the *_raw columns: they are NOT symmetric.  extractGluSNFR3
% has already applied its own boundary exclusion (opts.bound) and minimum event
% count internally, so its "unfiltered" count is itself pre-filtered, whereas
% GIAnT's is the full peak list.  Use the filtered columns as the primary
% result; the raw ones only to see how much the shared filter moved each method.
cfgRaw = cfg;
cfgRaw.minComponentArea_px = 0;
cfgRaw.maxComponentArea_px = Inf;
cfgRaw.boundaryMargin_px   = 0;
Sraw = standardizeFootprints(result, cfgRaw);
pairsRaw = matchToTruth(Sraw, truth, cfgRaw);
TPraw = size(pairsRaw,1);
FNraw = nTrue - TPraw;
FPraw = size(Sraw.footprints,3) - TPraw;

nRows = max(nTrue, 1);
rows = repmat(emptyScoreRow(), nRows, 1);
for i = 1:nRows
    rows(i).movie_id       = row.movie_id;
    rows(i).seed           = row.seed;
    rows(i).dataset        = row.dataset;
    rows(i).method         = string(methodName);
    rows(i).amplitude_dFF  = row.amplitude_dFF;
    rows(i).sigma_um       = row.sigma_um;
    rows(i).tauRise_s      = row.tauRise_s;
    rows(i).tauDecay_s     = row.tauDecay_s;
    rows(i).nSourcesTrue   = nTrue;
    rows(i).nDetected      = nDet;
    rows(i).TP = TP;  rows(i).FP = FP;  rows(i).FN = FN;
    rows(i).nDetected_raw = size(Sraw.footprints,3);
    rows(i).TP_raw = TPraw;  rows(i).FP_raw = FPraw;  rows(i).FN_raw = FNraw;
    rows(i).rejectedTooSmall  = S.rejectReason.tooSmall;
    rows(i).rejectedTooLarge  = S.rejectReason.tooLarge;
    rows(i).rejectedAtBoundary = S.rejectReason.atBoundary;
    rows(i).runtime_s      = result.runtime_s;
    rows(i).failed         = double(result.failed);
    rows(i).location_class = "none";
end

for i = 1:nTrue
    rows(i).location_class = truth.locationClass(i);
    rows(i).x_true      = truth.x(i);
    rows(i).y_true      = truth.y(i);
    rows(i).event_count = numel(truth.eventFrames{i});

    m = find(pairs(:,1) == i, 1);
    if isempty(m), continue; end
    j = pairs(m,2);

    rows(i).detected          = 1;
    rows(i).matched_component = j;
    rows(i).centroid_error_um = dpair(m);
    rows(i).footprint_iou     = ipair(m);
    rows(i).spatial_corr      = cpair(m);
    rows(i).area_ratio        = nnz(S.bw(:,:,j)) / max(nnz(truth.bw(:,:,i)), 1);

    % ---- temporal ------------------------------------------------------
    est   = nativeTrace(S, j);
    trueT = double(truth.trace(i,:)) * truth.amplitude_dFF(i);   % dF/F units

    [pc, ~, lag] = traceAgreement(est, trueT, cfg);
    rows(i).trace_corr        = pc;
    rows(i).peak_time_error_s = lag;

    estAmp  = robustPeakAmp(est, cfg);
    trueAmp = max(trueT);
    rows(i).amplitude_bias = (estAmp - trueAmp) / max(trueAmp, eps);
    rows(i).decay_bias     = decayBias(est, truth.eventFrames{i}, row.tauDecay_s, cfg);

    [rec, prec, nFalse]    = eventAgreement(est, truth.peakFrames{i}, cfg);
    rows(i).event_recall    = rec;
    rows(i).event_precision = prec;
    rows(i).false_events    = nFalse;

    rows(i).baseline_noise = robustSD(est);
    rows(i).output_snr     = traceSNR(est, cfg);
end

Tout = struct2table(rows, 'AsArray', true);
end

function r = emptyScoreRow()
r = struct('movie_id',NaN,'seed',NaN,'dataset',"",'method',"", ...
    'amplitude_dFF',NaN,'sigma_um',NaN,'tauRise_s',NaN,'tauDecay_s',NaN, ...
    'nSourcesTrue',NaN,'event_count',NaN,'location_class',"", ...
    'x_true',NaN,'y_true',NaN,'detected',0,'matched_component',NaN, ...
    'centroid_error_um',NaN,'footprint_iou',NaN,'spatial_corr',NaN,'area_ratio',NaN, ...
    'trace_corr',NaN,'peak_time_error_s',NaN,'amplitude_bias',NaN,'decay_bias',NaN, ...
    'event_recall',NaN,'event_precision',NaN,'false_events',NaN, ...
    'baseline_noise',NaN,'output_snr',NaN,'nDetected',NaN, ...
    'TP',NaN,'FP',NaN,'FN',NaN, ...
    'nDetected_raw',NaN,'TP_raw',NaN,'FP_raw',NaN,'FN_raw',NaN, ...
    'rejectedTooSmall',NaN,'rejectedTooLarge',NaN,'rejectedAtBoundary',NaN, ...
    'runtime_s',NaN,'failed',0);
end

function [pairs, dpair, ipair, cpair] = matchToTruth(S, truth, cfg)
%MATCHTOTRUTH  One-to-one Hungarian assignment of detected components to true
%sources.  A pair qualifies only if it passes BOTH preregistered criteria:
%centroid error <= cfg.maxCentroidDistance_um AND IoU >= cfg.minFootprintIoU.
nTrue = size(truth.footprint,3);
nDet  = size(S.footprints,3);
pairs = zeros(0,2);  dpair = [];  ipair = [];  cpair = [];
if nTrue == 0 || nDet == 0, return; end

truthStd = struct('bw', truth.bw, 'normFp', truth.footprint, ...
                  'centroid_xy', [truth.x(:) truth.y(:)]);
[D, IOU, SC] = pairwiseSpatial(truthStd, S, cfg);
cost    = 0.5*(D/cfg.maxCentroidDistance_um) + 0.5*(1-IOU);
allowed = D <= cfg.maxCentroidDistance_um & IOU >= cfg.minFootprintIoU;
cost(~allowed) = 1e6;
p = matchpairs(cost, 1.0);
if ~isempty(p)
    ok = false(size(p,1),1);
    for k = 1:size(p,1), ok(k) = allowed(p(k,1), p(k,2)); end
    pairs = p(ok,:);
end
if ~isempty(pairs)
    idx   = sub2ind([nTrue nDet], pairs(:,1), pairs(:,2));
    dpair = D(idx);  ipair = IOU(idx);  cpair = SC(idx);
end
end

function [rec, prec, nFalse] = eventAgreement(est, trueFrames, cfg)
detFrames = detectEvents(est, cfg);
tol = max(1, round(cfg.eventTol_s * cfg.frameRate_Hz));
if isempty(trueFrames)
    rec = NaN;  prec = NaN;  nFalse = numel(detFrames);  return
end
matchedTrue = false(size(trueFrames));
usedDet     = false(size(detFrames));
for a = 1:numel(trueFrames)
    if isempty(detFrames), break; end
    d = abs(detFrames - trueFrames(a));
    d(usedDet) = Inf;
    [mn, im] = min(d);
    if mn <= tol
        matchedTrue(a) = true;  usedDet(im) = true;
    end
end
rec    = mean(matchedTrue);
prec   = nnz(usedDet) / max(numel(detFrames),1);
nFalse = numel(detFrames) - nnz(usedDet);
end

function b = decayBias(est, trueFrames, tauDecayTrue, cfg)
%  Single-exponential fit to the mean event-triggered decay of the estimated
%  trace, compared with the injected tauDecay.
b = NaN;
if isempty(trueFrames), return; end
n = round(3 * tauDecayTrue * cfg.frameRate_Hz);
if n < 4, return; end
segs = [];
for f = trueFrames
    if f + n <= numel(est)
        segs(end+1,:) = est(f:f+n); %#ok<AGROW>
    end
end
if isempty(segs), return; end
m = mean(segs, 1);
[pk, ipk] = max(m);
if pk <= 0 || ipk >= numel(m)-3, return; end
y = m(ipk:end);  y = y - min(y);
t = (0:numel(y)-1) / cfg.frameRate_Hz;
good = y > 0.05*max(y);
if nnz(good) < 4, return; end
p = polyfit(t(good), log(y(good)), 1);
if p(1) >= 0, return; end
b = (-1/p(1) - tauDecayTrue) / tauDecayTrue;
end

% =========================================================================
%  PART C — end-to-end (inject before artifact regression)
% =========================================================================

function E = runEndToEndSubset(meta, win, bgInfo, bgStrata, cfg)
%RUNENDTOENDSUBSET  Inject into the RAW window, then run the whole pipeline.
%  Reported separately: it measures signal attenuation caused by the artifact
%  regression itself, not extractor quality.

E = struct('rows', table(), 'attenuation', table(), 'note', ...
    ['synthetic sources injected into raw camera counts before ' ...
     'regressArtifactWrapper; recovered dF/F compared with injected dF/F']);

amps = cfg.amplitudes_dFF([2 4]);
nRep = 3;
rowsOut = table();
attRows = [];

rawWin = readGluFrames(meta, win, cfg);      % same frozen window as Part A
qc0 = win.benchIdx;
y0 = bgInfo.y0;  x0 = bgInfo.x0;  sz = cfg.simSize;
yr = y0:y0+sz(1)-1;  xr = x0:x0+sz(2)-1;

% F0 in RAW units over the crop, used to convert dF/F into counts
F0raw = single(mean(rawWin(yr, xr, qc0), 3));

id = 0;
for a = amps
    for r = 1:nRep
        id = id + 1;
        seed = cfg.seed*100000 + id;

        row = table(id, seed, "end_to_end", a, cfg.sigma_um(2), cfg.tauRise_s(2), ...
                    cfg.tauDecay_s(2), 5, "mixed", NaN, NaN, ...
            'VariableNames', {'movie_id','seed','dataset','amplitude_dFF','sigma_um', ...
                              'tauRise_s','tauDecay_s','nSources','location_class', ...
                              'pairDistance_um','temporalOverlap'});

        rawCrop = single(rawWin(yr, xr, qc0));
        [injCrop, truth] = injectSyntheticSignals(rawCrop, F0raw, bgStrata, row, cfg);

        % Put the injected crop back and regress exactly as in Part A.  The
        % window is patched in place and restored afterwards; copying it would
        % duplicate several GB on every repeat.
        origCrop = rawWin(yr, xr, qc0);
        rawWin(yr, xr, qc0) = double(injCrop);
        cl = regressArtifactWrapper(rawWin, meta, cfg);
        movE2E = single(cl.F(yr, xr, qc0));
        rawWin(yr, xr, qc0) = origCrop;
        clear cl origCrop

        % Signal attenuation at the true footprints, before vs after regression
        for s = 1:size(truth.footprint,3)
            mask  = truth.bw(:,:,s);
            preT  = commonMaskDFF(injCrop, F0raw, mask, cfg);
            postT = commonMaskDFF(movE2E,  F0raw, mask, cfg);
            attRows(end+1,:) = [id, a, s, max(preT), max(postT), ...
                                max(postT)/max(max(preT),eps)]; %#ok<AGROW>
        end

        outE = tryRun(@() runExtractGluSNFR3(movE2E, cfg.extractParams, cfg, ...
                    sprintf('e2e%03d_extract', id)), size(movE2E));
        outG = tryRun(@() runGIAnT(movE2E, cfg.giantParams, cfg, ...
                    sprintf('e2e%03d_giant', id)), size(movE2E));

        rowsOut = [rowsOut; scoreAgainstTruth(outE, truth, row, cfg, 'extractGluSNFR3'); ...
                            scoreAgainstTruth(outG, truth, row, cfg, 'GIAnT')]; %#ok<AGROW>
    end
end

E.rows = rowsOut;
if ~isempty(attRows)
    E.attenuation = array2table(attRows, 'VariableNames', ...
        {'movie_id','amplitude_dFF','source','peakDFF_preRegression', ...
         'peakDFF_postRegression','retainedFraction'});
end
end

% =========================================================================
%  SUMMARY
% =========================================================================

function results = summarizeBenchmark(realMetrics, simRows, endToEnd, simFail, cfg)
results = struct();
results.cfg         = cfg;
results.realMetrics = realMetrics;
results.simRows     = simRows;
results.endToEnd    = endToEnd;
results.failures    = simFail;

methods = ["extractGluSNFR3","GIAnT"];
det = table();
for m = 1:numel(methods)
    sub = simRows(simRows.method == methods(m), :);
    if isempty(sub), continue; end
    perMovie = groupsummary(sub, 'movie_id', 'max', {'TP','FP','FN'});
    TP = sum(perMovie.max_TP);  FP = sum(perMovie.max_FP);  FN = sum(perMovie.max_FN);
    rec  = TP / max(TP+FN, 1);
    prec = TP / max(TP+FP, 1);

    emp = sub(sub.dataset == "empty_control", :);
    empMovies = unique(emp.movie_id);
    fpEmpty = 0;
    for k = 1:numel(empMovies)
        fpEmpty = fpEmpty + max(emp.FP(emp.movie_id == empMovies(k)));
    end
    nEmpty = max(numel(empMovies), 1);
    nPix   = prod(cfg.simSize(1:2));

    % Same numbers with the shared area/boundary filter switched off
    TPr = sum(perMovieRaw(perMovie, sub, 'TP_raw'));
    FPr = sum(perMovieRaw(perMovie, sub, 'FP_raw'));
    FNr = sum(perMovieRaw(perMovie, sub, 'FN_raw'));
    recR  = TPr / max(TPr+FNr, 1);
    precR = TPr / max(TPr+FPr, 1);

    keptFrac = mean(sub.nDetected ./ max(sub.nDetected_raw,1), 'omitnan');

    det = [det; table(methods(m), TP, FP, FN, rec, prec, ...
        2*rec*prec/max(rec+prec, eps), fpEmpty/nEmpty, ...
        fpEmpty/nEmpty * 40000/nPix, mean(sub.runtime_s,'omitnan'), ...
        recR, precR, 2*recR*precR/max(recR+precR, eps), keptFrac, ...
        'VariableNames', {'method','TP','FP','FN','recall','precision','F1', ...
                          'FP_per_empty_movie','FP_per_40000px','meanRuntime_s', ...
                          'recall_unfiltered','precision_unfiltered','F1_unfiltered', ...
                          'fractionKeptByFilter'})]; %#ok<AGROW>
end
results.detection = det;

% Fairness check: the shared area/boundary filter must not remove a very
% different share of each method's components (runbook section 11).
if height(det) == 2 && all(isfinite(det.fractionKeptByFilter))
    kf = det.fractionKeptByFilter;
    if max(kf)/max(min(kf), eps) > 2
        warning('Benchmark:FilterAsymmetry', ...
            ['The shared component filter kept %.0f%% of %s components but only ' ...
             '%.0f%% of %s components.  Report the *_unfiltered columns alongside ' ...
             'the primary numbers and check cfg.minComponentArea_um2 / ' ...
             'cfg.maxComponentArea_um2 against both methods'' footprint sizes.'], ...
            100*kf(1), det.method(1), 100*kf(2), det.method(2));
    end
end

% Detection probability vs dF/F: 50 % and 90 % thresholds with bootstrap CIs
results.logistic = table();
for m = 1:numel(methods)
    sub = simRows(simRows.method == methods(m) & simRows.nSourcesTrue > 0 & ...
                  simRows.amplitude_dFF > 0, :);
    L = logisticThresholds(sub);
    L.method = methods(m);
    results.logistic = [results.logistic; struct2table(L)]; %#ok<AGROW>
end

% Paired per-source comparison — both methods saw the same movie and source
results.paired = buildPairedTable(simRows);
P = results.paired;
if isempty(P)
    empt = pairedBootstrap([], []);
    results.pairedStats = struct('iou',empt,'centroid',empt,'traceCorr',empt,'ampBias',empt);
else
    results.pairedStats = struct( ...
        'iou',       pairedBootstrap(P.iou_extract,  P.iou_giant), ...
        'centroid',  pairedBootstrap(P.err_extract,  P.err_giant), ...
        'traceCorr', pairedBootstrap(P.corr_extract, P.corr_giant), ...
        'ampBias',   pairedBootstrap(abs(P.ampbias_extract), abs(P.ampbias_giant)));
end

if ~isempty(endToEnd.attenuation)
    results.attenuationSummary = groupsummary(endToEnd.attenuation, 'amplitude_dFF', ...
        {'mean','std'}, 'retainedFraction');
end
end

function v = perMovieRaw(~, sub, varName)
%  One value per simulated movie (the per-movie counts are repeated on every
%  source row, so take the max within each movie rather than summing rows).
g = groupsummary(sub, 'movie_id', 'max', varName);
v = g.(['max_' varName]);
end

function P = buildPairedTable(simRows)
P = table();
if isempty(simRows), return; end
key = strcat(string(simRows.movie_id), '_', string(simRows.x_true), '_', string(simRows.y_true));
A = simRows(simRows.method == "extractGluSNFR3", :);  keyA = key(simRows.method == "extractGluSNFR3");
B = simRows(simRows.method == "GIAnT", :);            keyB = key(simRows.method == "GIAnT");
[~, ia, ib] = intersect(keyA, keyB, 'stable');
if isempty(ia), return; end
P = table(keyA(ia), A.amplitude_dFF(ia), A.dataset(ia), ...
          A.detected(ia), B.detected(ib), ...
          A.footprint_iou(ia), B.footprint_iou(ib), ...
          A.centroid_error_um(ia), B.centroid_error_um(ib), ...
          A.trace_corr(ia), B.trace_corr(ib), ...
          A.amplitude_bias(ia), B.amplitude_bias(ib), ...
    'VariableNames', {'srcKey','amplitude_dFF','dataset','det_extract','det_giant', ...
                      'iou_extract','iou_giant','err_extract','err_giant', ...
                      'corr_extract','corr_giant','ampbias_extract','ampbias_giant'});
end

function L = logisticThresholds(sub)
L = struct('n',height(sub),'dFF50',NaN,'dFF90',NaN, ...
           'dFF50_lo',NaN,'dFF50_hi',NaN,'dFF90_lo',NaN,'dFF90_hi',NaN);
if isempty(sub), return; end
good = isfinite(sub.amplitude_dFF) & sub.amplitude_dFF > 0 & isfinite(sub.detected);
if nnz(good) < 20 || numel(unique(sub.detected(good))) < 2, return; end
x = log10(sub.amplitude_dFF(good));
y = sub.detected(good);
seeds = sub.seed(good);

[L.dFF50, L.dFF90] = fitOne(x, y);

% Bootstrap BY SIMULATION SEED (movie), never by frame or by source.
nb = 200;
b50 = nan(nb,1);  b90 = nan(nb,1);
us = unique(seeds);
for b = 1:nb
    pick = us(randi(numel(us), numel(us), 1));
    sel  = ismember(seeds, pick);
    if numel(unique(y(sel))) < 2, continue; end
    [b50(b), b90(b)] = fitOne(x(sel), y(sel));
end
L.dFF50_lo = prctile(b50, 2.5);  L.dFF50_hi = prctile(b50, 97.5);
L.dFF90_lo = prctile(b90, 2.5);  L.dFF90_hi = prctile(b90, 97.5);

    function [a50, a90] = fitOne(xx, yy)
        a50 = NaN;  a90 = NaN;
        try
            bfit = glmfit(xx, yy, 'binomial', 'link', 'logit');
            a50 = 10^(( 0      - bfit(1)) / bfit(2));
            a90 = 10^(( log(9) - bfit(1)) / bfit(2));
        catch
        end
    end
end

function S = pairedBootstrap(a, b)
%  Bootstrap the paired difference by resampling SOURCES (each from an
%  independent simulation seed) — never by frame, which would treat temporally
%  correlated samples as independent.
S = struct('n',0,'meanDiff',NaN,'ci',[NaN NaN],'p',NaN);
if isempty(a) || isempty(b), return; end
good = isfinite(a) & isfinite(b);
a = a(good);  b = b(good);
S.n = numel(a);
if S.n < 5, return; end
d = a - b;
S.meanDiff = mean(d);
nb = 2000;
bd = nan(nb,1);
for k = 1:nb
    bd(k) = mean(d(randi(numel(d), numel(d), 1)));
end
S.ci = [prctile(bd,2.5) prctile(bd,97.5)];
S.p  = 2*min(mean(bd <= 0), mean(bd >= 0));
end

% =========================================================================
%  FIGURES
% =========================================================================

function makeArtifactQCFigures(qc, cfg, figDir)
f = figure('Name','Artifact regression QC','Position',[50 50 1400 800]);
tiledlayout(2,3,'TileSpacing','compact','Padding','compact');

ax = nexttile; imagesc(qc.varBefore); axis image off; colorbar;
title('Pixelwise variance BEFORE regression');  cl = clim(ax);
nexttile; imagesc(qc.varAfter); axis image off; colorbar; clim(cl);
title('Pixelwise variance AFTER (same scale)');
nexttile; imagesc(qc.varAfter ./ max(qc.varBefore, eps)); axis image off; colorbar;
title('Variance ratio after / before');

nexttile;
t = (1:numel(qc.globalRaw))/cfg.frameRate_Hz;
plot(t, qc.globalRaw, 'k'); hold on; plot(t, qc.globalFit, 'r', 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Mean F'); legend('raw','bleach fit'); title('Global trace');

nexttile;
loglog(qc.psdF, qc.psdVessel, 'r'); hold on; loglog(qc.psdF, qc.psdTissue, 'b');
xlabel('Hz'); ylabel('PSD'); legend('vessel-rich','vessel-poor');
title('Residual power spectra');

nexttile;
bar(qc.residualRegressorCorr);
set(gca,'XTick',1:numel(qc.regressorNames),'XTickLabel',qc.regressorNames,'XTickLabelRotation',45);
ylabel('corr(residual global trace, regressor)'); title('Residual vs regressors');

sgtitle('Artifact regression QC');
safeExport(f, fullfile(figDir,'fig01_artifact_QC.png'));
end

function makeRealMovieFigures(movCleanF, A, B, M, R, cfg, figDir)
mImg = mean(double(movCleanF),3);

% (2) overview with all footprints
f = figure('Name','Real movie overview','Position',[50 50 1500 700]);
tiledlayout(1,2,'TileSpacing','compact');
for s = 1:2
    if s == 1, S = A; else, S = B; end
    nexttile; imagesc(mImg); colormap(gca,gray); axis image off; hold on;
    for k = 1:size(S.bw,3)
        cc = bwboundaries(S.bw(:,:,k));
        for q = 1:numel(cc), plot(cc{q}(:,2), cc{q}(:,1), 'r', 'LineWidth', 0.7); end
    end
    title(sprintf('%s : %d kept (%d before area/boundary filter)', ...
        S.name, S.nAfter, S.nBefore), 'Interpreter','none');
end
sgtitle('All footprints on the mean image');
safeExport(f, fullfile(figDir,'fig02_real_overview.png'));

% (3) matched-component gallery
n = min(4, size(M.pairs,1));
if n > 0
    f = figure('Name','Matched components','Position',[50 50 1500 220*n]);
    tiledlayout(n, 4, 'TileSpacing','compact','Padding','compact');
    pad = 15;
    for k = 1:n
        i = M.pairs(k,1);  j = M.pairs(k,2);
        c  = mean([A.centroid_xy(i,:); B.centroid_xy(j,:)], 1);
        xr = max(1,round(c(1))-pad) : min(size(mImg,2), round(c(1))+pad);
        yr = max(1,round(c(2))-pad) : min(size(mImg,1), round(c(2))+pad);

        nexttile; imagesc(mImg(yr,xr)); colormap(gca,gray); axis image off;
        title(sprintf('mean image   d = %.2f \\mum', M.dist_um(k)),'FontSize',8);

        nexttile; imagesc(A.normFp(yr,xr,i)); axis image off; colormap(gca,hot);
        title(sprintf('extractGluSNFR3 #%d', i),'FontSize',8);

        nexttile; imagesc(B.normFp(yr,xr,j)); axis image off; colormap(gca,hot);
        title(sprintf('GIAnT #%d   IoU = %.2f', j, M.iou(k)),'FontSize',8);

        nexttile;
        tv = (0:numel(R.nativeTraces{k,1})-1)/cfg.frameRate_Hz;
        plot(tv, zscore2(R.nativeTraces{k,1}), 'b'); hold on;
        plot(tv, zscore2(R.nativeTraces{k,2}), 'r');
        plot(tv, zscore2(R.commonTraces{k}) - 8, 'k');
        legend({'extract','GIAnT','common mask'},'FontSize',6,'Location','best');
        title(sprintf('native r = %.2f', R.per.pearsonNative(k)),'FontSize',8);
        xlabel('s','FontSize',7);
    end
    sgtitle('Matched components: footprints, native traces and common-mask trace');
    safeExport(f, fullfile(figDir,'fig03_matched_gallery.png'));
end

% (4) blinded unmatched gallery
if ~isempty(R.blindedUnmatched.key)
    key = R.blindedUnmatched.key;
    n = size(key,1);
    f = figure('Name','Unmatched (blinded)','Position',[50 50 1200 170*ceil(n/4)]);
    tiledlayout(ceil(n/4), 4, 'TileSpacing','compact');
    pad = 15;
    for k = 1:n
        if key(k,1) == 1, S = A; else, S = B; end
        idx = key(k,2);
        c  = S.centroid_xy(idx,:);
        xr = max(1,round(c(1))-pad) : min(size(mImg,2), round(c(1))+pad);
        yr = max(1,round(c(2))-pad) : min(size(mImg,1), round(c(2))+pad);
        nexttile; imagesc(S.normFp(yr,xr,idx)); axis image off; colormap(gca,hot);
        title(sprintf('candidate %d', k),'FontSize',8);
    end
    sgtitle('Unmatched components, blinded (key in realMetrics.blindedUnmatched)');
    safeExport(f, fullfile(figDir,'fig04_unmatched_blinded.png'));
end
end

function z = zscore2(x)
x = double(x(:))';
g = isfinite(x);
if ~any(g), z = x; return; end
s = std(x(g));
if s == 0, z = x - mean(x(g)); else, z = (x - mean(x(g)))/s; end
end

function makeSimulationFigures(simRows, results, cfg, figDir)
if isempty(simRows), return; end
methods = ["extractGluSNFR3","GIAnT"];
cols = [0 0.4 0.9; 0.9 0.3 0];

f = figure('Name','Simulation curves','Position',[50 50 1500 850]);
tiledlayout(2,3,'TileSpacing','compact','Padding','compact');

% (5) detection probability vs dF/F
nexttile; hold on;
for m = 1:2
    sub = simRows(simRows.method == methods(m) & simRows.nSourcesTrue > 0, :);
    if isempty(sub), continue; end
    g = groupsummary(sub, 'amplitude_dFF', 'mean', 'detected');
    se = sqrt(max(g.mean_detected.*(1-g.mean_detected),0) ./ max(g.GroupCount,1));
    errorbar(g.amplitude_dFF, g.mean_detected, se, '-o','Color',cols(m,:),'LineWidth',1.5);
end
set(gca,'XScale','log'); xlabel('injected \DeltaF/F'); ylabel('detection probability');
ylim([0 1]); legend(cellstr(methods),'Location','southeast','Interpreter','none');
title('Detection probability');

% (6) F1 vs dF/F
nexttile; hold on;
for m = 1:2
    sub = simRows(simRows.method == methods(m) & simRows.nSourcesTrue > 0, :);
    if isempty(sub), continue; end
    amps = unique(sub.amplitude_dFF);
    f1 = nan(size(amps));
    for a = 1:numel(amps)
        s2 = sub(sub.amplitude_dFF == amps(a), :);
        pm = groupsummary(s2,'movie_id','max',{'TP','FP','FN'});
        TP = sum(pm.max_TP);  FP = sum(pm.max_FP);  FN = sum(pm.max_FN);
        rc = TP/max(TP+FN,1);  pr = TP/max(TP+FP,1);
        f1(a) = 2*rc*pr/max(rc+pr, eps);
    end
    plot(amps, f1, '-o','Color',cols(m,:),'LineWidth',1.5);
end
set(gca,'XScale','log'); xlabel('injected \DeltaF/F'); ylabel('F1'); ylim([0 1]);
title('F1 vs amplitude');

% (7) footprint IoU, with individual paired points
nexttile; hold on;
for m = 1:2
    sub = simRows(simRows.method == methods(m) & simRows.detected == 1, :);
    if isempty(sub), continue; end
    scatter(sub.amplitude_dFF .* (1+0.06*randn(height(sub),1)), sub.footprint_iou, ...
            6, cols(m,:), 'filled', 'MarkerFaceAlpha', 0.15);
    g = groupsummary(sub,'amplitude_dFF','median','footprint_iou');
    plot(g.amplitude_dFF, g.median_footprint_iou, '-o','Color',cols(m,:),'LineWidth',1.5);
end
set(gca,'XScale','log'); xlabel('injected \DeltaF/F'); ylabel('footprint IoU');
title('Spatial recovery');

% localization error
nexttile; hold on;
for m = 1:2
    sub = simRows(simRows.method == methods(m) & simRows.detected == 1, :);
    if isempty(sub), continue; end
    g = groupsummary(sub,'amplitude_dFF','median','centroid_error_um');
    plot(g.amplitude_dFF, g.median_centroid_error_um, '-o','Color',cols(m,:),'LineWidth',1.5);
end
set(gca,'XScale','log'); xlabel('injected \DeltaF/F'); ylabel('centroid error (\mum)');
title('Localization error');

% (8) trace correlation and amplitude bias
nexttile; hold on;
for m = 1:2
    sub = simRows(simRows.method == methods(m) & simRows.detected == 1, :);
    if isempty(sub), continue; end
    scatter(sub.amplitude_dFF .* (1+0.06*randn(height(sub),1)), sub.trace_corr, ...
            6, cols(m,:), 'filled', 'MarkerFaceAlpha', 0.15);
    g = groupsummary(sub,'amplitude_dFF','median','trace_corr');
    plot(g.amplitude_dFF, g.median_trace_corr, '-o','Color',cols(m,:),'LineWidth',1.5);
end
set(gca,'XScale','log'); xlabel('injected \DeltaF/F'); ylabel('trace corr with truth');
title('Temporal recovery');

nexttile; hold on;
for m = 1:2
    sub = simRows(simRows.method == methods(m) & simRows.detected == 1, :);
    if isempty(sub), continue; end
    g = groupsummary(sub,'amplitude_dFF','median','amplitude_bias');
    plot(g.amplitude_dFF, g.median_amplitude_bias, '-o','Color',cols(m,:),'LineWidth',1.5);
end
yline(0,'k--'); set(gca,'XScale','log');
xlabel('injected \DeltaF/F'); ylabel('(est - true)/true'); title('Amplitude bias');
sgtitle('Simulation: detection, spatial and temporal recovery');
safeExport(f, fullfile(figDir,'fig05_simulation_curves.png'));

% (9) stratification by background class
sub = simRows(simRows.dataset == "background" & simRows.nSourcesTrue > 0, :);
if height(sub) > 0
    f = figure('Name','Background stratification','Position',[50 50 1200 500]);
    g = groupsummary(sub, {'location_class','method'}, 'mean', 'detected');
    classes = unique(g.location_class);
    v = nan(numel(classes),2);
    for c = 1:numel(classes)
        for m = 1:2
            sel = g.location_class == classes(c) & g.method == methods(m);
            if any(sel), v(c,m) = g.mean_detected(find(sel,1)); end
        end
    end
    bar(v);
    set(gca,'XTick',1:numel(classes),'XTickLabel',cellstr(classes),'XTickLabelRotation',45);
    ylabel('detection probability'); legend(cellstr(methods),'Interpreter','none');
    title('Detection by background class');
    safeExport(f, fullfile(figDir,'fig09_background_stratification.png'));
end

% (10) runtime and component counts
f = figure('Name','Computational cost','Position',[50 50 1400 450]);
tiledlayout(1,3,'TileSpacing','compact');
mcat = categorical(simRows.method);
nexttile; boxchart(mcat, simRows.runtime_s);
ylabel('runtime per movie (s)');
title({'Runtime','(GIAnT uses a thread pool; extractGluSNFR3 is serial)'});
nexttile; boxchart(mcat, simRows.nDetected);
ylabel('# components kept'); title('Components after the shared filter');
nexttile; boxchart(mcat, simRows.nDetected_raw);
ylabel('# components returned'); title('Components before the shared filter');
safeExport(f, fullfile(figDir,'fig10_cost.png'));

% (11) individual paired points
P = results.paired;
if height(P) > 0
    f = figure('Name','Paired per-source','Position',[50 50 1300 420]);
    tiledlayout(1,3,'TileSpacing','compact');
    nexttile; scatter(P.iou_giant, P.iou_extract, 12, 'filled','MarkerFaceAlpha',0.3);
    hold on; plot([0 1],[0 1],'k--'); axis square;
    xlabel('GIAnT IoU'); ylabel('extractGluSNFR3 IoU'); title('Footprint IoU');
    nexttile; scatter(P.corr_giant, P.corr_extract, 12, 'filled','MarkerFaceAlpha',0.3);
    hold on; plot([-1 1],[-1 1],'k--'); axis square;
    xlabel('GIAnT trace corr'); ylabel('extract trace corr'); title('Trace correlation');
    nexttile; scatter(abs(P.ampbias_giant), abs(P.ampbias_extract), 12, 'filled','MarkerFaceAlpha',0.3);
    hold on; mx = max([1; abs(P.ampbias_giant); abs(P.ampbias_extract)]);
    plot([0 mx],[0 mx],'k--'); axis square;
    xlabel('|GIAnT amp bias|'); ylabel('|extract amp bias|'); title('Amplitude bias');
    sgtitle('Paired per-source comparison (same movie, same injected source)');
    safeExport(f, fullfile(figDir,'fig11_paired_points.png'));
end
end

function safeExport(fh, fname)
% R2026a figures follow the desktop theme, so on a dark desktop every exported
% panel comes out with black axes and near-invisible labels.  Force light.
try, theme(fh, 'light'); catch, end
try
    exportgraphics(fh, fname, 'Resolution', 150);
catch
    try
        print(fh, fname, '-dpng', '-r150', '-painters');
    catch ME
        warning('Figure export failed for %s : %s', fname, ME.message);
    end
end
end

% =========================================================================
%  SAVING AND REPORTING
% =========================================================================

function saveBenchmark(results, resExtract, resGIAnT, artifactModel, artifactQC, cfg)
outRoot = cfg.outRoot;
save(fullfile(outRoot,'results','benchmark_results.mat'), ...
    'results','artifactModel','artifactQC','cfg','-v7.3');
save(fullfile(outRoot,'results','raw_method_outputs.mat'), ...
    'resExtract','resGIAnT','-v7.3');
writeIf(results.detection, fullfile(outRoot,'results','detection_summary.csv'));
writeIf(results.logistic,  fullfile(outRoot,'results','detection_thresholds.csv'));
writeIf(results.paired,    fullfile(outRoot,'results','paired_per_source.csv'));
writeIf(results.simRows,   fullfile(outRoot,'results','simulation_rows.csv'));
writeIf(results.endToEnd.rows, fullfile(outRoot,'results','end_to_end_rows.csv'));
writeIf(results.endToEnd.attenuation, fullfile(outRoot,'results','signal_attenuation.csv'));
if ~isempty(results.realMetrics.per)
    writeIf(results.realMetrics.per, fullfile(outRoot,'results','real_matched_components.csv'));
end
fprintf('Saved all results under %s\n', fullfile(outRoot,'results'));
end

function writeIf(T, fname)
if istable(T) && ~isempty(T), writetable(T, fname); end
end

function printDecisionSummary(results)
fprintf('\n================ BENCHMARK SUMMARY ================\n');
disp(results.detection);
fprintf('\nDetection thresholds (dF/F at 50%% and 90%% detection, bootstrap CI):\n');
disp(results.logistic);
fprintf('\nPaired differences (extractGluSNFR3 - GIAnT), bootstrapped by source:\n');
fn = fieldnames(results.pairedStats);
for k = 1:numel(fn)
    s = results.pairedStats.(fn{k});
    fprintf('  %-10s n = %4d   mean diff %+0.4f   95%% CI [%+0.4f %+0.4f]   p = %.3f\n', ...
        fn{k}, s.n, s.meanDiff, s.ci(1), s.ci(2), s.p);
end
if isfield(results,'attenuationSummary')
    fprintf('\nPart C — fraction of the injected dF/F surviving artifact regression:\n');
    disp(results.attenuationSummary);
end
if ~isempty(results.failures)
    fprintf('\n%d run(s) failed and were preserved:\n', numel(results.failures));
    for k = 1:numel(results.failures)
        fprintf('  movie %d  %-16s : %s\n', results.failures(k).movie_id, ...
            results.failures(k).method, results.failures(k).message);
    end
end
fprintf(['\nNo single winner is declared (runbook section 13): read the operating\n' ...
         'curves in figures/fig05 together with the empty-movie false-positive rate\n' ...
         'in detection_summary.csv, and choose per biological question.\n']);
fprintf('==================================================\n');
end

% -------------------------------------------------------------------------
function [out, msg] = tryRun(fcn, movSize)
%TRYRUN  Run an adapter, preserving rather than hiding a failure
%(runbook section 11: "Preserve failed runs and report why they failed").
msg = '';
try
    out = fcn();
catch ME
    msg = ME.message;
    out = struct('name','failed', ...
                 'footprints', zeros(movSize(1), movSize(2), 0), ...
                 'traces',     zeros(0, movSize(3)), ...
                 'tracesDFF',  zeros(0, movSize(3)), ...
                 'params', struct(), 'runtime_s', NaN, ...
                 'extra', struct('error', getReport(ME,'extended','hyperlinks','off')), ...
                 'failed', true, 'message', msg);
    fprintf(2,'Run failed and was recorded: %s\n', msg);
end
end
