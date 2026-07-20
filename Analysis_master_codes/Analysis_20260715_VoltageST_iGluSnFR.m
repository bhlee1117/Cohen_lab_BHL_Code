%% Section 0 : Setup
clear
clc;

addpath('D:\Labmember\Data\ByungHun\VU_analysiscode');
addpath(genpath('C:\Users\Lab Member\Documents\GitHub\Cohen_lab_BHL_Code'));

% Remove any shadow of the built-in 'graph' (old cvx/sdpt3 defines graph.m,
% which breaks the centroid clustering in extractVoltronST / extractGluSNFR3).
gshadow = which('graph','-all');
for gi = 1:numel(gshadow)
    if ~contains(gshadow{gi}, fullfile('toolbox','matlab'))   % keep MATLAB's built-in
        rmpath(fileparts(gshadow{gi}));
        fprintf('Removed shadowing graph.m folder from path: %s\n', fileparts(gshadow{gi}));
    end
end
iswindow=0;
try
    [~, ~, raw] = xlsread('/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/Prism_V2+Glu_Data_Arrangement.xlsx', 'Sheet1', 'C5:AA31');
catch
    raw = readcell(macToWindowsPath('/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/Prism_V2+Glu_Data_Arrangement.xlsx'),...
        'Sheet', 'Sheet1', 'Range', 'C5:AA31');
    iswindow=1;
end

fpath           = raw(:,1)';
VRfpath           = raw(:,8)';
VRbw= (raw(:,9));
fpath_valid=cell2mat(cellfun(@(x) any(~ismissing(x)),fpath,'UniformOutput',false));
VRfpath_valid=cell2mat(cellfun(@(x) any(~ismissing(x)),VRfpath,'UniformOutput',false));
if iswindow
    for f=find(fpath_valid>0)
        fpath{f}=macToWindowsPath(fpath{f});        
    end
    for f=find(VRfpath_valid>0)
        VRfpath{f}=macToWindowsPath(VRfpath{f});        
    end
end
V2moviemaxTime  = 15000;   % voltage movie chunk length (frames)
GlumoviemaxTime = 3500;    % glutamate movie chunk length (frames)
o_laser         = 200001;  % 607/488 modulation onset (DAQ sample index, ~2 s after acq start)
mod488          = 200001;
foi             = [15 16 18];   % files of interest
set(0,'DefaultFigureWindowStyle','docked')

% --- extra columns + transform used by the f=3/f=5 sections below (cross-platform) ---
StructureData = raw(:,6);            % NB: paths in here also need macToWindowsPath on Windows
tformFile = "/Users/bhlee1117/Documents/GitHub/Cohen_lab_BHL_Code/16X_transformationMatrix.mat";
if iswindow, tformFile = macToWindowsPath(tformFile); end
load(tformFile);                     % tformReg
%% Map Glu Footprint onto morphology
for f=foi;
GluResult=loadGlu(fullfile(fpath{f},"Glu_Result.mat"));   % skips F0_glu/T_glu (unused, ~3.7 GB)
VoltResult=importdata(fullfile(fpath{f},"Volt_Result.mat"));
load(fullfile(fpath{f},"output_data.mat"));
% swcfiles=dir(fullfile(fpath{f}, 'Tracing*.swc'));
% if isempty(swcfiles)
% swcfiles=dir(fullfile(fileparts(StructureData{f}), 'Tracing*.swc'));
% end
sz=size(GluResult.AvgGluImg);

cam1_vsyn = Device_Data{1,2}.Counter_Inputs(1,1).data;
start_idx = find(cam1_vsyn == max(cam1_vsyn), 1);
end_idx   = numel(Device_Data{1,2}.buffered_tasks(1,2).channels(1,1).data);
segment_size = 103 * floor(Device_Data{3}.exposuretime/0.001);
n_to_add = end_idx - start_idx + 1;
vals = repelem((cam1_vsyn(start_idx-1)+1 : cam1_vsyn(start_idx-1)+ceil(n_to_add/segment_size))', segment_size);
cam1_vsyn(start_idx:end_idx) = vals(1:n_to_add);

cam2_vsyn = Device_Data{1, 2}.Counter_Inputs(1, 2).data;
cam2_vsyn_trig = find (diff(cam2_vsyn) == 1)+1;
segment_size2 = cam2_vsyn_trig (10) - cam2_vsyn_trig(9);
start_idx = min(find(cam2_vsyn ==max (cam2_vsyn)));
end_idx = length(Device_Data{1, 2}.buffered_tasks(1, 2).channels(1, 2).data);
last_val = cam2_vsyn(start_idx - 1);
n_to_add = end_idx - start_idx + 1;
n_segments = ceil(n_to_add / segment_size2);
added_part = repelem((last_val + 1 : last_val + n_segments)', segment_size2);
added_part = added_part(1:n_to_add);
cam2_vsyn(start_idx:end_idx) = added_part;

alignedVolt = transformCamera_O2B(Device_Data, tformReg, VoltResult.ref_im, GluResult.AvgGluImg);
alignedVoltftprnt =[];
for ft=1:size(VoltResult.ftprnt,3)
    alignedVoltftprnt(:,:,ft) = transformCamera_O2B(Device_Data, tformReg, VoltResult.ftprnt(:,:,ft), GluResult.AvgGluImg);
end
GluCoord=getGluCoord(GluResult);
figure(f); clf;
show_footprnt(alignedVoltftprnt,GluResult.AvgGluImg);

%loadVR data
if ~isfile(fullfile(fpath{f},'VR_Result.mat'))
fid = fopen(VRfpath{f});
VRdata = fread(fid,[12 inf],'double');
WorldITrack=find(VRdata(2,:)==1); %world 1
VRdata(5,WorldITrack)=(VRdata(5,WorldITrack)+6)*115/121;
VRdata(10,:)=bwlabel(VRdata(10,:));
VRdata=VRdata(:,VRdata(10,:)==VRbw{f});
DAQ_rate=Device_Data{1, 2}.Counter_Inputs.rate;
CamTrigger=find(cam1_vsyn(2:end)-cam1_vsyn(1:end-1));

t_DAQ=CamTrigger/DAQ_rate;
t_VR = datetime(datetime(VRdata(1,:),'ConvertFrom','datenum'), 'InputFormat', 'yyyy-MM-dd HH:mm:ss.SSS');
t_VR= t_VR-t_VR(1);
t_VR= milliseconds(t_VR)/1000;
t_VR= t_VR*t_DAQ(end)/t_VR(end);
VRdata(1,:)=t_VR;
[Virmen_data_int vel_trace]=virmen_interpolate(VRdata,115,t_DAQ);
Virmen_data_int(end+1,:)=vel_trace;
Virmen_data_int=Virmen_data_int(:,cam1_vsyn(o_laser)+1:cam1_vsyn(o_laser)+size(VoltResult.traces,2));
VR_Result=Virmen_data_int;
save(fullfile(fpath{f},'VR_Result.mat'),'VR_Result','-v7.3')
else
    load(fullfile(fpath{f},'VR_Result.mat'));
end

disp(['file #' num2str(f) '; GluResult, VoltResult, and VR data have loaded.'])
end
%% see the traces
vTax=VoltResult.tax(1:size(VoltResult.traces,2));
interactive2chViewer(GluResult.dFF_glu, getS_glu(GluResult), GluResult.t_ax(1:end-1), VoltResult.traces, alignedVoltftprnt, vTax)
%% =====================================================================
%  SPIKE-TRIGGERED-AVERAGE + PLACE-TUNING PIPELINE  (files f = 14:19)
%  Soma voltage (Voltron-ST) + dendritic glutamate (iGluSNFR) + treadmill VR.
%
%  Per file, in order, run the cells below:
%    S1  detect simple (SS) + complex (CS) spikes on the soma voltage
%    S2  STA voltage trace  (SS- and CS-triggered)
%    S3  STA voltage movie  (streamed from mc*.bin voltage chunks)
%    S4  STA glutamate movie + per-ROI STA glutamate trace (voltage-spike-triggered)
%    S5  glutamate place tuning vs treadmill position (PlaceTrigger_average)
%
%  Inputs expected in each fpath{f} (produced by Analysis_20260713_SD_VoltronST_plus_Glu.m):
%    Volt_Result.mat  -> Result: .traces (nROI x T), .t_ax, .ref_im, .ftprnt,
%                                 .VmovTimesegments, .mc
%    Glu_Result.mat   -> fields : .S_glu, .dFF_glu (Nsyn x Tg), .t_ax,
%                                 .AvgGluImg, .bvMask, .RegressComponentTile
%      After S0b, .S_glu is sparse [H*W x Nsyn] (+ .sz_glu, .GluCoord); read
%      footprints/centroids via getS_glu(GluResult) / getGluCoord(GluResult).
%    output_data.mat  -> Device_Data
%  Voltage bins  : mc%02d.bin  + mcTrace%02d.mat   (orange cam, Device_Data{3})
%  Glutamate bins: mc2%02d.bin + mc2Trace%02d.mat  (blue   cam, Device_Data{4})
%
%  Results are collected in STA_Result and saved as STA_Glu_Volt_Result.mat per file.
% =====================================================================
foi_STA = [15 16 18];   % from Section 0 (14:19)

% -- shared parameters --------------------------------------------------
% NB: o_laser, mod488 and GlumoviemaxTime are defined in Section 0 (cross-platform
%     setup) and are used as-is here — do NOT redefine them, or the glutamate tile
%     boundaries reconstructed in S4/S5 will not match how the movies were chunked.
staTauV   = [100 200];   % STA voltage window  [pre post] in voltage frames (~ms)
staTauG   = [15 15];     % STA glutamate window [pre post] in glutamate frames
matchTol  = 1;           % max |dt| (glu frames) allowed when mapping a V-spike to a glu frame

SS_thr    = [4 2];       % find_spike_bh [height prominence] (z-scored high-pass)
CS_onthr  = [5 1.5];     % detect_transient2 [on off] thresholds for complex-spike bursts
CS_minSp  = 3;           % >= this many spikes in a burst to call it a complex spike
CS_maxISI = 20;          % mean of first ISIs (frames) below this to call it a complex spike

% -- S3/S4 spike-triggered movie windows (ms) & raw-stack saving --------
staWin_ms.CS   = [300 500];   % [pre post] ms  -> CS-triggered window (glu & volt)
staWin_ms.SS   = [50 20];     % [pre post] ms  -> SS-triggered window (glu & volt)
isoSS_ms       = 20;          % isolated leading SS: require NO other SS within the preceding this-many ms
staStackOffset = 0;           % offset added before writing raw stacked .bin (raw camera counts are >=0 -> 0 is lossless)
stackFlushGB   = 3;           % flush the in-memory raw stack to a new .bin part above this size (GB)

% -- place-tuning parameters (S5) ---------------------------------------
place_bin     = 150;     % # of position bins along the track
vel_thresh    = 0.004;   % running threshold, in VR track-units / voltage-frame (VR_Result velocity row)
% VRtrackLength = one lap length in the SAME units as VR_Result row 5 (track position).
% Must match the lap_dist passed to virmen_interpolate when VR_Result was built (115).
VRtrackLength = 115;
SI_nShuffle   = 1000;    % # circular shuffles for the spatial-information p-value (S5)

% %% S0b : one-time maintenance -- compact Glu_Result (sparsify S_glu) for fast loading
% %  S_glu is a mostly-zero H x W x Nsyn footprint stack; storing it dense makes
% %  Glu_Result.mat huge and slow to load. This converts it to a sparse [H*W x Nsyn]
% %  matrix (+ sz_glu, + precomputed GluCoord) and re-saves. Run ONCE per file; the
% %  getS_glu / getGluCoord accessors read either the old (dense) or new (sparse) form.
% for f = foi_STA
%     fn = fullfile(fpath{f},'Glu_Result.mat');
%     G  = load(fn);                                   % top-level fields (saved with -struct)
%     if issparse(G.S_glu)
%         fprintf('file #%d: Glu_Result already compact.\n', f); continue;
%     end
%     G.sz_glu   = [size(G.S_glu,1) size(G.S_glu,2)];
%     G.GluCoord = get_coord(G.S_glu);                 % precompute centroids while still dense
%     G.S_glu    = sparse(double(reshape(G.S_glu, [], size(G.S_glu,3))));  % [H*W x Nsyn], sparse
%     save(fn, '-struct', 'G', '-v7.3');
%     fprintf('file #%d: S_glu sparsified (%.2f%% nonzero); Glu_Result re-saved.\n', ...
%         f, 100*nnz(G.S_glu)/numel(G.S_glu));
% end

%% S1i : Interactive spike + complex-spike detection w/ bleach correction (single file)
%  Curate ONE file at a time (set f below).  For each ROI:
%    (1) interactive_spike_finder on the z-scored high-pass trace
%           -> rough spikes + the [thr1 thr2] you set
%    (2) fit an exponential to the rough-spike HEIGHTS (photobleaching envelope)
%    (3) divide the trace by the (normalised) envelope so spike height is flat in time
%    (4) re-run find_spike_bh on the corrected trace with the SAME [thr1 thr2]
%    (5) interactive_cs_finder on the corrected trace -> CS trace + CS params
%    (6) you are asked whether this neuron is good to keep
%  All thresholds / params / bleach envelope / good flag are written back into
%  Volt_Result.mat, with SpClass kept in the (2 x nT x nROI) format S2-S5 expect.
f = 16;            % <-- set which file to curate
fprintf('=== S1i interactive spike/CS detection, file #%d ===\n', f);
VoltResult = importdata(fullfile(fpath{f},'Volt_Result.mat'));
[nROI, nT] = size(VoltResult.traces);

spike        = zeros(nROI,nT);
CStrace      = zeros(nROI,nT);
SpClass      = zeros(2,nT,nROI);
SpikeThr     = nan(nROI,2);      % [thr1 thr2] from interactive_spike_finder
CSThr        = nan(nROI,2);      % [thr1 thr2] from interactive_cs_finder
CS_NSpike    = nan(nROI,1);
CS_MaxAUC    = nan(nROI,1);
CS_NSpike2ISI= nan(nROI,1);
BleachEnv    = ones(nROI,nT);    % normalised photobleaching envelope (per ROI)
GoodNeuron   = false(nROI,1);
stopband=[[49 51];[184.5 186.5]];

% --- Regress the motion trace (VoltResult.mc) out of the raw traces, per 15000-frame chunk ---
% The movie was chunked at V2moviemaxTime (=15000) frames, so the concatenated traces are
% contiguous per chunk.  Regress the xy-mean motion out WITHIN each chunk separately, so
% chunk-boundary discontinuities are not mixed across segments. (SeeResiduals is a linear
% projection, so re-running this is ~idempotent.)
VoltResult.traces_filt=VoltResult.traces;
if isfield(VoltResult,'mc') && ~isempty(VoltResult.mc)
    mcTr = VoltResult.mc;
    if size(mcTr,1) ~= nT, mcTr = mcTr'; end        % ensure [T x 2]
    segEdges = [1:V2moviemaxTime:nT, nT+1];         % chunk boundaries (last chunk may be shorter)
    for s = 1:numel(segEdges)-1
        idx = segEdges(s):segEdges(s+1)-1;
        seg = reshape(VoltResult.traces(:,idx), nROI, 1, numel(idx));  % [nROI x 1 x t]
        seg = SeeResiduals(seg, mcTr(idx,:));                          % regress out motion
        seg = SeeResiduals(seg, mcTr(idx,:).^2);                          % regress out motion
        VoltResult.traces_filt(:,idx) = reshape(seg, nROI, numel(idx));
    end
    fprintf('Regressed motion (VoltResult.mc) out of traces over %d chunk(s) of <=%d frames.\n', ...
        numel(segEdges)-1, V2moviemaxTime);
else
    warning('VoltResult.mc not found or empty; skipping motion regression.');
end

for sb=1:size(stopband,1)
 VoltResult.traces_filt=get_bandstop(VoltResult.traces_filt,1000,stopband(sb,:));
end

for n = 1:nROI
    fprintf('--- ROI %d / %d ---\n', n, nROI);

    % z-scored high-pass trace (spikes are positive-going in .traces)
    tr_ref = VoltResult.traces(n,:) ./ get_threshold(VoltResult.traces(n,:),1);

    % (3) divide the trace by the envelope -> spike height flat in time
    tr_corr = tr_ref;

    % (4) re-detect on the corrected trace with the SAME thresholds you set
    [sp, thr1, thr2] = interactive_spike_finder(tr_corr - movmedian(tr_corr,500,2));
    sp=ind2vec(size(VoltResult.traces,2),sp,1);

    % (5) interactive complex-spike detection on the corrected trace
    [CS_tr, cs_thres, cs_params] = interactive_cs_finder(sp, tr_corr);

    % (5b) manual curation: add/remove spikes & mark complex-spike spans.
    %      Edits sp AND CS_tr together, so the CS onsets computed below stay in sync.
    [sp, CS_tr] = interactive_spike_checker(tr_corr, sp, CS_tr);

    % --- CS onsets = first spike of each burst (matches S1 / lab convention) ---
    bwCS   = bwlabel(CS_tr);
    CS_spk = sp .* bwCS;
    [~, CS_spk_t] = unique(CS_spk); CS_spk_t(CS_spk_t<6) = [];

    % (6) neuron QC
    ansQC  = questdlg(sprintf('ROI %d: keep this neuron?', n),'Neuron QC','Good','Bad','Good');
    isGood = strcmp(ansQC,'Good');

    % --- store (SpClass kept compatible with S2-S5) ---
    spike(n,:)      = sp;
    CStrace(n,:)    = CS_tr;
    SpClass(1,:,n)  = sp .* (1-CS_tr);          % simple spikes (bAPs)
    SpClass(2,CS_spk_t,n) = 1;                  % complex-spike onsets
    SpikeThr(n,:)   = [thr1 thr2];
    CSThr(n,:)      = cs_thres;
    CS_NSpike(n)    = cs_params.N_Spike;
    CS_MaxAUC(n)    = cs_params.Max_AUC;
    CS_NSpike2ISI(n)= cs_params.N_Spike2ISI;
    GoodNeuron(n)   = isGood;

    qcstr = {'BAD','GOOD'};
    fprintf('  ROI %d: %d SS, %d CS | thr=[%.2f %.2f] CSthr=[%.2f %.2f] -> %s\n', ...
        n, sum(SpClass(1,:,n)), sum(SpClass(2,:,n)), thr1, thr2, ...
        cs_thres(1), cs_thres(2), qcstr{isGood+1});
end

VoltResult.spike         = spike;
VoltResult.CStrace       = CStrace;
VoltResult.SpClass       = SpClass;
VoltResult.SpikeThr      = SpikeThr;
VoltResult.CSThr         = CSThr;
VoltResult.CS_NSpike     = CS_NSpike;
VoltResult.CS_MaxAUC     = CS_MaxAUC;
VoltResult.CS_NSpike2ISI = CS_NSpike2ISI;
VoltResult.BleachEnv     = BleachEnv;
VoltResult.GoodNeuron    = GoodNeuron;
Result = VoltResult;  save(fullfile(fpath{f},'Volt_Result.mat'),'Result','-v7.3');
fprintf('Saved interactive detection params for file #%d (%d/%d ROIs marked good).\n', ...
    f, sum(GoodNeuron), nROI);

%% S2 : STA voltage trace  (SS- and CS-triggered, f = 14:19)
for f = foi_STA
    fprintf('=== S2 STA voltage trace, file #%d ===\n', f);
    VoltResult = importdata(fullfile(fpath{f},'Volt_Result.mat'));
    [nROI, nT] = size(VoltResult.traces);
    tr = VoltResult.traces - movmedian(VoltResult.traces,1000,2);
    tr = tr ./ get_threshold(tr,1);

    tV = (-staTauV(1):staTauV(2)) / 1000;           % s, assuming ~1 kHz voltage
    STA_V_SS = nan(nROI,sum(staTauV)+1);
    STA_V_CS = nan(nROI,sum(staTauV)+1);
    for n = find(VoltResult.GoodNeuron)'
        ssIdx = find(VoltResult.SpClass(1,:,n) > 0);
        csIdx = find(VoltResult.SpClass(2,:,n) > 0);
        s = get_STA(tr(n,:), ssIdx, staTauV(1), staTauV(2)); if ~isempty(s), STA_V_SS(n,:) = s; end
        c = get_STA(tr(n,:), csIdx, staTauV(1), staTauV(2)); if ~isempty(c), STA_V_CS(n,:) = c; end
    end

    STA_Result = struct();
    STA_Result.tV = tV;
    STA_Result.STA_V_SS = STA_V_SS;
    STA_Result.STA_V_CS = STA_V_CS;
    save(fullfile(fpath{f},'STA_Glu_Volt_Result.mat'),'-struct','STA_Result','-v7.3');

    figure(320+f); clf;
    goodROI = find(VoltResult.GoodNeuron)';
    if isempty(goodROI), warning('File %d: no good ROIs to plot.',f); else
    tiledlayout('flow','TileSpacing','compact','Padding','compact');
    for n = goodROI
        nexttile; hold on;
        plot(tV, STA_V_SS(n,:),'k'); plot(tV, STA_V_CS(n,:),'r');
        xlabel('Time from spike (s)'); ylabel('V (\sigma)');
        title(sprintf('ROI %d',n)); axis tight;
    end
    sgtitle(sprintf('File %d STA voltage: black=SS red=CS',f));
    end
    drawnow;
end

%% S3 : STA voltage movie (avg) + raw stacked event windows (per good neuron; SS & CS)
%  Per good neuron: SS = isolated leading simple spikes (no SS within preceding isoSS_ms),
%  CS = complex-spike onsets.  Windows: SS = staWin_ms.SS, CS = staWin_ms.CS (ms -> volt frames).
%  Saves STA_V_movie_Result.mat with:
%    STA_V_mov(k)  : .roi .cls .win .n .mov (H x W x winLen, averaged, spike-positive residual)
%    STA_V_tr(k)   : .roi .cls .win .trace (nROI x winLen, from VoltResult.traces)
%    STAinfo_V(k)  : raw-stack metadata (.centerFrames global volt frame of every stored window,
%                    .binFiles, .evPerFile, .win, .H .W .pixels .offset) -> use loadSTAevent() to reload.
for f = foi_STA
    fprintf('=== S3 STA voltage movie, file #%d ===\n', f);
    VoltResult = importdata(fullfile(fpath{f},'Volt_Result.mat'));
    load(fullfile(fpath{f},'output_data.mat'));                 % Device_Data
    sz   = double(Device_Data{3}.ROI([2 4]));                   % [nCol nRow]
    H = sz(2); W = sz(1);
    dtV  = median(diff(VoltResult.t_ax));                       % s/frame (voltage), from precise time axis
    readFrames = diff(VoltResult.VmovTimesegments);
    nSeg = numel(readFrames);
    flushBytes = stackFlushGB*1e9;

    goodROI = find(VoltResult.GoodNeuron(:))';
    if isempty(goodROI), warning('File %d: no good neurons; skipping S3.',f); continue; end

    % windows in voltage frames (ms -> frames via VoltResult.t_ax)
    preSS = ms2frames(staWin_ms.SS(1),dtV); postSS = ms2frames(staWin_ms.SS(2),dtV);
    preCS = ms2frames(staWin_ms.CS(1),dtV); postCS = ms2frames(staWin_ms.CS(2),dtV);
    isoV  = ms2frames(isoSS_ms,dtV);

    % build stackers: for each good neuron -> SS(isolated leading) + CS(onset)
    clear EV; e = 0;
    for n = goodROI
        ssAll = find(VoltResult.SpClass(1,:,n) > 0);
        if isempty(ssAll), ssLead = ssAll; else, ssLead = ssAll([true, diff(ssAll) > isoV]); end
        csOn  = find(VoltResult.SpClass(2,:,n) > 0);
        e=e+1; EV(e)=newStacker(n,1,[preSS postSS],ssLead,H,W);
        e=e+1; EV(e)=newStacker(n,2,[preCS postCS],csOn ,H,W);
    end

    % STA voltage TRACE (all ROIs) per neuron/class, with the new windows
    clear STA_V_tr;
    for kk=1:numel(EV)
        STA_V_tr(kk).roi=EV(kk).roi; STA_V_tr(kk).cls=clsName(EV(kk).cls); STA_V_tr(kk).win=EV(kk).win;
        STA_V_tr(kk).trace = get_STA(VoltResult.traces, EV(kk).ev, EV(kk).win(1), EV(kk).win(2));
    end

    % stream chunks once; feed every stacker
    off = 0;
    for j = 1:nSeg
        hasEv=false;
        for kk=1:numel(EV), if any(EV(kk).ev>off & EV(kk).ev<=off+readFrames(j)), hasEv=true; break; end, end
        if ~hasEv, off=off+readFrames(j); continue; end
        [movR, mov_mc] = loadVoltChunk(fpath{f}, j, sz, readFrames(j)); % processed (spike-pos), raw counts
        nT = size(movR,3); movVecP = tovec(movR); movVecR = tovec(mov_mc);
        for kk=1:numel(EV)
            locEv = EV(kk).ev(EV(kk).ev>off & EV(kk).ev<=off+nT) - off;
            if isempty(locEv), continue; end
            [~, AddMovP, AddT] = get_STA(movVecP, locEv, EV(kk).win(1), EV(kk).win(2));  % processed (for avg)
            if isempty(AddT), continue; end
            [~, AddMovR]       = get_STA(movVecR, AddT,  EV(kk).win(1), EV(kk).win(2));  % raw counts (for stack)
            EV(kk) = accumStack(EV(kk), AddMovP, AddMovR, AddT+off, fpath{f}, 'STA_V', staStackOffset, flushBytes);
        end
        off = off + nT;
        fprintf('  chunk %2d/%2d\n', j, nSeg);
    end

    % finalize: flush remaining buffers, average, collect metadata
    clear STAinfo_V STA_V_mov
    for kk=1:numel(EV)
        EV(kk) = flushStacker(EV(kk), fpath{f}, 'STA_V', staStackOffset);
        STA_V_mov(kk).roi=EV(kk).roi; STA_V_mov(kk).cls=clsName(EV(kk).cls);
        STA_V_mov(kk).win=EV(kk).win; STA_V_mov(kk).n=EV(kk).cnt;
        STA_V_mov(kk).mov = EV(kk).sum ./ max(EV(kk).cnt,1);
        STAinfo_V(kk) = stackInfo(EV(kk), staStackOffset);
    end

    save(fullfile(fpath{f},'STA_V_movie_Result.mat'),'STAinfo_V','STA_V_mov','STA_V_tr','staWin_ms','-v7.3');
    fprintf('  saved STA_V_movie_Result.mat: %d stackers, n=%s\n', numel(EV), mat2str([STA_V_mov.n]));

    figure(340+f); clf; tiledlayout('flow','TileSpacing','compact');
    for kk=1:numel(EV)
        mv = STA_V_mov(kk).mov; base = median(mv(:,:,1:max(1,round(EV(kk).win(1)/2))),3);
        nexttile; imshow2(max(mv-base,[],3),[]);
        title(sprintf('ROI%d %s n=%d',STA_V_mov(kk).roi,STA_V_mov(kk).cls,STA_V_mov(kk).n));
    end
    sgtitle(sprintf('File %d STA voltage movie (peak dF)',f)); drawnow;
end

%% S4 : STA glutamate movie (avg) + raw stacked windows + STA glu trace (per good neuron; SS & CS)
%  Triggers are the SAME soma V-spike sets as S3 (SS = isolated leading, CS = onset), mapped
%  to the nearest glutamate frame.  Windows in glu frames from staWin_ms.  Saves
%  STA_Glu_movie_Result.mat with STA_glu_mov / STA_glu_tr / STAinfo_Glu (see S3 header;
%  STAinfo_Glu.centerFrames are GLUTAMATE frames).  Reload single events with loadSTAevent().
for f = foi_STA
    fprintf('=== S4 STA glutamate, file #%d ===\n', f);
    VoltResult = importdata(fullfile(fpath{f},'Volt_Result.mat'));
    GluResult  = loadGlu(fullfile(fpath{f},'Glu_Result.mat'));   % skips F0_glu/T_glu (unused, ~3.7 GB)
    load(fullfile(fpath{f},'output_data.mat'));                 % Device_Data
    nCol2 = double(Device_Data{4}.ROI(2)); nRow2 = double(Device_Data{4}.ROI(4));
    dtV   = median(diff(VoltResult.t_ax));         % s/frame voltage   (SS isolation) - precise time axis
    dtG   = median(diff(GluResult.t_ax(1:size(GluResult.dFF_glu,2)))); % s/frame glutamate (windows)
    exp2  = Device_Data{4}.exposuretime;           % nominal glu exposure (only for the 5-s low-pass, to match extraction)
    flushBytes = stackFlushGB*1e9;

    goodROI = find(VoltResult.GoodNeuron(:))';
    if isempty(goodROI), warning('File %d: no good neurons; skipping S4.',f); continue; end

    % ---- reconstruct glutamate (blue) camera clock -> tile boundaries ----
    cam2_vsyn = Device_Data{1,2}.Counter_Inputs(1,2).data;
    cam2_vsyn_trig = find(diff(cam2_vsyn)==1)+1;
    segment_size2  = cam2_vsyn_trig(10) - cam2_vsyn_trig(9);
    start_idx = min(find(cam2_vsyn==max(cam2_vsyn)));
    end_idx   = length(Device_Data{1,2}.buffered_tasks(1,2).channels(1,2).data);
    last_val  = cam2_vsyn(start_idx-1); n_to_add = end_idx-start_idx+1;
    added = repelem((last_val+1 : last_val+ceil(n_to_add/segment_size2))', segment_size2);
    cam2_vsyn(start_idx:end_idx) = added(1:n_to_add);
    if cam2_vsyn(end) < GlumoviemaxTime
        GluSeg = [(cam2_vsyn(mod488)+2) cam2_vsyn(end)];
    else
        GluSeg = (cam2_vsyn(mod488)+2) : GlumoviemaxTime : cam2_vsyn(end);
        GluSeg(end+1) = cam2_vsyn(end);
    end
    Ntile   = numel(GluSeg)-1;
    tileLen = diff(GluSeg);
    nTg     = size(GluResult.dFF_glu,2);

    % windows in glutamate frames (via GluResult.t_ax); isolation window in voltage frames (via VoltResult.t_ax)
    preSS = ms2frames(staWin_ms.SS(1),dtG); postSS = ms2frames(staWin_ms.SS(2),dtG);
    preCS = ms2frames(staWin_ms.CS(1),dtG); postCS = ms2frames(staWin_ms.CS(2),dtG);
    isoV  = ms2frames(isoSS_ms,dtV);

    % map voltage spikes -> nearest glutamate frame (both t_ax in seconds; keep |gap| <= matchTol glu frames)
    gtax = GluResult.t_ax(1:nTg);
    map2glu = @(vf) mapV2Glu(VoltResult.t_ax(vf), gtax, matchTol*dtG);

    % build stackers per good neuron (events are glu frames)
    clear EV; e=0;
    for n = goodROI
        ssAll = find(VoltResult.SpClass(1,:,n) > 0);
        if isempty(ssAll), ssLead=ssAll; else, ssLead = ssAll([true, diff(ssAll) > isoV]); end
        csOn  = find(VoltResult.SpClass(2,:,n) > 0);
        e=e+1; EV(e)=newStacker(n,1,[preSS postSS],map2glu(ssLead),nRow2,nCol2);
        e=e+1; EV(e)=newStacker(n,2,[preCS postCS],map2glu(csOn) ,nRow2,nCol2);
    end

    % STA glutamate TRACE (per synapse) per neuron/class
    clear STA_glu_tr;
    for kk=1:numel(EV)
        STA_glu_tr(kk).roi=EV(kk).roi; STA_glu_tr(kk).cls=clsName(EV(kk).cls); STA_glu_tr(kk).win=EV(kk).win;
        STA_glu_tr(kk).t     = (-EV(kk).win(1):EV(kk).win(2))*dtG;
        STA_glu_tr(kk).trace = get_STA(GluResult.dFF_glu, EV(kk).ev, EV(kk).win(1), EV(kk).win(2));
    end

    % stream tiles once; reprocess (bleach fit, stored PCs, motion, bvMask); feed stackers
    bvMask = max(GluResult.bvMask,[],3)==0;
    off = 0;
    for k = 1:Ntile
        hasEv=false;
        for kk=1:numel(EV), if any(EV(kk).ev>off & EV(kk).ev<=off+tileLen(k)), hasEv=true; break; end, end
        if ~hasEv, off=off+tileLen(k); continue; end

        mov2 = double(readBinMov_times(fullfile(fpath{f}, ['mc2' num2str(k,'%02d') '.bin']), nRow2, nCol2, 1:tileLen(k)));
        load(fullfile(fpath{f}, ['mc2Trace' num2str(k,'%02d') '.mat']));   % mcTrace2
        meanF = squeeze(mean(mov2,[1 2]));
        y_fit = expfitDM_2((1:size(mov2,3))', meanF, (1:size(mov2,3))', [10000 100]);
        mov_f = mov2 ./ reshape(y_fit,1,1,[]);
        mov_lw= movmedian(mov_f, round(5/exp2), 3);
        mov_s = mov_f - mov_lw;
        pcatr2regress = tovec(mov_s)'*tovec(GluResult.RegressComponentFootprint);
        mov_s = SeeResiduals(mov_s, pcatr2regress);
        mov_s = SeeResiduals(mov_s, mcTrace2.xymean(1:tileLen(k),:));
        mov_s = SeeResiduals(mov_s, mcTrace2.xymean(1:tileLen(k),:).^2);
        mov_s = SeeResiduals(mov_s, mcTrace2.xymean(1:tileLen(k),1).*mcTrace2.xymean(1:tileLen(k),end));
        mov_s = mov_s .* bvMask;
        movVecP = tovec(mov_s);        % processed (for average)
        movVecR = tovec(mov2);         % raw counts (for stack)

        for kk=1:numel(EV)
            locEv = EV(kk).ev(EV(kk).ev>off & EV(kk).ev<=off+tileLen(k)) - off;
            if isempty(locEv), continue; end
            [~, AddMovP, AddT] = get_STA(movVecP, locEv, EV(kk).win(1), EV(kk).win(2));
            if isempty(AddT), continue; end
            [~, AddMovR]       = get_STA(movVecR, AddT,  EV(kk).win(1), EV(kk).win(2));
            EV(kk) = accumStack(EV(kk), AddMovP, AddMovR, AddT+off, fpath{f}, 'STA_Glu', staStackOffset, flushBytes);
        end
        off = off + tileLen(k);
        fprintf('  glu tile %2d/%2d\n', k, Ntile);
    end

    % finalize
    clear STAinfo_Glu STA_glu_mov
    for kk=1:numel(EV)
        EV(kk) = flushStacker(EV(kk), fpath{f}, 'STA_Glu', staStackOffset);
        STA_glu_mov(kk).roi=EV(kk).roi; STA_glu_mov(kk).cls=clsName(EV(kk).cls);
        STA_glu_mov(kk).win=EV(kk).win; STA_glu_mov(kk).n=EV(kk).cnt;
        STA_glu_mov(kk).mov = EV(kk).sum ./ max(EV(kk).cnt,1);
        STAinfo_Glu(kk) = stackInfo(EV(kk), staStackOffset);
    end

    save(fullfile(fpath{f},'STA_Glu_movie_Result.mat'),'STAinfo_Glu','STA_glu_mov','STA_glu_tr','staWin_ms','-v7.3');
    fprintf('  saved STA_Glu_movie_Result.mat: %d stackers, n=%s\n', numel(EV), mat2str([STA_glu_mov.n]));

    figure(360+f); clf; tiledlayout('flow','TileSpacing','compact');
    for kk=1:numel(EV)
        nexttile; imshow2(max(STA_glu_mov(kk).mov,[],3),[]);
        title(sprintf('ROI%d %s n=%d',STA_glu_mov(kk).roi,STA_glu_mov(kk).cls,STA_glu_mov(kk).n));
    end
    sgtitle(sprintf('File %d STA glutamate movie (peak)',f)); drawnow;
end

%% S4b : Voltage-spike-triggered glutamate STA (-100..100 ms) + pre/post-spike tuning significance
%  Per good neuron, for SS and CS: STA of every synapse's dF/F around the soma spikes, plus a
%  1000x circular-shuffle test of the glutamate in TWO windows (staGlu_prewin):
%    row 1 = -50..0 ms  PRE-spike  (is glu elevated just BEFORE the spike?)
%    row 2 =   0..50 ms POST-spike (is glu elevated just AFTER  the spike?)
%  score/pval are returned Nsyn x 2 (col 1 = pre, col 2 = post). Results are saved as VtrigGlu
%  in STA_Glu_Volt_Result.mat and reused by S7 (S7 does NOT recompute).
%  Triggers: SS = ISOLATED simple spike (NO other SS spike AND NO complex-spike (CS) burst
%  in the preceding VtrigIso_ms, so the pre-spike window is clean); CS = complex-spike onset.
VtrigWin_ms   = [100 200];       % STA window [pre post] ms
VtrigScore_ms = [-50 0; 0 50];   % scoring windows (ms): row1 = pre-spike, row2 = post-spike
VtrigIso_ms   = 50;             % SS isolation: require no other SS / CS-burst in the preceding this-many ms
for f = foi_STA
    fprintf('=== S4b V-spike-triggered glu tuning, file #%d ===\n', f);
    VoltResult = importdata(fullfile(fpath{f},'Volt_Result.mat'));
    GluResult  = loadGlu(fullfile(fpath{f},'Glu_Result.mat'));   % skips F0_glu/T_glu (unused, ~3.7 GB)
    nTg  = size(GluResult.dFF_glu,2);
    dtG  = median(diff(GluResult.t_ax(1:nTg)));
    dtV  = median(diff(VoltResult.t_ax));                           % s/frame voltage (for SS isolation)
    isoW = ms2frames(VtrigIso_ms, dtV);                            % isolation window (frames) before each SS
    gtax = GluResult.t_ax(1:nTg);
    GluCoord = getGluCoord(GluResult);  AvgImg = GluResult.AvgGluImg;  Nsyn = size(GluResult.dFF_glu,1);
    goodROI = find(VoltResult.GoodNeuron(:))';
    classes = {'SS',1; 'CS',2};

    % voltage footprints aligned into the glutamate image space (for the soma contour)
    load(fullfile(fpath{f},'output_data.mat'));                     % Device_Data
    alignedVoltftprnt = zeros([size(AvgImg,1) size(AvgImg,2) size(VoltResult.ftprnt,3)]);
    for ft = 1:size(VoltResult.ftprnt,3)
        alignedVoltftprnt(:,:,ft) = transformCamera_O2B(Device_Data, tformReg, VoltResult.ftprnt(:,:,ft), AvgImg);
    end

    clear VtrigGlu
    for gi = 1:numel(goodROI)
        n = goodROI(gi); VtrigGlu(gi).roi = n;
        for ci = 1:2
            vfr  = find(VoltResult.SpClass(classes{ci,2},:,n) > 0);           % SS or CS-onset frames
            if ci==1 && ~isempty(vfr)   % SS: keep only spikes with NO other SS and NO CS burst in the preceding isoW
                blocker = double(VoltResult.SpClass(1,:,n)>0 | VoltResult.CStrace(n,:)>0);  % 1 x nT : other SS or CS trace
                cb  = [0, cumsum(blocker)];                                   % prefix sum -> fast window counts
                lo  = max(1, vfr - isoW);                                     % window (s-isoW .. s-1), excludes s itself
                cnt = cb(vfr) - cb(lo);                                       % # SS/CS-trace frames strictly before s
                vfr = vfr(cnt==0);                                           % keep only fully isolated SS
            end
            trig = mapV2Glu(VoltResult.t_ax(vfr), gtax, matchTol*dtG);        % -> glu frames
            [STA, tMs, score, pval] = staGlu_prewin(GluResult.dFF_glu, trig, dtG, VtrigWin_ms, VtrigScore_ms, SI_nShuffle);
            VtrigGlu(gi).(classes{ci,1}) = struct('ntrig',numel(trig),'STA',STA,'tMs',tMs,'score',score,'pval',pval);
            fprintf('  ROI %d %s: %d triggers | pre-spike tuned %d/%d, post-spike tuned %d/%d (p<0.05)\n', ...
                n, classes{ci,1}, numel(trig), sum(pval(:,1)<0.05), Nsyn, sum(pval(:,2)<0.05), Nsyn);
        end
    end

    STA_Result = importdata(fullfile(fpath{f},'STA_Glu_Volt_Result.mat'));
    STA_Result.VtrigGlu          = VtrigGlu;
    STA_Result.VtrigGlu_win_ms   = VtrigWin_ms;
    STA_Result.VtrigGlu_score_ms = VtrigScore_ms;
    STA_Result.VtrigGlu_nShuffle = SI_nShuffle;
    save(fullfile(fpath{f},'STA_Glu_Volt_Result.mat'),'-struct','STA_Result','-v7.3');

    plot_VtrigGlu(VtrigGlu, GluCoord, AvgImg, VtrigScore_ms, 400+f, alignedVoltftprnt);   % SS/CS heatmap + pre/post-spike maps + soma contour
end

%% S5 : Glutamate place tuning vs treadmill position  (VR_Result-based, f = 14:19)
for f = foi_STA(3)
    fprintf('=== S5 glutamate place tuning, file #%d ===\n', f); tS5 = tic;
    GluResult  = loadGlu(fullfile(fpath{f},'Glu_Result.mat'));   % skips F0_glu/T_glu (unused, ~3.7 GB)
    VoltResult = importdata(fullfile(fpath{f},'Volt_Result.mat'));   % only for the voltage time axis
    VR_Result  = importdata(fullfile(fpath{f},'VR_Result.mat'));     % Virmen matrix on the VOLTAGE timebase
    nTg = size(GluResult.dFF_glu,2);

    % ---- resample VR (voltage timebase) onto the glutamate frames ----
    % VR_Result(:,i) is voltage frame i (time VoltResult.t_ax(i)); map each glutamate
    % frame to its nearest voltage frame in time (both t_ax in seconds).
    % Rows used by PlaceTrigger_average: 5 = track position, 8 = lap, end = velocity.
    vtax   = VoltResult.t_ax(1:size(VR_Result,2));
    vIdx   = match_nearest(GluResult.t_ax(1:nTg), vtax);            % nTg x 1 index into voltage frames
    Virmen = VR_Result(:, vIdx);                                    % VR rows resampled to glutamate frames

    % ---- glutamate events (per synapse) ----
    GluSpike = find_spike_bh(zscore(GluResult.dFF_glu,0,2), 3, 1);

    % ---- lap x position maps: event rate, and mean dF/F ----
    [Lap_FR, Lap_SpkN, Lap_Nvalid] = PlaceTrigger_average(GluSpike,          place_bin, Virmen, vel_thresh, VRtrackLength);
    [Lap_dFF]                      = PlaceTrigger_average(GluResult.dFF_glu,  place_bin, Virmen, vel_thresh, VRtrackLength);

    % ---- GLUTAMATE spatial information + circular-shuffle p-value (running frames only) ----
    binTr = ceil(Virmen(5,:)/(VRtrackLength/place_bin));
    binTr(Virmen(end,:) < vel_thresh) = NaN;                       % drop non-running (velocity = last row)
    Nsyn  = size(GluResult.dFF_glu,1);
    [SI, SI_p] = si_pvalue(GluSpike, binTr, GluResult.t_ax, SI_nShuffle, 0, sprintf('f%d glu (%d syn)',f,Nsyn));
    fprintf('  glutamate: %d/%d synapses spatially tuned (p<0.05)\n', sum(SI_p<0.05), Nsyn);

    % ---- VOLTAGE spatial information + p-value (good neurons; VR is already on the volt timebase) ----
    goodROI = find(VoltResult.GoodNeuron(:))';
    Tvr    = min(size(VR_Result,2), size(VoltResult.spike,2));      % align lengths 1:1
    VRv    = VR_Result(:,1:Tvr);
    spikeV = double(VoltResult.spike(goodROI,1:Tvr) > 0);          % nGood x Tvr
    binTrV = ceil(VRv(5,:)/(VRtrackLength/place_bin));
    binTrV(VRv(end,:) < vel_thresh) = NaN;
    [Lap_FR_V, Lap_SpkN_V, Lap_Nvalid_V] = PlaceTrigger_average(spikeV, place_bin, VRv, vel_thresh, VRtrackLength);
    [SI_V, SI_p_V] = si_pvalue(spikeV, binTrV, VoltResult.t_ax, SI_nShuffle, 0, sprintf('f%d volt (%d roi)',f,numel(goodROI)));
    fprintf('  voltage:   %d/%d good neurons spatially tuned (p<0.05)\n', sum(SI_p_V<0.05), numel(goodROI));

    STA_Result = importdata(fullfile(fpath{f},'STA_Glu_Volt_Result.mat'));
    STA_Result.Virmen     = Virmen;
    STA_Result.Lap_FR     = Lap_FR;      % nLap x place_bin x Nsyn (glu event rate)
    STA_Result.Lap_dFF    = Lap_dFF;     % nLap x place_bin x Nsyn (mean dF/F)
    STA_Result.Lap_SpkN   = Lap_SpkN;
    STA_Result.Lap_Nvalid = Lap_Nvalid;
    STA_Result.SpatialInfo   = SI;                   % glutamate, per synapse
    STA_Result.SpatialInfo_p = SI_p;                 % glutamate circular-shuffle p-value
    STA_Result.V_goodROI       = goodROI;            % voltage neurons (ROI indices) for rows below
    STA_Result.V_SpatialInfo   = SI_V;               % voltage, per good neuron
    STA_Result.V_SpatialInfo_p = SI_p_V;             % voltage circular-shuffle p-value
    STA_Result.V_Lap_FR        = Lap_FR_V;           % nLap x place_bin x nGood (voltage event rate)
    STA_Result.V_Lap_SpkN      = Lap_SpkN_V;
    STA_Result.V_Lap_Nvalid    = Lap_Nvalid_V;
    STA_Result.SI_nShuffle   = SI_nShuffle;
    STA_Result.place_bin  = place_bin;
    STA_Result.VRtrackLength = VRtrackLength; 
    save(fullfile(fpath{f},'STA_Glu_Volt_Result.mat'),'-struct','STA_Result','-v7.3');

    GluCoord=getGluCoord(GluResult);
    % population place map (mean over laps) sorted by peak position
    popMap = squeeze(mean(Lap_dFF(:,:,:),1,'omitnan'))';               % Nsyn x place_bin
    [~,pk] = max(ringmovMean(popMap,5),[],2); [~,ord] = sort(pk);
    pm  = popMap(ord,:);
    pmz = (pm - mean(pm,2,'omitnan')) ./ std(pm,0,2,'omitnan');    % per-synapse z-score across position
    
    figure(380+f); clf; tiledlayout(2,2);
    nexttile; imagesc(pmz,'AlphaData',~isnan(pmz)); axis tight; xlabel('Position bin'); ylabel('synapse (sorted)');
    title(sprintf('File %d glu place map (peak-sorted)',f)); colorbar;
    nexttile; hold on;
    histogram(SI(SI_p>=0.05),30,'FaceColor',[.6 .6 .6]);
    histogram(SI(SI_p<0.05), 30,'FaceColor',[.85 .2 .2]);
    xlabel('Spatial information (bits/event)'); ylabel('# synapses');
    title(sprintf('Glu tuning: %d/%d sig (p<0.05)',sum(SI_p<0.05),Nsyn));
    nexttile; hold on;                                         % voltage SI per good neuron
    sigV = SI_p_V<0.05;
    bar(find(~sigV), SI_V(~sigV), 'FaceColor',[.6 .6 .6]);
    bar(find(sigV),  SI_V(sigV),  'FaceColor',[.85 .2 .2]);
    set(gca,'XTick',1:numel(goodROI),'XTickLabel',goodROI);
    xlabel('Voltage neuron (ROI)'); ylabel('Spatial info (bits/spike)');
    title(sprintf('Voltage tuning: %d/%d sig',sum(sigV),numel(goodROI))); drawnow;
   
    nexttile;
    image(repmat(mat2gray(double(GluResult.AvgGluImg),[30 300]),1,1,3)); hold on; axis image off;   % avg image (truecolor bg)
    cvals = -log10(SI_p(:));                                        % -log10 p-value per synapse
    scatter(GluCoord(:,1),GluCoord(:,2),2,cvals,'filled');
    colormap(gca,'turbo');
    cmax = max([prctile(cvals(isfinite(cvals)),99), -log10(0.05)]); if ~isfinite(cmax)||cmax<=0, cmax=3; end
    caxis([0 cmax]);
    pt = [0.5 0.05 0.01 0.001]; tp = -log10(pt); keep = tp<=cmax+eps;   % show p-value on the colorbar
    cb = colorbar; cb.Ticks = tp(keep); cb.TickLabels = arrayfun(@(x) sprintf('%.3g',x), pt(keep), 'uni', 0);
    cb.Label.String = 'SI p-value (shuffle)';
    title('glu synapses colored by SI p-value');
    fprintf('=== S5 file #%d done in %.1f s ===\n', f, toc(tS5));
end

%% S6b : Interactive VOLTAGE + GLUTAMATE ROI browser (single file: set f)
%  Top: glu image + glu synapses (dots, -log10 SI_p) + voltage footprint CONTOURS (-log10 V_SI_p).
%  Click a synapse -> glu row (trace, dF/F map, GluSpike map); click a footprint -> voltage row
%  (trace, place map, mean tuning). Needs S5 outputs and tformReg (from Section 0).
f = foi_STA(3);            % <-- pick the file to browse
GluResult  = loadGlu(fullfile(fpath{f},'Glu_Result.mat'));
VoltResult = importdata(fullfile(fpath{f},'Volt_Result.mat'));
STA_Result = importdata(fullfile(fpath{f},'STA_Glu_Volt_Result.mat'));   % needs S5 outputs
load(fullfile(fpath{f},'output_data.mat'));                              % Device_Data
nTg = size(GluResult.dFF_glu,2);

% voltage footprints aligned into the glutamate image space (contours)
alignedVoltftprnt = [];
for ft = 1:size(VoltResult.ftprnt,3)
    alignedVoltftprnt(:,:,ft) = transformCamera_O2B(Device_Data, tformReg, VoltResult.ftprnt(:,:,ft), GluResult.AvgGluImg);
end
goodROI = VoltResult.GoodNeuron;                 % voltage neurons with place maps (S5 order)

% voltage-trace place map (mean V per lap x position) for the middle voltage panel
VR_Result = importdata(fullfile(fpath{f},'VR_Result.mat'));
Tvr = min(size(VR_Result,2), size(VoltResult.traces,2));
V_Lap_trace = PlaceTrigger_average(VoltResult.traces(goodROI,1:Tvr), place_bin, VR_Result(:,1:Tvr), vel_thresh, VRtrackLength);

G = struct('coord',getGluCoord(GluResult), 'cvec',-log10(STA_Result.SpatialInfo_p), ...
           'dFF',GluResult.dFF_glu, 'spike',find_spike_bh(zscore(GluResult.dFF_glu,0,2),3,1), ...
           't',GluResult.t_ax(1:nTg), 'Lap_dFF',STA_Result.Lap_dFF, 'Lap_FR',STA_Result.Lap_FR);
V = struct('ftprnt',alignedVoltftprnt(:,:,goodROI), 'cvec',-log10(STA_Result.V_SpatialInfo_p), ...
           'trace',VoltResult.traces(goodROI,:), 'spike',VoltResult.spike(goodROI,:), ...
           't',VoltResult.t_ax(1:size(VoltResult.traces,2)), ...
           'Lap_FR',STA_Result.V_Lap_FR, 'Lap_trace',V_Lap_trace);
interactive_voltglu_ROI(GluResult.AvgGluImg, G, V, {'Position (cm)', VRtrackLength});

%% S6c : Tuning-curve correlation browser (significant synapses only; set f)
%  Top    : significant synapses (SI_p<0.05) colored by PEAK position as an angle (circular, hsv).
%  Bottom : click a synapse -> synapses recolored by tuning-curve correlation to it (caxis [-1 1]).
f = foi_STA(3);            % <-- pick the file to browse
GluResult  = loadGlu(fullfile(fpath{f},'Glu_Result.mat'));
STA_Result = importdata(fullfile(fpath{f},'STA_Glu_Volt_Result.mat'));   % needs S5 outputs
coord   = getGluCoord(GluResult);
popMap  = squeeze(mean(STA_Result.Lap_FR,1,'omitnan'))';    % Nsyn x place_bin (mean place firing rate)
sig     = STA_Result.SpatialInfo_p < 0.05;                  % significant synapses (NaN -> excluded)
fprintf('S6c: %d/%d synapses significant (p<0.05)\n', sum(sig), numel(sig));
interactive_glu_tuningcorr(GluResult.AvgGluImg, coord(sig,:), popMap(sig,:));

%% S6d : Max-lag dF/F trace-correlation browser (significant synapses only; set f)
%  Top    : significant synapses colored by PEAK position (angle, hsv).
%  Bottom : click a synapse -> synapses recolored by MAX cross-correlation of their raw dF/F
%           trace with the clicked one over lags +/- 5 frames.
f = foi_STA(3);            % <-- pick the file to browse
GluResult  = loadGlu(fullfile(fpath{f},'Glu_Result.mat'));
STA_Result = importdata(fullfile(fpath{f},'STA_Glu_Volt_Result.mat'));
coord   = getGluCoord(GluResult);
popMap  = squeeze(mean(STA_Result.Lap_FR,1,'omitnan'))';    % Nsyn x place_bin (peak-position color)
sig     = STA_Result.SpatialInfo_p < 0.05;
fprintf('S6d: %d/%d synapses significant (p<0.05)\n', sum(sig), numel(sig));
interactive_glu_tracecorr(GluResult.AvgGluImg, coord(sig,:), popMap(sig,:), GluResult.dFF_glu(sig,:), 5);

%% S7 : Voltage + glutamate example (nominate file) : spike-event table
%  (A) Table of voltage spike EVENTS (SS singles + CS bursts) per good neuron:
%      ROI, Type, nSpike, voltage frame, glutamate frame, times, V->G gap, CS burst duration.
%  Then: S7b = peri-spike traces for a nominated event; S7c = SS/CS pre-spike glu tuning maps
%  (S7c REUSES VtrigGlu computed in S4b -> no recomputation here).
f = foi_STA(2);            % <-- nominate the file to view
fprintf('=== S7 volt+glu example, file #%d ===\n', f);
VoltResult = importdata(fullfile(fpath{f},'Volt_Result.mat'));
GluResult  = loadGlu(fullfile(fpath{f},'Glu_Result.mat'));
nTg  = size(GluResult.dFF_glu,2);
dtG  = median(diff(GluResult.t_ax(1:nTg)));
dtV  = median(diff(VoltResult.t_ax));
gtax = GluResult.t_ax(1:nTg);
nTv  = size(VoltResult.traces,2);
coord   = getGluCoord(GluResult);
AvgImg  = GluResult.AvgGluImg;
goodROI = find(VoltResult.GoodNeuron(:))';

% ---------- (A) spike-event table ----------
rows = {};
for n = goodROI
    ssV  = find(VoltResult.SpClass(1,:,n) > 0);              % simple spikes
    csV  = find(VoltResult.SpClass(2,:,n) > 0);              % complex-spike onsets
    bwCS = bwlabel(VoltResult.CStrace(n,:));
    for k = 1:numel(ssV)
        fv = ssV(k); gi = match_nearest(VoltResult.t_ax(fv), gtax);
        rows(end+1,:) = {n,"SS",1,fv,VoltResult.t_ax(fv),gi,gtax(gi),(VoltResult.t_ax(fv)-gtax(gi))*1000,0}; %#ok<SAGROW>
    end
    for k = 1:numel(csV)
        fo = csV(k); w = max(1,fo-2):min(nTv,fo+2); lbl = max(bwCS(w));
        if lbl>0, reg = find(bwCS==lbl); nspk = sum(VoltResult.spike(n,reg)); dur = (reg(end)-reg(1))*dtV*1000;
        else,     nspk = NaN;            dur = NaN; end
        gi = match_nearest(VoltResult.t_ax(fo), gtax);
        rows(end+1,:) = {n,"CS",nspk,fo,VoltResult.t_ax(fo),gi,gtax(gi),(VoltResult.t_ax(fo)-gtax(gi))*1000,dur}; %#ok<SAGROW>
    end
end
SpikeTable = cell2table(rows,'VariableNames', ...
    {'ROI','Type','nSpike','Vframe','Vtime_s','Gframe','Gtime_s','V_G_gap_ms','CSdur_ms'});
SpikeTable = sortrows(SpikeTable,{'ROI','Vframe'});
fprintf('  %d spike events (%d SS, %d CS) across %d good neurons\n', ...
    height(SpikeTable), sum(SpikeTable.Type=="SS"), sum(SpikeTable.Type=="CS"), numel(goodROI));
disp(SpikeTable(1:min(25,height(SpikeTable)),:));   % full table is in the SpikeTable variable

%% S7b : Peri-spike glutamate (tuned-sorted kymograph) + pre/post tuning maps (run S7 first; set evIdx)
%  Pick a row of SpikeTable; show (1) the ROI's voltage trace around the event,
%  (2) the peri-event glutamate kymograph with synapses sorted by THIS neuron's spike-tuning,
%  (3-5) avgImg + synapse scatter colored by THIS EVENT's glutamate amplitude (mean dF/F) in
%        the -50..0 ms (pre), 0..50 ms (post), and their difference windows. These now change
%        with evIdx (the STA/neuron-average version depended only on ROI+class, so it looked
%        identical across events of the same neuron+type).
evIdx = 100;                                   % <-- nominate an event = row index of SpikeTable
nbin=10;
ev = SpikeTable(evIdx,:);
n  = ev.ROI;  vf = ev.Vframe;  gf = ev.Gframe;  cls = char(ev.Type);   % 'SS' or 'CS'
winMs = 500;                                 % +/- window to display (ms)
vw = round(winMs/1000/dtV);  vI = max(1,vf-vw):min(nTv,vf+vw);  tVms = (vI-vf)*dtV*1000;
gw = round(winMs/1000/dtG);  gI = max(1,gf-gw):min(nTg,gf+gw);  tGms = (gI-gf)*dtG*1000;

% ---- this neuron's spike-triggered glutamate tuning (from S4b), for sorting + maps ----
STA_Result = importdata(fullfile(fpath{f},'STA_Glu_Volt_Result.mat'));
VG    = STA_Result.VtrigGlu;   gidx = find([VG.roi]==n,1);
C     = VG(gidx).(cls);                       % .score/.pval Nsyn x 2 (col1 = -50..0, col2 = 0..50 ms)
scPre = C.score(:,1);  scPost = C.score(:,2);
tunedSSCS = any(VG(gidx).SS.pval<0.05,2) | any(VG(gidx).CS.pval<0.05,2);   % tuned to this neuron's SS OR CS
zscoredGluTrace=zscore(GluResult.dFF_glu,0,2);

figure(680); clf; tiledlayout(6,2,'TileSpacing','compact','Padding','compact');
% (1) voltage trace + spikes in window
nexttile(1,[2 1]);
plot(tVms, VoltResult.traces(n,vI),'k'); hold on;
sp = vI(VoltResult.spike(n,vI)>0); plot((sp-vf)*dtV*1000, VoltResult.traces(n,sp),'r.','MarkerSize',10);
xline(0,'r'); axis tight; box off; ylabel('V');
title(sprintf('Event %d: ROI %d %s (nSpike=%d)   Vframe %d / Gframe %d', evIdx, n, string(ev.Type), ev.nSpike, vf, gf));

% (2) peri-event glutamate averaged in 5 camera-x bins (every 20%), tuned (top) vs untuned (bottom)
xE = linspace(min(coord(:,1)), max(coord(:,1)), nbin+1);          % 5 equal-width x bins
xC = round((xE(1:end-1)+xE(2:end))/2);  binOf = discretize(coord(:,1), xE);  binCol = jet(nbin);
grpMask = {tunedSSCS, 'tuned (SS/CS)'; ~tunedSSCS, 'untuned'};
tilePos = [5 2; 9 2];                                        % [firstTile rowspan] : tuned rows2-4, untuned rows5-6
for g = 1:2
    nexttile(tilePos(g,1), [tilePos(g,2) 1]); hold on;
    lg = gobjects(1,5); lgtxt = cell(1,5);
    for b = 1:nbin
        inB = grpMask{g,1} & (binOf==b);
        if ~any(inB), continue; end
        M = movmean(zscoredGluTrace(inB, gI),5,2);
        m = mean(M,1,'omitnan');  s = std(M,0,1,'omitnan') ./ sqrt(sum(~isnan(M),1));
        gd = ~isnan(m);
        fill([tGms(gd) fliplr(tGms(gd))], [m(gd)+s(gd) fliplr(m(gd)-s(gd))], binCol(b,:), 'FaceAlpha',0.2,'EdgeColor','none');
        lg(b) = plot(tGms(gd), m(gd), '-', 'Color', binCol(b,:), 'LineWidth',1.5);
        lgtxt{b} = sprintf('x~%d (n=%d)', xC(b), sum(inB));
    end
    xline(0,'r'); axis tight; box off; xlabel('Time from event (ms)'); ylabel('Z score');
    title(sprintf('%s synapses: glu by x-bin (mean\\pmSEM)', grpMask{g,2}));
    ok = isgraphics(lg); if any(ok), legend(lg(ok), lgtxt(ok), 'Location','eastoutside','Box','off','fontsize',10); end
end

% (3-4) avgImg + synapse scatter colored by THIS EVENT's glu amplitude (pre / post window)
imgRGB = repmat(mat2gray(double(AvgImg),[30 150]),1,1,3);
nc=256; hc=nc/2; cmapBWR=[[linspace(0,1,hc)';ones(hc,1)],[linspace(0,1,hc)';linspace(1,0,hc)'],[ones(hc,1);linspace(1,0,hc)']];
gPre  = gI(tGms>=-50 & tGms<0);   evPre  = mean(GluResult.dFF_glu(:,gPre), 2,'omitnan');   % this event, -50..0 ms
gPost = gI(tGms>0   & tGms<=50);  evPost = mean(GluResult.dFF_glu(:,gPost),2,'omitnan');   % this event,  0..50 ms
cl = prctile(abs([evPre(tunedSSCS);evPost(tunedSSCS)]),95); if ~isfinite(cl)||cl<=0, cl=1; end
mapDat = {evPre,'from -50 to 0 ms (pre)'; evPost,'from 0 to 50 ms (post)'; evPost-evPre,'Post - pre'};
% this neuron's voltage footprint aligned into the glutamate image space (soma contour)
load(fullfile(fpath{f},'output_data.mat'));                                % Device_Data
fpAligned = transformCamera_O2B(Device_Data, tformReg, VoltResult.ftprnt(:,:,n), AvgImg);
for pp = 1:3
    nexttile(2+(pp-1)*4,[2 1]); image(imgRGB); hold on; axis image off;
    for xb = 1:numel(xE), xline(xE(xb),'LineStyle',':','Color',[.85 .85 .85],'LineWidth',0.5); end   % x-bin edges
    mfp = max(fpAligned(:));
    if isfinite(mfp)&&mfp>0, contour(fpAligned/mfp, [0.2 0.5], 'c', 'LineWidth', 1); end   % soma footprint contour
     scatter(coord(tunedSSCS,1),coord(tunedSSCS,2),5,mapDat{pp,1}(tunedSSCS),'filled');   % SS/CS-tuned synapses only
    colormap(gca,cmapBWR); caxis([-cl cl]); cb=colorbar; cb.Label.String='dF/F';
     title(sprintf('%s  ev%d glu (%d SS/CS-tuned syn)',mapDat{pp,2},evIdx,sum(tunedSSCS)));
end
drawnow;
set_font('Arial'); set_fontsize(16); set_figsize(500,270);

%% S7b_2 : SS (isolated) & CS TRIGGERED-AVERAGE glutamate, S7b-style (run S7 first; set n)
%  Same layout idea as S7b but for the spike-triggered AVERAGE (VtrigGlu from S4b), per neuron:
%    SS = isolated leading simple spike (nothing in the preceding VtrigIso_ms=100 ms), CS = burst onset.
%  3x3 layout: col1 = STA voltage (SS/CS) + SS-avg + CS-avg glu by camera-x bin;
%              col2 = PRE (-50..0 ms) amplitude maps (SS | CS | CS-SS);
%              col3 = POST (0..50 ms) amplitude maps (SS | CS | CS-SS).
%             (baseline = STA over tMs < -50 ms; 3-frame smoothed)
%  Traces/maps are in Z-SCORE units (each synapse's STA divided by its dF/F std) so artifact
%  synapses with a huge dF/F are normalised to their own variability.
n = SpikeTable.ROI(evIdx);      % <-- neuron to view (defaults to the S7b event's ROI)
nbin2 = 10;                     % # camera-x bins for the trace panels
STA_Result = importdata(fullfile(fpath{f},'STA_Glu_Volt_Result.mat'));
VG   = STA_Result.VtrigGlu;  gidx = find([VG.roi]==n,1);
assert(~isempty(gidx), 'ROI %d not found in VtrigGlu - run S4b for file %d.', n, f);
Sss = VG(gidx).SS;  Scs = VG(gidx).CS;  tMs = Sss.tMs;
tunedSSCS = any(Sss.pval<0.05,2) | any(Scs.pval<0.05,2);              % tuned to SS or CS
baseEdge  = min(STA_Result.VtrigGlu_score_ms(:,1));                   % far-pre baseline edge (= -50 ms)
% z-score by each synapse's dF/F std (STA is a linear mean, so STA/std == STA of the z-scored
% trace). This normalises artifact synapses with a huge dF/F to their own variability.
gStd = std(GluResult.dFF_glu,0,2);  gStd(gStd==0) = NaN;
SSc = (Sss.STA - mean(Sss.STA(:, tMs<baseEdge),2,'omitnan')) ./ gStd;   % z-scored, baseline-subtracted SS STA
CSc = (Scs.STA - mean(Scs.STA(:, tMs<baseEdge),2,'omitnan')) ./ gStd;   % z-scored, baseline-subtracted CS STA
ssAmp    = Sss.score(:,2) ./ gStd;  csAmp    = Scs.score(:,2) ./ gStd;  % 0..50 ms POST-window amplitude (z-score)
ssAmpPre = Sss.score(:,1) ./ gStd;  csAmpPre = Scs.score(:,1) ./ gStd;  % -50..0 ms PRE-window amplitude (z-score)

xE = linspace(min(coord(:,1)), max(coord(:,1)), nbin2+1);            % camera-x bins
xC = round((xE(1:end-1)+xE(2:end))/2);  binOf = discretize(coord(:,1), xE);  binCol = jet(nbin2);
imgRGB = repmat(mat2gray(double(AvgImg),[30 150]),1,1,3);
nc=256; hc=nc/2; cmapBWR=[[linspace(0,1,hc)';ones(hc,1)],[linspace(0,1,hc)';linspace(1,0,hc)'],[ones(hc,1);linspace(1,0,hc)']];
load(fullfile(fpath{f},'output_data.mat'));                          % Device_Data
fpAligned = transformCamera_O2B(Device_Data, tformReg, VoltResult.ftprnt(:,:,n), AvgImg);

figure(681); clf; tiledlayout(3,3,'TileSpacing','compact','Padding','compact');
% col1 = traces (voltage, SS, CS) ; col2 = PRE (-50..0) maps ; col3 = POST (0..50) maps
% --- col1, top: STA voltage (SS/CS) ---
nexttile(1); hold on;
if all(isfield(STA_Result,{'STA_V_SS','STA_V_CS','tV'}))
    tVms = STA_Result.tV*1000;
    plot(tVms, STA_Result.STA_V_SS(n,:), 'Color',[.2 .5 .9],'LineWidth',1.5);
    plot(tVms, STA_Result.STA_V_CS(n,:), 'Color',[.9 .3 .2],'LineWidth',1.5);
    legend({'SS','CS'},'Box','off','Location','best'); ylabel('V (\sigma)');
end
xline(0,'r'); axis tight; box off; xlabel('Time from spike (ms)');
title(sprintf('ROI %d STA voltage (n_{SS}=%d, n_{CS}=%d)', n, Sss.ntrig, Scs.ntrig));

% --- col1: SS then CS triggered-average glu by camera-x bin (tuned synapses) ---
staDat = {SSc,'SS',4; CSc,'CS',7};
for q = 1:2
    nexttile(staDat{q,3}); hold on;  lg = gobjects(1,nbin2); lgtxt = cell(1,nbin2);
    for b = 1:nbin2
        inB = tunedSSCS & (binOf==b);
        if ~any(inB), continue; end
        M = movmean(staDat{q,1}(inB,:), 3, 2);
        m = mean(M,1,'omitnan');  s = std(M,0,1,'omitnan') ./ sqrt(sum(~isnan(M),1));
        gd = ~isnan(m);
        fill([tMs(gd) fliplr(tMs(gd))], [m(gd)+s(gd) fliplr(m(gd)-s(gd))], binCol(b,:), 'FaceAlpha',0.2,'EdgeColor','none');
        lg(b) = plot(tMs(gd), m(gd), '-', 'Color', binCol(b,:), 'LineWidth',1.5);
        lgtxt{b} = sprintf('x~%d (n=%d)', xC(b), sum(inB));
    end
    xline(0,'r'); axis tight; box off; xlabel('Time from spike (ms)'); ylabel('z-score');
    title(sprintf('%s-triggered avg glu by x-bin', staDat{q,2}));
    ok = isgraphics(lg); if any(ok), legend(lg(ok), lgtxt(ok), 'Location','eastoutside','Box','off','FontSize',9); end
end

% --- col2 = PRE (-50..0 ms), col3 = POST (0..50 ms) amplitude maps (SS / CS / CS-SS) ---
clA = prctile(abs([ssAmpPre(tunedSSCS);csAmpPre(tunedSSCS);ssAmp(tunedSSCS);csAmp(tunedSSCS)]),90);
clD = prctile(abs([csAmpPre(tunedSSCS)-ssAmpPre(tunedSSCS); csAmp(tunedSSCS)-ssAmp(tunedSSCS)]),90);
if ~isfinite(clA)||clA<=0, clA=1; end
if ~isfinite(clD)||clD<=0, clD=1; end
mapDat = {ssAmpPre,'SS pre',2,clA; csAmpPre,'CS pre',5,clA; csAmpPre-ssAmpPre,'CS-SS pre',8,clD; ...
          ssAmp,   'SS post',3,clA; csAmp,   'CS post',6,clA; csAmp-ssAmp,      'CS-SS post',9,clD};
for q = 1:size(mapDat,1)
    nexttile(mapDat{q,3}); image(imgRGB); hold on; axis image off;
    for xb = 1:numel(xE), xline(xE(xb),'LineStyle',':','Color',[.85 .85 .85],'LineWidth',0.5); end
    mfp = max(fpAligned(:)); if isfinite(mfp)&&mfp>0, contour(fpAligned/mfp, [0.2 0.5], 'c', 'LineWidth', 1); end
    scatter(coord(tunedSSCS,1),coord(tunedSSCS,2),4,mapDat{q,1}(tunedSSCS),'filled');
    colormap(gca,cmapBWR); caxis([-mapDat{q,4} mapDat{q,4}]); cb=colorbar; cb.Label.String='z-score';
    title(sprintf('%s glu (%d tuned syn)', mapDat{q,2}, sum(tunedSSCS)));
end
sgtitle(sprintf('File %d ROI %d : SS (isolated) & CS triggered-average glutamate', f, n));
drawnow;

%% S7c : SS/CS pre-spike glutamate tuning maps (REUSE S4b results; run S4b first; set f)
%  Loads VtrigGlu from STA_Glu_Volt_Result.mat (computed in S4b) and plots per neuron:
%  SS/CS STA-glutamate heatmaps + maps of synapses tuned in the -50..0 ms pre-spike window.
f = foi_STA(2);            % <-- file to view
GluResult  = loadGlu(fullfile(fpath{f},'Glu_Result.mat'));
VoltResult = importdata(fullfile(fpath{f},'Volt_Result.mat'));           % for footprint contour
STA_Result = importdata(fullfile(fpath{f},'STA_Glu_Volt_Result.mat'));   % contains VtrigGlu (from S4b)
if ~isfield(STA_Result,'VtrigGlu'), error('VtrigGlu not found - run S4b for file %d first.', f); end

% place-tuned synapse mask (from S5). Restrict the display to these synapses.
if isfield(STA_Result,'SpatialInfo_p')
    placeSel = STA_Result.SpatialInfo_p(:) < 0.05;
    fprintf('S7c/S7d: restricting to %d/%d place-tuned synapses (SI_p<0.05)\n', ...
        sum(placeSel), numel(placeSel));
else
    warning('SpatialInfo_p not found - run S5 for file %d to get place tuning. Showing ALL synapses.', f);
    placeSel = true(size(GluResult.dFF_glu,1),1);
end

% voltage footprints aligned into the glutamate image space (soma contour)
load(fullfile(fpath{f},'output_data.mat'));                             % Device_Data
alignedVoltftprnt = zeros([size(GluResult.AvgGluImg,1) size(GluResult.AvgGluImg,2) size(VoltResult.ftprnt,3)]);
for ft = 1:size(VoltResult.ftprnt,3)
    alignedVoltftprnt(:,:,ft) = transformCamera_O2B(Device_Data, tformReg, VoltResult.ftprnt(:,:,ft), GluResult.AvgGluImg);
end
% 
% plot_VtrigGlu(STA_Result.VtrigGlu, getGluCoord(GluResult), GluResult.AvgGluImg, ...
%               STA_Result.VtrigGlu_score_ms, 700, alignedVoltftprnt, placeSel);

% S7d : SS / CS / (CS-SS) triggered glutamate MOVIE over time (per neuron; run S7c first)
%  Per good neuron, write a 3-panel movie (SS-triggered glu | CS-triggered glu | CS-SS
%  difference) animating AvgGluImg + the spike-tuned synapses colored by baseline-subtracted
%  dF/F across the STA time axis (each STA baseline-subtracted over its far-pre window,
%  tMs < -50 ms). SS and CS share a color scale; the difference panel has its own. One .mp4
%  per neuron is saved into fpath{f} (Motion-JPEG .avi fallback if the MPEG-4 codec is
%  unavailable).  Synapses = tuned to SS OR CS in ANY window, per neuron.
%  Row 2 (if S2 was run) shows the STA voltage (SS | CS | overlay) with a sweeping dotted
%  line at the current movie time. Panel/super titles report the # of SS and CS averaged.
%  NB: SS = simple spikes OUTSIDE complex-spike bursts (all of them); CS = complex-spike
%      burst ONSET (first spike of the burst).
VtrigMov_ms = [-100 200];                                    % movie display window (ms from spike)
VtrigMovFPS = 10;                                          % playback frame rate
VtrigGlu = STA_Result.VtrigGlu;
baseEdge = min(STA_Result.VtrigGlu_score_ms(:,1));         % far-pre baseline edge (= -50 ms)
coordAll = getGluCoord(GluResult);                         % all synapse coordinates
imgRGB   = repmat(mat2gray(double(GluResult.AvgGluImg),[30 150]),1,1,3);   % avg image (windowed [30 150]) as truecolor bg
nc = 256; hc = nc/2;                                       % blue-white-red diverging colormap
cmapBWR = [ [linspace(0,1,hc)'; ones(hc,1)], ...
            [linspace(0,1,hc)'; linspace(1,0,hc)'], ...
            [ones(hc,1); linspace(1,0,hc)'] ];
% STA voltage (from S2) for the 2nd movie row, if available
hasV = all(isfield(STA_Result,{'STA_V_SS','STA_V_CS','tV'}));
if hasV, tVms = STA_Result.tV*1000; else, warning('STA voltage (STA_V_SS/CS) not found - run S2; movie will have no voltage row.'); end
for gi = 1:numel(VtrigGlu)
    roi = VtrigGlu(gi).roi;
    Sss = VtrigGlu(gi).SS;  Scs = VtrigGlu(gi).CS;  tMs = Sss.tMs;
    nSS = Sss.ntrig;  nCS = Scs.ntrig;                     % # spikes averaged (glu STA triggers)
    tunedSel = any(Sss.pval < 0.05, 2) | any(Scs.pval < 0.05, 2);   % tuned to SS or CS, any window
    if ~any(tunedSel), fprintf('  ROI %d: no spike-tuned synapse; skipping movie.\n', roi); continue; end
    coordP = coordAll(tunedSel,:);
    SSsta = movmean(Sss.STA(tunedSel,:), 3, 2);  CSsta = movmean(Scs.STA(tunedSel,:), 3, 2);   % spike-tuned + 3-frame temporal smoothing
    bSS = mean(SSsta(:, tMs < baseEdge), 2, 'omitnan');   % Nsyn x 1 baseline (SS)
    bCS = mean(CSsta(:, tMs < baseEdge), 2, 'omitnan');   % Nsyn x 1 baseline (CS)
    SSc = SSsta - bSS;                                     % baseline-subtracted SS-triggered glu
    CSc = CSsta - bCS;                                     % baseline-subtracted CS-triggered glu
    D   = CSc - SSc;                                       % CS - SS difference

    frIdx = find(tMs >= VtrigMov_ms(1) & tMs <= VtrigMov_ms(2));   % movie frames (native STA resolution)
    clAmp = prctile(abs([SSc(:,frIdx) CSc(:,frIdx)]),95,'all');    % shared SS/CS scale
    clDif = prctile(abs(D(:,frIdx)),80,'all');                     % separate difference scale
    if ~isfinite(clAmp)||clAmp<=0, clAmp = 1; end
    if ~isfinite(clDif)||clDif<=0, clDif = 1; end
    panels = { sprintf('SS glu (n=%d)',nSS), SSc, clAmp; ...
               sprintf('CS glu (n=%d)',nCS), CSc, clAmp; ...
               'CS - SS glu',                D,   clDif };

    % --- set up a fixed-size, undocked figure so every captured frame is identical ---
    %  LEFT column = 3 glu scatter maps ([2 1] each); RIGHT column = 2 STA voltage traces ([3 1] each).
    nCol = 1 + hasV;                                       % col 1 = scatter maps, col 2 = voltage traces
    fig = figure(760+roi); clf;
    set(fig,'WindowStyle','normal','Color','k','Position',[80 40 350*4 250*3],'theme','dark');
    tl  = tiledlayout(fig, 6, nCol,'Padding','compact');
    hSc = gobjects(1,3);
    for p = 1:3                                            % --- LEFT column: glutamate spatial maps ---
        ax = nexttile(tl, 2*(p-1)*nCol + 1, [2 1]); image(ax, imgRGB); hold(ax,'on'); axis(ax,'image','off');  % avg image bg
        if roi<=size(alignedVoltftprnt,3)                 % soma footprint contour
            fp = alignedVoltftprnt(:,:,roi); mfp = max(fp(:));
            if isfinite(mfp)&&mfp>0, contour(ax, fp/mfp, [0.2 0.5], 'c', 'LineWidth', 1); end
        end
        hSc(p) = scatter(ax, coordP(:,1), coordP(:,2), 3, panels{p,2}(:,frIdx(1)), 'filled');  % synapses over image
        colormap(ax, cmapBWR); caxis(ax, [-panels{p,3} panels{p,3}]);
        cb = colorbar(ax); cb.Label.String = '\DeltaF/F'; title(ax, panels{p,1});
        cb.Ticks = [-panels{p,3} panels{p,3}];  cb.TickLabels = num2str([-panels{p,3} panels{p,3}]',2);
    end
    hVL = gobjects(1,2);
    if hasV                                               % --- RIGHT column: STA voltage traces ---
        vSS = STA_Result.STA_V_SS(roi,:);  vCS = STA_Result.STA_V_CS(roi,:);
        % trace 1 (rows 1-3, col 2): SS & CS STA voltage overlaid, different colors
        ax = nexttile(tl, 2, [3 1]); hold(ax,'on');
        plot(ax,tVms,vSS,'Color',[.3 .6 1] ,'LineWidth',1.5);     % SS = blue
        plot(ax,tVms,vCS,'Color',[1 .35 .35],'LineWidth',1.5);    % CS = red
        legend(ax,{'SS','CS'},'TextColor','w','Color','none','Box','off','Location','best');
        hVL(1) = xline(ax, tMs(frIdx(1)), 'y:', 'LineWidth', 1.5);
        xlim(ax, VtrigMov_ms); xlabel(ax,'Time from spike (ms)'); ylabel(ax,'V (\sigma)');
        title(ax,'SS (blue) & CS (red) STA voltage'); box(ax,'off');
        % trace 2 (rows 4-6, col 2): CS - SS voltage
        ax = nexttile(tl, 3*nCol+2, [3 1]); hold(ax,'on');
        plot(ax,tVms,vCS-vSS,'Color',[.4 1 .4],'LineWidth',1.5);
        hVL(2) = xline(ax, tMs(frIdx(1)), 'y:', 'LineWidth', 1.5);
        xlim(ax, VtrigMov_ms); xlabel(ax,'Time from spike (ms)'); ylabel(ax,'V (\sigma)');
        title(ax,'CS - SS voltage'); box(ax,'off');
    end
    st = sgtitle(tl, '', 'Color','w');   % layout title (reserves space -> no overlap with tile titles)

    % --- open the movie writer (MPEG-4, fall back to Motion JPEG AVI) ---
    movBase = fullfile(fpath{f}, sprintf('S7d_SS_CS_diff_glu_ROI%02d', roi));
    try
        vw = VideoWriter([movBase '.mp4'], 'MPEG-4');   movFile = [movBase '.mp4'];
    catch
        vw = VideoWriter([movBase '.avi'], 'Motion JPEG AVI'); movFile = [movBase '.avi'];
    end
    vw.FrameRate = VtrigMovFPS; open(vw);

    for it = 1:numel(frIdx)
        tNow = tMs(frIdx(it));
        for p = 1:3, set(hSc(p), 'CData', panels{p,2}(:,frIdx(it))); end
        if hasV, for p = 1:2, set(hVL(p),'Value',tNow); end, end
        set(st, 'String', sprintf('File %d ROI %d : SS(n=%d) | CS(n=%d) | CS-SS   %+d ms  (%d tuned syn)', ...
            f, roi, nSS, nCS, round(tNow), sum(tunedSel)), 'Color','w');
        set_font('Arial'); set_fontsize(16);
        drawnow;
        im = frame2im(getframe(fig));
        if it==1, movSz = 2*floor([size(im,1) size(im,2)]/2); end        % lock even frame size on 1st frame
        im = im(1:min(size(im,1),movSz(1)), 1:min(size(im,2),movSz(2)), :);          % crop to locked size
        if size(im,1)<movSz(1) || size(im,2)<movSz(2), im(movSz(1),movSz(2),size(im,3)) = 0; end  % pad if getframe drifted small
        writeVideo(vw, im);
    end
    close(vw);
    fprintf('  ROI %d: saved %d-frame movie -> %s\n', roi, numel(frIdx), movFile);
end

%% =====================================================================
%  SECTION 8 : SPATIAL STATISTICS  (pairwise correlation vs distance; set f)
% =====================================================================
% S8 : pairwise correlation vs Euclidean distance + tuned-count profile (set f; needs S5)
%  Between every pair of synapses:
%    - dF/F trace cross-correlation: max over +/- S8_maxlag frames of the z-scored dF/F
%    - tuning-curve correlation   : Pearson corr of the place tuning curves (Lap_dFF), which
%      are first smoothed with a 5-bin CIRCULAR moving average (VR position wraps).
%  are plotted vs their Euclidean distance on the camera (um, via pixelsize) as mean +/- SEM
%  (errorbar_shade). Figure = 2 subplots: LEFT dF/F xcorr, RIGHT tuning corr; within each,
%  one curve per camera-x bin of the reference synapse (S8_nXbin equal-width bins).
%  A second figure shows the number of place-tuned synapses (SpatialInfo_p<0.05) vs camera x.
f = foi_STA(3);            % <-- file to analyse (must have S5 outputs)
GluResult  = loadGlu(fullfile(fpath{f},'Glu_Result.mat'));
STA_Result = importdata(fullfile(fpath{f},'STA_Glu_Volt_Result.mat'));
assert(isfield(STA_Result,'Lap_dFF') && isfield(STA_Result,'SpatialInfo_p'), ...
    'S8 needs S5 outputs (Lap_dFF, SpatialInfo_p) for file %d - run S5 first.', f);

coordAll = getGluCoord(GluResult);  NsynAll = size(coordAll,1);   % ALL synapses (for the count profile)
W = size(GluResult.AvgGluImg,2);                             % camera frame width (x extent, px)

% -- parameters --
S8_maxlag = 5;      % dF/F cross-correlation max lag (glu frames)
S8_nXbin  = 5;      % # camera-x bins for reference synapses (= # subplots)
S8_nDbin  = 15;     % # Euclidean-distance bins
S8_nXcnt  = 20;     % # camera-x bins for the tuned-count profile
pixelsize= 6.5/16*180/100; %µm

% -- tuning curves (place), significance, dF/F traces --
%  VR position is circular (teleports to the start), so the tuning curve is smoothed with a
%  5-bin CIRCULAR moving average (ringmovMean wraps end->start); the correlation is then a
%  plain Pearson between the smoothed curves (position-specific).
popMap   = squeeze(mean(ringmovMean(STA_Result.Lap_dFF,5),1,'omitnan'))';  % NsynAll x place_bin, 5-bin circular movmean
placeSel = STA_Result.SpatialInfo_p(:) < 0.05;              % place-tuned synapses
dFF      = GluResult.dFF_glu;

% -- restrict the correlation analysis to PLACE-TUNED synapses (SI_p<0.05) --
coord  = coordAll(placeSel,:);  Nsyn = size(coord,1);
popMap = popMap(placeSel,:);    dFF  = dFF(placeSel,:);
fprintf('S8: correlation analysis on %d/%d place-tuned synapses (SI_p<0.05)\n', Nsyn, NsynAll);

% -- pairwise tuning-curve correlation (Pearson on the circularly-smoothed curves, NaN-safe) --
Rtune = corr(popMap','rows','pairwise');                    % Nsyn x Nsyn

% -- pairwise max-lag dF/F cross-correlation (over +/- S8_maxlag) --
%  This matrix form is numerically IDENTICAL to xcorr(Z',S8_maxlag,'coeff') (verified) but is
%  vectorised over all pairs instead of looping. NB: taking the MAX over the 2*maxlag+1 lags
%  biases the value UP - even independent traces give E[max] > 0 (~0.008 at T=4e4, and more for
%  autocorrelated traces), which is why the curve floors above 0 rather than decaying to 0.
%  S8_deBias subtracts that chance level, estimated with the SAME contiguous lag window centred
%  at random circular offsets, so uncorrelated (far) pairs sit at ~0.
S8_deBias = true;   % subtract the max-over-lag chance level
S8_nullSh = 8;      % # circular shifts used to estimate it
Z = zscore(dFF,0,2);  Tn = size(Z,2);  offd = ~logical(eye(Nsyn));
Rdff = -inf(Nsyn,Nsyn);
for L = 0:S8_maxlag
    C = (Z(:,1:Tn-L) * Z(:,1+L:Tn)') / (Tn-L);              % corr at lag +L (z-scored traces)
    Rdff = max(Rdff, max(C, C'));                           % C' covers lag -L -> +/- lag max
end
if S8_deBias
    rng(0); sh = randi([round(0.1*Tn) round(0.9*Tn)], 1, S8_nullSh);   ch = nan(1,S8_nullSh);
    for k = 1:S8_nullSh
        Zs = circshift(Z, -sh(k), 2);                       % decorrelate: Zs_j(t) = Z_j(t+s)
        Rn = -inf(Nsyn,Nsyn);
        for L = -S8_maxlag:S8_maxlag                        % same contiguous lag window, centred at s
            if L>=0, C = (Z(:,1:Tn-L)*Zs(:,1+L:Tn)')/(Tn-L);
            else,    C = (Z(:,1-L:Tn)*Zs(:,1:Tn+L)')/(Tn+L); end
            Rn = max(Rn, C);
        end
        ch(k) = mean(Rn(offd),'omitnan');
    end
    S8_chance = mean(ch);  Rdff = Rdff - S8_chance;
    fprintf('S8: max-lag chance level = %.4f (sd %.4f over %d shifts) subtracted from dF/F xcorr\n', ...
        S8_chance, std(ch), S8_nullSh);
end

% -- Euclidean distance on the camera (um) --
Dist = pdist2(coord, coord) * pixelsize;                    % um
selfPair = logical(eye(Nsyn));

% -- correlation vs distance: LEFT = dF/F xcorr, RIGHT = tuning corr; one curve per camera-x bin --
xEdges = linspace(1, W, S8_nXbin+1);
dEdges = linspace(0, prctile(Dist(~selfPair),99), S8_nDbin+1);
dCent  = (dEdges(1:end-1)+dEdges(2:end))/2;
binCol = nebula(S8_nXbin);                                   % color per reference-x bin
metrics = {Rdff, sprintf('dF/F xcorr (max over \\pm%d lag)',S8_maxlag); ...
           Rtune,'tuning-curve corr (Pearson, circ-smoothed)'};

% -- x-bin layout on the average image (colors match the curves below) --
figure('Color','w','Name','S8 x-bin layout');
image(repmat(mat2gray(double(GluResult.AvgGluImg),[30 200]),1,1,3)); hold on; axis image off;
scatter(coordAll(~placeSel,1), coordAll(~placeSel,2), 4, [.6 .6 .6], 'filled');   % excluded (not place-tuned)
for xb = 1:S8_nXbin
    s = coord(:,1)>=xEdges(xb) & coord(:,1)<xEdges(xb+1);   % same mask as the analysis below
    if any(s), scatter(coord(s,1), coord(s,2), 5, binCol(xb,:), 'filled'); end
end
for xb = 1:numel(xEdges), xline(xEdges(xb),'LineStyle',':','Color',[1 1 .3],'LineWidth',1); end   % x-bin edges
title(sprintf('File %d: camera-x bins - colored = %d place-tuned (used), grey = %d excluded', ...
    f, Nsyn, NsynAll-Nsyn));

figure('Color','w','Name','S8 corr vs distance'); tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
for mi = 1:2
    nexttile; hold on;
    lg = gobjects(1,S8_nXbin); lgtxt = cell(1,S8_nXbin);
    for xb = 1:S8_nXbin
        ri = coord(:,1)>=xEdges(xb) & coord(:,1)<xEdges(xb+1);      % reference synapses in this x-bin
        if ~any(ri), continue; end
        mask = false(Nsyn); mask(ri,:) = true; mask(selfPair) = false;   % their pairs (drop self)
        bi = discretize(Dist(mask), dEdges);
        [m,s] = binMeanSem(metrics{mi,1}(mask), bi, S8_nDbin);
        lg(xb)  = errorbar_shade(dCent, m, s, binCol(xb,:));
        lgtxt{xb} = sprintf('x %d-%d px (%d syn)', round(xEdges(xb)), round(xEdges(xb+1)), sum(ri));
    end
    yline(0,'k:'); box off; xlabel('Euclidean distance (\mum)'); ylabel('correlation');
    title(metrics{mi,2});
    ok = isgraphics(lg); if any(ok), legend(lg(ok), lgtxt(ok),'Location','northeast','Box','off'); end
end
sgtitle(sprintf('File %d: pairwise correlation vs distance (colored by camera-x bin of reference synapse)', f));

% -- # and FRACTION of place-tuned synapses vs camera x (uses ALL synapses) --
xEc = linspace(1, W, S8_nXcnt+1);  xCc = (xEc(1:end-1)+xEc(2:end))/2;
nTot = histcounts(coordAll(:,1),          xEc);      % all synapses per x bin
nTun = histcounts(coordAll(placeSel,1),   xEc);      % place-tuned per x bin
frac = nTun ./ nTot;  frac(nTot==0) = NaN;           % fraction place-tuned
figure('Color','w','Name','S8 tuned count vs x');
yyaxis left;
bTot = bar(xCc*pixelsize, nTot, 1, 'FaceColor',[.75 .75 .75],'EdgeColor','none'); hold on;   % all synapses
bTun = bar(xCc*pixelsize, nTun, 1, 'FaceColor',[.85 .3 .3], 'EdgeColor','none');              % place-tuned (overlaid)
ylabel('# synapses');
yyaxis right; pF = plot(xCc*pixelsize, frac, 'k-o','LineWidth',1); ylabel('fraction place-tuned'); ylim([0 1]);
legend([bTot bTun pF], {'all','place-tuned','fraction'}, 'Location','northeast','Box','off');
xlabel('Larminar displacement (µm)'); box off;
title(sprintf('File %d: place-tuned synapses vs camera x (%d/%d tuned)', f, sum(placeSel), NsynAll));

%% S8b : 2D (|dx|,|dy|) displacement map of pairwise correlation, per camera-x bin (run S8 first)
%  For every synapse pair: the ABSOLUTE x/y displacement (um) and their correlation. Reference
%  synapses are split into S8_nXbin camera-x bins; within each bin the pairs are 2D-binned on
%  (|dx|,|dy|) and the correlation averaged -> one heatmap per x-bin (S8_nXbin subplots).
%  One figure per metric (dF/F xcorr, tuning corr). Reuses S8: metrics, coord, xEdges, selfPair.
S8_nDxBin = 12;    % # |dx| bins for the 2D displacement map
S8_nDyBin = 12;    % # |dy| bins

DX = abs(coord(:,1) - coord(:,1)') * pixelsize;             % Nsyn x Nsyn, |dx| (um)
DY = abs(coord(:,2) - coord(:,2)') * pixelsize;             % Nsyn x Nsyn, |dy| (um)
dxEdges = linspace(0, prctile(DX(~selfPair),99), S8_nDxBin+1);
dyEdges = linspace(0, prctile(DY(~selfPair),99), S8_nDyBin+1);
dxC = (dxEdges(1:end-1)+dxEdges(2:end))/2;  dyC = (dyEdges(1:end-1)+dyEdges(2:end))/2;
nc=256; hc=nc/2; cmapBWR=[[linspace(0,1,hc)';ones(hc,1)],[linspace(0,1,hc)';linspace(1,0,hc)'],[ones(hc,1);linspace(1,0,hc)']];

for mi = 1:2
    R = metrics{mi,1};
    Mst = nan(S8_nDyBin, S8_nDxBin, S8_nXbin);              % mean corr per (dy,dx) per x-bin
    nRef = zeros(1,S8_nXbin);
    for xb = 1:S8_nXbin
        ri = coord(:,1)>=xEdges(xb) & coord(:,1)<xEdges(xb+1);       % reference synapses in this x-bin
        nRef(xb) = sum(ri); if ~any(ri), continue; end
        mask = false(Nsyn); mask(ri,:) = true; mask(selfPair) = false;   % their pairs (drop self)
        bx = discretize(DX(mask), dxEdges);  by = discretize(DY(mask), dyEdges);  rv = R(mask);
        ok = ~isnan(bx) & ~isnan(by) & ~isnan(rv);
        if ~any(ok), continue; end
        Mst(:,:,xb) = accumarray([by(ok) bx(ok)], rv(ok), [S8_nDyBin S8_nDxBin], @mean, NaN);
    end
    v = Mst(isfinite(Mst)); cl = prctile(abs(v),98);
    if isempty(cl)||~isfinite(cl)||cl<=0, cl = 1; end       % shared color scale across x-bins

    figure('Color','w','Name',sprintf('S8b displacement map: %s',metrics{mi,2}));
    tiledlayout(1,S8_nXbin,'TileSpacing','compact','Padding','compact');
    for xb = 1:S8_nXbin
        nexttile; imagesc(dxC, dyC, Mst(:,:,xb), 'AlphaData', isfinite(Mst(:,:,xb)));
        set(gca,'YDir','normal'); colormap(gca,cmapBWR); caxis([-cl cl]); axis square; box off;
        xlabel('|\Deltax| (\mum)'); ylabel('|\Deltay| (\mum)');
        title(sprintf('x %d-%d px (%d syn)', round(xEdges(xb)), round(xEdges(xb+1)), nRef(xb)));
        if xb==S8_nXbin, cb = colorbar; cb.Label.String = 'mean correlation'; end
    end
    sgtitle(sprintf('File %d: %s vs (|\\Deltax|,|\\Deltay|) by reference camera-x bin', f, metrics{mi,2}));
end

%% S8c : cluster the correlation matrices + map clusters spatially (run S8 first)
%  Hierarchical clustering (average linkage on 1-corr distance) of the pairwise correlation
%  matrices, for dF/F xcorr and tuning-curve corr, on BOTH the place-tuned set and ALL synapses.
%  ONE figure, 4 rows (set x metric) x 3 cols:
%    col 1 = cluster-ordered correlation matrix
%    col 2 = avgImg + synapses colored by cluster
%    col 3 = per-cluster laminar-position (camera-x, um) histogram
S8c_K = 5;                                                  % # clusters
% -- full correlation matrices over ALL synapses (submatrix gives any subset) --
popMapAll = squeeze(mean(ringmovMean(STA_Result.Lap_dFF,5),1,'omitnan'))';   % NsynAll x place_bin
Rtune_all = corr(popMapAll','rows','pairwise');
Zall = zscore(GluResult.dFF_glu,0,2);  Tn = size(Zall,2);
Rdff_all = -inf(NsynAll,NsynAll);
for L = 0:S8_maxlag, C = (Zall(:,1:Tn-L)*Zall(:,1+L:Tn)')/(Tn-L); Rdff_all = max(Rdff_all, max(C,C')); end

setsC = {true(NsynAll,1),'all'; placeSel,'place-tuned'};
metsC = {Rdff_all,'dF/F xcorr'; Rtune_all,'tuning corr'};
imgRGBc = repmat(mat2gray(double(GluResult.AvgGluImg),[30 300]),1,1,3);
xumAll  = coordAll(:,1)*pixelsize;                          % laminar position (um)
xEdgesH = linspace(min(xumAll), max(xumAll), 16);           % shared histogram bins
clustID = cell(2,2);                                        % keep cluster labels (rows=set, cols=metric)

figure('Color','w','Name','S8c clusters'); tiledlayout(4,3,'TileSpacing','compact','Padding','compact');
for si = 1:2
    sel = find(setsC{si,1});
    for mi = 1:2
        R = metsC{mi,1}(sel,sel);
        good = ~all(isnan(R),2);  idx = sel(good);  Rg = R(good,good);   % drop all-NaN synapses
        Rg(isnan(Rg)) = 0;                                  % residual NaN pairs -> 0 corr (distance 1)
        Dd = 1 - Rg;  Dd = max((Dd+Dd')/2, 0);  Dd(1:size(Dd,1)+1:end) = 0;   % symmetric nonneg zero-diag
        clab = cluster(linkage(squareform(Dd),'average'),'maxclust',S8c_K);
        cvec = nan(NsynAll,1); cvec(idx) = clab;  clustID{si,mi} = cvec;
        colc = lines(S8c_K);

        % --- col 1: cluster-ordered correlation matrix ---
        [~,so] = sort(clab);
        nexttile; imagesc(Rg(so,so), [-1 1]); axis image off; colormap(gca,cmapBWR);
        ylabel(sprintf('%s | %s', setsC{si,2}, metsC{mi,2}));
        title(sprintf('corr matrix (K=%d)',S8c_K)); cb = colorbar; cb.Label.String='corr';

        % --- col 2: spatial cluster map ---
        nexttile; image(imgRGBc); hold on; axis image off;
        out = true(NsynAll,1); out(idx) = false;
        scatter(coordAll(out,1), coordAll(out,2), 3, [.6 .6 .6], 'filled');   % not in this set
        for c = 1:S8c_K, m = idx(clab==c); scatter(coordAll(m,1), coordAll(m,2), 8, colc(c,:), 'filled'); end
        title(sprintf('%d syn colored by cluster', numel(idx)));

        % --- col 3: per-cluster laminar-position histogram ---
        nexttile; hold on;
        for c = 1:S8c_K
            histogram(xumAll(idx(clab==c)), xEdgesH, 'DisplayStyle','stairs', 'EdgeColor',colc(c,:), 'LineWidth',1.5);
        end
        box off; xlabel('laminar position (\mum)'); ylabel('# synapses');
        legend(arrayfun(@(c) sprintf('c%d (n=%d)',c,sum(clab==c)),1:S8c_K,'uni',0),'Box','off','Location','best');
        title('laminar distribution');
        fprintf('S8c %s | %s: cluster sizes = %s\n', setsC{si,2}, metsC{mi,2}, mat2str(accumarray(clab,1)'));
    end
end
sgtitle(sprintf('File %d: correlation clusters (rows: set x metric)', f));

%% ===================== Helpers (STA pipeline) =====================
% The S3/S4 helpers are now STANDALONE files (so they work when you Run Section):
%   ms2frames.m  clsName.m  newStacker.m  accumStack.m  flushStacker.m
%   stackInfo.m  loadSTAevent.m  mapV2Glu.m  (V-spike -> glu frame via match_nearest)
%   spatialInfo_mat.m  si_pvalue.m  interactive_glu_ROI.m  interactive_voltglu_ROI.m
%   interactive_glu_tuningcorr.m (S6c peak-angle + tuning-corr)  interactive_glu_tracecorr.m (S6d max-lag dFF corr)
%   sta_binmean.m  staGlu_prewin.m (S4b pre-spike glu tuning)  plot_VtrigGlu.m (S4b/S7c plots)
%   getS_glu.m  getGluCoord.m (compact S_glu accessors)  loadGlu.m
%   binMeanSem.m (S8 mean/SEM per bin)