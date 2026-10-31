%% Slice patch BTSP: EFS-evoked EPSPs before vs after the plateau potential (PP)
%  ChangHo Kim, whole-cell patch in slice, 20261030. One folder per DAQ acquisition, named
%  HHMMSS<tag> (acquisition start), e.g. 155910Cell2_EFS_15sec_probe_6V_I0_Probe_postBTSP
%  (cell 1 has no 'Cell' tag).
%  DAQ = Device_Data{3}.buffered_tasks, 100 kHz (found by deviceType: some early files hold 2 devices)
%    aif Dev2/ai0          amplifier output. Current clamp / I=0: 100 mV/V (rest ~ -0.6 V)
%    aof Dev2/ao0          current command, 400 pA/V (only in the PP and current-ramp trials)
%    dof Dev2/port0/line6  EFS trigger (0.5-1 ms TTL). Stimulus strength (V) is in the folder name,
%                          '5p5mV' = 5.5 V ('mV' in the names means V)
%    dof Dev2/Port0/Line15 camera trigger, not used
%  Cell 2 protocols
%    probe : EFS at 2.05, 8.05, 14.05 s, single pulse or paired pulse (50 ms), in I=0 or VC
%    PP    : 10 EFS pulses at 10 Hz (3.1-4.0 s), then 200 pA x 300 ms at 4.0 s, in IC
%    Iinj  : PulseRamp_CurrentInjection, no EFS: 10 increasing 100 ms current steps every 0.5 s
%            from 1 s, then a 5 s current ramp from 6.6 s
%  Trials are classified from the commands, not from the folder tag (154936..._IC_onePulse has
%  the PP command, so it counts as a PP trial). Clamp mode comes from the _I0/_IC/_VC tag; untagged
%  folders are VC when the amplifier output sits near 0 V.
%  Time zero is the current injection of the first PP trial (first BTSP); all times are relative to it.
%  Only the first pulse of a paired-pulse probe is shown and analysed.

clear; clc;
root     = '/Volumes/cohen_lab/Lab/Labmembers/ChangHo Kim/Data/20261030';
cellID   = 2;              % folders without a 'Cell#' tag are cell 1
ccScale  = 100;            % mV/V, amplifier output in current clamp / I=0
vcScale  = 400;  vcUnit = 'pA';  % pA/V, amplifier output in voltage clamp
aoScale  = 400;            % pA/V, current command
win      = [-20 150];      % ms, window around the first EFS pulse of each probe
blWin    = [-15 -1];       % ms, baseline window (before the stimulus artifact)
maxRest  = -50;            % mV, current-clamp sweeps with a baseline above this are left out
pkWin    = [1.5 30];       % ms, EPSC/EPSP peak search window (stimulus artifact ends by ~1 ms)
pkSmooth = 1;              % ms, gaussian smoothing before the peak search
spkWin   = [1 40];         % ms, window for an EFS-evoked spike
spikeThr = 0;              % mV, a current-clamp stimulation counts as a spike when Vm crosses this
spkProm  = 20;             % mV, current injection: min. prominence of an action potential
spkMaxW  = 10;             % ms, current injection: max. width at half prominence (broader = plateau)
dsFull   = 10;             % downsampling of the full-length traces for display (100 kHz -> 10 kHz)
ppXlim   = [2.8 5];        % s, zoom on the PP trials (EFS train + current step)
colPre   = [0 0.45 0.74];  colPP = [0.49 0.18 0.56];  colPost = [0.85 0.33 0.1];
colRec   = [0.47 0.67 0.19; 0.64 0.08 0.18; 0.93 0.69 0.13; 0.30 0.75 0.93];   % one per current-injection recording

%% Load the recordings of this cell and classify them from the DAQ commands
d = dir(root);
d = d([d.isdir] & ~cellfun(@isempty, regexp({d.name}, '^\d{6}', 'once')));
R = struct([]);
for i = 1:numel(d)
    nm  = d(i).name;  tag = nm(7:end);
    c = regexp(tag, 'Cell(\d+)', 'tokens', 'once');
    if isempty(c), c = 1; else, c = str2double(c{1}); end
    if c ~= cellID, continue; end
    f = fullfile(root, nm, 'output_data.mat');
    if ~isfile(f), fprintf('Skipping %s: no output_data.mat\n', nm); continue; end

    ws = warning('off', 'all');   % function handles of the rig (autofocus) cannot be resolved here
    S  = load(f, 'Device_Data');
    warning(ws);
    isDAQ = cellfun(@(x) isstruct(x) && isfield(x,'deviceType') && strcmp(x.deviceType,'DAQ'), S.Device_Data);
    bt  = S.Device_Data{find(isDAQ,1)}.buffered_tasks;
    ai  = getChan(bt, 'aif', 'ai0');
    [ao, aoCh] = getChan(bt, 'aof', 'ao0');
    efs = getChan(bt, 'dof', 'port0/line6');
    fs  = bt(strcmp({bt.task_type}, 'aif')).rate;

    if     contains(tag, '_VC'),  mode = 'VC';
    elseif contains(tag, '_IC'),  mode = 'IC';
    elseif contains(tag, '_I0'),  mode = 'I0';
    elseif abs(median(ai)) < 0.2, mode = 'VC';   % holding current ~0 V vs. rest ~ -0.6 V
    else,                         mode = 'CC';
    end
    isVC = strcmp(mode, 'VC');

    v = regexp(tag, '(\d+(?:p\d+)?)m?V(?=_|$|EFS)', 'tokens', 'once');
    if isempty(v), efsV = NaN; else, efsV = str2double(strrep(v{1}, 'p', '.')); end

    if isempty(efs), efs = zeros(size(ai)); end
    efsOn   = (find(diff([0; efs > 0.5]) == 1) - 1)/fs;   % s, pulse onsets
    trainOn = efsOn(diff([-Inf; efsOn]) > 1);              % s, first pulse of each probe / train
    hasAO   = ~isempty(ao) && any(abs(ao) > 1e-3);
    if     hasAO && ~isempty(efsOn), type = 'PP';
    elseif ~isempty(efsOn),          type = 'probe';
    elseif hasAO,                    type = 'Iinj';    % current injection only
    else,                            type = 'other';
    end

    % time of the first event (current injection, else first EFS pulse), s of the day
    tday = [3600 60 1]*str2double({nm(1:2); nm(3:4); nm(5:6)});
    if hasAO,                tFirst = tday + (find(abs(ao) > 1e-3, 1) - 1)/fs;
    elseif ~isempty(trainOn), tFirst = tday + trainOn(1);
    else,                    tFirst = tday;
    end

    y = ai*ccScale;  if isVC, y = ai*vcScale; end
    if isempty(trainOn), base = median(y(1:fs)); else, base = median(y(1:round(trainOn(1)*fs))); end

    r = struct('name', nm, 'tday', tday, 'tFirst', tFirst, 'mode', mode, 'isVC', isVC, ...
        'efsV', efsV, 'type', type, 'fs', fs, 'y', single(y), 'ao', single(ao*aoScale), 'aoCh', aoCh, ...
        'efs', efs > 0.5, 'efsOn', efsOn, 'trainOn', trainOn, ...
        'nPerTrain', numel(efsOn)/max(numel(trainOn),1), 'base', base);
    if ~hasAO, r.ao = []; end
    R = [R r];
end
[~, o] = sort([R.tday]);  R = R(o);

% pre / post = probes before the first and after the last PP trial; t = 0 at the first BTSP
isPP  = strcmp({R.type}, 'PP');
isPrb = strcmp({R.type}, 'probe');
if ~any(isPP), error('No PP trial found for cell %d', cellID); end
tPP   = [R(isPP).tday];
phase = repmat({''}, 1, numel(R));
phase(isPP) = {'PP'};
phase(isPrb & [R.tday] < min(tPP)) = {'pre'};
phase(isPrb & [R.tday] > max(tPP)) = {'post'};
phase(strcmp({R.type}, 'Iinj')) = {'Iinj'};
[R.phase] = phase{:};

% current commands sent to another channel than in the PP trials may not have reached the cell
ppCh = R(find(isPP, 1)).aoCh;
note = repmat({''}, 1, numel(R));
for k = find(~cellfun(@isempty, {R.ao}) & ~strcmpi({R.aoCh}, ppCh))
    note{k} = sprintf(' (cmd on %s)', R(k).aoCh);
    fprintf('%s: current command on %s, not %s as in the PP trials\n', R(k).name, R(k).aoCh, ppCh);
end
[R.note] = note{:};
t0    = R(find(isPP, 1)).tFirst;                       % s of the day, first BTSP current injection
tRel  = num2cell(([R.tFirst] - t0)/60);  [R.tRel] = tRel{:};   % min, first event of each recording
tPPend = R(find(isPP, 1, 'last')).tRel;                % min, last BTSP current injection

T = table(round([R.tRel]',2), string({R.phase})', string({R.type})', string({R.mode})', [R.efsV]', ...
    [R.nPerTrain]', round([R.base]',3), string({R.name})', 'VariableNames', ...
    {'Min_from_BTSP','Phase','Type','Mode','EFS_V','PulsesPerProbe','Baseline','Folder'});
disp(T);

%% Section 1: traces before vs after PP, with the command triggers
% Figures 1 (current clamp) and 2 (voltage clamp): full-length recordings, stacked in time order
% (earliest on top), with the EFS trigger (height = stimulus strength) and the current command below.
grp  = {'pre', 'PP', 'post'};
gTtl = {'Pre-PP probes', 'PP induction', 'Post-PP probes'};
gCol = {colPre, colPP, colPost};
clampSel = {~[R.isVC], [R.isVC]};
clampTtl = {'Current clamp (I=0 / IC)', 'Voltage clamp'};
unitOf   = {'mV', vcUnit};
for f = 1:2
    gUse = find(cellfun(@(g) any(strcmp({R.phase}, g) & clampSel{f}), grp));
    if isempty(gUse), continue; end
    gap = 1.1*max(arrayfun(@(r) diff(prctile(r.y, [0.05 99.95])), R(ismember({R.phase}, grp) & clampSel{f})));
    dockFig(f);
    tl = tiledlayout(4, numel(gUse), 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, sprintf('Cell %d, %s: full recordings (t = 0 at the first BTSP)', cellID, clampTtl{f}));
    for j = 1:numel(gUse)
        g  = gUse(j);
        Rs = R(strcmp({R.phase}, grp{g}) & clampSel{f});
        xl = [0 max(arrayfun(@(r) numel(r.y)/r.fs, Rs))];
        if strcmp(grp{g}, 'PP'), xl = ppXlim; end
        axT = nexttile(j, [3 1]);
        axC = nexttile(3*numel(gUse) + j);
        plotFull(axT, axC, Rs, shades(gCol{g}, numel(Rs)), xl, dsFull, gap, unitOf{f});
        title(axT, sprintf('%s (n = %d)', gTtl{g}, numel(Rs)));
    end
end

% Figure 3: current-injection trials (steps + ramp) with the current command, the F-I curve from
% the steps and the rheobase from the ramp, one color per recording. Spikes (marked on the traces)
% are peaks above spikeThr with a prominence > spkProm and a width at half prominence < spkMaxW,
% so the broad plateaus that follow the first spike of a step are not counted. Firing rate =
% spikes peaking within a step / step duration. Rheobase = injected current at the first spike on
% the ramp. Trials whose command went to another channel than in the PP trials (gray, no response)
% are left out of the F-I curve and rheobase.
Rs = R(strcmp({R.type}, 'Iinj'));
if ~isempty(Rs)
    nI   = numel(Rs);
    use  = cellfun(@isempty, {Rs.note});
    cols = repmat([0.6 0.6 0.6], nI, 1);
    cols(use,:) = colRec(mod(0:sum(use)-1, size(colRec,1)) + 1, :);

    FI = struct('I', {}, 'rate', {}, 'nSpk', {}, 'iSpk', {}, 'vSpk', {}, 'rheo', {}, 'tRheo', {}, 'slope', {});
    for k = 1:nI
        r = Rs(k);  V = double(r.y);  I = double(r.ao);
        [vSpk, iSpk] = findpeaks(V, 'MinPeakHeight', spikeThr, 'MinPeakProminence', spkProm, ...
            'MaxPeakWidth', spkMaxW/1e3*r.fs, 'WidthReference', 'halfprom', ...
            'MinPeakDistance', round(2e-3*r.fs));          % one peak per spike (flat tops give two)
        m  = abs(I) > 0.5;                                  % pA
        on = find(diff([0; m]) == 1);  off = find(diff([m; 0]) == -1);
        isRamp = arrayfun(@(s) max(I(on(s):off(s))) - min(I(on(s):off(s))) > 1, 1:numel(on));   % steps are flat
        stp  = find(~isRamp);  rmp = find(isRamp, 1);
        nSpk = arrayfun(@(s) sum(iSpk >= on(s) & iSpk <= off(s)), stp);
        FI(k).I     = arrayfun(@(s) median(I(on(s):off(s))), stp);
        FI(k).nSpk  = nSpk;
        FI(k).rate  = nSpk ./ ((off(stp) - on(stp) + 1)'/r.fs);
        FI(k).iSpk  = iSpk;  FI(k).vSpk = vSpk;
        FI(k).rheo  = NaN;   FI(k).tRheo = NaN;  FI(k).slope = NaN;
        if ~isempty(rmp)
            FI(k).slope = (I(off(rmp)) - I(on(rmp)))/((off(rmp) - on(rmp))/r.fs);   % pA/s
            j = iSpk(find(iSpk >= on(rmp) & iSpk <= off(rmp), 1));
            if ~isempty(j), FI(k).rheo = I(j);  FI(k).tRheo = (j-1)/r.fs; end
        end
    end
    Ti = table(round([Rs.tRel]',2), round([FI.slope]',1), round([FI.rheo]',1), ...
        string(arrayfun(@(f) mat2str(f.nSpk), FI, 'UniformOutput', false))', string({Rs.name})', ...
        'VariableNames', {'Min_from_BTSP','Ramp_pA_per_s','Rheobase_pA','SpikesPerStep','Folder'});
    disp(Ti);

    dockFig(3);
    tl = tiledlayout(4, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, sprintf('Cell %d: current injection, steps + ramp (t = 0 at the first BTSP)', cellID));
    axT = nexttile(1, [3 2]);  axC = nexttile(10, [1 2]);
    gap = 1.1*max(arrayfun(@(r) diff(prctile(r.y, [0.05 99.95])), Rs));
    plotFull(axT, axC, Rs, cols, [0 max(arrayfun(@(r) numel(r.y)/r.fs, Rs))], dsFull, gap, 'mV');
    title(axT, sprintf('Current injection (n = %d), \\nabla = detected spike', nI));
    for k = 1:nI
        plot(axT, (FI(k).iSpk - 1)/Rs(k).fs, FI(k).vSpk - (k-1)*gap + 0.04*gap, 'v', 'Color', cols(k,:), ...
            'MarkerFaceColor', cols(k,:), 'MarkerSize', 3);
        if use(k) && ~isnan(FI(k).rheo)
            plot(axC, FI(k).tRheo, FI(k).rheo, 'o', 'Color', cols(k,:), 'MarkerFaceColor', cols(k,:));
        end
    end

    nexttile(3, [2 1]);  hold on;
    for k = find(use)
        plot(FI(k).I, FI(k).rate, '-o', 'Color', cols(k,:), 'MarkerFaceColor', cols(k,:), ...
            'DisplayName', sprintf('%+.1f min', Rs(k).tRel));
        xline(FI(k).rheo, '--', 'Color', cols(k,:), 'HandleVisibility', 'off');
    end
    xlabel('Step current (pA)');  ylabel('Firing rate (Hz)');
    title({'F-I, 100 ms steps', 'dashed: rheobase from the ramp'});
    legend('Location', 'southeast');

    nexttile(9, [2 1]);  hold on;
    u = find(use);
    b = bar(1:numel(u), [FI(u).rheo], 'FaceColor', 'flat');
    b.CData = cols(u,:);
    for q = 1:numel(u)
        text(q, FI(u(q)).rheo, {sprintf('%.1f pA', FI(u(q)).rheo), sprintf('ramp %.0f pA/s', FI(u(q)).slope)}, ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    end
    set(gca, 'XTick', 1:numel(u), 'XTickLabel', compose('%+.1f min', [Rs(u).tRel]));
    ylabel('Rheobase (pA)');  title('Rheobase, first spike on the ramp');
    ylim([0 1.35*max([FI(u).rheo 1])]);
end

% Figures 4, 5, ...: probe sweeps aligned to the first EFS pulse, one figure per clamp mode x EFS
% strength recorded both before and after PP. Sweeps are baseline-subtracted and paired-pulse sweeps stop
% at the second pulse; current-clamp sweeps with a baseline above maxRest (depolarized cell) are
% drawn dotted gray and left out of the mean.
key   = strcat(string({R.mode}), " ", string([R.efsV]), " V");
key(~[R.isVC]) = strcat("CC ", string([R(~[R.isVC]).efsV]), " V");
conds = intersect(key(strcmp({R.phase}, 'pre')), key(strcmp({R.phase}, 'post')), 'stable');
for c = 1:numel(conds)
    dockFig(3 + c);
    tl = tiledlayout(4, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, sprintf('Cell %d, %s: EFS-evoked responses before vs after PP', cellID, conds(c)));
    vc = startsWith(conds(c), "VC");
    if vc, ylab = sprintf('\\DeltaI (%s)', vcUnit); else, ylab = '\DeltaV_m (mV)'; end
    axT = gobjects(1,3);  axC = gobjects(1,3);
    m = cell(1,2);  se = cell(1,2);  nOk = zeros(1,2);  Wj = cell(1,2);
    for j = 1:2
        Rs = R(key == conds(c) & strcmp({R.phase}, grp{2*j-1}));
        [W, tt] = getSweeps(Rs, win, blWin, maxRest);
        Wj{j} = W;
        cols = shades(gCol{2*j-1}, numel(Rs));
        axT(j) = nexttile(j, [3 1]);  hold on;
        for s = 1:numel(W)
            if W(s).ok, plot(tt, W(s).y, 'Color', cols(W(s).rec,:), 'LineWidth', 0.5);
            else,       plot(tt, W(s).y, ':', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.5); end
        end
        Ym = vertcat(W([W.ok]).y);
        nOk(j) = size(Ym,1);
        m{j}  = mean(Ym, 1, 'omitnan');
        se{j} = std(Ym, 0, 1, 'omitnan')./sqrt(sum(~isnan(Ym), 1));
        plot(tt, m{j}, 'Color', 0.6*gCol{2*j-1}, 'LineWidth', 2);
        ylabel(ylab);
        title({gTtl{2*j-1}, sprintf('%d rec, %+.1f to %+.1f min, %d/%d sweeps', numel(Rs), ...
            min([Rs.tRel]), max([Rs.tRel]), nOk(j), numel(W))});

        axC(j) = nexttile(9 + j);  hold on;
        for s = 1:numel(W), plot(tt, W(s).efs, 'Color', cols(W(s).rec,:)); end
        ylabel('EFS (V)');  xlabel('Time from EFS (ms)');
    end

    axT(3) = nexttile(3, [3 1]);  hold on;
    h = gobjects(1,2);
    for j = 1:2
        v = find(~isnan(m{j}));    % paired-pulse means stop at the 2nd pulse
        if isempty(v), continue; end
        fill([tt(v) fliplr(tt(v))], [m{j}(v)+se{j}(v) fliplr(m{j}(v)-se{j}(v))], gCol{2*j-1}, ...
            'FaceAlpha', 0.25, 'EdgeColor', 'none');
        h(j) = plot(tt, m{j}, 'Color', gCol{2*j-1}, 'LineWidth', 1.5);
    end
    legend(h, {sprintf('pre (n = %d)', nOk(1)), sprintf('post (n = %d)', nOk(2))}, 'Location', 'northeast');
    ylabel(ylab);  title({'Pre vs post', 'mean \pm SEM'});

    axC(3) = nexttile(12);  hold on;
    for j = 1:2
        for s = 1:numel(Wj{j}), plot(tt, Wj{j}(s).efs, 'Color', gCol{2*j-1}); end
    end
    ylabel('EFS (V)');  xlabel('Time from EFS (ms)');
    linkaxes([axT axC], 'x');  xlim(axT(1), win);
    Yok = [vertcat(Wj{1}([Wj{1}.ok]).y); vertcat(Wj{2}([Wj{2}.ok]).y)];   % excluded sweeps do not set the range
    set(axT, 'YLim', [min(Yok(:)) max(Yok(:))] + [-0.05 0.05]*(max(Yok(:)) - min(Yok(:))));
    vMax = max([R(key == conds(c)).efsV]);
    set(axC, 'YLim', [-0.5 1.2*vMax], 'YTick', [0 vMax]);
end

%% Section 2: EPSC / EPSP peak and spike probability vs time from the first BTSP
% First EFS pulse of every probe. Peak = extreme of the smoothed (pkSmooth), baseline-subtracted
% response in pkWin: min for the VC EPSC (reported as amplitude = -peak), max for the CC EPSP.
% A CC stimulation counts as a spike when Vm crosses spikeThr within spkWin; spiking stimulations
% are left out of the EPSP peak. CC sweeps from a depolarized cell (baseline > maxRest) are left
% out of both.
P = R(ismember({R.phase}, {'pre', 'post'}));
[W, tt] = getSweeps(P, win, blWin, maxRest);
Y    = vertcat(W.y);
Ys   = smoothdata(Y, 2, 'gaussian', round(pkSmooth/1e3*P(1).fs));
inPk = tt >= pkWin(1) & tt <= pkWin(2);
rec  = [W.rec]';
vcW  = [P(rec).isVC]';
ok   = [W.ok]';
tMin = ([W.tPulse]' - t0)/60;
spk  = ~vcW & max(Y(:, tt >= spkWin(1) & tt <= spkWin(2)) + [W.bl]', [], 2) > spikeThr;
amp  = max(Ys(:, inPk), [], 2);
amp(vcW) = -min(Ys(vcW, inPk), [], 2);
amp(~ok | spk) = NaN;

% per recording
nP = numel(P);
pSpk = NaN(nP,1);  ampM = NaN(nP,1);  ampS = NaN(nP,1);  nAmp = zeros(nP,1);  nOk = zeros(nP,1);
for k = 1:nP
    sel = rec == k & ok;
    nOk(k)  = sum(sel);
    nAmp(k) = sum(sel & ~isnan(amp));
    ampM(k) = mean(amp(sel), 'omitnan');
    ampS(k) = std(amp(sel), 'omitnan');
    if ~P(k).isVC, pSpk(k) = mean(spk(sel)); end
end
tRec = arrayfun(@(k) mean(tMin(rec == k)), (1:nP)');
Tm = table(round(tRec,2), string({P.phase})', string({P.mode})', [P.efsV]', nOk, nAmp, ...
    round(ampM,4), round(ampS,4), pSpk, string({P.name})', 'VariableNames', ...
    {'Min_from_BTSP','Phase','Mode','EFS_V','nSweeps','nPeak','Peak','PeakSD','pSpike','Folder'});
disp(Tm);

dockFig(4 + numel(conds));
tl = tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf('Cell %d: response to the first EFS pulse vs time from the first BTSP', cellID));
efsList = unique([P.efsV]);
mk  = {'^', 'o', 's', 'd', 'v', 'p'};
pnl = {vcW, ~vcW, ~vcW};
ylb = {sprintf('EPSC (%s)', vcUnit), 'EPSP (mV)', 'P(spike)'};
ttl = {'Voltage clamp: EPSC peak', 'Current clamp: EPSP peak (stimulations that fired a spike omitted)', ...
       'Current clamp: spike probability'};
ax = gobjects(1,3);
for p = 1:3
    ax(p) = nexttile;  hold on;
    xregion(0, max(tPPend, 0.05), 'FaceColor', colPP, 'FaceAlpha', 0.15, 'EdgeColor', 'none');
    xline(0, '--', 'Color', colPP);
    if p < 3, val = amp; else, val = double(spk); end
    for k = find(arrayfun(@(r) r.isVC == (p == 1), P))
        col = colPre;  if strcmp(P(k).phase, 'post'), col = colPost; end
        sym = mk{find(efsList == P(k).efsV, 1)};
        sel = rec == k & ok;
        plot(tMin(sel), val(sel), sym, 'Color', col, 'MarkerSize', 4);
        if p < 3
            if nAmp(k) > 0
                errorbar(tRec(k), ampM(k), ampS(k), sym, 'Color', col, 'MarkerFaceColor', col, ...
                    'MarkerSize', 8, 'LineWidth', 1);
            end
        else
            plot(tRec(k), pSpk(k), sym, 'Color', col, 'MarkerFaceColor', col, 'MarkerSize', 8);
        end
    end
    if p == 3, ylim([-0.1 1.1]); end
    if p == 2   % recordings with no EPSP left because every stimulation fired a spike
        yl = ylim;
        for k = find(~[P.isVC] & nOk' > 0 & nAmp' == 0)
            text(tRec(k), yl(1) + 0.05*diff(yl), 'all spiked', 'Rotation', 90, 'FontSize', 8, 'Color', colPost);
        end
    end
    ylabel(ylb{p});  title(ttl{p});
end
hl = gobjects(1, numel(efsList));
for e = 1:numel(efsList)
    hl(e) = plot(ax(3), NaN, NaN, mk{e}, 'Color', [0.4 0.4 0.4], 'MarkerFaceColor', [0.4 0.4 0.4]);
end
legend(hl, compose('EFS %g V', efsList), 'Location', 'northwest');
linkaxes(ax, 'x');
xlabel(ax(3), 'Time from the first BTSP (min)');

%% Local functions
function dockFig(k)
% Figure k, cleared, docked and in the light theme
figure(k);  clf;
set(gcf, 'WindowStyle', 'docked');
theme(gcf, 'light');
end

function [x, ch] = getChan(bt, type, port)
% Data and name of the first channel of a task type whose physical channel contains port ([] if absent)
x = [];  ch = '';
for k = find(strcmp({bt.task_type}, type))
    for c = 1:numel(bt(k).channels)
        if contains(bt(k).channels(c).phys_channel, port, 'IgnoreCase', true)
            x  = double(bt(k).channels(c).data(:));
            ch = strtrim(bt(k).channels(c).phys_channel);  return
        end
    end
end
end

function cols = shades(base, n)
% n shades of base, from light (earliest) to the full color (latest)
cols = base + (1 - base).*linspace(0.5, 0, n)';
end

function [W, tt] = getSweeps(Rs, win, blWin, maxRest)
% One baseline-subtracted sweep per probe, aligned to its first EFS pulse. Sweep and EFS command
% are cut 0.5 ms before the next pulse within the window (paired pulse: first pulse only).
W = struct('rec', {}, 'y', {}, 'efs', {}, 'bl', {}, 'tPulse', {}, 'ok', {});
tt = [];
for k = 1:numel(Rs)
    r  = Rs(k);
    w  = round(win(1)/1e3*r.fs):round(win(2)/1e3*r.fs);
    tt = w/r.fs*1e3;
    for p = 1:numel(r.trainOn)
        i0 = round(r.trainOn(p)*r.fs) + 1;
        if i0 + w(1) < 1 || i0 + w(end) > numel(r.y), continue; end
        y  = double(r.y(i0 + w))';    % row
        e  = double(r.efs(i0 + w))'*r.efsV;
        bl = median(y(tt >= blWin(1) & tt <= blWin(2)));
        nxt = r.efsOn(r.efsOn > r.trainOn(p) + 1e-4) - r.trainOn(p);
        if ~isempty(nxt)
            cut = tt >= nxt(1)*1e3 - 0.5;
            y(cut) = NaN;  e(cut) = NaN;
        end
        W(end+1) = struct('rec', k, 'y', y - bl, 'efs', e, 'bl', bl, 'tPulse', r.tday + r.trainOn(p), ...
            'ok', r.isVC || bl < maxRest); %#ok<AGROW>
    end
end
end

function plotFull(axT, axC, Rs, cols, xl, ds, gap, unit)
% Traces stacked gap apart (earliest on top) with a scale bar, and below them the commands:
% EFS trigger x stimulus strength (left axis, V) and current command (right axis, pA).
hold(axT, 'on');
lab = cell(1, numel(Rs));
for k = 1:numel(Rs)
    r = Rs(k);  idx = 1:ds:numel(r.y);  t = (idx-1)/r.fs;
    plot(axT, t, r.y(idx) - (k-1)*gap, 'Color', cols(k,:), 'LineWidth', 0.5);
    if any(r.efs), stim = sprintf('%gV\\times%g', r.efsV, r.nPerTrain);
    else,          stim = sprintf('max %.0f pA', max(r.ao)); end
    lab{k} = sprintf('\\color[rgb]{%.2f,%.2f,%.2f}%+.1f min  %s %s%s', 0.75*cols(k,:), r.tRel, r.mode, stim, r.note);
end
% each trace is labelled at its baseline on the y axis: time from the first BTSP, clamp mode,
% EFS strength x pulses per probe (or the largest injected current)
yb = [Rs.base] - (0:numel(Rs)-1)*gap;
[yb, o] = sort(yb);
set(axT, 'YTick', yb, 'YTickLabel', lab(o), 'TickLabelInterpreter', 'tex', 'TickLength', [0 0]);
axT.YAxis.FontSize = 8;  axT.XTickLabel = [];
xlim(axT, xl);  axis(axT, 'tight');  xlim(axT, xl);
yl = ylim(axT);
cand = [1 2 5]'*10.^(-4:3);  sb = max(cand(cand <= gap/3));
ylim(axT, [yl(1) - 1.3*sb, yl(2) + 0.03*diff(yl)]);    % scale bar below the lowest trace
plot(axT, xl(2)*[1 1] - 0.01*diff(xl), yl(1) - [1.2 0.2]*sb, 'Color', axT.XColor, 'LineWidth', 2);
text(axT, xl(2) - 0.02*diff(xl), yl(1) - 0.7*sb, sprintf('%g %s', sb, unit), 'HorizontalAlignment', 'right');

hold(axC, 'on');
anyAO  = any(arrayfun(@(r) ~isempty(r.ao), Rs));
anyEFS = any(arrayfun(@(r) any(r.efs), Rs));
if anyEFS
    if anyAO, yyaxis(axC, 'left'); end
    for k = 1:numel(Rs)
        r = Rs(k);  idx = 1:ds:numel(r.y);
        plot(axC, (idx-1)/r.fs, double(r.efs(idx))*max(r.efsV, 0), '-', 'Color', cols(k,:));
    end
    ylabel(axC, 'EFS (V)');  ylim(axC, [-0.5 1.2*max([Rs.efsV 1])]);
end
if anyAO
    if anyEFS, yyaxis(axC, 'right'); end
    for k = find(arrayfun(@(r) ~isempty(r.ao), Rs))
        idx = 1:ds:numel(Rs(k).ao);
        plot(axC, (idx-1)/Rs(k).fs, Rs(k).ao(idx), '-', 'Color', cols(k,:), 'LineWidth', 1);
    end
    ylabel(axC, 'I_{cmd} (pA)');
    if anyEFS, axC.YAxis(1).Color = axC.XColor;  axC.YAxis(2).Color = axC.XColor; end
end
xlabel(axC, 'Time in recording (s)');
linkaxes([axT axC], 'x');  xlim(axT, xl);
end
