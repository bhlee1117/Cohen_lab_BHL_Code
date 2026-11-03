%% Slice patch inhibitory BTSP: SST -> CA1 pyramidal IPSC before vs after the plateau potential (PP)
%  ChangHo Kim, 20261101. ChR2 in SST interneurons, whole-cell patch of a CA1 pyramidal neuron.
%  One folder per DAQ acquisition, named HHMMSS<tag> (acquisition start), e.g.
%  174200Cell3_VC_IPSC_Probe_burst1. 'burst#' is the repetition number.
%  DAQ = Device_Data{3}.buffered_tasks, 100 kHz (found by deviceType)
%    aif Dev2/ai0           amplifier output. VC: 400 pA/V, IC: 100 mV/V
%    aof Dev2/ao0           amplifier command. VC: holding step, 20 mV/V; IC: current, 400 pA/V
%    aof Dev2/ao2           488 nm laser modulation (V)
%    dof Dev2/port0/line11  shutter488;  dof Dev2/port0/line8  enable488
%                           light = mod488 x shutter488 x enable488
%    dof Dev2/port0/line6   EFS trigger, ignored: there was no EFS electrode in this experiment
%    dof Dev2/Port0/Line15  camera trigger, not used
%  Cell 3 protocols
%    probe     (VC): +20 mV holding step for 2 s at 1.5, 11.5, 21.5 s with a 10 ms light pulse 1 s
%                    into each step, where the SST IPSC is outward
%    PP        (IC): 300 ms current step at 4.0 s, alone ('PP_induction', 200 or 400 pA) or after
%                    10 x 1 ms light pulses at 10 Hz (3.1-4.0 s, 'wSSTstimulation', 200 pA)
%    titration (VC): before the probes, a light pulse every 1 s alternating between the resting holding
%                    potential and a 1 s depolarizing step: constant light with steps of +10 to +26 mV,
%                    or +20 mV with increasing light (Section 3)
%  During the +20 mV step the cell fires unclamped action currents (inward, ~ -200 pA), more often
%  after PP. They are detected, and a probe with one near the light pulse (within win) is omitted.
%  Time zero is the current injection of the first PP trial (first BTSP); all times refer to it.

clear; clc;
root     = '/Volumes/cohen_lab/Lab/Labmembers/ChangHo Kim/Data/20261101';
cellID   = 3;
ccScale  = 100;            % mV/V, amplifier output in current clamp
vcScale  = 400;            % pA/V, amplifier output in voltage clamp
icCmd    = 400;            % pA/V, ao0 current command in current clamp
vcCmd    = 20;             % mV/V, ao0 holding command in voltage clamp
lpCut    = 2000;           % Hz, low-pass of the voltage-clamp current (display and measurement)
win      = [-50 60];       % ms, window around each light pulse; a probe with an action current in it is omitted
blWin    = [-20 -1];       % ms, IPSC baseline
ipscWin  = [2 50];         % ms, IPSC peak search window (outward at +20 mV)
ivWin    = [5 20];         % ms, titration IPSC = mean current in this window (unbiased at small amplitudes)
acProm   = 80;             % pA, action current = inward peak with at least this prominence
dsFull   = 10;             % downsampling of the full-length traces for display (100 kHz -> 10 kHz)
ppXlim   = [2.8 5];        % s, zoom on the PP trials (light pulses + current step)
colPre   = [0 0.45 0.74];  colPost = [0.85 0.33 0.1];
colPP    = [0.49 0.18 0.56];  colSST = [0.47 0.67 0.19];   % PP alone / paired with SST light
colLight = [0.2 0.55 1];   % light command

%% Load the recordings of this cell and classify them from the DAQ commands
d = dir(root);
d = d([d.isdir] & ~cellfun(@isempty, regexp({d.name}, sprintf('^\\d{6}Cell%d_', cellID), 'once')));
R = struct([]);
for i = 1:numel(d)
    nm  = d(i).name;  tag = nm(7:end);
    f = fullfile(root, nm, 'output_data.mat');
    if ~isfile(f), fprintf('Skipping %s: no output_data.mat\n', nm); continue; end

    ws = warning('off', 'all');   % function handles of the rig (autofocus) cannot be resolved here
    S  = load(f, 'Device_Data');
    warning(ws);
    isDAQ = cellfun(@(x) isstruct(x) && isfield(x,'deviceType') && strcmp(x.deviceType,'DAQ'), S.Device_Data);
    bt  = S.Device_Data{find(isDAQ,1)}.buffered_tasks;
    ai  = getChan(bt, 'aif', 'ai0');
    ao0 = getChan(bt, 'aof', 'ao0');
    mod = getChan(bt, 'aof', 'ao2');
    sh  = getChan(bt, 'dof', 'line11');
    en  = getChan(bt, 'dof', 'line8');
    fs  = bt(strcmp({bt.task_type}, 'aif')).rate;
    n   = numel(ai);
    if isempty(ao0), ao0 = zeros(n,1); end
    light = zeros(n,1);
    if ~isempty(mod)
        light = mod;
        if ~isempty(sh), light = light.*sh; end
        if ~isempty(en), light = light.*en; end
    end

    isVC = contains(tag, '_VC');
    if isVC, y = lowpass(ai*vcScale, lpCut, fs);  cmd = ao0*vcCmd;   % pA, mV
    else,    y = ai*ccScale;                      cmd = ao0*icCmd;   % mV, pA
    end
    lightOn = (find(diff([0; light > 0.5]) == 1) - 1)/fs;   % s, light pulse onsets

    if ~isVC && any(cmd > 1)
        if ~isempty(lightOn), type = 'PP-SST'; else, type = 'PP'; end
        label = sprintf('%s %.0f pA', type, max(cmd));
    elseif isVC && contains(tag, 'titration'), type = 'titration';  label = type;
    elseif isVC && ~isempty(lightOn),          type = 'probe';      label = 'probe';
    else,                                      type = 'other';      label = type;
    end

    % time of the first event (current injection, else first light pulse), s of the day
    tday = [3600 60 1]*str2double({nm(1:2); nm(3:4); nm(5:6)});
    if startsWith(type, 'PP'),    tFirst = tday + (find(cmd > 1, 1) - 1)/fs;
    elseif ~isempty(lightOn),     tFirst = tday + lightOn(1);
    else,                         tFirst = tday;
    end

    r = struct('name', nm, 'tday', tday, 'tFirst', tFirst, 'isVC', isVC, 'type', type, 'label', label, ...
        'fs', fs, 'y', single(y), 'cmd', single(cmd), 'light', single(light), ...
        'lightOn', lightOn, 'base', median(y(1:fs)));
    R = [R r];
end
[~, o] = sort([R.tday]);  R = R(o);

% pre / post = probes before the first and after the last PP trial; t = 0 at the first BTSP
isPP  = startsWith({R.type}, 'PP');
isPrb = strcmp({R.type}, 'probe');
if ~any(isPP), error('No PP trial found for cell %d', cellID); end
tPP   = [R(isPP).tday];
phase = {R.type};
phase(isPrb & [R.tday] < min(tPP)) = {'pre'};
phase(isPrb & [R.tday] > max(tPP)) = {'post'};
[R.phase] = phase{:};
t0   = R(find(isPP, 1)).tFirst;                        % s of the day, first BTSP current injection
tRel = num2cell(([R.tFirst] - t0)/60);  [R.tRel] = tRel{:};    % min, first event of each recording
epPP  = [min([R(strcmp({R.type}, 'PP')).tRel]) max([R(strcmp({R.type}, 'PP')).tRel])];
epSST = [min([R(strcmp({R.type}, 'PP-SST')).tRel]) max([R(strcmp({R.type}, 'PP-SST')).tRel])];

T = table(round([R.tRel]',2), string({R.phase})', string({R.label})', round([R.base]',1), ...
    string({R.name})', 'VariableNames', {'Min_from_BTSP','Phase','Protocol','Baseline','Folder'});
disp(T);

%% Section 1: traces before vs after PP, with the commands
% Figure 1: VC probes, full recordings (pre | post), stacked in time order with the holding step,
% light commands below. Figure 2: PP induction trials (PP alone | paired with SST light).
grpF = {{'pre', 'post'}, {'PP', 'PP-SST'}};
ttlF = {{'Pre-PP probes', 'Post-PP probes'}, {'PP alone', 'PP paired with SST light'}};
for f = 1:2
    g   = grpF{f};
    sel = ismember({R.phase}, g);
    if ~any(sel), continue; end
    vc  = R(find(sel, 1)).isVC;
    if vc, unit = 'pA'; cmdLab = '\DeltaV_{hold} (mV)'; else, unit = 'mV'; cmdLab = 'I_{cmd} (pA)'; end
    gap = 1.1*max(arrayfun(@(r) diff(prctile(r.y, [0.05 99.95])), R(sel)));
    dockFig(f);
    tl = tiledlayout(4, numel(g), 'TileSpacing', 'compact', 'Padding', 'compact');
    if vc, ttl = 'VC probes, full recordings'; else, ttl = 'PP induction trials'; end
    title(tl, sprintf('Cell %d: %s (t = 0 at the first BTSP)', cellID, ttl));
    for j = 1:numel(g)
        Rs = R(strcmp({R.phase}, g{j}));
        if isempty(Rs), continue; end
        xl = [0 max(arrayfun(@(r) numel(r.y)/r.fs, Rs))];
        if ~vc, xl = ppXlim; end
        axT = nexttile(j, [3 1]);
        axC = nexttile(3*numel(g) + j);
        plotFull(axT, axC, Rs, recCols(numel(Rs)), xl, dsFull, gap, unit, cmdLab, colLight);
        title(axT, sprintf('%s (n = %d)', ttlF{f}{j}, numel(Rs)));
    end
end

% Light-evoked sweeps of all probes. A probe with an action current within win of the light pulse
% is omitted from the IPSC analysis (Figures 3 and 5); Figure 4 shows which ones and why.
P = R(ismember({R.phase}, {'pre', 'post'}));
[Wall, tt] = getLightSweeps(P, win, blWin, ipscWin, ivWin, acProm);
postAll = strcmp({P([Wall.rec]).phase}, 'post');
W       = Wall([Wall.clean]);
isPost  = postAll([Wall.clean]);
nAll    = [sum(~postAll) sum(postAll)];
fprintf('Probes kept: %d/%d pre, %d/%d post\n', sum(~isPost), nAll(1), sum(isPost), nAll(2));

% Figure 3: kept probes aligned to the light pulse (baseline-subtracted, low-passed at lpCut), one
% color per recording and the mean in black; right: pre vs post mean +- SEM.
dockFig(3);
tl = tiledlayout(4, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf('Cell %d: SST light-evoked IPSC at +20 mV, before vs after PP', cellID));
axT = gobjects(1,3);  axC = gobjects(1,3);
m = cell(1,2);  se = cell(1,2);  nOk = zeros(1,2);
phs = {'pre', 'post'};  phCol = {colPre, colPost};  phTtl = {'Pre-PP probes', 'Post-PP probes'};
for j = 1:2
    sw   = find(isPost == (j == 2));
    recs = unique([W(sw).rec]);
    cols = recCols(numel(recs));
    nOk(j) = numel(sw);
    axT(j) = nexttile(j, [3 1]);  hold on;
    for s = sw, plot(tt, W(s).y, 'Color', cols(recs == W(s).rec,:), 'LineWidth', 0.75); end
    if nOk(j) > 0
        Ym = vertcat(W(sw).y);
        m{j}  = mean(Ym, 1);
        se{j} = std(Ym, 0, 1)/sqrt(nOk(j));
        plot(tt, m{j}, 'k', 'LineWidth', 2);
        title({phTtl{j}, sprintf('%d/%d probes kept, from %d rec (%+.1f to %+.1f min)', nOk(j), nAll(j), ...
            numel(recs), min([P(recs).tRel]), max([P(recs).tRel]))});
    else
        title({phTtl{j}, sprintf('0/%d probes kept', nAll(j))});
    end
    ylabel('\DeltaI (pA)');
    axC(j) = nexttile(9 + j);  hold on;
    for s = sw, plot(tt, W(s).light, 'Color', colLight); end
    ylabel('Light (V)');  xlabel('Time from light (ms)');
end
axT(3) = nexttile(3, [3 1]);  hold on;
h = gobjects(0);  lg = {};
for j = find(nOk > 0)
    fill([tt fliplr(tt)], [m{j}+se{j} fliplr(m{j}-se{j})], phCol{j}, 'FaceAlpha', 0.25, 'EdgeColor', 'none');
    h(end+1) = plot(tt, m{j}, 'Color', phCol{j}, 'LineWidth', 1.5);  %#ok<SAGROW>
    lg{end+1} = sprintf('%s (n = %d)', phs{j}, nOk(j));               %#ok<SAGROW>
end
legend(h, lg, 'Location', 'northeast');
ylabel('\DeltaI (pA)');  title({'Pre vs post', 'mean \pm SEM'});
axC(3) = nexttile(12);  hold on;
for s = 1:numel(W), plot(tt, W(s).light, 'Color', colLight); end
ylabel('Light (V)');  xlabel('Time from light (ms)');
linkaxes([axT axC], 'x');  xlim(axT(1), win);
Yok = vertcat(W.y);
set(axT, 'YLim', [min(Yok(:)) max(Yok(:))] + [-0.1 0.1]*(max(Yok(:)) - min(Yok(:))));
set(axC, 'YLim', [-0.5 1.2*max([Wall.lightMax 1])]);

% Figure 4: which probes were omitted. Action currents during the +20 mV step, one row per probe
% in time order; the window that decides is shaded, omitted probes are gray.
dockFig(4);
tl = tiledlayout(1, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf('Cell %d: unclamped action currents during the +20 mV step', cellID));
ax1 = nexttile;  hold on;
xregion(win(1), win(2), 'FaceColor', [0.5 0.5 0.5], 'FaceAlpha', 0.15);
for s = 1:numel(Wall)
    c = colPre;  if postAll(s), c = colPost; end
    if ~Wall(s).clean, c = [0.65 0.65 0.65]; end
    plot([Wall(s).acT; Wall(s).acT], s + [-0.4; 0.4], 'Color', c, 'LineWidth', 1);
end
xline(0, '-', 'Color', colLight, 'LineWidth', 1.5);
lab = compose("%+.1f min", ([Wall.tLight] - t0)/60);
lab(~[Wall.clean]) = strcat(lab(~[Wall.clean]), " (omitted)");
set(ax1, 'YDir', 'reverse', 'YTick', 1:numel(Wall), 'YTickLabel', lab, 'FontSize', 7);
ylim([0.5 numel(Wall) + 0.5]);  xlim([-1000 1000]);
xlabel('Time from light (ms)');  ylabel('Probe (time from the first BTSP)');
title(sprintf('Kept %d/%d pre, %d/%d post; light at 0 ms (blue), omission window %d to %d ms (shaded)', ...
    sum(~isPost), nAll(1), sum(isPost), nAll(2), win));

%% Section 2: IPSC amplitude vs time from the first BTSP
% IPSC = peak outward current in ipscWin minus the baseline, kept probes only.
tMin  = ([W.tLight] - t0)/60;
amp   = [W.amp];
rec   = [W.rec];
nP    = numel(P);
tRec  = arrayfun(@(k) mean(([Wall([Wall.rec] == k).tLight] - t0)/60), 1:nP);
nKept = arrayfun(@(k) sum(rec == k), 1:nP);
nTot  = arrayfun(@(k) sum([Wall.rec] == k), 1:nP);
ampM  = arrayfun(@(k) mean(amp(rec == k)), 1:nP);
ampSE = arrayfun(@(k) std(amp(rec == k))/sqrt(sum(rec == k)), 1:nP);

Tm = table(round(tRec',2), string({P.phase})', compose('%d/%d', [nKept' nTot']), round(ampM',1), ...
    round(ampSE',1), string({P.name})', 'VariableNames', ...
    {'Min_from_BTSP','Phase','Kept','IPSC_pA','IPSC_SEM','Folder'});
disp(Tm);
for j = 1:2
    a = amp(isPost == (j == 2));
    fprintf('%-4s: IPSC %.1f +- %.1f pA (mean +- SD, n = %d probes)\n', phs{j}, mean(a), std(a), numel(a));
end

dockFig(5);
tl = tiledlayout(1, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf('Cell %d: SST IPSC at +20 mV vs time from the first BTSP', cellID));
ax = nexttile;  hold on;
xregion(epPP(1), epPP(2), 'FaceColor', colPP, 'FaceAlpha', 0.15);
xregion(epSST(1), epSST(2), 'FaceColor', colSST, 'FaceAlpha', 0.2);
for k = 1:nP    % small dots: probes; large markers: recording mean +- SEM
    c = colPre;  if strcmp(P(k).phase, 'post'), c = colPost; end
    plot(tMin(rec == k), amp(rec == k), '.', 'Color', 0.4*c + 0.6, 'MarkerSize', 10);
    errorbar(tRec(k), ampM(k), ampSE(k), 'o', 'Color', c, 'MarkerFaceColor', c, 'MarkerSize', 7, 'LineWidth', 1);
end
ylabel('IPSC (pA)');  xlabel('Time from the first BTSP (min)');
title('IPSC peak, kept probes (numbers: probes kept / recorded)');
ylim([0 1.3*max([amp 1])]);  yl = ylim;
for k = 1:nP
    c = colPre;  if strcmp(P(k).phase, 'post'), c = colPost; end
    text(tRec(k), yl(1), sprintf(' %d/%d', nKept(k), nTot(k)), 'Color', c, 'FontSize', 8, 'Rotation', 90);
end
text(mean(epPP), yl(2), 'PP', 'Color', colPP, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
text(mean(epSST), yl(2), 'PP + SST light', 'Color', colSST, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');

%% Section 3: IPSC vs holding step (titration)
% IPSC = mean current in ivWin minus the baseline (the peak is biased upward by noise at these small
% amplitudes); pulses with an action current within win are omitted (gray dotted), pulses at the
% resting holding potential are dashed. A line fitted to IPSC vs step in the constant-light recordings
% gives the conductance (slope, pA/mV = nS, with the p-value of slope = 0) and the reversal potential
% relative to the resting holding potential (x-intercept); the IPSC grows with the step if its driving
% force does. R_step = step / change of the holding current, a check of the access to the cell.
Q = R(strcmp({R.type}, 'titration'));
if ~isempty(Q)
    V    = getLightSweeps(Q, win, blWin, ipscWin, ivWin, acProm);
    nQ   = numel(Q);
    vrec = [V.rec];  vh = [V.hold];  vl = [V.lightV];  va = [V.ampMean];  vok = [V.clean];
    gS = NaN(nQ,1);  pS = NaN(nQ,1);  Erev = NaN(nQ,1);  Rstep = NaN(nQ,1);
    fitIV = false(nQ,1);  byLight = false(nQ,1);  lightDur = NaN(nQ,1);
    for k = 1:nQ
        s = vrec == k & vok;
        byLight(k) = max(vl(vrec == k)) - min(vl(vrec == k)) > 0.05;
        fitIV(k)   = ~byLight(k) && numel(unique(vh(s))) >= 3;
        if numel(unique(vh(s))) >= 2, Rstep(k) = 1e3/linFit(vh(s), [V(s).bl]); end   % mV/pA = GOhm
        if fitIV(k)
            [gS(k), b, pS(k)] = linFit(vh(s), va(s));
            if pS(k) < 0.05, Erev(k) = -b/gS(k); end
        end
        L = Q(k).light > 0.5;
        lightDur(k) = (find(diff([L; 0]) == -1, 1) - find(diff([0; L]) == 1, 1) + 1)/Q(k).fs*1e3;
    end
    Tq = table(round([Q.tRel]',2), compose('%.1f-%.1f V, %.0f ms', [arrayfun(@(k) min(vl(vrec == k)), 1:nQ)' ...
        arrayfun(@(k) max(vl(vrec == k)), 1:nQ)' lightDur]), ...
        string(arrayfun(@(k) mat2str(unique(vh(vrec == k))), (1:nQ)', 'UniformOutput', false)), ...
        compose('%d/%d', [arrayfun(@(k) sum(vrec == k & vok), 1:nQ)' arrayfun(@(k) sum(vrec == k), 1:nQ)']), ...
        round(Rstep), round(gS,2), round(pS,3), round(Erev,1), string({Q.name})', 'VariableNames', ...
        {'Min_from_BTSP','Light','Steps_mV','Kept','Rstep_MOhm','g_nS','p_slope','Erev_rel_mV','Folder'});
    disp(Tq);

    % Figure 6: titration sweeps aligned to the light pulse, one tile per recording, colored by the
    % holding step (constant light) or by the light intensity (light series); light command below.
    dockFig(6);
    tl = tiledlayout(4, nQ, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, sprintf('Cell %d: titration, light-evoked IPSC at different holding steps', cellID));
    cmHold  = interp1([0 1], [0.30 0.45 0.80; 0.80 0.10 0.10], linspace(0, 1, 256));   % resting -> largest step
    cmLight = interp1([0 1], [0.70 0.80 1.00; 0.05 0.20 0.55], linspace(0, 1, 256));   % dim -> bright
    axT = gobjects(1,nQ);  axC = gobjects(1,nQ);
    for k = 1:nQ
        s = find(vrec == k);
        if byLight(k), cm = cmLight; z = vl(s); zl = 'Light (V)';
        else,          cm = cmHold;  z = vh(s); zl = '\DeltaV_{hold} (mV)'; end
        zr  = [min(z) max(z) + eps];
        col = interp1(linspace(zr(1), zr(2), 256), cm, z);
        axT(k) = nexttile(k, [3 1]);  hold on;
        for q = 1:numel(s)
            ls = '-';  if vh(s(q)) < 1, ls = '--'; end
            if vok(s(q)), plot(tt, V(s(q)).y, ls, 'Color', col(q,:), 'LineWidth', 1);
            else,         plot(tt, V(s(q)).y, ':', 'Color', [0.6 0.6 0.6]); end
        end
        colormap(axT(k), cm);  clim(axT(k), zr);
        cb = colorbar(axT(k));  cb.Label.String = zl;
        ylabel('\DeltaI (pA)');
        title({sprintf('%s, %+.1f min', Q(k).name(1:6), Q(k).tRel), sprintf('light %s', Tq.Light{k}), ...
            sprintf('%d/%d kept, R_{step} %.0f M\\Omega', sum(vok(s)), numel(s), Rstep(k))}, 'FontSize', 9);
        axC(k) = nexttile(3*nQ + k);  hold on;
        for q = 1:numel(s), plot(tt, V(s(q)).light, 'Color', col(q,:)); end
        ylabel('Light (V)');  xlabel('Time from light (ms)');
    end
    linkaxes([axT axC], 'x');  xlim(axT(1), win);
    Yk = vertcat(V(vok).y);
    set(axT, 'YLim', [min(Yk(:)) max(Yk(:))] + [-0.1 0.1]*(max(Yk(:)) - min(Yk(:))));
    set(axC, 'YLim', [-0.5 1.2*max([vl 1])]);

    % Figure 7: I-V at constant light (with the probes at +20 mV, same measure, for reference) and,
    % for the light series, IPSC vs light intensity at each holding step.
    dockFig(7);
    tl = tiledlayout(1, 2 + any(byLight), 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, sprintf('Cell %d: IPSC (mean %d-%d ms) vs holding step and light', cellID, ivWin));
    cq = recCols(nQ);
    nexttile([1 2]);  hold on;
    h = gobjects(0);  lg = {};
    for k = find(fitIV')
        s = vrec == k & vok;
        h(end+1) = plot(vh(s), va(s), 'o', 'Color', cq(k,:), 'MarkerFaceColor', cq(k,:));   %#ok<SAGROW>
        xf = [min(vh(s)) max(vh(s))];
        plot(xf, gS(k)*xf + mean(va(s)) - gS(k)*mean(vh(s)), '-', 'Color', cq(k,:));
        if isnan(Erev(k)), es = 'n.d.'; else, es = sprintf('%+.0f mV', Erev(k)); end
        lg{end+1} = sprintf('%s (%+.1f min): %.2f nS, p = %.2g, E_{rev} %s', Q(k).name(1:6), Q(k).tRel, ...
            gS(k), pS(k), es);                                                              %#ok<SAGROW>
    end
    for j = 1:2   % kept probes, same measure
        sw = isPost == (j == 2);
        if ~any(sw), continue; end
        a = [W(sw).ampMean];
        h(end+1) = errorbar(mean([W(sw).hold]) + (j - 1.5)*0.6, mean(a), std(a)/sqrt(numel(a)), 'd', ...
            'Color', phCol{j}, 'MarkerFaceColor', phCol{j}, 'MarkerSize', 7, 'LineWidth', 1);      %#ok<SAGROW>
        lg{end+1} = sprintf('%s-PP probes (n = %d), mean \\pm SEM', phs{j}, numel(a));            %#ok<SAGROW>
    end
    xline(0, ':');  yline(0, ':');
    legend(h, lg, 'Location', 'southoutside');
    xlabel('\DeltaV_{hold} (mV, from the resting holding potential)');  ylabel('IPSC (pA)');
    title('I-V at constant light');
    for j = 1:2
        a = [W(isPost == (j == 2)).ampMean];
        fprintf('%-4s-PP probes at +20 mV, mean %d-%d ms: %.1f +- %.1f pA (mean +- SEM, n = %d)\n', ...
            phs{j}, ivWin, mean(a), std(a)/sqrt(numel(a)), numel(a));
    end
    if any(byLight)
        nexttile;  hold on;
        hl = gobjects(0);  ll = {};
        for k = find(byLight')
            for hv = unique(vh(vrec == k))
                s = vrec == k & vok & vh == hv;
                c = interp1([0 1], cmHold([1 end],:), hv/max([vh 1]));
                hl(end+1) = plot(vl(s), va(s), '-o', 'Color', c, 'MarkerFaceColor', c);    %#ok<SAGROW>
                ll{end+1} = sprintf('%s, %+.0f mV', Q(k).name(1:6), hv);                    %#ok<SAGROW>
            end
        end
        yline(0, ':');
        legend(hl, ll, 'Location', 'southoutside');
        xlabel('Light (V)');  ylabel('IPSC (pA)');  title({'IPSC vs light intensity', '(omitted pulses not shown)'});
    end
end

%% Local functions
function dockFig(k)
% Figure k, cleared, docked and in the light theme
figure(k);  clf;
set(gcf, 'WindowStyle', 'docked');
theme(gcf, 'light');
end

function x = getChan(bt, type, port)
% Data of the first channel of a task type whose physical channel contains port ([] if absent)
x = [];
for k = find(strcmp({bt.task_type}, type))
    for c = 1:numel(bt(k).channels)
        if contains(bt(k).channels(c).phys_channel, port, 'IgnoreCase', true)
            x = double(bt(k).channels(c).data(:));  return
        end
    end
end
end

function cols = recCols(n)
% n distinct colors, one per recording (hue wheel at a saturation / value readable on white)
cols = hsv2rgb([linspace(0, 0.8, n)' 0.8*ones(n,1) 0.75*ones(n,1)]);
end

function [W, tt] = getLightSweeps(Rs, win, blWin, ipscWin, ivWin, acProm)
% One sweep per light pulse: baseline-subtracted current (row), light command, holding step (mV),
% IPSC peak (in ipscWin) and mean (in ivWin), and the action currents (inward peaks with prominence
% > acProm) as times from the light onset (ms), searched over the holding step that contains the
% pulse, or over win for a pulse at the resting holding potential. clean = no action current in win.
W = struct('rec', {}, 'y', {}, 'light', {}, 'lightMax', {}, 'lightV', {}, 'hold', {}, 'bl', {}, ...
    'tLight', {}, 'acT', {}, 'clean', {}, 'amp', {}, 'ampMean', {});
tt = [];
for k = 1:numel(Rs)
    r  = Rs(k);
    w  = round(win(1)/1e3*r.fs):round(win(2)/1e3*r.fs);
    tt = w/r.fs*1e3;
    step = r.cmd > 1;                                   % mV, depolarizing holding step
    sOn  = find(diff([0; step]) == 1);  sOff = find(diff([step; 0]) == -1);
    for p = 1:numel(r.lightOn)
        i0 = round(r.lightOn(p)*r.fs) + 1;
        if i0 + w(1) < 1 || i0 + w(end) > numel(r.y), continue; end
        y  = double(r.y(i0 + w))';
        bl = median(y(tt >= blWin(1) & tt <= blWin(2)));
        q  = find(sOn <= i0 & sOff >= i0, 1);
        if isempty(q), a = i0 + w(1);  b = i0 + w(end);                        % at the resting holding
        else,          a = sOn(q) + round(5e-3*r.fs);  b = sOff(q) - round(2e-3*r.fs);   % skip the capacitive transients
        end
        [~, iAC] = findpeaks(-double(r.y(a:b)), 'MinPeakProminence', acProm, 'MinPeakDistance', round(2e-3*r.fs));
        acT = (a + iAC - 1 - i0)'/r.fs*1e3;
        W(end+1) = struct('rec', k, 'y', y - bl, 'light', double(r.light(i0 + w))', ...
            'lightMax', max(r.light(i0 + w)), 'lightV', double(r.light(i0 + 10)), 'hold', double(r.cmd(i0)), ...
            'bl', bl, 'tLight', r.tday + r.lightOn(p), 'acT', acT, 'clean', ~any(acT >= win(1) & acT <= win(2)), ...
            'amp', max(y(tt >= ipscWin(1) & tt <= ipscWin(2))) - bl, ...
            'ampMean', mean(y(tt >= ivWin(1) & tt <= ivWin(2))) - bl); %#ok<AGROW>
    end
end
end

function [g, b, p] = linFit(x, y)
% Least-squares line y = g*x + b, with the two-sided p-value of g = 0 (t-test, n - 2 dof)
x = x(:);  y = y(:);  n = numel(x);
c = [x ones(n,1)] \ y;  g = c(1);  b = c(2);
df = n - 2;  p = NaN;
if df > 0
    se = sqrt(sum((y - [x ones(n,1)]*c).^2)/df / sum((x - mean(x)).^2));
    t  = g/se;  p = betainc(df/(df + t^2), df/2, 0.5);
end
end

function plotFull(axT, axC, Rs, cols, xl, ds, gap, unit, cmdLab, colLight)
% Traces stacked gap apart (earliest on top) with a scale bar, and below them the commands:
% amplifier command (left axis, recording colors) and light (right axis, V).
hold(axT, 'on');
lab = cell(1, numel(Rs));
for k = 1:numel(Rs)
    r = Rs(k);  idx = 1:ds:numel(r.y);  t = (idx-1)/r.fs;
    plot(axT, t, r.y(idx) - (k-1)*gap, 'Color', cols(k,:), 'LineWidth', 0.5);
    lab{k} = sprintf('\\color[rgb]{%.2f,%.2f,%.2f}%+.1f min  %s', 0.75*cols(k,:), r.tRel, r.label);
end
% each trace is labelled at its baseline on the y axis: time from the first BTSP and protocol
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
yyaxis(axC, 'left');
for k = 1:numel(Rs)
    idx = 1:ds:numel(Rs(k).cmd);
    plot(axC, (idx-1)/Rs(k).fs, Rs(k).cmd(idx), '-', 'Color', cols(k,:), 'LineWidth', 1);
end
ylabel(axC, cmdLab);
yyaxis(axC, 'right');
for k = 1:numel(Rs)
    idx = 1:ds:numel(Rs(k).cmd);  t = (idx-1)/Rs(k).fs;
    if any(Rs(k).light), plot(axC, t, Rs(k).light(idx), '-', 'Color', colLight); end
end
ylabel(axC, 'Light (V)');  ylim(axC, [-0.2 1.2*max([arrayfun(@(r) max(r.light), Rs) 1])]);
axC.YAxis(1).Color = axC.XColor;  axC.YAxis(2).Color = axC.XColor;
xlabel(axC, 'Time in recording (s)');
linkaxes([axT axC], 'x');  xlim(axT, xl);
end
