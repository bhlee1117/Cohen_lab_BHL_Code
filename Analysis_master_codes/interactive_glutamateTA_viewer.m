function figG = interactive_glutamateTA_viewer(GluResult, VoltResult, VoltTrace_Norm, ...
                                               GluSpike, AlignedVolt_ftprnt, VoltRefIm, ...
                                               SkelDend, Device_Data, tformReg, varargin)
% INTERACTIVE_GLUTAMATETA_VIEWER  Interactively explore glutamate-triggered averages.
%
%   figG = interactive_glutamateTA_viewer(GluResult, VoltResult, VoltTrace_Norm, ...
%             GluSpike, AlignedVolt_ftprnt, VoltRefIm, SkelDend, ...
%             Device_Data, tformReg, Name, Value, ...)
%
%   Layout: the large merged reference image (top-left; glutamate reference in
%   red, voltage reference in light blue) is clickable and marks the glutamate
%   ROIs with 'o' markers colored by their number of events. The selected synapse
%   is circled in red and its nearest neighbours in the same colors as their
%   overlaid traces. Panels:
%     top-right  : Glu-triggered average voltage kymograph (voltage ROIs sorted
%                  by geodesic distance), with vertical guides at the glutamate
%                  frame interval (the temporal tolerance of the event timing)
%     bottom row : voltage footprints | Glu-triggered average dF/F trace of the
%                  selected ROI (+ NumNeighbors nearest ROIs) | glutamate trial x
%                  peri-spike matrix | voltage trial x peri-spike matrix
%
%   Because glutamate is imaged at a lower frame rate than voltage, each spike's
%   peak time is refined to sub-frame resolution with a local spline fit; that
%   refined time is used as the trigger for both the glutamate and voltage
%   averages, and all traces are resampled (interp1) onto a common time grid
%   built directly from GluResult.t_ax / VoltResult.t_ax.
%
%   REQUIRED INPUTS
%     GluResult          struct with fields .dFF_glu (nGluROI x nGluFrame),
%                        .t_ax (1 x nGluFrame), .S_glu (H x W x nGluROI),
%                        .AvgGluImg (H x W)
%     VoltResult         struct with field .t_ax (1 x nVoltFrame)
%     VoltTrace_Norm     normalized voltage traces (nVoltROI x nVoltFrame)
%     GluSpike           logical/double glutamate spike raster (nGluROI x nGluFrame)
%     AlignedVolt_ftprnt voltage footprints TRANSFORMED into glu-camera space
%                        (H x W x nVoltROI). Voltage ROIs must be co-registered
%                        with the glutamate data for the distance sorting to be
%                        meaningful (e.g. transformCamera_O2B(...VoltResult.ftprnt...)).
%     VoltRefIm          voltage reference image in glu-camera space (H x W),
%                        used only as the background for the footprint outlines
%     SkelDend           dendrite skeleton / cell mask in VOLTAGE-camera space
%                        (same space as VoltResult.ftprnt / VoltResult.ref_im),
%                        used for geodesic distance (e.g. DendriteMask_N1)
%     Device_Data        device metadata cell array (for camera ROI geometry)
%     tformReg           Volt<->Glu registration transform (affine2d)
%
%   Geodesic distances are computed in VOLTAGE space: the native voltage
%   footprints (VoltResult.ftprnt) give the ROI positions, and the selected
%   glutamate synapse is mapped from glu space into voltage space via
%   Device_Data / tformReg, so SkelDend, the ROI coordinates, and the synapse
%   all live in the same frame as DendriteMask_N1.
%
%   NAME-VALUE OPTIONS
%     'PixelSize'     micrometers per pixel                    (default 0.4680)
%     'VoltLag'       voltage peri-spike half-window in ms     (default 100)
%     'GluWin_ms'     glutamate peri-spike half-window in ms   (default 200)
%     'CLimVolt'      color limits for the voltage STA panel   (default [-0.005 0.01])
%     'InitROI'       glutamate ROI index to display on open   (default 1)
%     'PeakWinFrames' half-window (glu frames) searched for the sub-frame
%                     spline peak around each detected spike   (default 2)
%     'SplineUp'      spline upsampling factor for the peak search (default 20)
%     'NumNeighbors'  # of nearest glutamate ROIs overlaid on panel 1 (default 5)
%
%   OUTPUT
%     figG         handle to the created figure
%
%   Example:
%     interactive_glutamateTA_viewer(GluResult, VoltResult, VoltTrace_Norm, ...
%         GluSpike, AlignedVolt_ftprnt, VoltRefIm, DendriteMask_N1, ...
%         Device_Data, tformReg, 'PixelSize', 0.468, 'InitROI', 7);

% ---------------------- Parse options ----------------------
p = inputParser;
addParameter(p, 'PixelSize', 0.4680, @(x) isscalar(x) && x > 0);
addParameter(p, 'VoltLag',   100,    @(x) isscalar(x) && x > 0);
addParameter(p, 'GluWin_ms', 200,    @(x) isscalar(x) && x > 0);
addParameter(p, 'CLimVolt',  [-0.005 0.01], @(x) numel(x) == 2);
addParameter(p, 'InitROI',   1,      @(x) isscalar(x) && x >= 1);
addParameter(p, 'PeakWinFrames', 2,  @(x) isscalar(x) && x >= 1);
addParameter(p, 'SplineUp', 20,      @(x) isscalar(x) && x >= 1);
addParameter(p, 'NumNeighbors', 5,   @(x) isscalar(x) && x >= 0);
parse(p, varargin{:});
opt = p.Results;

% ---------------------- Assemble state ----------------------
% Trim traces and their time axes to a common length (t_ax and the trace
% columns can differ by a frame), and read the native glutamate frame period
% straight from GluResult.t_ax.
Ng = min(numel(GluResult.t_ax), size(GluResult.dFF_glu, 2));
Nv = min(numel(VoltResult.t_ax), size(VoltTrace_Norm, 2));
dt_glu_ms = median(diff(GluResult.t_ax(1:Ng))) * 1000;   % native glu frame period (ms)

IP = struct();
IP.VoltLag   = opt.VoltLag;                             % voltage half-window (ms, 1 ms/frame)
IP.GluWin_ms = opt.GluWin_ms;                           % glutamate half-window (ms)
IP.PeakWin   = opt.PeakWinFrames;                       % half-window (glu frames) for peak spline
IP.SplineUp  = opt.SplineUp;                            % spline upsampling factor for peak search
IP.CLimVolt  = opt.CLimVolt;
IP.pixelsize = opt.PixelSize;

% Common relative time axes (ms), referenced to the refined glutamate peak time.
IP.tax_volt   = -IP.VoltLag : IP.VoltLag;                        % 1 ms grid (voltage)
IP.GluStep_ms = max(dt_glu_ms/2, eps);                           % ~2x upsample of native glu rate
IP.tax_glu    = -IP.GluWin_ms : IP.GluStep_ms : IP.GluWin_ms;    % ms grid (glutamate)
IP.dt_glu_ms  = dt_glu_ms;                                       % native glu frame interval (ms)

IP.dFF_glu        = GluResult.dFF_glu(:, 1:Ng);
IP.GluResult_t_ax = GluResult.t_ax(1:Ng);
IP.Volt_t_ax      = VoltResult.t_ax(1:Nv);
IP.VoltTrace_Norm = VoltTrace_Norm(:, 1:Nv);
IP.GluSpike       = GluSpike(:, 1:min(size(GluSpike,2), Ng));
IP.Glu_coord      = get_coord(GluResult.S_glu);        % nGluROI x 2 (glu space, for display/click)
IP.VoltCoord      = get_coord(AlignedVolt_ftprnt);     % nVoltROI x 2 (glu space, for display)
IP.nVoltROI       = size(AlignedVolt_ftprnt, 3);
IP.nGluROI        = size(GluResult.S_glu, 3);
IP.cmapRank       = parula(IP.nVoltROI);               % color per distance rank (near->far)
IP.distCache      = nan(IP.nVoltROI, IP.nGluROI);      % geodesic dist cache (volt x glu)

% ---- Geodesic distance is computed in VOLTAGE space (same frame as SkelDend) ----
IP.SkelDend    = SkelDend;                             % voltage-space dendrite mask
IP.VoltCoord_v = get_coord(VoltResult.ftprnt);         % nVoltROI x 2 (voltage space)
% Referencing objects to map the selected glu synapse -> voltage space, using the
% same pipeline as the main analysis (intrinsicToWorld/transformPointsForward).
IP.tformReg = tformReg;
IP.Rvolt    = refFromROI(size(VoltResult.ref_im),   double(Device_Data{1,3}.ROI));
IP.Rglu     = refFromROI(size(GluResult.AvgGluImg), double(Device_Data{1,4}.ROI));

% Glutamate ROI coordinates in voltage space, for finding the nearest neighbours
% of a selected synapse (geodesic distance along SkelDend, cached per synapse).
IP.GluCoord_v = zeros(IP.nGluROI, 2);
for gi = 1:IP.nGluROI
    IP.GluCoord_v(gi,:) = glu2volt(IP, IP.Glu_coord(gi,:));
end
IP.gluDistCache = nan(IP.nGluROI, IP.nGluROI);   % glu-glu geodesic dist cache
IP.nNeighbors   = opt.NumNeighbors;              % # of neighbour ROIs to overlay

% ---------------------- Build figure ----------------------
% Grid: 3 rows x 12 cols. Row 1-2 hold the big merged reference image (left 2/3)
% and the Glu-triggered average voltage kymo (right 1/3); row 3 holds the voltage
% footprint, the Glu-triggered average trace, and the two trial peri-spike matrices.
figG = figure('Name','Glu-triggered explorer','Color','w', ...
              'Units','normalized','Position',[0.05 0.1 0.9 0.8]);
tlG = tiledlayout(figG, 3, 12, 'TileSpacing','compact','Padding','compact');
set(figG,'theme','light');
% ---- Merged reference image (top-left, large, clickable) ----
% Glutamate reference in RED + voltage reference in LIGHT BLUE. Glutamate ROI
% locations are marked with 'o' markers colored by their number of events. The
% selected synapse (red circle) and its nearest neighbours (colored like their
% STA traces) are added on each click.
IP.lightblue = [0.40 0.70 1.00];
nrm   = @(im) mat2gray(double(im), double(prctile(double(im(:)), [1 99])));
gluN  = nrm(GluResult.AvgGluImg);
voltN = nrm(VoltRefIm);
mergedRef = cat(3, gluN*1 + voltN*IP.lightblue(1), ...
                   gluN*0 + voltN*IP.lightblue(2), ...
                   gluN*0 + voltN*IP.lightblue(3));
mergedRef = min(mergedRef, 1);

IP.axGlu = nexttile(tlG, [2 8]);
imGlu = imagesc(IP.axGlu, mergedRef); axis(IP.axGlu,'image','off');
hold(IP.axGlu,'on');
nGluEvents = sum(IP.GluSpike > 0, 2);                    % # of events per glutamate ROI
gluCols    = vec2cmap(nGluEvents, turbo(256));           % color each ROI by event count
scatter(IP.axGlu, IP.Glu_coord(:,1), IP.Glu_coord(:,2), 40, gluCols, 'o', 'filled', ...
        'HitTest','off', 'PickableParts','none');
title(IP.axGlu, 'Merged reference (glu = red, volt = light blue); ROI color = # events');

% ---- Glutamate-triggered average voltage kymograph (top-right) ----
IP.ax2 = nexttile(tlG, [2 4]);   % voltage STA, ROIs sorted by distance

% ---- Voltage footprint image (bottom row) ----
% Outline of every (transformed) voltage footprint over the reference image;
% boundaries recolored on each click by distance rank from the selected synapse.
IP.axVolt = nexttile(tlG, [1 3]);
imagesc(IP.axVolt, VoltRefIm); axis(IP.axVolt,'image','off');
colormap(IP.axVolt, gray); hold(IP.axVolt,'on');
IP.hVoltBounds = gobjects(IP.nVoltROI, 1);
for i = 1:IP.nVoltROI
    B = bwboundaries(AlignedVolt_ftprnt(:,:,i) > 0);
    if isempty(B)
        IP.hVoltBounds(i) = plot(IP.axVolt, nan, nan);   % placeholder keeps handle array aligned
        continue;
    end
    [~, bi] = max(cellfun(@(b) size(b,1), B));           % largest boundary for this ROI
    IP.hVoltBounds(i) = plot(IP.axVolt, B{bi}(:,2), B{bi}(:,1), '-', ...
        'LineWidth', 1.5, 'Color', [0.5 0.5 0.5], ...
        'HitTest','off', 'PickableParts','none');
end
title(IP.axVolt, 'Voltage footprints (color = distance rank)');

% ---- Bottom row: Glu-triggered average trace + trial peri-spike matrices ----
IP.ax1 = nexttile(tlG, [1 3]);   % glu-triggered average trace (+ neighbours)
IP.ax3 = nexttile(tlG, [1 3]);   % glutamate trial x peri-spike matrix
IP.ax4 = nexttile(tlG, [1 3]);   % voltage  trial x peri-spike matrix

IP.hSelGlu  = gobjects(1);       % selected-synapse marker (glu image)
IP.hNbGlu   = gobjects(0);       % neighbour markers (glu image)
IP.hSelVolt = gobjects(1);       % nearest-volt-ROI marker (volt image)

setappdata(figG, 'GTA', IP);
set([imGlu, IP.axGlu], 'ButtonDownFcn', @GTA_onClickGlu);

% Draw an initial synapse
if opt.InitROI <= IP.nGluROI
    GTA_update(figG, opt.InitROI);
end
fprintf('Interactive explorer ready: click a glutamate synapse to update panels.\n');
end

% ======================= CALLBACKS (local functions) =======================
function GTA_onClickGlu(src, ~)
% Click handler: find the glutamate ROI nearest the clicked point and refresh.
    fig = ancestor(src, 'figure');
    ax  = ancestor(src, 'axes');
    cp  = get(ax, 'CurrentPoint');
    xy  = cp(1, 1:2);
    IP  = getappdata(fig, 'GTA');
    [~, si] = min(sum((IP.Glu_coord - xy).^2, 2));
    GTA_update(fig, si);
end

function v = glu2volt(IP, g)
% Map a point g = [x y] from glutamate-camera intrinsic coordinates into
% voltage-camera intrinsic coordinates (same pipeline as the main analysis).
    [xw, yw] = intrinsicToWorld(IP.Rglu, g(1), g(2));
    [xw, yw] = transformPointsForward(IP.tformReg, xw, yw);
    [vx, vy] = worldToIntrinsic(IP.Rvolt, xw, yw);
    v = [vx, vy];
end

function GTA_update(fig, si)
% Recompute and redraw all four panels for the selected glutamate ROI si.
    IP = getappdata(fig, 'GTA');

    gluSpkFrames = find(IP.GluSpike(si, :) > 0);   % glutamate spike frames (selected ROI)
    if isempty(gluSpkFrames)
        cla(IP.ax1); cla(IP.ax2); cla(IP.ax3); cla(IP.ax4);
        title(IP.ax1, sprintf('Synapse %d: no glutamate spikes detected', si));
        return;
    end

    % ---------- Sub-frame glutamate peak time (spline) ----------
    % Glutamate is sampled slower than voltage, so refine each spike's peak time
    % to sub-frame resolution: upsample a short window around the detected spike
    % with a spline and take the location of its maximum as the trigger time.
    Ng   = numel(IP.GluResult_t_ax);
    spkF = gluSpkFrames(gluSpkFrames >= 1 & gluSpkFrames <= Ng);
    tpk  = nan(numel(spkF), 1);                            % refined peak time (s) per trial
    for k = 1:numel(spkF)
        wf = max(1, spkF(k)-IP.PeakWin) : min(Ng, spkF(k)+IP.PeakWin);
        tl = IP.GluResult_t_ax(wf);
        yl = IP.dFF_glu(si, wf);
        if numel(tl) < 3
            tpk(k) = IP.GluResult_t_ax(spkF(k));          % window too short: keep frame time
            continue;
        end
        tq = linspace(tl(1), tl(end), (numel(tl)-1)*IP.SplineUp + 1);
        yq = interp1(tl, yl, tq, 'spline');
        [~, im] = max(yq);
        tpk(k) = tq(im);
    end
    nTrial = numel(tpk);

    % ---------- Glutamate peri-spike matrix + STA (selected ROI) ----------
    % Resample the glutamate trace on the common ms grid, referenced to each
    % refined peak time (spline interpolation, sub-frame aligned).
    gluPeri = nan(nTrial, numel(IP.tax_glu));
    for k = 1:nTrial
        tq = tpk(k) + IP.tax_glu/1000;                    % query times (s)
        gluPeri(k,:) = interp1(IP.GluResult_t_ax, IP.dFF_glu(si,:), tq, 'spline', NaN);
    end
    gluPeri  = gluPeri - median(gluPeri, 2, 'omitnan');
    gluTrace = mean(gluPeri, 1, 'omitnan');

    % ---------- Voltage peri-spike matrices + STA (all ROIs) ----------
    % Same refined trigger times, resampled onto the 1 ms voltage grid.
    Mat_volt = nan(IP.nVoltROI, nTrial, numel(IP.tax_volt));
    for k = 1:nTrial
        tq = tpk(k) + IP.tax_volt/1000;                   % query times (s)
        vq = interp1(IP.Volt_t_ax, IP.VoltTrace_Norm.', tq, 'linear', NaN);  % time x ROI
        Mat_volt(:, k, :) = vq.';
    end
    Mat_volt = Mat_volt - median(Mat_volt, 3, 'omitnan');  % baseline per ROI/trial
    TA_volt  = squeeze(mean(Mat_volt, 2, 'omitnan'));      % voltROI x time

    % ---------- Distance of each voltage ROI from selected synapse ----------
    % Computed in VOLTAGE space: geodesic distance along SkelDend (DendriteMask_N1)
    % from each native voltage footprint to the selected synapse, after mapping the
    % synapse from glu space into voltage space. Cached per synapse.
    if any(isnan(IP.distCache(:, si)))
        synV = glu2volt(IP, IP.Glu_coord(si, :));      % selected synapse in voltage space
        d = nan(IP.nVoltROI, 1);
        for i = 1:IP.nVoltROI
            try
                d(i) = geodesic_distance(IP.SkelDend, IP.VoltCoord_v(i,:), synV);
            catch
                d(i) = norm(IP.VoltCoord_v(i,:) - synV);   % fallback: euclidean
            end
        end
        IP.distCache(:, si) = d;
        setappdata(fig, 'GTA', IP);
    end
    [distSorted, order] = sort(IP.distCache(:, si), 'ascend');
    nearestVolt = order(1);

    % Distance rank of each voltage ROI (1 = nearest to the selected synapse).
    % Recolor the footprint boundaries so their color matches the row order of
    % the voltage STA panel below.
    rankOfROI = zeros(IP.nVoltROI, 1);
    rankOfROI(order) = 1:IP.nVoltROI;
    for i = 1:IP.nVoltROI
        if isgraphics(IP.hVoltBounds(i))
            set(IP.hVoltBounds(i), 'Color', IP.cmapRank(rankOfROI(i), :));
        end
    end

    voltPeri = squeeze(Mat_volt(nearestVolt, :, :));       % trials x time (already baselined)

    % ---------- Nearest neighbour glutamate ROIs ----------
    % Pre-filter by euclidean distance (fast), then rank the candidates by
    % geodesic distance along SkelDend (voltage space). Same trigger times as
    % the selected ROI, so the overlaid traces show how nearby synapses
    % co-activate when the selected synapse fires.
    eucl = sqrt(sum((IP.GluCoord_v - IP.GluCoord_v(si,:)).^2, 2));
    [~, euOrder] = sort(eucl, 'ascend');
    cand = euOrder(2 : min(1 + 3*IP.nNeighbors, IP.nGluROI));   % skip self
    need = cand(isnan(IP.gluDistCache(cand, si)));
    for j = need(:)'
        try
            IP.gluDistCache(j, si) = geodesic_distance(IP.SkelDend, IP.GluCoord_v(j,:), IP.GluCoord_v(si,:));
        catch
            IP.gluDistCache(j, si) = norm(IP.GluCoord_v(j,:) - IP.GluCoord_v(si,:));
        end
    end
    setappdata(fig, 'GTA', IP);
    [gdSorted, gord] = sort(IP.gluDistCache(cand, si), 'ascend');
    nbN    = min(IP.nNeighbors, numel(cand));
    nbROIs = cand(gord(1:nbN));
    nbDist = gdSorted(1:nbN);

    % STA trace of each neighbour ROI, triggered by the SELECTED ROI's peaks.
    nbTrace = nan(nbN, numel(IP.tax_glu));
    for m = 1:nbN
        pm = nan(nTrial, numel(IP.tax_glu));
        for k = 1:nTrial
            tq = tpk(k) + IP.tax_glu/1000;
            pm(k,:) = interp1(IP.GluResult_t_ax, IP.dFF_glu(nbROIs(m),:), tq, 'spline', NaN);
        end
        pm = pm - median(pm, 2, 'omitnan');
        nbTrace(m,:) = mean(pm, 1, 'omitnan');
    end

    % ---------- (1) Glu-triggered average trace (+ nearest neighbours) ----------
    cla(IP.ax1); hold(IP.ax1, 'on');
    cmapNb = cool(max(nbN,1));
    hNb = gobjects(nbN, 1);
    for m = 1:nbN
        hNb(m) = plot(IP.ax1, IP.tax_glu, nbTrace(m,:), 'Color', cmapNb(m,:), 'LineWidth', 1);
    end
    hSel = plot(IP.ax1, IP.tax_glu, gluTrace, 'r', 'LineWidth', 2);   % selected (matches red circle)
    hold(IP.ax1, 'off');
    xlabel(IP.ax1,'Time from glu spike (ms)'); ylabel(IP.ax1,'\DeltaF/F');
    title(IP.ax1, sprintf('Glu-triggered average (spline-aligned) — synapse %d (n = %d)', ...
          si, nTrial));
    xlim(IP.ax1, IP.tax_glu([1 end])); box(IP.ax1,'off');
    legStr = [{sprintf('Selected ROI %d', si)}, ...
              arrayfun(@(m) sprintf('ROI %d (%.0f \\mum)', nbROIs(m), nbDist(m)*IP.pixelsize), ...
                       1:nbN, 'UniformOutput', false)];
    legend([hSel; hNb], legStr, 'Location','best', 'FontSize',7, 'Box','off');

    % ---------- (2) Glu-triggered average voltage kymograph (sorted by distance) ----------
    cla(IP.ax2);
    imagesc(IP.ax2, IP.tax_volt, 1:IP.nVoltROI, TA_volt(order, :));
    colormap(IP.ax2, turbo); caxis(IP.ax2, IP.CLimVolt);
    axis(IP.ax2,'tight'); hold(IP.ax2,'on');
    % Vertical guides at the glutamate frame interval: the kymo is at 1 ms but the
    % glutamate event can only be localized to within one glu frame (dt_glu_ms).
    kmax = floor(IP.tax_volt(end) / IP.dt_glu_ms);
    for kk = -kmax:kmax
        xline(IP.ax2, kk*IP.dt_glu_ms, ':', 'Color', [1 1 1]*0.85, 'LineWidth', 0.5);
    end
    % Left color strip: each row is colored by the SAME distance-rank colormap
    % used for the footprint boundaries (row 1 = nearest ROI, at the top).
    xStrip = IP.tax_volt(1) - 0.04*(IP.tax_volt(end) - IP.tax_volt(1));
    scatter(IP.ax2, xStrip*ones(IP.nVoltROI,1), 1:IP.nVoltROI, 40, IP.cmapRank, ...
            'filled', 's', 'Clipping','off');
    xlim(IP.ax2, [xStrip IP.tax_volt(end)]);
    hold(IP.ax2,'off');
    xlabel(IP.ax2,'Time from glu spike (ms)');
    ylabel(IP.ax2,'Voltage ROI (near \rightarrow far)');
    title(IP.ax2, 'Glu-triggered average voltage kymograph (sorted by distance)');
    cb = colorbar(IP.ax2); cb.Label.String = '-\DeltaF/F';

    % ---------- (3) trial x peri-spike matrix of glutamate ----------
    cla(IP.ax3);
    imagesc(IP.ax3, IP.tax_glu, 1:size(gluPeri,1), gluPeri);
    colormap(IP.ax3, turbo); axis(IP.ax3,'tight');
    xlabel(IP.ax3,'Time from glu spike (ms)'); ylabel(IP.ax3,'Trial #');
    title(IP.ax3, sprintf('Glutamate peri-spike matrix (synapse %d)', si));

    % ---------- (4) trial x peri-spike matrix of voltage ----------
    cla(IP.ax4);
    imagesc(IP.ax4, IP.tax_volt, 1:size(voltPeri,1), voltPeri);
    colormap(IP.ax4, turbo); axis(IP.ax4,'tight');
    xlabel(IP.ax4,'Time from glu spike (ms)'); ylabel(IP.ax4,'Trial #');
    title(IP.ax4, sprintf('Voltage peri-spike matrix (ROI %d, %.1f µm away)', ...
          nearestVolt, distSorted(1)*IP.pixelsize));

    % ---------- Update footprint markers ----------
    if isgraphics(IP.hSelGlu),  delete(IP.hSelGlu);  end
    if isgraphics(IP.hSelVolt), delete(IP.hSelVolt); end
    if ~isempty(IP.hNbGlu),     delete(IP.hNbGlu(isgraphics(IP.hNbGlu))); end
    % Neighbours on the glutamate footprint, circled in their trace colors
    IP.hNbGlu = gobjects(nbN, 1);
    for m = 1:nbN
        IP.hNbGlu(m) = scatter(IP.axGlu, IP.Glu_coord(nbROIs(m),1), IP.Glu_coord(nbROIs(m),2), ...
                       70, cmapNb(m,:), 'LineWidth', 2, 'HitTest','off', 'PickableParts','none');
    end
    % Selected synapse: red circle
    IP.hSelGlu  = scatter(IP.axGlu, IP.Glu_coord(si,1), IP.Glu_coord(si,2), ...
                          120, 'r', 'LineWidth', 2.5, 'HitTest','off', 'PickableParts','none');
    % Nearest voltage ROI on the voltage footprint
    IP.hSelVolt = scatter(IP.axVolt, IP.VoltCoord(nearestVolt,1), IP.VoltCoord(nearestVolt,2), ...
                          90, 'g', 'LineWidth', 2, 'HitTest','off', 'PickableParts','none');
    setappdata(fig, 'GTA', IP);
end
