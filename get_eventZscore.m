function eventTable = get_eventZscore(v_dF, spikeVec, ftprint, shotVar)
%% Compute event-wise shot-noise-limited z-score / d' for detected voltage events
%
% Required inputs:
%   v_dF      : N x T voltage trace in dF
%   spikeVec  : N x T binary or event-count matrix
%   ftprint   : X x Y x N footprint variable
%   shotVar   : X x Y empirical shot-noise variance map
%
% Interpretation:
%   event z-score / d' = event amplitude in dF / shot-noise std of ROI trace
%
% Assumptions:
%   - Events are positive-going
%   - v_dF was extracted using footprint-weighted average
%   - shot noise is independent across pixels
%   - shotVar is in the same intensity units as the movie used to compute v_dF

%% Parameters

eventWin = 0;          % frames after detected event to search peak
baselineWin = 0;       % frames before event for local baseline subtraction

normalizeFootprint = false;  
% true  if v_dF(n,t) = sum(ftprint(:,:,n).*mov(:,:,t),'all') / sum(ftprint(:,:,n),'all')
% false if v_dF(n,t) = sum(ftprint(:,:,n).*mov(:,:,t),'all')

makePlots = true;

%% Basic checks

[N, T] = size(spikeVec);

assert(size(v_dF,1) == N && size(v_dF,2) == T, ...
    'v_dF must be N x T');

assert(size(ftprint,3) == N, ...
    'ftprint must be X x Y x N');

assert(all(size(shotVar) == size(ftprint(:,:,1))), ...
    'shotVar must be X x Y');

%% Compute ROI-level shot-noise std from footprints

roiShotVar = NaN(N,1);
roiShotStd = NaN(N,1);
roiWeightSum = NaN(N,1);
roiNPixels = NaN(N,1);

for n = 1:N
    w = ftprint(:,:,n);
    w(~isfinite(w)) = 0;
    w(w < 0) = 0;

    roiWeightSum(n) = sum(w(:));
    roiNPixels(n) = nnz(w > 0);

    if roiWeightSum(n) <= 0
        continue;
    end

    if normalizeFootprint
        wEff = w ./ roiWeightSum(n);
    else
        wEff = w;
    end

    roiShotVar(n) = sum((wEff(:).^2) .* shotVar(:), 'omitnan');
    roiShotStd(n) = sqrt(roiShotVar(n));
end

%% Extract detected events and compute event-wise z-score / d'

[eventROI, eventT] = find(spikeVec > 0);
nEvents = numel(eventROI);

eventAmp_dF = NaN(nEvents,1);
eventZ_dF = NaN(nEvents,1);
eventPeakT = NaN(nEvents,1);

for e = 1:nEvents
    n = eventROI(e);
    t0 = eventT(e);

    % Search window for positive event peak
    t1 = t0;
    t2 = min(T, t0 + eventWin);
    idx = t1:t2;

    % Local baseline window
    b1 = max(1, t0 - baselineWin);
    b2 = max(1, t0 - 1);
    bidx = b1:b2;

    if baselineWin > 0 && ~isempty(bidx)
        base_dF = mean(v_dF(n,bidx), 'omitnan');
    else
        base_dF = 0;
    end

    seg_dF = v_dF(n,idx) - base_dF;

    % Positive-going event amplitude
    [amp_dF, relIdx] = max(seg_dF);

    eventAmp_dF(e) = amp_dF;
    eventPeakT(e) = idx(relIdx);

    % Shot-noise-limited event z-score / d'
    eventZ_dF(e) = amp_dF ./ roiShotStd(n);
end

%% Store results in table

eventTable = table;
eventTable.roi = eventROI;
eventTable.t_event = eventT;
eventTable.t_peak = eventPeakT;
eventTable.amp_dF = eventAmp_dF;
eventTable.roiShotStd_dF = roiShotStd(eventROI);
eventTable.eventZ_dF = eventZ_dF;
eventTable.roiShotVar_dF = roiShotVar(eventROI);
eventTable.roiWeightSum = roiWeightSum(eventROI);
eventTable.roiNPixels = roiNPixels(eventROI);

%% Quick visualization

if makePlots
    figure;
    nexttile([1 1]);
    histogram(eventTable.eventZ_dF, 50);
    xlabel('Event shot-noise z-score / d''');
    ylabel('Number of events');
    title('Detected event shot-noise-limited z-score distribution');

    nexttile([1 1]);
    scatter(eventTable.amp_dF, eventTable.eventZ_dF, 15, 'filled');
    xlabel('Event amplitude \DeltaF');
    ylabel('Event shot-noise z-score / d''');
    title('Event amplitude vs shot-noise-limited z-score');
    grid on;
end

end