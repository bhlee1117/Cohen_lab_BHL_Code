function cfg = flow_defaults()
% cfg = flow_defaults()
%
% Every parameter of the vascular-flow pipeline, in one place, grouped by the stage that
% consumes it. The stages are:
%
%   flow_bandpass    spatial band-pass of the motion-residual movie      (was cell 7c)
%   flow_track       RBC detection and linking into tracklets            (was cell 7d)
%   flow_segments    vessel segments from the detection density          (was cell 7e)
%   flow_regressors  per-segment oscillator -> the design matrix         (was cell 7f)
%   flow_apply       fit it to the whole recording                       (was cell 7g)
%   flow_compare     what it does to the extracted traces                (was cell 7h)
%
% WHY A CONFIG STRUCT. The script version carried 131 loose workspace variables across
% six cells. Nothing recorded which stage owned which knob, changing one meant knowing by
% heart which cells had to be re-run, and a stale variable from an earlier run produced a
% wrong answer rather than an error. Here each stage reads only cfg.<stage>, and
% flow_pipeline hashes that sub-struct to decide what actually needs recomputing.
%
% Only the knobs live here. Anything MEASURED from the data (pxUm from Pixelsize, the
% achieved passband, the per-segment frequency) belongs in the result, not the config.
%
% 2026.08.28, Byung Hun Lee

% ===================== flow_bandpass  (7c) =========================================
c = struct();
c.voltRef     = 'ring'; % WHICH VOLTAGE ESTIMATE IS REMOVED BEFORE THE SPATIAL FILTER.
                        % Something has to be: the dendrite is a thin bright structure, so
                        % the spatial band-pass responds to its EDGES and only partly
                        % suppresses Vm. Removing it at source leaves something much closer
                        % to "blood only".
                        %
                        % 'ring'   from a ring of footprint just OUTSIDE the vessels, where
                        %          blood cannot be. See get_ringRef. This is the default
                        %          because the alternative below is contaminated at source.
                        % 'traces' data.voltTr, the per-ROI voltage traces. Those are read
                        %          off the WHOLE footprint, which OVERLAPS the vessels, so
                        %          the "voltage" being removed already contains blood:
                        %          regressing it out takes vascular signal away before
                        %          anything has measured it. Kept for comparison, and for
                        %          recordings with no usable footprint.
                        % 'none'   use mov_res unchanged and let the spatial filter do all
                        %          the work. Honest baseline; expect Vm at the edges.
% -- the 'ring' reference (ignored otherwise) -- see get_ringRef for the reasoning
c.ringIn      = 1;      % px, gap between the vessel mask and the inner edge of the ring
c.ringOut     = 6;      % px, outer edge of the ring
c.ringThr     = 0.2;    % fraction of an ROI's peak below which its footprint is absent
c.ringMinPx   = 50;     % px, smallest CONNECTED patch of overlap that becomes a trace.
                        % Each ROI's overlap with the ring is split into its connected
                        % components first: a vessel is usually longer than the stretch of
                        % dendrite carrying one ROI, so the ring clips that footprint in
                        % several disconnected places, and those sit at different points
                        % on the neurite. Averaging them gives a reference matching none.
c.ringSmoothMs = 5;     % ms, moving average on each reference trace. Keeps the part that
                        % is voltage and drops shot noise. It does NOT confine the trace
                        % to the sub-fMin band -- it does not have to, because the ring
                        % contains no vessel pixel and so no RBC train to pass.
c.ringDedupR  = 0.95;   % |corr| above which two reference traces count as the same trace
c.ringMax     = 40;     % cap on reference traces: each one is a degree of freedom taken
                        % out of the movie, and past a few dozen the removal starts to be
                        % measurable on the dendrite for its own sake.
c.loUm        = 6;      % um, fine end of the spatial passband
c.hiUm        = 30;     % um, coarse end
c.method      = 'dog';  % 'dog'     Difference of Gaussians: no ringing, but it cannot
                        %           make a band narrower than ~3.4:1 whatever the sigmas
                        % 'fourier' radial Butterworth: hits the requested band, at the
                        %           cost of a -0.24 side lobe
c.order       = 4;      % Butterworth order, 'fourier' only
c.chunk       = 200;    % frames filtered per pass. PURELY a memory/speed knob -- the
                        % answer is identical whatever it is set to.
cfg.bandpass = c;

% ===================== flow_track  (7d) ============================================
c = struct();
c.frames     = [];      % frames to track. [] = 1:min(2000, nT). Tracking is the slow
                        % stage and the flow is stationary over seconds, so a window is
                        % enough to MEASURE the geometry; it is applied to everything.
c.thrSD      = 4;       % detection threshold, in SD of S over the tracked frames
c.featPx     = [];      % odd feature size for pkfnd. [] = derived from the passband
c.maxDisp    = 6;       % px, largest displacement linked between consecutive frames
c.mem        = 0;       % frames a particle may vanish for and still be linked. 0 is
                        % deliberate: a longer memory almost always links an RBC to the
                        % plasma gap beside it rather than to itself.
c.minLen     = 4;       % frames, shortest tracklet track() will return. This is the
                        % LINKING minimum, not a quality cut -- the quality cuts live in
                        % cfg.segments, because deciding which tracklets are trustworthy
                        % is part of defining the segments, not part of finding them.
c.polarity   = 'both';  % 'bright' | 'dark' | 'both'. RBCs are dark on plasma but the
                        % filter inverts contrast at some depths, so 'both' runs the
                        % detector twice and labels which population each tracklet is.
c.binPx      = 12;      % px per bin of the quiver flow map
cfg.track = c;

% ===================== flow_segments  (7e) =========================================
c = struct();
% -- step 1 of this stage: SORT THE TRACKLETS. Only the ones that survive build the
%    density map, so these thresholds decide the segmentation more than anything below.
c.minNetUm     = 3;     % um. NET DISPLACEMENT: how far the cell actually GOT, from the
                        % first point of the tracklet to the last. The most direct
                        % statement of "this thing travelled", and the cut to reach for
                        % first: it is in physical units, so it can be set from what a
                        % red cell is expected to do rather than from a distribution.
                        % A tracklet that wandered in place fails it however long it is.
c.minLenKeep   = 6;     % frames. Minimum DURATION of one tracklet -- how many frames it
                        % was followed for. Not to be confused with minTrk below, which
                        % counts tracklets per segment: this cut judges an individual
                        % tracklet, that one judges a vessel patch.
c.minAmpPrc    = 20;    % drop tracklets dimmer than this percentile of |amplitude|
c.minPersist   = 0.6;   % mean cos(angle between consecutive steps). 0 = random walk
                        % REGARDLESS of length, 1 = straight travel. Length-fair (unlike
                        % straightness, whose chance level rises as tracklets shorten),
                        % but secondary to minNetUm: a cell can travel a long way while
                        % jittering, and distance is the thing that matters.
% -- step 2: the density map and the patches cut out of it
c.densSm       = 3;     % px, gaussian smoothing of the detection-density map
c.densPrc      = 70;    % percentile of the non-empty density map that becomes the mask
c.minPx        = 10;    % minimum pixels for a connected patch to be a segment
c.minTrk       = 8;     % how many tracklets must fall inside a connected patch before
                        % it is accepted as a segment. Counts TRACKLETS PER PATCH --
                        % minLenKeep above counts FRAMES PER TRACKLET. A patch crossed
                        % once by one long cell is not a vessel you can measure a
                        % frequency from; this is what rejects it.
cfg.segments = c;

% ===================== flow_regressors  (7f) =======================================
c = struct();
% -- which segments qualify --
c.minPxPCA    = 20;     % px. Segments at least this big get PCA (they can hold a phase
                        % gradient); smaller ones use their own trace, because their
                        % leading PC would just be the mean with extra noise.
c.binPix      = 12;     % px per mode SCANNED: a segment is searched with
                        % ceil(nPixels/binPix) modes, so a bigger segment gets more.
% -- the decomposition --
c.cpca        = true;   % COMPLEX PCA (band-limit, analytic signal, then SVD) instead of
                        % ordinary PCA. A travelling wave is rank 2 in the real domain and
                        % rank 1 in the complex one: measured on a synthetic 60 Hz wave,
                        % real PCA split it 49%/49% across two modes while complex PCA put
                        % 96% in one, recovering the phase gradient at corr 1.000.
c.maxPC       = 10;     % at most this many modes kept per segment
c.scoreMin    = 0.01;   % get_oscScore threshold. Calibrated on synthetic traces --
                        % oscillatory ones score >0.2, drift and noise <0.06.
c.dropBelowScore = true;% drop design rows whose re-scored value is below scoreMin,
                        % instead of regressing them out anyway
% -- where in frequency to look --
c.search      = [1 3];  % the spectral peak is searched in fTrk .* [1/search(2), search(2)]
c.bandFrac    = 0.6;    % passband half-width as a fraction of the centre frequency
c.fMin        = 40;     % Hz, HARD FLOOR on every frequency this stage uses. Dendritic Vm
                        % lives below it; RBC passage (v/d with v ~ 0.5-2 mm/s, d ~ 5-15
                        % um) is 40-400 Hz, so the floor costs no vascular signal. It
                        % clamps the search band, the passband AFTER the +/-bandFrac
                        % widening, and the fallback -- all three matter, because a legal
                        % fPeak of 50 Hz would otherwise open the band down to 20 Hz.
c.dropPinned  = true;   % drop a segment whose peak lands ON the floor: its real peak is
                        % below fMin, so the band was imposed rather than measured.
c.smoothMs    = 50;     % ms, smoothing of the DIAGNOSTIC-ONLY detection-count trace
% -- the voltage-map criterion --
c.voltPen     = true;   % VETO a regressor whose SPATIAL MAP resembles the neuron
c.voltRef     = 'auto'; % 'auto' = 'sta' when there are more than voltMinSpk spikes,
                        %          else 'ftprnt'. The STA is the direct spatial signature
                        %          of voltage and blood cannot contaminate it (the RBC
                        %          train is not spike-locked); the anatomical footprint is
                        %          a proxy that also flags a vessel running along the
                        %          dendrite. But an STA over too few spikes is noise, and
                        %          a noisy reference correlates with nothing -- which
                        %          switches the criterion off silently. Hence the fallback.
c.voltMinSpk  = 20;     % spikes required before the STA is trusted
c.voltCut     = 0.5;    % |r| at or above this and the row is DROPPED, whatever it scored.
                        % A hard cut, not a ramp: a regressor whose map looks like the
                        % neuron is not a weak vascular regressor to be down-weighted, it
                        % is the wrong thing entirely, and subtracting it removes signal.
                        % There is no |r| at which eating some voltage becomes acceptable,
                        % so there is nothing for a graded multiplier to express.
                        % The cut is on |r|: a map that is the NEGATIVE of the footprint is
                        % the same pattern with the regressor's arbitrary sign flipped.
                        % A vessel running ALONG the dendrite overlaps the footprint and
                        % can be vetoed honestly -- res.rVoltAlt carries the correlation
                        % against the other reference so you can check that case by hand.
% -- design matrix --
c.addQuad     = false; % give every design row a 90-degree-shifted twin, imag(hilbert(x)),
                       % where it does not already have one.
                       % READ THIS BEFORE TURNING IT ON: cos and sin are ALREADY 90 apart,
                       % and under complex PCA so are the real and imaginary parts of every
                       % extra mode. A quadrature pair spans EVERY phase -- a cos b sin
                       % reaches any A cos(phi - theta) -- so the per-pixel fit already
                       % matches any phase offset between a segment and a pixel, and
                       % copying those rows adds nothing but collinearity. This therefore
                       % only fills rows that genuinely lack a partner: extra modes under
                       % REAL pca, and the derivative rows. With the defaults (cpca = true,
                       % useDeriv = 0) it adds zero rows, and it says so.
                       % If vessel signal survives the correction, that zero is itself the
                       % answer: phase coverage is not what is missing.
c.useDeriv    = 0;      % 0 = cos/sin only, 1 = also d/dt, 2 = also d2/dt2. Mostly
                        % redundant: d/dt(A cos) is the sin term weighted by the
                        % instantaneous frequency (corr -1.000 at constant frequency,
                        % -0.964 at 40% FM), so it only adds a direction when omega swings.
cfg.regressors = c;

% ===================== flow_apply  (7g) ============================================
c = struct();
c.fitOn      = 'mov_res';% WHICH MOVIE THE CORRECTION IS FITTED TO AND SUBTRACTED FROM.
                        % 'mov_res' = the RAW motion residual, with the neuron's voltage
                        %             still in it. This is the honest final test: what
                        %             does the regressor actually do to the data you
                        %             started from? The cost is that voltage is present
                        %             for the regressor to eat if the two correlate,
                        %             which is exactly what the dendrite-vs-vessel
                        %             numbers and the voltage-map penalty are there to
                        %             catch.
                        % 'mov_nv'  = the voltage-regressed movie the band-pass built.
                        %             Safer -- the regressor cannot take voltage that is
                        %             already gone -- but then the reported "variance
                        %             removed" is not variance removed from your data,
                        %             and flow_compare is comparing against a movie you
                        %             would not otherwise use.
                        % Either way data.mov_res is never modified, so flow_compare can
                        % always put the corrected traces next to the ORIGINAL ones.
c.bandMode   = 'refit'; % 'refit' = re-find the PSD peak on the full-length trace
                        % 'reuse' = keep the value flow_regressors measured
c.refitPCA   = true;    % re-run the per-segment PCA on the whole recording and re-select
                        % the modes, reusing only the segment masks and the v/d band.
                        % The masks are geometry and do not change with time; the modes
                        % are a STATISTICAL estimate from a short window and are the one
                        % reused quantity worth redoing on more data.
c.fMin       = [];      % Hz floor for the refit. [] = inherit cfg.regressors.fMin
c.fitChunk   = 500;     % frames per pass when accumulating the normal equations. PURELY a
                        % memory/speed knob -- the answer is identical whatever it is set
                        % to; lower it if the fit runs out of memory.
c.blockS     = 1;       % s, block length for the drift / stability diagnostic
% -- the time-varying-coefficient test --
c.split      = true;    % refit the coefficients in overlapping windows and compare
c.splitWin   = 1000;    % frames per window
c.splitOvl   = 200;     % frames of overlap (linear cross-fade, weights sum to 1)
c.cvChunk    = 50;      % frames per held-out block. WHOLE BLOCKS, not single frames: a
                        % lone held-out frame is predictable from its neighbours, which
                        % would hand the more flexible model a free pass.
c.cvFrac     = 0.2;     % fraction held out. The split fit has ~nWin times the parameters
                        % so it wins IN-SAMPLE always; only the held-out score decides.
c.splitMinR  = 5;       % training frames per regressor a window needs, else it falls back
                        % to the global coefficients
cfg.apply = c;

% ===================== flow_compare  (7h) ==========================================
c = struct();
c.bvFrac    = 0.2;      % a sub-ROI is "on a vessel" if this fraction of it is in bvMask
c.baseWin   = 1001;     % ms, display baseline (running median)
c.psdWin    = 1024;     % samples per PSD window (50% overlap, Hann)
c.staWin    = [-30 30]; % ms around each spike for the STA
c.staBase   = [-30 -10];% ms, STA baseline window
c.staPeak   = [0 5];    % ms, where the spike amplitude is measured
c.normMode  = 'bap';    % 'bap'  divide by the bAP peak measured on each trace's own STA
                        %        -> 1.0 = one spike, and ROIs of different depth become
                        %        comparable. Measuring the bAP from the trace itself is
                        %        what keeps this independent of how the trace was
                        %        extracted: these project onto the UNNORMALISED
                        %        footprints, so an externally fitted scale factor would
                        %        be off by an unknown gain.
                        % 'own'  each trace by its own fast-band noise SD (no spikes
                        %        needed; 'bap' falls back to this when there are too few)
cfg.compare = c;

% ===================== shared ======================================================
c = struct();
% Everything the GUI needs for display it works out for itself; only quantities that
% change a RESULT live here.
c.fs        = 1000;     % Hz
c.smoothMs  = 20;       % ms, the paper's subthreshold window (the "fast" split point)
cfg.common = c;
end
