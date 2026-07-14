function volt_frames = gluPeaks2VoltFrames(spikeFrames, gluTrace, glu_tax, volt_tax, peakWin)
% Refine glutamate spike times to voltage-frame triggers via spline interpolation.
% For each glutamate spike, the peri-spike glutamate signal (+/- peakWin glu
% frames) is spline-interpolated onto the (faster) voltage time axis, and the
% voltage frame index at the interpolated peak is returned as the trigger. This
% removes the ~half-glu-frame jitter of nearest-frame matching, so downstream
% voltage triggered-averages are better aligned. Falls back to nearest-frame
% matching for windows too short or off the voltage axis.
nGluFrame   = numel(glu_tax);
volt_frames = zeros(1, numel(spikeFrames));
for k = 1:numel(spikeFrames)
    w  = max(1, spikeFrames(k)-peakWin) : min(nGluFrame, spikeFrames(k)+peakWin);
    tl = glu_tax(w);          % glu times in peri-spike window (s)
    yl = gluTrace(w);         % glu signal in window
    vw = find(volt_tax >= tl(1) & volt_tax <= tl(end));   % voltage frames spanned
    if numel(tl) < 3 || isempty(vw)
        volt_frames(k) = match_nearest(glu_tax(spikeFrames(k)), volt_tax);   % fallback
    else
        yq = interp1(tl, yl, volt_tax(vw), 'spline');     % glu signal on voltage grid
        [~, im] = max(yq);
        volt_frames(k) = vw(im);                          % refined trigger (voltage frame)
    end
end
end
