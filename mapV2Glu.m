function gluIdx = mapV2Glu(spikeTimes, gluAxis, tol_s)
% MAPV2GLU  Map voltage spike times to the nearest glutamate-frame index.
%   gluIdx = mapV2Glu(spikeTimes, gluAxis, tol_s)
%   spikeTimes, gluAxis : times in SECONDS (e.g. VoltResult.t_ax(spk), GluResult.t_ax)
%   tol_s : keep only matches whose |time gap| <= tol_s (seconds)
% Returns a row vector of glutamate-frame indices (into gluAxis).
if isempty(spikeTimes), gluIdx = []; return; end
[idx, d] = match_nearest(spikeTimes, gluAxis);
gluIdx = idx(abs(d) <= tol_s);
gluIdx = gluIdx(:)';
end
