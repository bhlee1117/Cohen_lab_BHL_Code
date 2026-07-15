function [Lap_FiringRate, Lap_SpikeNumber, Lap_NframeValid]=PlaceTrigger_average(spike,place_bin,Virmen_data,vel_thresh,lap_dist)
% Bin a spike train (or per-ROI spike matrix) by lap x track-position.
%
%   spike        nROI x nT   spike counts (0/1 or rate); NaN samples are ignored
%   place_bin    number of position bins along the track
%   Virmen_data  VR matrix: row 5 = track position, row 8 = lap, last row = velocity
%   vel_thresh   running threshold on velocity (frames below it are "not running")
%   lap_dist     track length (same units as position)
%
%   Lap_FiringRate   nLaps x place_bin x nROI   spikes / valid-running-frames
%   Lap_SpikeNumber  nLaps x place_bin x nROI   summed spikes (running only)
%   Lap_NframeValid  nLaps x place_bin x nROI   running & ~NaN frame counts
%
% Vectorized with accumarray. The previous implementation built a dense
% nLaps x place_bin x nT one-hot array (tens of GB on real recordings) inside a
% place_bin x nLaps double loop; this computes the same binning in one O(nT) pass.

lap_trace = Virmen_data(8,:);
vel_trace = Virmen_data(end,:);
nLaps     = max(lap_trace);
nROI      = size(spike,1);

% position bin (1..place_bin) and running mask
bin_dist    = ceil(Virmen_data(5,:)/(lap_dist/place_bin));
not_running = vel_trace < vel_thresh;

% running-gated spikes with NaN -> 0 (matches original spike_run)
spike_run = double(spike);
spike_run(:,not_running)   = 0;
spike_run(isnan(spike_run)) = 0;

% frames that fall in a valid (integer lap in 1..nLaps, bin in 1..place_bin)
valid = bin_dist>=1 & bin_dist<=place_bin & lap_trace>=1 & lap_trace<=nLaps & lap_trace==round(lap_trace);
K  = nLaps*place_bin;
c  = lap_trace + (bin_dist-1)*nLaps;   % linear (lap,bin) index, column-major over [nLaps place_bin]
cv = c(valid)';

Lap_SpikeNumber = zeros(nLaps,place_bin,nROI);
Lap_NframeValid = zeros(nLaps,place_bin,nROI);
runNoNaN = double(~not_running);
for r=1:nROI
    sr = spike_run(r,:);
    Lap_SpikeNumber(:,:,r) = reshape(accumarray(cv,sr(valid)',[K 1]),nLaps,place_bin);
    w  = runNoNaN .* double(~isnan(spike(r,:)));   % running & ~NaN, per ROI
    Lap_NframeValid(:,:,r) = reshape(accumarray(cv,w(valid)',[K 1]),nLaps,place_bin);
end

Lap_FiringRate = Lap_SpikeNumber ./ Lap_NframeValid;
end
