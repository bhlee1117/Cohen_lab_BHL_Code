function [MaxLocs, EdgeLocs] = spikefindconv(dat, ker, thresh)
% Function PeakLocs = spikefindconv(dat, ker, thresh);
% Inputs: 
% dat - data
% ker - convolution kernel, 
% thresh - faction of maximum to make threshold (between 0 and 1)

% Outputs:
% PeakLocs - indices of the peaks

% Convolves the data with the kernel, which should be an approximation for 
% peak shape.  Pulls out the locations where the trace rises above and 
% later drops below the threshold, and finds the index of the raw data on
% that range.  Will not find peaks within (length(ker)-1)/2 of the data edge. 
% DRH and AEC 24 Feb. 2011
% Kit Sept 2013

C = conv(dat,fliplr(ker),'valid');
C = C/max(C);
ConvOffset = ceil((length(ker)-1)/2);

% figure(111); clf;
% plot([zeros(1,ConvOffset) C'])
% hold all
% plot(dat)
% plot(ker)
% plot(thresh*ones(1,length(C)))

spikeon = find(C(2:end) > thresh & C(1:end-1) < thresh);
spikeoff = find(C(2:end) < thresh & C(1:end-1) > thresh);
% correct for convolution edges
spikeon = spikeon + ConvOffset;
spikeoff = spikeoff + ConvOffset;

if isempty(spikeon) || isempty (spikeoff)
    MaxLocs = [];
    return
end;
if spikeoff(1) < spikeon(1);
    spikeoff(1) = [];
end;

if spikeon(end) > spikeoff(end);
    spikeon(end) = [];
end;

npeaks = length(spikeon);
MaxLocs = zeros(1,npeaks);
EdgeLocs = zeros(1,npeaks);
for j = 1:npeaks;
    [~, indx] = max(dat(spikeon(j):spikeoff(j)));
    MaxLocs(j) = indx + spikeon(j)-1;
    [~, indx] = max(C((spikeon(j):spikeoff(j))-ConvOffset));
    EdgeLocs(j) = indx + spikeon(j)-1;
end;

