function out = spikefind(dat, thresh);
% Function out = spikefind(dat, thresh);
% for each time that dat goes above thresh, returns the index of the local
% maximum in dat
% DRH and AEC 24 Feb. 2011

spikeon = find(dat(2:end) > thresh & dat(1:end-1) < thresh);
spikeoff = find(dat(2:end) < thresh & dat(1:end-1) > thresh);


if isempty(spikeon) | isempty (spikeoff)
    out = [];
    return
end;

if spikeoff(1) < spikeon(1);
    spikeoff(1) = [];
end;
if spikeon(end) > spikeoff(end);
    spikeon(end) = [];
end;

for j = 1:length(spikeon);
    [y, indx] = max(dat(spikeon(j):spikeoff(j)));
    out(j) = indx + spikeon(j)-1;
end;

