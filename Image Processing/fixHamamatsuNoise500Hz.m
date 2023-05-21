function [out, correctionTrace] = fixHamamatsuNoise500Hz(in);

% function [out, correctionTrace] = fixHamamatsuNoise500Hz(in);
% Assumes that the Hamamatsu noise is a multiplicative error approximately
% every fourth frame; with phase shifts approximately every 120 frames.
% To apply the correction trace to another intensity trace run: 
%
% newIntens = rawIntens./correctionTrace
% 
% Eli Weinstein and AEC 5/1/2015

nt = length(in);
inSub = in - medfilt1(in, 4);
nsub = ceil(nt/4);
inReshape = zeros(4, nsub);  % input intensity keeping background
subDat = zeros(4, nsub);  % input intensity without background
% Reshape the data into four phase-shifted sets, sampled every fourth point
for j = 1:4;
    inReshape(j,1:length(j:4:nt)) = in(j:4:nt);
    subDat(j,1:length(j:4:nt)) = inSub(j:4:nt);
end;

% plot(subDat')

dSubDat = subDat - repmat(mean(subDat), [4 1]);  % Subtract off the mean trace
% plot(dSubDat')
dSubDatS = zeros(size(dSubDat));
for j = 1:4;
    dSubDatS(j,:) = medfilt1(dSubDat(j,:), 10);  % median filter to identify transitions.
end;
% plot(dSubDatS')

[maxVal, maxIdx] = max(dSubDatS);  % At each timepoint find the phase of the extra-amplified signal
% hold all;
% plot(maxVal, 'r*')
% hold off;

% Check the transition times and adjust for errors
nChange = 1;
nIter = 0;
while nIter < 10 && nChange > 0;
    nChange = 0;
    phaseIdx = find(diff(maxIdx) ~= 0);  % find the transition times
    for j = phaseIdx;
        if dSubDat(maxIdx(j+1),j) > dSubDat(maxIdx(j),j)  % If the transition happened too early
            maxIdx(j) = maxIdx(j+1);  % shift it back by one.
            nChange = nChange + 1;
        end;
    end;
    nIter = nIter + 1;
    [nIter nChange]
end;

% plot(maxIdx)
maxTrace = zeros(nsub, 1);  % The trace of all the extra-amplified signals
minTrace = zeros(nsub, 1);  % The trace of all the correct signals
for j = 1:nsub;
    maxTrace(j) = inReshape(maxIdx(j), j);
    minTrace(j) = mean(inReshape(setdiff(1:4, maxIdx(j)), j));
end;

% plot(inReshape'); hold all
% plot(1:nsub, maxTrace, 'r*', 1:nsub, minTrace, 'k-')
% plot(phaseIdx, maxTrace(phaseIdx), 'co') 
% hold off

scaleFact = median(maxTrace./minTrace);  % Scale factor

correctionMat = ones(4, nsub);
for j = 1:nsub;
    correctionMat(maxIdx(j),j) = scaleFact;
end;

fixDat = inReshape./correctionMat;

out = reshape(fixDat, [1, nsub*4]);
correctionTrace = reshape(correctionMat, [1, nsub*4]);
out = out(1:nt);
correctionTrace = correctionTrace(1:nt);
