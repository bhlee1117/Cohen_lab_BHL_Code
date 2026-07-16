function [spike, pks, prom]=find_spike_bh(volt_hi,threshold,prom_th,min_interval)

% Function written by Byung Hun Lee, 2022.06, Cohen Lab. (sped up 2026)
% Given voltage trace, N cell X T frames, finding a peak that the amplitude is higher than
% threshold & the local amplitude is higher than prom_th.
% For typical z-scored voltage trace, threshold of ~ 6, prom_th of ~ 4, works ok.
%
% pks, prom (optional) correspond to the LAST processed row, as before.

if nargin<3 || isempty(prom_th)
    prom_th=0;
end
if nargin<4
    min_interval=[];
end

[nR, nT] = size(volt_hi);
spike = zeros(nR, nT);
% 7-pt high-pass used as the local-amplitude ("prominence") proxy; vectorized once
volt_hi2 = volt_hi - movmedian(volt_hi, 7, 2);
pks = []; prom = [];

for i=1:nR
    % Ask findpeaks for peaks/locations only (skips its costly prominence & width),
    % and let it prune on height during detection.
    if ~isempty(min_interval)
        [pk, s_tmp] = findpeaks(volt_hi(i,:),'MinPeakHeight',threshold,'MinPeakDistance',min_interval);
    else
        [pk, s_tmp] = findpeaks(volt_hi(i,:),'MinPeakHeight',threshold);
    end
    pr = volt_hi2(i,s_tmp);
    keep = pk>threshold & pr>prom_th;
    s_tmp = s_tmp(keep);
    spike(i,s_tmp) = 1;
    pks = pk(keep); prom = pr(keep);   % last-row outputs (unchanged contract)
end
end
