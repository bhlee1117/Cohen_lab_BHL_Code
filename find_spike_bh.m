function [spike pks prom]=find_spike_bh(volt_hi,threshold,prom_th)

% Function written by Byung Hun Lee, 2022.06, Cohen Lab.
% Given voltage trace, N cell X T frames, finding a peak that the amplitude is higher than
% threshold & the local amplitude is higher than prom_th. 
% For typical z-scored voltage trace, threshold of ~ 6, prom_th of ~ 4, works ok.

if nargin<3
    prom_th=0;
end
for i=1:size(volt_hi,1)
    spike(i,:)=zeros(1,size(volt_hi,2));
    [pks s_tmp width prom]=findpeaks(volt_hi(i,:));
    volt_hi2(i,:)=volt_hi(i,:)-movmedian(volt_hi(i,:),7,2);
    prom=volt_hi2(i,s_tmp);
    S_available=find(pks>threshold & prom>prom_th);
    s_tmp=s_tmp(S_available);
    pks=pks(S_available);
    prom=prom(S_available);
    spike(i,s_tmp)=1;
end
end