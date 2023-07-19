function [spike pks prom]=find_spike_bh(volt_hi,threshold,prom_th)
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