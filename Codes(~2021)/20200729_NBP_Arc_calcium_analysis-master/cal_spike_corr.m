function binned_spike=cal_spike_corr(peaks,locs,bins,t,p_thres)
location=locs(peaks>p_thres);
binned_spike=zeros(ceil(t/bins),2);
for i=1:size(location,1)
    binned_spike(floor(location(i,1)/bins)+1,2)=...
    [binned_spike(floor(location(i,1)/bins)+1,2)+peaks(i,1)]; % bin time, sum_peak
end
binned_spike(:,1)=[0:bins:bins*(ceil(t/bins)-1)];
end