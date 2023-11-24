function reg_spike = remove_closeSpike(spike,tr,threshold)
reg_spike=zeros(1,length(spike));
sp_time=find(spike);
ISI=sp_time(2:end)-sp_time(1:end-1);
while ~isempty(find(ISI<threshold))
rmv_spike=[];
for s=find(ISI<threshold)

closeSpikes=[s s+1];
if tr(s)>tr(s+1)
rmv_spike=[rmv_spike s+1];
else
rmv_spike=[rmv_spike s];
end

end
sp_time(rmv_spike)=[];
ISI=sp_time(2:end)-sp_time(1:end-1);
end

reg_spike(sp_time)=1;
end