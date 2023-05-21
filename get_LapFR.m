function [Lap_FR lap_spot_stay Fr_NaN_pc]=get_LapFR(spike,place_bin,Ard_data,vel_thresh,noi)

lap_length=max(Ard_data(:,2));
Ard_data(end-1:end,2:4)=repmat(Ard_data(end-2,2:4),2,1);
lap_end=[0; find(abs(Ard_data(2:end,2)-Ard_data(1:end-1,2))>6500); size(Ard_data,1)];
laps=[lap_end(1:end-1)+1 lap_end(2:end)];

for l=1:size(laps,1)
lap_trace(laps(l,1):laps(l,2))=l; 
cum_trace(laps(l,1):laps(l,2))=Ard_data(laps(l,1):laps(l,2),2)+lap_length*l; end

cum_trace=movmean(cum_trace,6);
vel_trace=(cum_trace(2:end)-cum_trace(1:end-1))/1.25*1000;
vel_trace(end+1)=vel_trace(end);
Ard_data(:,2)=Ard_data(:,2)-min(Ard_data(:,2))+1;
reward_spot=mean(Ard_data(Ard_data(:,3)==1,2));

lap_dist=max(Ard_data(:,2));
bin_dist=ceil(Ard_data(:,2)/((lap_dist)/place_bin));
cmap=jet(place_bin);
spike_run=spike; spike_run(:,vel_trace<vel_thresh)=NaN;
for p=1:place_bin
spot_stay=find(bin_dist==p); 
lap_spot_stay(p,:)=(bin_dist==p)'.*lap_trace;
end

Lap_FR=NaN(max(lap_spot_stay(:)),place_bin,length(noi));

for p=1:place_bin
    for l=1:max(lap_spot_stay(p,:))
        if ~isempty(find(lap_spot_stay(p,:)==l)) %only the laps the mouse went to the place bin

    %Lap_FR(l,p,:)=mean(spike_run(noi,find(lap_spot_stay(p,:)==l)),2,'omitnan')/(1.25*1e-3);
    lg(p,l)=length(find(lap_spot_stay(p,:)==l & vel_trace>vel_thresh)); %Divided by time spent at bin + running
    if lg(p,l)==0
    Lap_FR(l,p,:)=NaN;    
    else
    Lap_FR(l,p,:)=sum(spike_run(noi,find(lap_spot_stay(p,:)==l)),2,'omitnan')/(lg(p,l)*1.25*1e-3);
    %Lap_FR(l,p,:)=mean(spike_run(noi,find(lap_spot_stay(p,:)==l)),2,'omitnan')/(1.25*1e-3);
    end
    %Lap_FR(l,p,:)=sum(spike_run(noi,find(lap_spot_stay(p,:)==l)),2,'omitnan'); %total spikes
    %Lap_FR(l,p,:)=sum(spike_run(noi,find(lap_spot_stay(p,:)==l)),2,'omitnan')/(sum(~isnan(spike_run(1,find(lap_spot_stay(p,:)==l))))*1.25*1e-3);
        end
    end
end
%figure; imagesc(lg)
end
