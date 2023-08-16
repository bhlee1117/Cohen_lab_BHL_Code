function [Lap_FR lap_spot_stay]=get_LapFR_Virmen(spike,place_bin,lap_length,Virmen_data,vel_thresh  )
noi=[1:size(spike,1)];
Virmen_data(:,5)=Virmen_data(:,5)-min(Virmen_data(:,5))+1;
%Ard_data(end-1:end,2:4)=repmat(Ard_data(end-2,2:4),2,1);
Virmen_data(:,8)=round(Virmen_data(:,8));
lap_end=[0; find(abs(Virmen_data(2:end,8)-Virmen_data(1:end-1,8))>0); size(Virmen_data,1)];
laps=[lap_end(1:end-1)+1 lap_end(2:end)];
rate=median(Virmen_data(2:end,1)-Virmen_data(1:end-1,1));

for l=1:size(laps,1)
lap_trace(laps(l,1):laps(l,2))=l; end

cum_trace=movmean(Virmen_data(:,13),50);
vel_trace=(cum_trace(2:end)-cum_trace(1:end-1))';
vel_trace(end+1)=vel_trace(end);

lap_dist=lap_length;
bin_dist=ceil(Virmen_data(:,5)/((lap_dist)/place_bin));
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
    Lap_FR(l,p,:)=sum(spike_run(noi,find(lap_spot_stay(p,:)==l)),2,'omitnan')/(lg(p,l)*rate);
    %Lap_FR(l,p,:)=mean(spike_run(noi,find(lap_spot_stay(p,:)==l)),2,'omitnan')/(1.25*1e-3);
    end
    %Lap_FR(l,p,:)=sum(spike_run(noi,find(lap_spot_stay(p,:)==l)),2,'omitnan'); %total spikes
    %Lap_FR(l,p,:)=sum(spike_run(noi,find(lap_spot_stay(p,:)==l)),2,'omitnan')/(sum(~isnan(spike_run(1,find(lap_spot_stay(p,:)==l))))*1.25*1e-3);
        end
    end
end
%figure; imagesc(lg)
end
