function [Lap_FR lap_spot_stay]=get_LapFR_Virmen(spike,place_bin,Virmen_data,vel_thresh,noi)

Virmen_data(:,3)=Virmen_data(:,3)-min(Virmen_data(:,3))+1;
lap_length=max(Virmen_data(:,3));
%Ard_data(end-1:end,2:4)=repmat(Ard_data(end-2,2:4),2,1);
lap_end=[0; find(abs(Virmen_data(2:end,4)-Virmen_data(1:end-1,4))>0); size(Virmen_data,1)];
laps=[lap_end(1:end-1)+1 lap_end(2:end)];

for l=1:size(laps,1)
lap_trace(laps(l,1):laps(l,2))=l; end

cum_trace=Virmen_data(:,5);
cum_trace=cum_trace([1:300:length(cum_trace)]);
cum_trace=interp1([1:300:size(Virmen_data,1)],cum_trace,[1:size(Virmen_data,1)],'linear');
vel_trace=(cum_trace(2:end)-cum_trace(1:end-1));
vel_trace(end+1)=vel_trace(end);

bw=bwlabel(Virmen_data(:,2));
for b=1:max(bw)
    tmp=find(bw==b);
    R(b)=tmp(1);
end
reward_spot=mean(Virmen_data(R,3));

lap_dist=max(Virmen_data(:,3));
bin_dist=ceil(Virmen_data(:,3)/((lap_dist)/place_bin));
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
