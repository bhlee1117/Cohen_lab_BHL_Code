function [Lap_FR Lap_V]=get_LapFR_VU(spike,spike_time,place_bin,Virmen_data,vel_thresh,noi,lap_dist)

%take VR time
t_VR = datetime(datetime(Virmen_data(1,:),'ConvertFrom','datenum'), 'InputFormat', 'yyyy-MM-dd HH:mm:ss.SSS');
t_VR= t_VR-t_VR(1);
t_VR= milliseconds(t_VR)/1000;
spike_time=spike_time-spike_time(1);

% Calculate converting bin time
[t_VR arg]=unique(t_VR);
Virmen_data=Virmen_data(:,arg);
%spike=spike(arg);
rate=(t_VR(end)-t_VR(1))/(size(spike,2)-1);
disp(['Calculated rate is ',num2str(rate) ,'sec/frame'])
new_bin=[t_VR(1):rate:t_VR(end)];

% Lap change
lap_end=[0; find(abs(Virmen_data(8,2:end)-Virmen_data(8,1:end-1))>0)'; size(Virmen_data,2)];
laps=[lap_end(1:end-1)+1 lap_end(2:end)];

%Cumulative track
cumTrack=[];
Virmen_data(5,:)=Virmen_data(5,:)-min(Virmen_data(5,:));
cumTrack=[cumTrack Virmen_data(5,laps(1,1):laps(1,2))];
for l2=2:size(laps,1)
    %cumTrack=[cumTrack Virmen_data(5,laps(l2,1):laps(l2,2))+cumTrack(laps(l2,1)-1)];
    cumTrack=[cumTrack Virmen_data(5,laps(l2,1):laps(l2,2))+(l2-1)*lap_dist];
end
Virmen_data(end+1,:)=cumTrack;

% interpolate
Virmen_data_int(1,:)=new_bin;
for i=2:size(Virmen_data,1)
    Virmen_data_int(i,:)=interp1(t_VR,Virmen_data(i,:),new_bin,'linear');
end
Virmen_data_int(5,:)=Virmen_data_int(end,:);

[spike_time arg]=unique(spike_time);
spike=spike(:,arg);
spike_int=zeros(1,length(new_bin));
if length(unique(spike))<1
distances = pdist2(spike_time', new_bin'); 
[~, nearest]=min(distances,[],1);
spike_int=spike(nearest);    
else
spike_int=interp1(spike_time,spike,new_bin,'linear');
end

ll=1; laps_int=[1];
lap_trace=zeros(1,size(Virmen_data_int,2));
% calculate track back from cumulative track
while ~isempty(find(Virmen_data_int(5,:)>lap_dist))
sub_ind=find(Virmen_data_int(5,:)>lap_dist);
laps_int(ll,2)=sub_ind(1);
lap_trace(laps_int(ll,1):laps_int(ll,2))=ll;
Virmen_data_int(5,sub_ind)=Virmen_data_int(5,sub_ind)-lap_dist;
ll=ll+1;
laps_int(ll,1)=sub_ind(1)+1;
end
laps_int(ll,2)=size(Virmen_data_int,2);
lap_trace(laps_int(ll-1,2)+1:end)=ll;

Virmen_data_int(8,:)=round(Virmen_data_int(8,:));

% Calculate speed
cum_trace=Virmen_data_int(end,:);%movmean(Virmen_data_int(end,:),200);
vel_trace=(cum_trace(2:end)-cum_trace(1:end-1));
vel_trace(end+1)=vel_trace(end);

% Calculate mean value at each position
bin_dist=ceil(Virmen_data_int(5,:)/((lap_dist)/place_bin));
cmap=jet(place_bin);
not_running=vel_trace<vel_thresh;
spike_run=spike_int; spike_run(:,not_running)=NaN;
for p=1:place_bin
    spot_stay=find(bin_dist==p);
    lap_spot_stay(p,:)=(bin_dist==p).*lap_trace;
end
%lap_spot_stay(:,not_running)=NaN;
Lap_FR=NaN(max(lap_spot_stay(:)),place_bin,length(noi));
Lap_V=NaN(max(lap_spot_stay(:)),place_bin,length(noi));

for p=1:place_bin
    for l=1:max(lap_spot_stay(p,:))
        if ~isempty(find(lap_spot_stay(p,:)==l)) %only the laps the mouse went to the place bin

            %Lap_FR(l,p,:)=mean(spike_run(noi,find(lap_spot_stay(p,:)==l)),2,'omitnan')/(1.25*1e-3);
            lg(p,l)=length(find(lap_spot_stay(p,:)==l & ~not_running)); % number of frames stayed in l th lap, p position bin
            lgV(p,l)=length(find(lap_spot_stay(p,:)==l));
            if lg(p,l)==0
                Lap_FR(l,p,:)=NaN;
            else
                Lap_FR(l,p,:)=sum(spike_run(noi,find(lap_spot_stay(p,:)==l)),2,'omitnan')/(lg(p,l)*rate); %number of spike divided by stayed time
            end

            if lgV(p,l)==0
                p
                l
                Lap_V(l,p,:)=NaN;
            else
                %Lap_V(l,p,:)=sum(vel_trace(noi,find(lap_spot_stay(p,:)==l)),2,'omitnan')/(lgV(p,l)*rate); %number of spike divided by stayed time
                Lap_V(l,p,:)=mean(vel_trace(noi,find(lap_spot_stay(p,:)==l)),2,'omitnan');%/(lgV(p,l)*rate); %number of spike divided by stayed time
            end
        end
    end
end
%figure; imagesc(lg)
end
