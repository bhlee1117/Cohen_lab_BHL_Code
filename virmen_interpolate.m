function [Virmen_data_int vel_trace]=virmen_interpolate(Virmen_data,lap_dist,interpolateBin)

% Calculate converting bin time
t_VR=Virmen_data(1,:);
[t_VR arg]=unique(t_VR);
Virmen_data=Virmen_data(:,arg);
new_bin=interpolateBin;

% Lap change
lap_end=[0; find(abs(Virmen_data(8,2:end)-Virmen_data(8,1:end-1))>0)'; size(Virmen_data,2)];
Virmen_data(5,lap_end(1:end-1)+1)=0;
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



Virmen_data_int(8,:)=cumsum([0 (Virmen_data_int(5,2:end)-Virmen_data_int(5,1:end-1))<-lap_dist*0.9])+1;
%Virmen_data_int(8,:)=round(Virmen_data_int(8,:));
Virmen_data_int(11,:)=round(Virmen_data_int(11,:));

% Calculate speed
cum_trace=Virmen_data_int(end,:);%movmean(Virmen_data_int(end,:),200);
vel_trace=(cum_trace(2:end)-cum_trace(1:end-1));
vel_trace(end+1)=vel_trace(end);
end