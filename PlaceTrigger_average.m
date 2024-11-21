function [Lap_FiringRate Lap_SpikeNumber Lap_NframeValid]=PlaceTrigger_average(spike,place_bin,Virmen_data,vel_thresh,lap_dist)

%take VR time
t_VR = Virmen_data(1,:);
rate=t_VR(3)-t_VR(2);
% Lap change
lap_end=[0; find(abs(Virmen_data(8,2:end)-Virmen_data(8,1:end-1))>0)'; size(Virmen_data,2)];
laps=[lap_end(1:end-1)+1 lap_end(2:end)];
disp('Calculating Place trigger average..')
% 
% %Cumulative track
% cumTrack=[];
% Virmen_data(5,:)=Virmen_data(5,:)-min(Virmen_data(5,:));
% cumTrack=[cumTrack Virmen_data(5,laps(1,1):laps(1,2))];
% for l2=2:size(laps,1)
%     %cumTrack=[cumTrack Virmen_data(5,laps(l2,1):laps(l2,2))+cumTrack(laps(l2,1)-1)];
%     cumTrack=[cumTrack Virmen_data(5,laps(l2,1):laps(l2,2))+(l2-1)*lap_dist];
% end
% Virmen_data(end+1,:)=cumTrack;
% 
% % interpolate
% Virmen_data_int(1,:)=new_bin;
% for i=2:size(Virmen_data,1)
%     Virmen_data_int(i,:)=interp1(t_VR,Virmen_data(i,:),new_bin,'linear');
% end
% Virmen_data_int(5,:)=Virmen_data_int(end,:);
% 
% [spike_time arg]=unique(spike_time);
% spike=spike(:,arg);
% spike_int=zeros(1,length(new_bin));
% if length(unique(spike))<1
% distances = pdist2(spike_time', new_bin'); 
% [~, nearest]=min(distances,[],1);
% spike_int=spike(nearest);    
% else
% spike_int=interp1(spike_time,spike,new_bin,'linear');
% end
% 
% ll=1; laps_int=[1];
% lap_trace=zeros(1,size(Virmen_data_int,2));
% calculate track back from cumulative track
lap_trace=Virmen_data(8,:);
% 
% Virmen_data_int(8,:)=round(Virmen_data_int(8,:));
% 
% % Calculate speed
% cum_trace=Virmen_data_int(end,:);%movmean(Virmen_data_int(end,:),200);
% vel_trace=(cum_trace(2:end)-cum_trace(1:end-1));
%vel_trace(end+1)=vel_trace(end);
vel_trace=Virmen_data(end,:);

% Calculate mean value at each position
bin_dist=ceil(Virmen_data(5,:)/((lap_dist)/place_bin));
cmap=jet(place_bin);
not_running=vel_trace<vel_thresh;
spike_run=double(spike); spike_run(:,not_running)=0; spike_run(isnan(spike_run))=0;
LapPosMat=zeros(max(lap_trace),place_bin,length(t_VR));

for p=1:place_bin
    %spot_stay=find(bin_dist==p);
    %lap_spot_stay(p,:)=(bin_dist==p).*lap_trace;
    for l=1:max(lap_trace)
    LapPosMat(l,p,find((bin_dist==p) & lap_trace ==l))=1;
    end
end
LapPosMat=tovec(LapPosMat);
%lap_spot_stay(:,not_running)=NaN;
%Lap_FR=NaN(max(lap_spot_stay(:)),place_bin);
%Lap_V=NaN(max(lap_spot_stay(:)),place_bin);

Lap_SpikeNumber=permute(reshape(spike_run*LapPosMat',[],max(lap_trace),place_bin),[2 3 1]);
Lap_NframeValid=permute(reshape((double(~not_running) .* double(~isnan(spike)))*LapPosMat',[],max(lap_trace),place_bin),[2 3 1]); %run & ~nan
Lap_Nframe=reshape(ones(1,length(t_VR))*LapPosMat',[],place_bin); %run & ~nan
Lap_FiringRate=Lap_SpikeNumber./Lap_NframeValid;
Lap_velocity=reshape(vel_trace*LapPosMat',[],place_bin)./Lap_Nframe;
% for p=1:place_bin
%     for l=1:max(lap_spot_stay(p,:))
% 
%         AtPosL=find(lap_spot_stay(p,:)==l);
%         AtPosL_Run=find(lap_spot_stay(p,:)==l & ~not_running);
% 
%         if ~isempty(find(lap_spot_stay(p,:)==l)) %only the laps the mouse went to the place bin
% 
%             %Lap_FR(l,p,:)=mean(spike_run(noi,find(lap_spot_stay(p,:)==l)),2,'omitnan')/(1.25*1e-3);
%             lg(p,l)=length(AtPosL_Run); % number of frames stayed in l th lap, p position bin
%             lgV(p,l)=length(AtPosL);
%             if lg(p,l)==0   % mouse ran so fast that hard to capture at the position
%                 Lap_FR(l,p,:)=NaN;
%             else
%                 switch rate_switch
%                     case 'rate'
%                 %Lap_FR(l,p,:)=sum(spike_run(1,AtPosL),2,'omitnan')/(lg(p,l)*rate); %number of spike divided by stayed time
%                 Lap_FR(l,p,:)=mean(spike_run(1,AtPosL),2,'omitnan'); %number of spike divided by stayed time
%                     case 'sum'
%                 Lap_FR(l,p,:)=sum(spike_run(1,AtPosL),2,'omitnan'); %number of spike divided by stayed time
%                 end
%             end
% 
%             if lgV(p,l)==0
%                 p
%                 l
%                 Lap_V(l,p,:)=NaN;
%             else
%                 %Lap_V(l,p,:)=sum(vel_trace(noi,find(lap_spot_stay(p,:)==l)),2,'omitnan')/(lgV(p,l)*rate); %number of spike divided by stayed time
%                 Lap_V(l,p,:)=mean(vel_trace(1,AtPosL),2,'omitnan');%/(lgV(p,l)*rate); %number of spike divided by stayed time
%             end
%         end
%     end
% end
%figure; imagesc(lg)
end
