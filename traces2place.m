function [trace_place SI bootstrap_SI]=traces2place(traces,place_bin,Virmen_data,vel_thresh,N_it)

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
disp(['Reward spot: ' num2str(ceil(reward_spot/((lap_dist)/place_bin))) 'th bin, Total laps: ' num2str(size(laps,1))])
traces_run=traces;
%traces_run(zscore(traces,0,2)<-4)=NaN;
traces_run(:,vel_trace<vel_thresh)=NaN;
trace_place=NaN(size(traces,1),place_bin);
for p=1:place_bin
spot_stay=find(bin_dist==p & (vel_trace>vel_thresh)'); 
lap_spot_stay(p,:)=(bin_dist==p)'.*lap_trace;
trace_place(:,p)=sum(traces_run(:,spot_stay),2,'omitnan')./(length(spot_stay)*1.25*1e-3); %divided by time resting in the place (even not walking)
%trace_place(:,p)=sum(traces_run(:,spot_stay),2,'omitnan')./(1.25*1e-3); 
%trace_place(:,p)=mean(traces_run(:,spot_stay),2,'omitnan')./(1.25*1e-3);
%spatial information
%pi(p)=sum(bin_dist_run==p)./sum(~isnan(bin_dist_run));
pi(p)=sum(bin_dist==p & (vel_trace>vel_thresh)')/sum((vel_trace>vel_thresh));
end

lambda=sum(traces_run,2,'omitnan')/sum(vel_trace>vel_thresh)*800;
SI=sum(pi.*trace_place.*log(trace_place./lambda)/log(2),2,'omitnan');

nbytes = fprintf('processing iteration %0 of %d', N_it);
for it=1:N_it

    valid_ind=find(~isnan(traces_run(1,:)));
    shuffle_ind=valid_ind(randperm(length(valid_ind),length(valid_ind)));
    shuffle_trace(:,valid_ind)=traces_run(:,shuffle_ind);
    shuffle_trace(:,setdiff([1:size(traces,2)],valid_ind))=NaN;
   
for p=1:place_bin
spot_stay=find(bin_dist==p & (vel_trace>vel_thresh)'); 
lap_spot_stay(p,:)=(bin_dist==p)'.*lap_trace;
shuffle_trace_place(:,p)=sum(shuffle_trace(:,spot_stay),2,'omitnan')./(length(spot_stay)*1.25*1e-3); %divided by time resting in the place (even not walking)
end
fprintf(repmat('\b',1,nbytes))
nbytes = fprintf('processing iteration %d of %d', it, N_it);
%shuffle=(randperm([place_bin],size(traces,1),place_bin)-1)*size(traces,1)+[1:size(traces,1)]';
bootstrap_SI(:,it)=sum(pi.*shuffle_trace_place.*log(shuffle_trace_place./lambda)/log(2),2,'omitnan');
end


%percentile_95=prctile(bootstrap_SI',perc)';
%sort_bSI=sort(bootstrap_SI,2,'ascend');
%percentile_95s=sort_bSI(:,950);
%put_pc=find((SI-percentile_95)>0);
%trace_place_norm=(trace_place-min(trace_place,[],2));
%trace_place_norm=trace_place_norm./max(trace_place_norm,[],2);

end


