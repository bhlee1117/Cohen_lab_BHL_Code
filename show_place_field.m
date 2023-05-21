function [place_volt order] = show_place_field(volt,trajectory,run,place_bin,lap_dist,reward_pos,im_on)
if nargin<7
    im_on=false;
end
bin_dist=ceil(trajectory/(lap_dist/place_bin)); 
bin_dist(~run)=NaN;
clear place_volt place_volt_lp place_spike
volt=movmean(volt,30,2);
for i=1:size(volt,1)
    for p=1:place_bin
        tmp=find(bin_dist==p);
        place_volt(i,p)=mean(volt(i,tmp));
    end
end
    %place_volt=movmean(place_volt,3,2);
    place_volt=(place_volt-min(place_volt,[],2));
    place_volt=place_volt./range(place_volt,2);

% sort by spike #
[a b]=sort(place_volt,2,'ascend');
put_pc=find(max(place_volt,[],2)>std(a(:,1:round(size(a,2)*0.3)),0,2)*2);
%calculate centroid of place field
%[a c]=sort(place_spike(put_pc,:)*[1:place_bin]'./sum(place_spike(put_pc,:),2),'ascend'); %calculate centroid
%calculate maximum of place field
[a arg]=max(place_volt(put_pc,:),[],2);
[a c]=sort(arg);
order=[c ;setdiff([1:size(volt,1)],put_pc)'];
if im_on
imagesc(place_volt(order,:)); axis tight; hold on;
line(repmat(ceil(reward_pos/(lap_dist/place_bin)),1,2),[0 length(order)+0.5],'color','r');
end
end