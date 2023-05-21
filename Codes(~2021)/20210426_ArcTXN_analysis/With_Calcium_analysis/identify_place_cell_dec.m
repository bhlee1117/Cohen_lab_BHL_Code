function [place_fr field is_pc running]=identify_place_cell_dec(C_df,bin,dec,position,peak_th,in_out_th)

frm_end=min(size(position,1),size(C_df,2));
min_track=-2000;
% calculate Mean Cdf(x)
place_cdf=[position(1:frm_end,2) C_df(:,1:frm_end)'];
place_cdf(1,1)=min_track;
% find the region the speed exceeds 7 cm/s and track longer than 40 cm

step=place_cdf(2:end,1)-place_cdf(1:end-1,1);
l=find(abs(step)>abs(min_track));
step(l,1)=(step(l-1,1)+step(l+1,1))/2;
cum_track=imgaussfilt([0;cumsum(step)],5);
runn=bwlabel(step*30>15);
running=zeros(frm_end,1);
for i=1:max(runn)
    seg=find(runn==i);
    if cum_track(seg(end,1))-cum_track(seg(1,1))>50
        running(seg,1)=1;
    end
end
place_cdf_thr=place_cdf;
place_cdf_thr(~running,1:2)=NaN;

bin_place=round((place_cdf_thr(:,1)-min_track)/bin);
clear place_fr conditions
place_fr=zeros(max(round((place_cdf(:,1)-min_track)/bin)),1);
for i=1:max(bin_place) %place
    place_fr(i,1)=mean(place_cdf_thr(find(bin_place==i),2),'omitnan');
end
place_fr=movmean(place_fr,3);
place_fr=place_fr-min(place_fr);

%s=sort(place_fr,'ascend');
%%%% Check the conditions %%%%
%baseline=mean(s(1:20,1));
field=zeros(size(place_fr,1),1);
is_pc=0;
if ~isempty(Peak)
pot_pf=bwlabel(place_fr>median(place_fr));
for i=1:max(pot_pf)  % potential place fields
    seg=find(pot_pf==i);
    nonseg=setdiff([1:1:size(place_fr,1)]',seg);
    if size(seg,1)*bin>200 % 1st condition
        conditions(1,i)=1; else conditions(1,i)=0;
    end
    if ~isempty(find(place_fr(seg,1)>mean(place_cdf_thr(:,2),'omitnan')*0.15)) % 2nd condition
        conditions(2,i)=1; else conditions(2,i)=0;
    end
    n_Place_field=[];
    stay_field=bwlabel(place_cdf(:,1)-min_track>seg(1,1)*bin & place_cdf(:,1)-min_track<seg(end,1)*bin); %when does the mouse stayed in the potential field
    for fd=1:max(stay_field)
        nth_stay=find(stay_field==fd);
        n_Place_field(fd,1)=size(nth_stay,1)*1/30;
        if ~isempty(find(Peak(:,3)>nth_stay(1,1)*1/30 & Peak(:,3)<nth_stay(end,1)*1/30)) % if the peak is in the n th place field.
            n_Place_field(fd,2)=1;  % Yes there was a peak.
        else
            n_Place_field(fd,2)=0;
        end
    end
    if ~isempty(n_Place_field)
        %sum(n_Place_field(find(n_Place_field(:,2)==1),1))/sum(n_Place_field(:,1))
    if sum(n_Place_field(find(n_Place_field(:,2)==1),1))/sum(n_Place_field(:,1))>peak_th % 4th condition 20% of time there is a peak
        conditions(3,i)=1; else conditions(3,i)=0;  
    end
    else
        conditions(3,i)=0;  
    end
    
    if sum(conditions(:,i))==3
    else
     pot_pf(seg,1)=0;  
      field=zeros(size(place_fr,1),1);
      is_pc=0;
    end
end
g=1;
pp=unique(pot_pf);
if size(pp,1)>1
for i=pp(2:end)'
     seg=find(pot_pf==i);
     nonseg=find(pot_pf==0 & place_fr~=0);
    if mean(place_fr(seg,1))>mean(place_fr(nonseg,1))*in_out_th && ~isempty(nonseg) % 3rd condition
        conditions(4,i)=1;
        field=zeros(size(place_fr,1),1);
        field(seg)=g;
        is_pc=1;
        g=g+1;
    else conditions(4,i)=0;
        field=zeros(size(place_fr,1),1);
        is_pc=0;
    end
    
end
end
else
    is_pc=0;
    field=zeros(size(place_fr,1),1);
end
end