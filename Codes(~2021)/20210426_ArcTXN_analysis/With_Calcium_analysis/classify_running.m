function running=classify_running(position)


frm_end=size(position,1);
min_track=-2000;
% calculate Mean Cdf(x)
place_cdf=[position(1:frm_end,2)];
place_cdf(1,1)=min_track;
% find the region the speed exceeds 7 cm/s and track longer than 40 cm

step=place_cdf(2:end,1)-place_cdf(1:end-1,1);
l=find(abs(step)>abs(min_track));
step(l,1)=(step(l-1,1)+step(l+1,1))/2;
cum_track=imgaussfilt([0;cumsum(step)],5);
runn=bwlabel(step*30>5); %mm
running=zeros(frm_end,1);
for i=1:max(runn)
    seg=find(runn==i);
    if cum_track(seg(end,1))-cum_track(seg(1,1))>5 %mm
        running(seg,1)=1;
    end
end


end