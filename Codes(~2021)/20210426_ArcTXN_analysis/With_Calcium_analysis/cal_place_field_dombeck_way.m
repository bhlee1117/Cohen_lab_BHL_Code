function cal_place_field_dombeck_way(cal_trace,position,min_track,bin) % in  mm 
    position(1,1)=min_track;

    step=position(2:end,1)-position(1:end-1,1);
    l=find(abs(step)>abs(min_track));
    cont_position=position;
m=min(size(cal_trace,2),size(position,1));
    for lap=1:size(l,1)
        if lap<size(l,1)
    cont_position(l(lap,1)+1:l(lap+1,1),1)=cont_position(l(lap,1)+1:l(lap+1,1),1)+abs(cont_position(l(lap,1)+1)-cont_position(l(lap,1)));
    P{lap,1}=[position(l(lap,1)+1:l(lap+1,1),1) cal_trace(1,l(lap,1)+1:l(lap+1,1))'];
        else
      cont_position(l(lap,1)+1:m,1)=cont_position(l(lap,1)+1:m,1)+abs(cont_position(l(lap,1)+1)-cont_position(l(lap,1)));
    P{lap,1}=[position(l(lap,1)+1:m,1) cal_trace(1,l(lap,1)+1:m)'];      
    end
    end

    
long_period=find_long_per(cont_position,vel_thr);

bin_place=round((position-min_track)/bin);
clear place_fr place_fr_gauss
for i=1:max(bin_place) %place
    place_fr(i,1)=full(mean(cal_trace(1,find(bin_place==i))));
end
f=filter((1/3)*ones(1,3),1,place_fr);
place_fr=[place_fr(1,1);f(3:end,1); place_fr(end,1)];
f=sort(place_fr,'ascend');   baseline=mean(f(1:20,1));
poten_field=bwlabel(place_fr>(max(place_fr)-baseline)*0.25+baseline);
for i=1:max(poten_field)
    if size(find(poten_field==i),1)>200/bin && sum(place_fr(find(poten_field==i),1)>mean(place_fr)*1.1)>0 && ...
       mean(place_fr(find(poten_field==i),1))>mean(place_fr(find(poten_field~=i),1))*4 && 
    end
end

for i=1:size(place_cdf,2)-1 %cell
    place_fr_gauss(i,:)=imgaussfilt(place_fr(:,i),1)';
end
place_fr_gauss=(place_fr_gauss-min(place_fr_gauss,[],2))./(max(place_fr_gauss,[],2)-min(place_fr_gauss,[],2));
end
function long_period=find_long_per(cont_position,vel_thr,min_distance)
vel=(cont_position(2:end)-cont_position(1:end-1))/(1/30);
labeledArray = bwlabel((vel>vel_thr)==1);
g=1;
for i=1:max(labeledArray)
    if max(cont_position(find(labeledArray==i),1))-min(cont_position(find(labeledArray==i),1))>min_distance
     long_period(find(labeledArray==i),1)=g;
     g=g+1;
    end
end
end