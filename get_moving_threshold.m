function thres=get_moving_threshold(volt_hp,ratio,bin)
seg=[1:bin:size(volt_hp,2)];
seg=[seg(1:end-1)' seg(2:end)'];

for i=1:size(seg,1)
th(i,:)=[mean(seg(i,:)) get_threshold(volt_hp(:,seg(i,1):seg(i,2)),ratio)'];
end

for n=2:size(th,2)
thres(n-1,:)=interp1(th(:,1),th(:,n),[1:size(volt_hp,2)],'linear');
end
end