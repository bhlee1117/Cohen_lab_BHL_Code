function pos_time=find_PosTime(VR_data,pos)

ls=unique(VR_data(8,:));
loi=[1:length(ls)];

pos_time=zeros(length(ls),length(pos));
for p=1:length(pos)
for l=1:length(ls)
[dist frm]=min(abs(VR_data(5,VR_data(8,:)==ls(l))-pos(p)));
if dist<1
pos_time(l,p)=frm+find(VR_data(8,:)==ls(l),1,'first')-1;
else
pos_time(l,p)=NaN;    
end
end
end