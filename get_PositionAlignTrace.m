function [PosTriggerTrace]=get_PositionAlignTrace(traces_F,pos,Window,VR_data)

ls=unique(VR_data(8,:));

for l=1:length(ls)
[dist frm]=min(abs(VR_data(5,VR_data(8,:)==ls(l))-pos));
if dist<1
pos_time(l)=frm+find(VR_data(8,:)==ls(l),1,'first')-1;
else
pos_time(l)=NaN;    
end
end

PosTriggerTrace=[];   
for n=1:size(traces_F,1)
    for l=1:length(ls)
        if ~isnan(pos_time(l))
    try
    plot_tr=traces_F(n,pos_time(l)+Window);
    PosTriggerTrace(l,:,n)=plot_tr;
    end
        end
    end
end
end
