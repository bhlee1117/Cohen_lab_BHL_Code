function show_traces_align_Position(traces_F,traces_sp,pos,Window,Ard_data,noi)

scale=-15;
vel_thresh=0;

lap_length=max(Ard_data(:,2));
Ard_data(end-1:end,2:4)=repmat(Ard_data(end-2,2:4),2,1);
lap_end=[0; find(abs(Ard_data(2:end,2)-Ard_data(1:end-1,2))>6500); size(Ard_data,1)];
laps=[lap_end(1:end-1)+1 lap_end(2:end)];

for l=1:size(laps,1)
lap_trace(laps(l,1):laps(l,2))=l; 
cum_trace(laps(l,1):laps(l,2))=Ard_data(laps(l,1):laps(l,2),2)+lap_length*l; end

cum_trace=movmean(cum_trace,6);
vel_trace=(cum_trace(2:end)-cum_trace(1:end-1))/1.25*1000;
vel_trace(end+1)=vel_trace(end);

t=[1:size(traces_F,2)];

figure;
ax1=[];
cmap=turbo(size(laps,1));

% find the frames that passing the position
b=bwlabel(Ard_data(:,2)>pos);
pos_reach_time=NaN(size(laps,1),1);

for l=1:size(laps,1)
[dist tmp]=min(abs(Ard_data(laps(l,1)+25:laps(l,2)-25,2)-pos));
tmp=tmp+laps(l,1)-1+25;

if dist>10 %not nearby
pos_reach_time(l)=NaN;    
else
pos_reach_time(l)=tmp(1);
end
end


for n=noi
      nexttile([1 1])
    for l=1:size(laps,1)
        if ~isnan(pos_reach_time(l))
    try
        plot_tr=traces_F(n,pos_reach_time(l)+Window);
    plot_sp=traces_sp(n,pos_reach_time(l)+Window);
    plot(Window,plot_tr+l*scale,'color',cmap(l,:))
    hold all
    plot(Window(find(plot_sp)),plot_tr(find(plot_sp))+l*scale,'r.')
    plot(Window,zeros(length(Window),1)+l*scale,'color',[0.7 0.7 0.7])
    end
        end
    end
end
line([0 0],[0 scale*size(laps,1)],'color',[0 0.2 1],'linewidth',2)
axis tight
xlabel('Frame')
ylabel('Laps')
set(gca,'ytick',[size(laps,1):-1:1]*scale,'yticklabel',[size(laps,1):-1:1])
end

