function show_traces_align_Position_wCS(traces_F,traces_CS,pos,Window,VR_data,loi)


figure;
%traces_F=rescale(traces_F);
scale=20;
ls=unique(VR_data(8,:));
if nargin<6
loi=[1:length(ls)];
end
cmap=[0.1 0.1 0.1];
for l=1:length(ls)
[dist frm]=min(abs(VR_data(5,VR_data(8,:)==ls(l))-pos));
if dist<1
pos_time(l)=frm+find(VR_data(8,:)==ls(l),1,'first')-1;
else
pos_time(l)=NaN;    
end
end

      nexttile([1 1])
      g=1;
    for l=loi
        if ~isnan(pos_time(l))
    try
    plot_tr=traces_F(1,pos_time(l)+Window);
    CS=traces_CS(1,pos_time(l)+Window);
    plot_CS=plot_tr; plot_CS(CS==0)=NaN;
    plot(Window,plot_tr-g*scale,'color',cmap)

    hold all
    plot(Window,plot_CS-g*scale,'color',[0.8 0 0.2])
    plot(Window,zeros(length(Window),1)-g*scale,'color',[0.7 0.7 0.7])
    g=g+1;
    end
        end
    end

%line([0 0],[-length(ls) 0]*scale,'color',[0 0.2 1],'linewidth',2)
axis tight
xlabel('Time (ms)')
ylabel('Laps')
set(gca,'ytick',[-g:3:-1]*scale,'yticklabel',loi([length(loi):-3:1]))
end
