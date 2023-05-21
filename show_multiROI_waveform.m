function [dmd_time Blue]=show_multiROI_waveform(DAQ_waves,n,frm_rate,figure_on)

    cmap=distinguishable_colors(n);

a=find(DAQ_waves.amplitude(1,[1:500]));
interval=a(2)-a(1);

dmd_time=mod(cumsum(DAQ_waves.amplitude(3,:)),n);
dmd_time(dmd_time==0)=n;
Blue=[];
plot(dmd_time)
for i=2:n
    t=find(dmd_time==i);
    offset=mean(find(bwlabel(dmd_time==i)==1));
    if offset<interval/2
    offset=round(mean(find(bwlabel(dmd_time==i)==2)));
    end
    
    Blue=[Blue; DAQ_waves.amplitude(4,round([offset:interval:size(DAQ_waves.amplitude,2)]))];
    if figure_on
    plot(t*1e-5*(1/frm_rate*1000)/(1/frm_rate*1000+0.02),DAQ_waves.amplitude(4,t),'color',cmap(i,:))
    hold all
    end
end
if figure_on
axis tight
end
end

