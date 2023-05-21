function plot_spike(spike)
for i=1:size(spike,1)
    plot(find(spike(i,:)),ones(sum((spike(i,:))),1)*i,'k.','markersize',4)
hold all
end
end