function ax1=plot_label(t,dat,label,cmap)
if nargin<4
    cmap=distinguishable_colors(1);
else
    label=double(label);
end

if ~(length(dat)==length(label))
error('check the size of data and label')
end

if max(label)==1

cmap=[cmap; 0.8 0 1]; % magenta
else
cmap=[cmap; jet(max(label))];
end

plot(t,dat,'color',cmap(1,:))
tr_label=NaN(max(label),size(dat,2));
for l=1:max(label)
tr_label(l,l==label)=dat(l==label);
end
hold all
if max(label)>0
lines=plot(t,tr_label);
arrayfun(@(l,c) set(l,'Color',c{:}),lines,num2cell(cmap(2:end,:),2))
end

end