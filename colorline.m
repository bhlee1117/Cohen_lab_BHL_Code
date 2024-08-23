function colorline(X,Y,cmap,markers)
if nargin<4
markers=[];
end
Xdata=[X(1:end-1) X(2:end)];
Ydata=[Y(1:end-1) Y(2:end)];
switch cmap
    case 'turbo'
cc=turbo(size(Xdata,1));
    case 'jet'
cc=jet(size(Xdata,1));        
    case 'winter'
cc=winter(size(Xdata,1));        
    case 'summer'
cc=summer(size(Xdata,1));        
case 'hot'
cc=hot(size(Xdata,1));        
case 'spring'
cc=spring(size(Xdata,1));        
end
if isempty(markers)
l=plot(Xdata',Ydata');
else
l=plot(Xdata',Ydata','Marker',markers);    
end
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(cc,2))
end