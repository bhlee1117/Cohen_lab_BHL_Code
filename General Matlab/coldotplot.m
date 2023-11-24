function coldotplot(x,y,s0,Ad)
% Color Scatter Plot for random data point visualization. It mimics a
% continuous 2D probability distribution.
% coldotplot(x,y,s0,Ad) creates at scatterplot with dots of sizes that
% correspond to their density in the swarm of points. The larger dots will
% also have a more "hot" color in the dense particle region.
% The data x and y are vectors of the same size, s0 is a parameter of local
% radii around each datapoint (defalult = 0.5). Ad is a visualization
% parameter for the area of the weighted dots (default = 1).
% Warning: May be slow for very large sizes of x and y.
%
% % Example:
% N=1000;
% x=randn(1,N);
% y=10*randn(1,N);
% s0=0.5;
% Ad=0.2;
% % Plot the density color plot
% coldotplot(x,y,s0,Ad)
%
% Created by: Per Sundqvist, ABB Corporate Research (Swe), 2010-12-07
radiix=s0*std(x);
radiiy=s0*std(y);

phi=linspace(0,2*pi,8);
xv0=radiix*cos(phi);
yv0=radiiy*sin(phi);

colval=zeros(size(x));
for j=1:length(x)
    xv=xv0+x(j);
    yv=yv0+y(j);
    in = inpolygon(x,y,xv,yv);
    colval(j)=sum(in);
end

[dummy,ixcol]=sort(colval);
figure, scatter(x(ixcol),y(ixcol),colval(ixcol)/Ad,colval(ixcol),'filled')
colormap hot
axis tight;
set(gca,'Color',[0.6784 0.9216 1],'FontSize',20)


