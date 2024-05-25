function h= errorbar_shade(x,y,error,cmap)
if nargin<4
    cmap=distinguishable_colors(1);
end
if size(error,1)>size(error,2)
    error=error';
end
if size(y,1)>size(y,2)
    y=y';
end
curve1 = y + error;
curve2 = y - error;
x2 = tovec([x, fliplr(x)]);
inBetween = tovec([curve1, fliplr(curve2)]);
fill(x2, inBetween, cmap , 'FaceAlpha', 0.3,'LineStyle','none');
hold on;
h= plot(x, y, 'color', cmap , 'LineWidth', 2);

end
