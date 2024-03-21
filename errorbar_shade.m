function errorbar_shade(x,y,error,cmap)

curve1 = y + error;
curve2 = y - error;
x2 = [x, fliplr(x)];
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, cmap , 'FaceAlpha', 0.1,'LineStyle','none');
hold on;
plot(x, y, 'color', cmap , 'LineWidth', 2);

end
