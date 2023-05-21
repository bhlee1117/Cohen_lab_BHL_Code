function inpoly=roi2mask(roi,r,c)
[x, y] = meshgrid(1:c, 1:r);
inpoly = inpolygon(x,y,roi(:,1),roi(:,2));
end