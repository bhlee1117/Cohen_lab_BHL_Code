function mask=roi2mask(roi,r,c)
mask=zeros(r,c);
[x, y] = meshgrid(1:c, 1:r);
for i=1:length(roi) 
inpoly = inpolygon(x,y,roi{i}(:,1),roi{i}(:,2));
mask=mask+inpoly;
end
mask=mask>0;
end