function [dist skeleton_path]=geodesic_distance(bwimg, point1, point2)

point1=round(point1); %reference
point2=round(point2);
bwss=bwlabel(bwimg);
for b=1:max(bwss(:))
    count(b)=sum(bwss(:)==b);
end
[~, max_b]=max(count);

if bwimg(point1(2),point1(1))==0

r=1;
    while ~(bwss(point1(2),point1(1))==max_b)
addimg=zeros(size(bwimg,1),size(bwimg,2));
addimg(point1(2),point1(1))=1;
SE = strel("disk",r);
addimg=imdilate(addimg,SE);
bwimg=max(cat(3,bwimg,addimg),[],3);
bwss=bwlabel(bwimg);
for b=1:max(bwss(:))
    count(b)=sum(bwss(:)==b);
end
[~, max_b]=max(count);
r=r+1;
    end
end

if bwimg(point2(2),point2(1))==0

r=1;
    while ~(bwss(point2(2),point2(1))==max_b)
addimg=zeros(size(bwimg,1),size(bwimg,2));
addimg(point2(2),point2(1))=1;
SE = strel("disk",r);
addimg=imdilate(addimg,SE);
bwimg=max(cat(3,bwimg,addimg),[],3);
bwss=bwlabel(bwimg);
for b=1:max(bwss(:))
    count(b)=sum(bwss(:)==b);
end
[~, max_b]=max(count);
r=r+1;
    end
end


if abs(bwss(point1(2),point1(1)) - bwss(point2(2),point2(1)))>0 %if the points are in different bwlabel

r=1;
    while abs(bwss(point1(2),point1(1)) - bwss(point2(2),point2(1)))>0
addimg=zeros(size(bwimg,1),size(bwimg,2));
addimg(point2(2),point2(1))=1;
SE = strel("disk",r);
addimg=imdilate(addimg,SE);
bwimg=max(cat(3,bwimg,addimg),[],3);
bwss=bwlabel(bwimg);
r=r+1;
    end
end
bwimg=bwimg>0;
D1 = bwdistgeodesic(bwimg, point1(1), point1(2), 'quasi-euclidean');
D2 = bwdistgeodesic(bwimg, point2(1), point2(2), 'quasi-euclidean');
D = D1 + D2;
D = round(D * 10) / 10;
D(isnan(D)) = inf;
skeleton_path = imregionalmin(D);
skeleton_path = bwmorph(skeleton_path, 'fill', Inf);
skeleton_path = bwmorph(skeleton_path, 'close', Inf);
skeleton_path = bwmorph(skeleton_path, 'skel', Inf);
dist=median(D(skeleton_path),'omitnan');
end