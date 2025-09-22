function [dist, skeleton_path] = geodesic_distance(bwimg, point1, point2)
% Computes geodesic distance and path between two points on a binary mask
% using quasi-euclidean distance transform.

    point1 = round(point1);
    point2 = round(point2);

    bwimg = imdilate(bwimg, strel('disk', 1));
    %bwimg=imfill(bwimg,'holes');
    bwimg(point1(2), point1(1))=1;
    bwimg(point2(2), point2(1))=1;

    % Ensure both points are inside the main component
    bwss = bwlabel(bwimg);
    stats = regionprops(bwss, 'Area');
    [~, max_b] = max([stats.Area]);

    bwimg = grow_until_in_component(bwimg, point1, max_b);
    bwimg = grow_until_in_component(bwimg, point2, max_b);

    % Ensure points are in the same connected component
    bwss = bwlabel(bwimg); r=3;
    while bwss(point1(2), point1(1)) ~= bwss(point2(2), point2(1))
        bwimg = grow_point(bwimg, point2, r); % coarser steps to reduce loops
        bwss = bwlabel(bwimg);
        r=r+3;
    end

    % Final geodesic calculation
    bwimg = bwimg > 0;
    D1 = bwdistgeodesic(bwimg, point1(1), point1(2), 'quasi-euclidean');
    D2 = bwdistgeodesic(bwimg, point2(1), point2(2), 'quasi-euclidean');
    D = D1 + D2;
    D = round(D * 10) / 10;
    D(isnan(D)) = inf;

    skeleton_path = imregionalmin(D);
    skeleton_path = bwmorph(skeleton_path, 'fill', Inf);
    skeleton_path = bwmorph(skeleton_path, 'close', Inf);
    skeleton_path = bwmorph(skeleton_path, 'skel', Inf);

    dist = median(D(skeleton_path), 'omitnan');
end

function bwimg = grow_until_in_component(bwimg, point, target_label)
    % Grow the point until it joins the largest connected component
    r = 1;
    while true
        bwss = bwlabel(bwimg);
        if bwss(point(2), point(1)) == target_label
            break;
        end
        bwimg = grow_point(bwimg, point, r);
        r = r + 5;
    end
end

function bwimg = grow_point(bwimg, point, radius)
    % Dilate a point and merge it with the existing image
    addimg = false(size(bwimg));
    addimg(point(2), point(1)) = true;
    addimg = imdilate(addimg, strel('disk', radius));
    bwimg = bwimg | addimg;
end

% function [dist skeleton_path]=geodesic_distance(bwimg, point1, point2)
% 
% point1=round(point1); %reference
% point2=round(point2);
% bwss=bwlabel(bwimg);
% for b=1:max(bwss(:))
%     count(b)=sum(bwss(:)==b);
% end
% [~, max_b]=max(count);
% 
% if bwimg(point1(2),point1(1))==0
% 
% r=1;
%     while ~(bwss(point1(2),point1(1))==max_b)
% addimg=zeros(size(bwimg,1),size(bwimg,2));
% addimg(point1(2),point1(1))=1;
% SE = strel("disk",r);
% addimg=imdilate(addimg,SE);
% bwimg=max(cat(3,bwimg,addimg),[],3);
% bwss=bwlabel(bwimg);
% for b=1:max(bwss(:))
%     count(b)=sum(bwss(:)==b);
% end
% [~, max_b]=max(count);
% r=r+1;
%     end
% end
% 
% if bwimg(point2(2),point2(1))==0
% 
% r=1;
%     while ~(bwss(point2(2),point2(1))==max_b) %if the point is not in the boundary, expand the bwimg
% addimg=zeros(size(bwimg,1),size(bwimg,2));
% addimg(point2(2),point2(1))=1;
% SE = strel("disk",r);
% addimg=imdilate(addimg,SE);
% bwimg=max(cat(3,bwimg,addimg),[],3);
% bwss=bwlabel(bwimg);
% for b=1:max(bwss(:))
%     count(b)=sum(bwss(:)==b);
% end
% [~, max_b]=max(count);
% r=r+1;
%     end
% end
% 
% 
% if abs(bwss(point1(2),point1(1)) - bwss(point2(2),point2(1)))>0 %if the points are in different bwlabel
% 
% r=1;
%     while abs(bwss(point1(2),point1(1)) - bwss(point2(2),point2(1)))>0
% addimg=zeros(size(bwimg,1),size(bwimg,2));
% addimg(point2(2),point2(1))=1;
% SE = strel("disk",r);
% addimg=imdilate(addimg,SE);
% bwimg=max(cat(3,bwimg,addimg),[],3);
% bwss=bwlabel(bwimg);
% r=r+3;
%     end
% end
% bwimg=bwimg>0;
% D1 = bwdistgeodesic(bwimg, point1(1), point1(2), 'quasi-euclidean');
% D2 = bwdistgeodesic(bwimg, point2(1), point2(2), 'quasi-euclidean');
% D = D1 + D2;
% D = round(D * 10) / 10;
% D(isnan(D)) = inf;
% skeleton_path = imregionalmin(D);
% skeleton_path = bwmorph(skeleton_path, 'fill', Inf);
% skeleton_path = bwmorph(skeleton_path, 'close', Inf);
% skeleton_path = bwmorph(skeleton_path, 'skel', Inf);
% dist=median(D(skeleton_path),'omitnan');
% end