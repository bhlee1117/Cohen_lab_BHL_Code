function [locmask] = findlocs(img,mindistance,maxpeaks)
%findrois   Find 2D peaks and calculate rois around them. 
%   Works well with low noise images of peaks with no background.
% 
%   2016 Vicente Parot
%   Cohen Lab - Harvard University
%
    for hpower = 4:-1:-1
        locmask = findpeaklocations(img,10^hpower,mindistance,maxpeaks);
        if nnz(locmask) >= maxpeaks
            break
        end
    end
end

function locmask = findpeaklocations(img,hthresholdfactor,mindistance,maxpeaks)
%% FindPeaks
% im3 = im2;
th0 = 0;
% im2(im2>prctile(im2(:),99)) = max(im2(:));
d1 = -conv2(img,[1 0 0; 0 -1 0; 0 0 0],'same');
d2 = -conv2(img,[0 1 0; 0 -1 0; 0 0 0],'same');
d3 = -conv2(img,[0 0 1; 0 -1 0; 0 0 0],'same');
d4 = -conv2(img,[0 0 0; 1 -1 0; 0 0 0],'same');
d6 = -conv2(img,[0 0 0; 0 -1 1; 0 0 0],'same');
d7 = -conv2(img,[0 0 0; 0 -1 0; 1 0 0],'same');
d8 = -conv2(img,[0 0 0; 0 -1 0; 0 1 0],'same');
d9 = -conv2(img,[0 0 0; 0 -1 0; 0 0 1],'same');
% p4 = d2>=th0 & d4>=th0 & d6>=th0 & d8>=th0; % will do 4-conn intead of 8-conn
% p4(im3==satlow) = false;
p8 = d1>=th0 & d2>=th0 & d3>=th0 & d4>=th0 & d6>=th0 & d7>=th0 & d8>=th0 & d9>=th0;
% p8(img < prctile(img(:),95)) = false;
% p8(img < max(img(:))/30) = false;
%% h-index Thresholding
s = sort(img(:),'descend');
x = linspace(0,max(img(:)),numel(img))';
[~, ind] = min(abs(s-x*hthresholdfactor));
% pr = ind/numel(img);
% [s(ind) prctile(img(:),100-pr*100) max(img(:))*(pr)]
p8(img < s(ind)) = false;
% nnz(p8)
% imshow(img > s(ind),[])
%% MinPeakDist
% nnz(p8)
% mindist = 7;
tp8 = p8;
it = 0;
tfp8 = find(tp8);
[tfp8r, tfp8c] = find(tp8);
[~, tix] = sort(img(tp8),'descend');
while true
    it = it + 1;
    if numel(tix) < it
        break
    end
    tooclose = find(sqrt((tfp8r(tix(it)) - tfp8r).^2 + (tfp8c(tix(it)) - tfp8c).^2) < mindistance);
    if numel(tooclose) > 1
        sd1 = setdiff(tooclose,tix(it));
        tp8(tfp8(sd1)) = false;
        sd2 = setdiff(1:numel(tfp8),sd1);
        tfp8 = tfp8(sd2);
        tfp8r = tfp8r(sd2); 
        tfp8c = tfp8c(sd2); 
        [~, tix] = sort(img(tp8),'descend');
    end
end
% imshow(tp8)
p8 = tp8;
%% NPeaks
% maxpeaks = 30;
fp8 = find(p8);
[~, ix] = sort(img(p8),'descend');
p8(fp8(ix(maxpeaks+1:end))) = false;
%% Intensity sorted bwlabel
% moved outside 
% fp8 = find(p8);
% [~, ix] = sort(img(p8),'descend');
% sp8 = zeros(size(p8));
% sp8(fp8(ix)) = 1:numel(fp8);
% % nnz(p8)
% % [fp8r, fp8c] = find(p8);
% %% calculate ROIs
% % roiradius = 3;
locmask = p8;
% dp8 = imdilate(sp8,strel('disk',roiradius));
% % dp8cc = bwconncomp(dp8,4);
% rois = regionprops(dp8,'ConvexHull');
% rois = {rois.ConvexHull};
% %% sort rois
% [~, ix] = sort(img(p8),'descend');


% clf
% imshow(im2,[])
% hold on
% plot(fp8c,fp8r,'o')
%%
% imshow(im3,[])
% imshow(dr1(2:end-1,2:end-1),[])
% imshow(dr2(2:end-1,2:end-1),[])
% imshow(dc1(2:end-1,2:end-1),[])
% imshow(dc2(2:end-1,2:end-1),[])
% imshow(p1(2:end-1,2:end-1),[])

%%
% [r, c] = find(p8);
% [v, cv] = voronoin([r c]);
% clf
% imshow(im2,[])
% hold on
% plot(v(:,2),v(:,1),'o')
%%
% forbid = unique([
%     find(v(:,1)<1)
%     find(v(:,2)<1)
%     find(v(:,1)>size(im3,1))
%     find(v(:,2)>size(im3,2))
%     ]);
% fp8 = find(p8);
% [fp8r, fp8c] = find(p8);
% [~, ix] = sort(im2(p8),'descend');
% select = ix(1:100);
% cs = cv;
% for it = 1:numel(forbid)
%     cs = cs(~cellfun(@(el)any(el==forbid(it)),cs)); % remove any objects in the border
% end
%%
% ic = ix(205);
% idc = cv{ic};
% clf
% imshow(im2,[])
% hold on
% plot(v(idc([1:end 1]),2),v(idc([1:end 1]),1),'-')
end
