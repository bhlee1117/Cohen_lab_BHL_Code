function [ROIpoly ROImask]=get_ROI(avgImg,appendmask,val)
if nargin<2
    appendmask=[];
    val=0;
end
figure(333); clf;
    
    if ~isempty(appendmask)
    appendcontour=bwboundaries(max(appendmask,[],3));
    appendcontour=cell2mat(appendcontour);
    avgImg(max(appendmask,[],3)>0)=0;
    imshow2(avgImg,[]); hold all
    plot(appendcontour(:,2),appendcontour(:,1),'w.')
    else
        imshow2(avgImg,[]); hold all
    end

    g=1; ROIpoly=[];
    while g
        h = drawpolygon('Color','r');
        if size(h.Position,1)==1 %no more ROI
            g=0;
        else
            ROIpoly=[ROIpoly; {h.Position}];
            hold all
            plot(h.Position(:,1),h.Position(:,2))
        end
    end
sz=size(avgImg);
for i=1:length(ROIpoly)

    ROImask(:,:,i)=poly2mask(ROIpoly{i}(:,1),ROIpoly{i}(:,2),sz(1),sz(2));
end
ROImask=cat(3,appendmask,ROImask);
close(figure(333));
end