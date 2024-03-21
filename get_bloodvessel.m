function ROIpoly = get_bloodvessel(mov_res_filt,tau,avgImg_mask,dilate_size)


mov_lag=mov_res_filt(:,:,1:end-tau).*mov_res_filt(:,:,1+tau:end);

bvImg=mean(mov_lag,3);
bvImg(avgImg_mask)=0;
f1 = figure(33); clf;
imshow2(bvImg,[0 prctile(bvImg(:),90)])

g=1; ROIpoly=[];
while g
    h = drawpolyline;
    if size(h.Position,1)==1 %no more ROI
        g=0;
    else

inD = 200;
pt = interparc(inD, h.Position(:,1), h.Position(:,2), 'spline');
        ROIpoly=[ROIpoly; {pt}];
    hold all
    h=[];
    plot(pt(:,1),pt(:,2))
    end
end
close(figure(33));

vesselmask=zeros(size(mov_res_filt,1),size(mov_res_filt,2));
vessel_points=round(cell2mat(ROIpoly));
vesselmask(sub2ind(size(vesselmask),vessel_points(:,2),vessel_points(:,1)))=1;

se = strel('sphere',dilate_size);
vesselmask=imdilate(vesselmask,se);

figure(33); clf;
imshow2(imfuse(mean(mov_lag,3),vesselmask),[])
end