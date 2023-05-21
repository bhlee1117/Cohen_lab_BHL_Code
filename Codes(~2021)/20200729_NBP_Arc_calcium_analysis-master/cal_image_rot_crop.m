function Y=cal_image_rot_crop(Y,angle,roi)
Y=imrotate(Y,angle,'crop');
if ceil(roi(1,2))+floor(roi(1,4))>size(Y,1) || ceil(roi(1,1))+floor(roi(1,3))>size(Y,2) 
Y(end:ceil(roi(1,2))+floor(roi(1,4)),end:ceil(roi(1,1))+floor(roi(1,3)),:)=0;
end
Y=Y(ceil(roi(1,2)):ceil(roi(1,2))+floor(roi(1,4)),ceil(roi(1,1)):ceil(roi(1,1))+floor(roi(1,3)),:);
end