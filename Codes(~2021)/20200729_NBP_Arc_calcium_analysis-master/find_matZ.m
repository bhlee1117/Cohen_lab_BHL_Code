function [angle_final argz roi]=find_matZ(MG,im_Arc_stack,z_interest)
figure1 = figure('InvertHardcopy','off','PaperUnits','centimeters',...
    'Color',[1 1 1],...
    'Renderer','painters','position',[100 100 500 700]);

marg=0;
manual_crop=false;
im_Arc=im_Arc_stack(:,:,z_interest);
MG=imresize(MG,[size(MG,1)*2 size(MG,2)*2]);

[crop_Arc roi_arc]=imcrop(im_Arc/30); crop_Arc=crop_Arc*30;

[ match_data, ~ ] = matchPattern( MG, crop_Arc, 0.01,2);  pattRows = size(crop_Arc,1);  pattCols = size(crop_Arc,2);
if size(match_data,1)>1
    [m argm]=max(match_data(:,3));
    match_data=match_data(argm,:);
end
match_data(:,1) = match_data(:,1)+ceil(pattRows/2);  match_data(:,2) = match_data(:,2)+ceil(pattCols/2);  candidatePos = match_data(:,[2,1]);
crop_G_mean=imcrop(MG,[match_data(:,2)-pattCols/2 match_data(:,1)-pattRows/2 pattCols-1 pattRows-1]);
angle=[-10:1:10];
cropz_Arc=crop_Arc;
close all
    for i=1:size(angle,2)
        for z=z_interest-20:2:z_interest+20
     c(i,z)=cal_corr(im_Arc_stack,MG,angle(1,i),z,roi_arc);
        end
    end
[maxc maxz]=find(c==max(c,[],'all'));
angle_final=angle(1,maxc); argz=maxz; %find the maximum correlated angle and z plane.

cropz_Arc=imcrop(im_Arc_stack(:,:,maxz),roi_arc);
[ match_data, ~ ] = matchPattern(imrotate(MG,angle_final,'crop'),cropz_Arc, 0.15,2);  
if size(match_data,1)>1
    [m argm]=max(match_data(:,3));
    match_data=match_data(argm,:);
end
pattRows = size(cropz_Arc,1);  pattCols = size(cropz_Arc,2);
match_data(:,1) = match_data(:,1)+ceil(pattRows/2);  match_data(:,2) = match_data(:,2)+ceil(pattCols/2);  candidatePos = match_data(:,[2,1]);
roi=[candidatePos-round(roi_arc(1,1:2))-round(roi_arc(1,3:4)/2) fliplr(size(im_Arc)-1)];

 MG_crop=imcrop(imrotate(MG,angle_final,'crop'),roi);

reg(:,:,1)=im_Arc_stack(:,:,maxz); 
reg(1:size(MG_crop,1),1:size(MG_crop,2),2)=MG_crop;
merge_image({reg(:,:,1),reg(:,:,2)},[1 0 0;1 1 1],{'Arc','Calcium'})
roi=roi/2;

end
function c=cal_corr(im_Arc_stack,MG,angle,z,roi_arc)
 cropz_Arc=imcrop(im_Arc_stack(:,:,z),roi_arc);
    rot_crop_G_mean=imrotate(MG,angle,'crop');
    [ match_data, ~ ] = matchPattern( rot_crop_G_mean, cropz_Arc, 0.01,2);
    c=max(match_data(:,3));
end