function c_image=depth_im(Seg_im,dir,rat)

cmap=hsv(size(Seg_im,dir));  
switch dir
            case 1
    c_image=zeros(size(Seg_im,2),size(Seg_im,3),3);
    case 2
    c_image=zeros(size(Seg_im,1),size(Seg_im,3),3);
    case 3
 c_image=zeros(size(Seg_im,1),size(Seg_im,2),3);
end
        
for z=1:size(Seg_im,dir)
    for c=1:3
        switch dir
            case 1
                [a aa]=size(squeeze(double(cmap(z,c)*adjust_dynamic_range(Seg_im(z,:,:)))))
    c_image(:,:,c)=c_image(:,:,c)+squeeze(double(cmap(z,c)*adjust_dynamic_range(Seg_im(z,:,:))));
    case 2
    c_image(:,:,c)=c_image(:,:,c)+squeeze(double(cmap(z,c)*adjust_dynamic_range(Seg_im(:,z,:))));
    case 3
    c_image(:,:,c)=c_image(:,:,c)+double(cmap(z,c)*adjust_dynamic_range(Seg_im(:,:,z)));
        end
    end
end
c_image=c_image-min(c_image(:));
c_image=(c_image)/median(max(sum(c_image,3)))*rat;
end
function ad_im=adjust_dynamic_range(im)
im=double(im);
ad_im=(im-min(im(:)))/max(im(:));
end