function crop_im=imcrop_3d(img,roi)
crop_im=[];
for i=1:size(img,3)
crop_im(:,:,i)=imcrop(img(:,:,i),roi);
end
end