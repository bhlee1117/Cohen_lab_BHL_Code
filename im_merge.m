function M_im=im_merge(im,cmap)

if size(im,3)~=size(cmap,1)
    error('Check the size of images')
else
M_im=zeros(size(im,1),size(im,2),3);
for i=1:size(im,3)
im_rescale(:,:,i)=rescale(im(:,:,i));
end
im_rescale(isnan(im_rescale))=0;
M_im=squeeze(sum(im_rescale.*reshape(cmap,1,1,[],3),3));
end

end