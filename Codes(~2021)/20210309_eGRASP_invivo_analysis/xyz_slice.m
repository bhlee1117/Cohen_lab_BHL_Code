function [c_image_xy]=xyz_slice(Im)

if length(size(Im))==4
    I=Im(:,:,:,3);
else
    I=Im;
end

cmap=jet(size(I,3));
c_image_xy=zeros(size(I,1),size(I,2),3);
for z=1:size(I,3)
    for c=1:3
    c_image_xy(:,:,c)=c_image_xy(:,:,c)+double(cmap(z,c)*adjust_dynamic_range(I(:,:,z)));
    end
end
c_image_xy=c_image_xy-min(c_image_xy(:));
c_image_xy=(c_image_xy)/median(max(sum(c_image_xy,3)))*3;


% cmap=jet(size(Im(:,:,:,3),1));
% c_image_yz=zeros(size(Im(:,:,:,3),2),size(Im(:,:,:,3),3),3);
% for x=1:size(Im(:,:,:,3),1)
%     for c=1:3
%     c_image_yz(:,:,c)=c_image_yz(:,:,c)+squeeze(double(cmap(x,c)*adjust_dynamic_range(Im(x,:,:,3))));
%     end
% end
% c_image_yz=c_image_yz-min(c_image_yz(:));
% c_image_yz=(c_image_yz)/median(max(sum(c_image_yz,3)))*2;
% 
% cmap=jet(size(Im(:,:,:,3),2));
% c_image_xz=zeros(size(Im(:,:,:,3),1),size(Im(:,:,:,3),3),3);
% for y=1:size(Im(:,:,:,3),2)
%     for c=1:3
%     c_image_xz(:,:,c)=c_image_xz(:,:,c)+squeeze(double(cmap(y,c)*adjust_dynamic_range(Im(:,y,:,3))));
%     end
% end
% c_image_xz=c_image_xz-min(c_image_xz(:));
% c_image_xz=(c_image_xz)/median(max(sum(c_image_xz,3)))*2;
end

function ad_im=adjust_dynamic_range(im)
im=double(im);
ad_im=(im-min(im(:)))/max(im(:));
end