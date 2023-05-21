function cir_im=point_image(data,sz,pt_sz)
cir_im=zeros(sz(1,1),sz(1,2));
data=round(data);
cir_im(createCirclesMask([sz(1,1) sz(1,2)],...
data(:,1:2,1),zeros(size(data,1),1)+pt_sz))=1;
% for i=1:size(data,1)
%     pt_im(data(i,2),data(i,1),data(i,3))=1;
% end
% h=ones(pt_sz,pt_sz,pt_sz)/pt_sz^3;
% pt_im=imgaussfilt3(pt_im,pt_sz);
% pt_im=pt_im>0.0005/pt_sz;
end