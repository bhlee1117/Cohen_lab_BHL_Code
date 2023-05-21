function im2=crx_reg(im)
clear drift drt xx

[t crop_pos]=imcrop(im(:,:,1),[0 1000]);
n=size(t);

for i=1:size(im,3)
    xx(:,:,i)=xcorr2(imcrop(mean(im(:,:,1),3),crop_pos),imcrop(im(:,:,i),crop_pos));
   [m m_arg]=max(reshape(xx(:,:,i),(2*n(1,1)-1)*(2*n(1,2)-1),1));
   drift(i,:)=[ceil(m_arg/(2*n(1,1)-1))-n(1,2) mod(m_arg,(2*n(1,1)-1))-n(1,1)];

end
drt=(-drift);
clear im2
sz=size(im);
im2=zeros(sz(1,1)+abs(min(drt(:,2)))*(1-heaviside(min(drt(:,2))))+abs(max(drt(:,2))),...
          sz(1,2)+abs(min(drt(:,1)))*(1-heaviside(min(drt(:,1))))+abs(max(drt(:,1))));
im3=zeros(sz(1,1)+abs(min(drt(:,2)))*(1-heaviside(min(drt(:,2))))+abs(max(drt(:,2))),...
          sz(1,2)+abs(min(drt(:,1)))*(1-heaviside(min(drt(:,1))))+abs(max(drt(:,1))));      
for j=1:size(im,3)
im2(abs(max(drt(:,2)))-drt(j,2)+1:abs(max(drt(:,2)))-drt(j,2)+sz(1,1),...
   abs(max(drt(:,1)))-drt(j,1)+1:abs(max(drt(:,1)))-drt(j,1)+sz(1,2),j)=im(:,:,j);
im3(abs(max(drt(:,2)))-drt(j,2)+1:abs(max(drt(:,2)))-drt(j,2)+sz(1,1),...
   abs(max(drt(:,1)))-drt(j,1)+1:abs(max(drt(:,1)))-drt(j,1)+sz(1,2),j)=1;
%     imwrite(uint16(imcrop(im,cro_pos)),[pth 'lip_reg_' char(fnm{1,i})],'Writemode','Append','Compression','none');
%imwrite(uint16(im2(:,:,j)),['H:\Image_data\OlympusTPM\20191205_jRCaMP1a\jRCAmp\' 'reg2_test.tif'],'Writemode','Append','Compression','none');
end

rng=[ceil([min(find(sum(im3==0,3)==0)) max(find(sum(im3==0,3)==0))]./(sz(1,1)+abs(min(drt(:,2)))*(1-heaviside(min(drt(:,2))))+abs(max(drt(:,2)))));...
      mod([min(find(sum(im3==0,3)==0)) max(find(sum(im3==0,3)==0))],(sz(1,1)+abs(min(drt(:,2)))*(1-heaviside(min(drt(:,2))))+abs(max(drt(:,2)))))];
im2=im2(rng(1,1):rng(1,2),rng(2,1):rng(2,2),:);
end