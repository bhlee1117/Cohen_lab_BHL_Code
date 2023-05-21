%%
clear
load('NF_map_.mat')
[fnm_ca pth_ca]=uigetfile('.tif','Load the registered & cropped calcium image');
fnm_ca={fnm_ca};
imginf=imfinfo([pth_ca fnm_ca{1,1}]);
numstack=size(imginf,1);
for i=1:numstack
    im_ca_crop(:,:,i)=imread([pth_ca fnm_ca{1,1}],i);
end
%% Set the ROI
clear imexp
cell_plane=NF_map.cell;
sz=200;
    im_ca_crop_ave=mean(im_ca_crop,3);
      imexp=zeros(size(im_ca_crop_ave,1)+2*floor(sz)/2,size(im_ca_crop_ave,2)+2*floor(sz/2),numstack);
      for i=1:numstack
    imexp(sz/2+1:sz/2+size(im_ca_crop_ave,1),sz/2+1:sz/2+size(im_ca_crop_ave,2),i)=im_ca_crop(:,:,i);
      end
      imexp_ave=mean(imexp,3);
    for k=1:size(cell_plane,2)
        cell_plane_shift(:,k)=cell_plane(:,k)+sz/2;
    end
    %%
          
for j=find(NF_map.list(:,4)~=3 & NF_map.list(:,5)~=3)'
    

clear crop_im
crop_im=imexp(round(cell_plane_shift(j,2))-sz/2:round(cell_plane_shift(j,2))+sz/2,round(cell_plane_shift(j,1))-sz/2:round(cell_plane_shift(j,1))+sz/2,:);
crop_im_std=std(crop_im,0,3);
crop_im_ave=mean(crop_im,3);
figure(1)
% imshow(mont,[0 max(max(mont))])

%mask(:,:,j)=crop_im>median(reshape(crop_im,size(crop_im,1)*size(crop_im,2),1));
subplot(2,2,1)
imagesc(crop_im_ave)
colormap('gray')
axis equal
hold all
plot(sz/2,sz/2,'ro')
subplot(2,2,2)
imagesc(max(crop_im,[],3))
hold all
colormap('gray')
axis equal
hold all
plot(sz/2,sz/2,'ro')
subplot(2,2,3)

df=sum((crop_im-median(crop_im,3))./median(crop_im,3)>0.2 & crop_im>900,3);
df_ave_RGB=zeros(size(crop_im,1),size(crop_im,2),3);
df_ave_RGB(:,:,1)=df/max(max(df));
df_ave_RGB(:,:,2)=crop_im_ave/max(max(crop_im_ave));

imagesc(df_ave_RGB)
subplot(2,2,4)
% sw=input(['( ' num2str(j) '/' num2str(size(cell_plane,1)) ' )th cell. Is the cell has TXN site?\n']);
% if sw==2
    %mask(:,:,j)=roipoly(crop_im/max(max(crop_im)));
    mask(:,:,j)=roipoly(df/max(max(df)));
    colormap('jet')
%end
end
  close(figure(1))
  
NF_map.mask=mask;
save([pth_ca,'NF_map_.mat'],'NF_map');

%%
ma_im=zeros(size(imexp,1),size(imexp,2),size(mask,3));
for j=1:size(mask,3)
    ma_im(round(cell_plane_shift(j,2))-sz/2:round(cell_plane_shift(j,2))+sz/2,round(cell_plane_shift(j,1))-sz/2:round(cell_plane_shift(j,1))+sz/2,j)=mask(:,:,j);
end
mask_im=max(ma_im,[],3);
imshowpair(imexp_ave, mask_im)
%%
clear imexpca
imexpca=zeros(size(im_ca_crop_ave,1)+2*floor(sz)/2,size(im_ca_crop_ave,2)+2*floor(sz/2));
for i=find(NF_map.list(:,4)~=3 & NF_map.list(:,5)~=3)' %¼¼Æ÷
  
    for j=1:size(im_ca_crop,3)
        clear tmp_im
    imexpca(sz/2+1:sz/2+size(im_ca_crop_ave,1),sz/2+1:sz/2+size(im_ca_crop_ave,2))=im_ca_crop(:,:,j);
    imexpca(ma_im(:,:,i)==0)=0;
    ca_trace(i,j)=sum(sum(imexpca));
    end
end

NF_map.Calcium=ca_trace;
save([pth_ca,'NF_map_.mat'],'NF_map');
%% Calculate dF/F

for i=1:size(ca_trace,1)
    dff(i,:)=(ca_trace(i,:)-median(ca_trace(i,:)))/median(ca_trace(i,:));
plot(dff)
hold all
end


%%

load([pth_ca,'NF_map_.mat'])

im