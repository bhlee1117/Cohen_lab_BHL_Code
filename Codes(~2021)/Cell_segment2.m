function Cell_cent=Cell_segment2(I,Cell_size,kernel,vac_ratio)
%% Control Tower
if nargin<4
Cell_size=200; %pixels
kernel=3;
vac_ratio=0.6;
end
%I=imread('20180606_Arc_oe1_HC_.tif',244);

%%

imfilt=imgaussfilt(I, kernel);
[BW,maskedImage] = segmentImage_BHL(imfilt,vac_ratio);
%%
% maskedImage2=imhmin(maskedImage,1);
D = bwdist(~BW);
D = -D;
D(~BW) = Inf;

%%
wat_filt=watershed(D,8);
wat_filt(~BW)=0;
%imagesc(wat_filt)

%%
% wat_bi=heaviside(double(wat_I)-0.5);
% wat_I(logical(wat_bi))=1;
raw_conn = bwconncomp(wat_filt);
S = regionprops(raw_conn,'Centroid');
g=1;
% imagesc(imfilt)
% colormap('gray')
% hold all
for i=1:size(raw_conn.PixelIdxList,2)
if size(raw_conn.PixelIdxList{i},1)>Cell_size
    Cell_cent(g,1:2)=S(i).Centroid;
    g=g+1;
%     plot(S(i).Centroid(1,1),S(i).Centroid(1,2),'ro')
%     text(S(i).Centroid(1,1),S(i).Centroid(1,2)+10,num2str(g),'color','r',...
% 'fontweight','bold')
end
end


