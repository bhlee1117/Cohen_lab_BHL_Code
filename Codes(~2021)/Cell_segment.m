%function Cell_cent=Cell_segment(I,Cell_size,kernel,vac_ratio)
%% Control Tower

%if nargin<4
Cell_size=250; %pixels
kernel=3;
vac_ratio=0.6;
%end

%imagesc(I)
%%

imfilt=imgaussfilt(I, kernel);

[his bin]=imhist(imfilt);
g=1;
while bin(2,1)>10
[his bin]=imhist(imfilt,g*100);
g=g+1;
end
gg=1;
cum_his=cumsum(his);
while cum_his(gg,1)<vac_ratio*size(imfilt,1)*size(imfilt,2)
    gg=gg+1;
    threshold=bin(gg,1);
end

bi_I=imfilt>threshold;
D = bwdist(~bi_I);
D = -D;
D(~BW) = Inf;
wat_filt=watershed(D);
wat_filt=wat_filt==0;
im2=imfilt;
im2(wat_filt)=0;
imagesc(im2)
%%

wat_filt(~BW) = 0;
wat_filt= wat_filt==0;
wat_filt=bwareaopen(wat_filt,50);
%imfilt(wat_filt)=0;
%imagesc(imfilt)
%imagesc(uint8([wat_filt;imfilt;zeros]))
%colormap('gray')
%%

mask(wat_filt)=0;
imagesc(mask);
%%



bi_I=logical(heaviside(double(imfilt)-threshold));
D = bwdist(~bi_I);

% imfilt(~bi_I)=0;

D = -D;
D(~bi_I) = Inf;
%%
wat_I=watershed(D);
wat_I(~bi_I) = 0;
wat_bi=heaviside(double(wat_I)-0.5);
wat_I(logical(wat_bi))=1;
raw_conn = bwconncomp(wat_I);
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


