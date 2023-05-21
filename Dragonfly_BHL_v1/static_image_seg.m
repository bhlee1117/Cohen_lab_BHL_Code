function [x_out,y_out] = static_image_seg(im,thres)
if ~exist('thres','var')
    thres = 50;
end
%%
im([1:3 end-3:end],:)=min(im(:));
[nr,nc] = size(im);
% bin = 2;
im_bin = double(im)+rand(size(im));%squeeze(sum(cat(3,im(1:bin:end,1:bin:end),im(2:bin:end,2:bin:end)),3));
% im_bin([1 2 end-2:end],:)=median(im_bin(:));

% figure(11);clf;
% imshow2(im,[])
% figure(12);clf
% imshow2(im_bin,[])
%%
[nr,nc] = size(im_bin);
fft_im = fftshift(fft2(im_bin));
band_width = 20;
f_ker = ones(size(im_bin));
f_ker(1:band_width,1:band_width)=0;
f_ker(end-band_width:end,1:band_width)=0;
f_ker(1:band_width,end-band_width:end)=0;
f_ker(end-band_width:end,end-band_width:end)=0;
f_ker = fftshift(f_ker);

im_f = abs(ifft2(ifftshift(fft_im.*imgaussfilt(f_ker,10)))).*im_bin;
% figure(13);clf
% imshow2(im_f,[]);%set(gca,'colorscale','log')
%%
im_f([1:2 end-2:end],:)=0;
im_fb = imgaussfilt(im_f,2);
% figure(14);clf
% imshow2(im_fb,[])
%%
% im_fb([1:6 end-6 end],:)=0;
L = watershed(-im_fb); 

L_bd = bwboundaries((L>0));
% figure(15);clf
% imshow2(im_fb,[]) 
% % set(gca,'colorscale','log')
% hold on
% cellfun(@(x) line(x(:,2),x(:,1),'color','r'),L_bd)
%%
% im_f([1 2 end-1 end],:) = 0;
stats = regionprops(L, im_fb, 'Area', 'PixelValues');
im_med = median(im_fb(:));
areas = [stats(:).Area];
pxl_max_avg_norm = arrayfun(@(x) max(x.PixelValues./im_med),stats);
pxl_max = arrayfun(@(x) max(x.PixelValues),stats);
pxl_min = arrayfun(@(x) min(x.PixelValues),stats);
pxl_bg_cnt = arrayfun(@(x) length(find(x.PixelValues>im_med)),stats);
pxl_kurt = arrayfun(@(x) length(kurtosis(x.PixelValues)),stats);
% pxl_diff = arrayfun(@(x) std(x.PixelValues),stats);


x = log(pxl_bg_cnt);
y = (pxl_max_avg_norm-1).*log(pxl_max);
% y = log(pxl_max);
% plot(x,y,'.')

[use_idx] = find(x>0);
b = robustfit(x(use_idx)+rand(size(x(use_idx)))*.1,y(use_idx),'bisquare',.4685,'off');

figure(50);clf
subplot(1,6,1)
% plot([prctile(im_f(:),95),prctile(im_f(:),95)],[0 30])
plot(x,y-(b*x),'.')

goodROI = (y-(b*x)>thres)& x>3 & (areas')<1000;
% goodROI = (y)>3.5;% & x<5;
sum(goodROI)

L2 = L;  % Make a new label matrix with only the good ROIs
for j = 1:max(L(:))
    if goodROI(j);
        L2(L2(:) == j) = sum(goodROI(1:j));
    else
        L2(L2(:) == j)=0;
    end;
end;
pxl = regionprops(L2,im_bin,'Pixelvalues');
 
ctr_val = arrayfun(@(x) max(x.PixelValues),pxl);
[ctr_r,ctr_c] = ind2sub(size(im_bin),find(ismember(im_bin(:),ctr_val)));

subplot(1,6,2:6)
% L2_bd = bwboundaries((L2>0));
imshow2(im_bin,[prctile(im_bin(:),1) prctile(im_bin(:),99.99)]) 
hold on
scatter(ctr_c,ctr_r,'*')

x_out = ctr_c;
y_out = ctr_r;