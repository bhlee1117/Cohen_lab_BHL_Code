function [x_out,y_out] = static_image_seg(im,thres)
if ~exist('thres','var')
    thres = 50;
end

%%
im_bin = double(im);

% gSig = 7; 
% gSiz = 7; 
% psf = fspecial('gaussian', round(gSiz), gSig);
% ind_nonzero = (psf(:)>=max(psf(:,1)));
% psf = psf-mean(psf(ind_nonzero));
% psf(~ind_nonzero) = 0;   % only use pixels within the center disk
% im_bin = imfilter(im_bin,psf,'same');
% % 
% im_bin([1:3 end-3:end],:) = 0;
% im_bin(:,[1:3 end-3:end]) = 0;
%%
im_f = im_high_pass_filt(im_bin,20);

gSig = 7; 
gSiz = 7; 
psf = fspecial('gaussian', round(gSiz), gSig);
ind_nonzero = (psf(:)>=max(psf(:,1)));
psf = psf-mean(psf(ind_nonzero));
psf(~ind_nonzero) = 0;   % only use pixels within the center disk
im_f = imfilter(im_f,psf,'same');

im_f([1:3 end-3:end],:) = 0;
im_f(:,[1:3 end-3:end]) = 0;

upsampl = 1;
im_fb = imresize(imgaussfilt(im_f,1.5),upsampl);
% im_fb = imresize(im_f,upsampl);
% im_fb = imgaussfilt(im_fb,3);

L = watershed(-im_fb); 
L_bd = bwboundaries((L>0));
%%
stats = regionprops(L, im_fb, 'Area', 'PixelValues');
im_med = median(im_fb(:));
areas = [stats(:).Area];

pxl_max = arrayfun(@(x) max(x.PixelValues),stats);

pxl_bg_cnt = arrayfun(@(x) length(find(x.PixelValues>im_med)),stats);


x = pxl_bg_cnt;
y = pxl_max;

figure(50);clf
subplot(1,6,1)
plot(x,y,'.')
%%
goodROI = (y>thres);
sum(goodROI)

L2 = L;  % Make a new label matrix with only the good ROIs
for j = 1:max(L(:))
    if goodROI(j);
        L2(L2(:) == j) = sum(goodROI(1:j));
    else
        L2(L2(:) == j)=0;
    end;
end;
% pxl = regionprops(L2,im_bin,'Pixelvalues');
%  
% ctr_val = arrayfun(@(x) max(x.PixelValues),pxl);
% [ctr_r,ctr_c] = ind2sub(size(im_bin),find(ismember(im_bin(:),ctr_val)));

pxl = regionprops(L2,im_fb,'Pixelvalues');
 
ctr_val = arrayfun(@(x) max(x.PixelValues),pxl);
[ctr_r,ctr_c] = ind2sub(size(im_fb),find(ismember(im_fb(:),ctr_val)));

im_fb = imresize(im_fb,1/upsampl);
subplot(1,6,2:6)
% imshow2(im_f,[prctile(im_f(:),1) prctile(im_f(:),99.5)]) 
% imshow2(im_bin,[prctile(im_bin(:),1) prctile(im_bin(:),99.9)]) 
% imshow2(im,[prctile(im(:),1) prctile(im(:),99.95)]) 
im_show = im_high_pass_filt(double(im),15);
imshow2(im_show,[prctile(im_show(:),1) prctile(im_show(:),99.95)]) 
hold on
scatter(ctr_c/upsampl,ctr_r/upsampl,'r*')

x_out = ctr_c;
y_out = ctr_r;

%%
function im_f = im_high_pass_filt(im,bw)

fft_im = fftshift(fft2(im));
band_width = bw;
f_ker = ones(size(im));
f_ker(1:band_width,1:band_width)=0;
f_ker(end-band_width:end,1:band_width)=0; 
f_ker(1:band_width,end-band_width:end)=0;
f_ker(end-band_width:end,end-band_width:end)=0;
f_ker = fftshift(f_ker);
im_f = abs(ifft2(ifftshift(fft_im.*imgaussfilt(f_ker,10)))).*im;

end
end