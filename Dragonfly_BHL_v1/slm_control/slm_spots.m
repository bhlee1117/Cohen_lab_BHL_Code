
x = (linspace(-1920/2,1920/2,1920))*9.2e-3;
y = (linspace(-1152/2,1152/2,1152))*9.2e-3;
[X,Y] = meshgrid(x,y);
theta = 45;
Gx = 1/(4/1920*max(X(:)*2));
Gy = 1/(4/1152*max(Y(:)*2));
lambda = .594e-3;
k=2*pi/lambda;
f=3100; %mm
% figure(101);clf;
% surf(abs(fftshift(fft2(fftshift(abs(input_intensity) .* exp(1i*pic_fft))))),'edgecolor','none');shading interp
% xlim([1200 1600]);ylim([200 400])

% figure(100);clf;imshow(mod(pic_fft_out,2*pi),[])
% figure(100);clf;imshow(mod(pic_fft_out,2*pi),[])
% pic_fft_out = mod(pic_fft_out+ 2*pi*Gx*X-2*pi*Gy*Y,2*pi); % blaze
pic_fft_out = mod(2*pi*Gx*X-2*pi*Gy*Y + k/2/f*(X.^2+Y.^2),2*pi); % blaze + lens
% pic_fft_out = mod((2*pi*Gx*X-2*pi*Gy*Y)/2 + (2*pi*Gx*.8*X-2*pi*Gy*.8*Y)/2 + k/2/f*(X.^2+Y.^2),2*pi); % 2x blaze + lens

% imshow(abs(fftshift(fft2(fftshift(abs(input_intensity_r) .* exp(1i*angle(A)))))),[]);
% % imshow(input_intensity,[])
% imshow(abs(fftshift(fft2(fftshift(abs(input_intensity) .* exp(1i*pic_fft_out))))),[]);


pic_fft_out = uint8(mat2gray(pic_fft_out)*255);
% figure(100);clf;imshow(pic_fft_uint8,[])
imwrite(pic_fft_out,'C:\Program Files\Meadowlark Optics\Blink 1920 HDMI\Image Files\single_spot.bmp')