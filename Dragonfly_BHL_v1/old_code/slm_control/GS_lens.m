x = (linspace(-1920/2,1920/2,1920))*9.2e-3;
y = (linspace(-1152/2,1152/2,1152))*9.2e-3;
[X,Y] = meshgrid(x,y);
theta = 45;
Gx = 1/(4/1920*max(X(:)*2));
Gy = 1/(4/1152*max(Y(:)*2));
lambda = .594e-3;
k=2*pi/lambda;
f=1000; %mm

pic_fft_out = mod(k/2/f*(X.^2+Y.^2),2*pi); % lens

% figure(101);clf
% imshow(abs(C),[])
% imshow(abs(fftshift(fft2(fftshift(abs(input_intensity_r) .* exp(1i*angle(A)))))),[]);
% % imshow(input_intensity,[])
% % imshow(abs(fftshift(fft2(fftshift(abs(input_intensity) .* exp(1i*pic_fft_out))))),[]);


pic_fft_out = uint8(mat2gray(pic_fft_out)*255);
% figure(100);clf;imshow(pic_fft_uint8,[])
imwrite(pic_fft_out,'C:\Program Files\Meadowlark Optics\Blink 1920 HDMI\Image Files\GS_lens.bmp')

pic_fft_out = uint8(mat2gray(ones(1152,1920))*255);
imwrite(pic_fft_out,'C:\Program Files\Meadowlark Optics\Blink 1920 HDMI\Image Files\white.bmp')