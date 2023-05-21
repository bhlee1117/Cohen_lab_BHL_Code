

    slm_mask = false(1152, 1920);
    slm_mask(1152/4+1:1152*3/4,1920/4+1:1920*3/4)=true;
    
% %-------------------------------------------------------
% calculate input beam or input intensity,
% in this example a gaussian beam was selected, 
% clear all; close all;


	x = linspace(-15,15,1920);
	y = linspace(-10,10,1152);
	[X,Y] = meshgrid(x,y);
	x0 = 0;     		% center
	y0 = 0;     		% center
	sigma = 2; 			% beam waist
	A = 1;      		% peak of the beam 
	res = ((X-x0).^2 + (Y-y0).^2)./(2*sigma^2);
	input_intensity = A  * exp(-res).*slm_mask;
    input_intensity = input_intensity/sqrt(sum(input_intensity.^2,[1 2]));
    power = sum(fft2(input_intensity).^2,[1 2]);
% 	surf(input_intensity);
% 	shading interp 
%---------------------------------------------------------
% 	Target=double(rgb2gray(imread('ThorLabsDog.png')));
% 	Target=1-double(Target);
	 Target = ones(1000,1000);
    m=2;
    n=4;
    Target = padarray(imresize(Target,floor([1152/n 1920]/m)),[ceil(1152*(1-1/m*1/n)/2) ceil(1920*(1-1/m)/2)],'both');
    Target = imresize(Target,[1152 1920]);
    Target = circshift(Target,[0 0]);
    Target = Target/sum(Target,[1 2])*power;
    
    s_mask = Target>0;
    f_mask = Target==0;
    m = 0.65;
    Target_sqrt = sqrt(Target);
    
    x = (linspace(-1920/2,1920/2,1920));
    y = (linspace(-1152/2,1152/2,1152));
    [X,Y] = meshgrid(x,y);
    
    R = 1.4e-3; a = .75;
    A = 4*R*(a*X.^2+(1-a)*Y.^2);
    A = mod(A,2*pi).*slm_mask;
    
    B = abs(input_intensity) .* exp(1i*A);
    C = fftshift(fft2(B));
    C_abs = abs(C);

    D = (m*abs(Target_sqrt.*s_mask) - (1-m)* C_abs.*f_mask).*exp(1i*angle(C));
    
     figure(122);subplot(2,1,1);imshowpair(A,abs(input_intensity));subplot(2,1,2);imshowpair(abs(C),abs(Target))
    
    %   
    %   D = abs(Target_sqrt) .* exp(1i*angle(C));
    A = ifft2(fftshift(D));
    

	error = [];
A = gpuArray(A);
Target = gpuArray(Target);
Target_sqrt = gpuArray(Target_sqrt);
input_intensity = gpuArray(input_intensity);
s_mask = gpuArray(s_mask);
f_mask = gpuArray(f_mask);
%%
	iteration_num = 200;
    tic;
for i=1:iteration_num
  B = abs(input_intensity) .* (A)./(abs(A));%
%   B = abs(input_intensity) .* exp(1i*angle(A));
  C = fftshift(fft2(B));

  D = (m*abs(Target_sqrt.*s_mask) - (1-m)* abs(C).*f_mask).*(C)./(abs(C));%
% D = (m*abs(Target_sqrt.*s_mask) - (1-m)* abs(C).*f_mask).* exp(1i*angle(C));
  A = ifft2(fftshift(D));
  error = [error; sum(abs(C) - abs(Target_sqrt),[1 2])];
%   i 
%   figure(10);clf;imshow(angle(D),[]);pause
end
toc;
A = gather(A);
C = gather(C);
figure(123);clf;plot(error)
 %%
x = (linspace(-1920/2,1920/2,1920))*9.2e-3;
y = (linspace(-1152/2,1152/2,1152))*9.2e-3;
[X,Y] = meshgrid(x,y);
theta = 45;
Gx = 1/(4/1920*max(X(:)*2));
Gy = 1/(4/1152*max(Y(:)*2));
lambda = .594e-3;
k=2*pi/lambda;
f=3050; %mm
% figure(101);clf;
% surf(abs(fftshift(fft2(fftshift(abs(input_intensity) .* exp(1i*pic_fft))))),'edgecolor','none');shading interp
% xlim([1200 1600]);ylim([200 400])
pic_fft_out = angle(A);
pic_fft_out = reshape(pic_fft_out(slm_mask),1152/2,1920/2);

pic_fft_out = imresize(pic_fft_out,[1152 1920],'nearest');
% figure(100);clf;imshow(mod(pic_fft_out,2*pi),[])
% figure(100);clf;imshow(mod(pic_fft_out,2*pi),[])
% pic_fft_out = mod(pic_fft_out+ 2*pi*Gx*X-2*pi*Gy*Y,2*pi); % blaze
pic_fft_out = mod(pic_fft_out+ 2*pi*Gx*X-2*pi*Gy*Y + k/2/f*(X.^2+Y.^2),2*pi); % blaze + lens


figure(101);clf
imshow(abs(C),[],'initialmagnification','fit')
% imshow(abs(fftshift(fft2(fftshift(abs(input_intensity_r) .* exp(1i*angle(A)))))),[]);
% % imshow(input_intensity,[])
% % imshow(abs(fftshift(fft2(fftshift(abs(input_intensity) .* exp(1i*pic_fft_out))))),[]);


pic_fft_out = uint8(mat2gray(pic_fft_out)*255);
% figure(100);clf;imshow(pic_fft_uint8,[])
imwrite(pic_fft_out,'C:\Program Files\Meadowlark Optics\Blink 1920 HDMI\Image Files\square_400_strip_MRAF_test.bmp')