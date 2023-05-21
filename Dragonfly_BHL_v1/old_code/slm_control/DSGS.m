% Double constrain Gerchberg Saxton algorithm

%--------------------------------------------------------
% pseudo code of GS algorithm
% Gerchberg–Saxton Algorithm(Source, Target, Retrieved_Phase)
%  A = IFT(Target)
%  while error criterion is not satisfied
%    B = Amplitude(Source) * exp(i*Phase(A))
%    C = FT(B)
%    D = Amplitude(Target) * exp(i*Phase(C))
%    A = IFT(D)
%  end while
%  Retrieved_Phase = Phase(A)
%---------------------------------------------------------

% %-------------------------------------------------------
% calculate input beam or input intensity,
% in this example a gaussian beam was selected, 
% clear all; close all;


	x = linspace(-15,15,1920);
	y = linspace(-10,10,1152);
	[X,Y] = meshgrid(x,y);
	x0 = 0;     		% center
	y0 = 0;     		% center
	sigma = 1.5; 			% beam waist
	A = 1;      		% peak of the beam 
	res = ((X-x0).^2 + (Y-y0).^2)./(2*sigma^2);
	input_intensity = A  * exp(-res);
    input_intensity = input_intensity/sqrt(sum(input_intensity.^2,[1 2]));
    power = sum(fft2(input_intensity).^2,[1 2]);
% 	surf(input_intensity);
% 	shading interp 
%---------------------------------------------------------

% 	Target=rgb2gray(imread('ThorLabsDog.png'));
% 	Target=1-double(Target);
    Target = ones(1000,1000);
    m=10;
    n=2;
    Target = padarray(imresize(Target,floor([1152/n 1920]/m)),[ceil(1152*(1-1/m*1/n)/2) ceil(1920*(1-1/m)/2)],'both');
    Target = imresize(Target,[1152 1920]);
    Target = circshift(Target,[0 0]);
    Target = Target/sum(Target,[1 2])*power;
    
    s_mask = Target>0;
    f_mask = Target==0;
    
    slm_mask = false(1152, 1920);
    slm_mask(1152/4+1:1152*3/4,1920/4+1:1920*3/4)=true;

    Target_sqrt = sqrt(Target);
    
	A = fftshift(ifft2(fftshift(Target)));
    A = A.*slm_mask;
	error = [];

%%
	iteration_num = 100;
	%hologram = |objectWave + referenceWave|.^2
tic;

alpha = .5; beta = 0; gamma = 1;

for i=1:iteration_num
  B = abs(input_intensity) .* exp(1i*angle(A));
  C = fftshift(fft2(fftshift(B)));
  C_abs = abs(C);
  C_ang = angle(C);
  D = C;
  D(s_mask) = (2*alpha*abs(Target_sqrt(s_mask)) - beta* C_abs(s_mask))*exp(1i*0);
  D(f_mask) = gamma*C_abs(f_mask).*exp(1i*C_ang(f_mask));
%   
%   D = abs(Target_sqrt) .* exp(1i*angle(C));
  A = fftshift(ifft2(fftshift(D)));
  error = [error; sum(sum(abs(1.32*abs(C) - abs(Target))))];
  i
%   figure(10);clf;imshow(angle(D),[]);pause
end
	toc;
 %%
x = (linspace(-1920/2,1920/2,1920))*9.2e-3;
y = (linspace(-1152/2,1152/2,1152))*9.2e-3;
[X,Y] = meshgrid(x,y);
theta = 45;
Gx = 1/(4/1920*max(X(:)*2));
Gy = 1/(4/1152*max(Y(:)*2));
lambda = .594e-3;
k=2*pi/lambda;
f=1; %mm
pic_fft = mod(angle(A)+ 2*pi*Gx*X-2*pi*Gy*Y,2*pi);
% pic_fft = angle(A);
% pic_fft = mod(2*pi*Gx*X-2*pi*Gy*Y,2*pi).*0;
figure(444);clf;imshow(abs(fftshift(fft2(fftshift(abs(input_intensity) .* exp(1i*pic_fft))))),[]);
% figure(102);clf;imshow(abs(C),[])
% pic_fft = angle(A);
% pic_fft_uint8 = uint8(mat2gray(pic_fft)*255);
% figure(100);clf;imshow(pic_fft_uint8,[])
% imwrite(pic_fft_uint8,'C:\Program Files\Meadowlark Optics\Blink 1920 HDMI\Image Files\square_5s_GS.bmp')