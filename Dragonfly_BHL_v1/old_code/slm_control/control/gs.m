function phase = gs(Target,it)


%%
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


    
%     A = fftshift(ifft2(fftshift(Target)));
%     A = A.*slm_mask;
% clear A B C D    
% Target = register_pattern_test(1152,1920,[650 450]);
% target_gen
% Target = circshift(Target,470,2);
% Target = circshift(Target,-300,1);
Target = Target/sum(Target,[1 2])*power;
Target_sqrt = sqrt(Target);


    x = (linspace(-1920/2,1920/2,1920));
    y = (linspace(-1152/2,1152/2,1152));
    [X,Y] = meshgrid(x,y);
    R = 2e-3; a = .6; %b = .8; u = -pi/3*1.2; 
%     R = 2e-3; a = .6; b = 9.2e-3*2*pi*27.2; u = -pi/2; 
    A = 4*R*(a*X.^2+(1-a)*Y.^2);% + b*(X*cos(u)+Y*sin(u));
%     A = 4*R*(a*X.^2+(1-a)*Y.^2) + b*(X-Y);
    A = mod(A,2*pi).*slm_mask;
    
    B = abs(input_intensity) .* exp(1i*A);
    C = fftshift(fft2(fftshift(B)));
    D = abs(Target) .* exp(1i*angle(C));

%      figure(122);subplot(2,1,1);imshowpair(A,abs(input_intensity));subplot(2,1,2);imshowpair(abs(C),abs(Target))
    
    A = fftshift(ifft2(fftshift(D)));
    

%%
A = gpuArray(A);
Target = gpuArray(Target);
Target_sqrt = gpuArray(Target_sqrt);
input_intensity = gpuArray(input_intensity);
    error = [];
	iteration_num = it;
    tic;
for i=1:iteration_num
  B = abs(input_intensity) .* (A+eps)./(abs(A)+eps);%exp(1i*angle(A));
  C = fftshift(fft2(B));
  D = abs(Target_sqrt) .* (C+eps)./(abs(C)+eps);%exp(1i*angle(C));
  A = ifft2(fftshift(D));
  error = [error; sum(sum(abs(1.32*abs(C) - abs(Target))))];
%   i
end
	toc;
    phase_uc = angle(gather(A));

    phase = phase_uc;
    %% wavefront correction
    x = (linspace(-1920/2,1920/2,1920))*9.2e-3;
y = (linspace(-1152/2,1152/2,1152))*9.2e-3;
[X,Y] = meshgrid(x,y);
theta = 45;
Gx = 1/(4/1920*max(X(:)*2));
Gy = 1/(4/1152*max(Y(:)*2));
lambda = .594e-3;
k=2*pi/lambda;
f=3050; %mm

[xx, yy] = ndgrid(1:1152,1:1920);
xx = (xx-1920/2)/1920;
yy = (yy-1152/2)/1152;
rr = sqrt(xx.^2+yy.^2);
th = atan2(yy,xx);

zpolys = vm(cat(3,...
2*yy,... % tilt
2*xx,... % tip
sqrt(6)*rr.^2.*sin(2*th),... % oblique astigmatism
sqrt(3)*(2*rr.^2-1),...      % defocus
sqrt(6)*rr.^2.*cos(2*th),... % vertical astigmatism
sqrt(8)*rr.^3.*sin(3*th),...        % vertical trefoil
sqrt(8)*(3*rr.^3-2*rr).*sin(th),... % vertical coma
sqrt(8)*(3*rr.^3-2*rr).*cos(th),... % horizontal coma
sqrt(8)*rr.^3.*cos(3*th),...        % oblique trefoil
sqrt(10)*rr.^4.*sin(4*th),...             % Oblique quadrafoil
sqrt(10)*(4*rr.^4-3*rr.^2).*sin(2*th),... % Oblique secondary astigmatism
sqrt(5)*(6*rr.^4-6*rr.^2+1),...           % Primary spherical
sqrt(10)*(4*rr.^4-3*rr.^2).*cos(2*th),... % Vertical secondary astigmatism
sqrt(10)*rr.^4.*cos(4*th)));              % Vertical quadrafoil

commonz = [
0 % Tilt (Y; vertical tilt)
0 % Tip (X; horizontal tilt)
0 % Oblique astigmatism
0 % Defocus
0 % Vertical astigmatism
0 % Vertical trefoil
0 % Vertical coma
-50 % Horizontal coma
0 % Oblique trefoil
0 % Oblique quadrafoil
0 % Oblique secondary astigmatism
0 % Primary spherical
0 % Vertical secondary astigmatism
0 % Vertical quadrafoil
];

% figure(101);clf;
% surf(abs(fftshift(fft2(fftshift(abs(input_intensity) .* exp(1i*pic_fft))))),'edgecolor','none');shading interp
% xlim([1200 1600]);ylim([200 400])
pic_fft_out = phase_uc;
pic_fft_out = reshape(pic_fft_out(slm_mask),1152/2,1920/2);

pic_fft_out = imresize(pic_fft_out,[1152 1920],'nearest');

% figure(100);imshow(mod(pic_fft_out,2*pi),[])
% figure(100);clf;imshow(mod(pic_fft_out,2*pi),[])
% pic_fft_out = mod(pic_fft_out+ 2*pi*Gx*X-2*pi*Gy*Y,2*pi); % blaze


% figure(101);clf
% imshow(abs(C),[])
% imshow(abs(fftshift(fft2(fftshift(abs(input_intensity_r) .* exp(1i*angle(A)))))),[]);
% % imshow(input_intensity,[])
% imshow(abs(fftshift(fft2(fftshift(abs(input_intensity) .* exp(1i*pic_fft_out))))),[]);

phase = mod(pic_fft_out+ 2*pi*Gx*X-2*pi*Gy*Y + k/2/f*(X.^2+Y.^2),2*pi); % blaze + lens
% phase = mod(pic_fft_out+ ...
%                 2*pi*Gx*X-2*pi*Gy*Y + ...
%                 k/2/f*(X.^2+Y.^2) +...
%                 squeeze(zpolys.data(:,:,9))*10 ...
%                 ,2*pi); % blaze + lens + coma

% phase = mod(pic_fft_out+ ...
%                 k/2/f*(X.^2+Y.^2) +...
%                 squeeze(zpolys.data(:,:,8))*35 ...
%                 ,2*pi); % lens + coma
% pic_fft_out = uint8(mat2gray(pic_fft_out)*255);
%     C = gather(C);
end