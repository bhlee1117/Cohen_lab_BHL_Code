f = 1; % meter focal length
d = 9.2e-6; % SLM pixel size
lambda = .594e-6;
n = 16;
slm_mask = true(1152,1920);
slm_mask([1:1152/4 round(1152*3/4):1152],:)=false;
slm_mask(:,[1:1920/4 1920*3/4:1920])=false;
%%
[x,y] = meshgrid(-1:1,-1:1);
x = reshape(x*10e-4,1,1,[]);
y = reshape(y*10e-4,1,1,[]);
z = reshape(linspace(-50,50,length(x))*0,1,1,[]);
%%
[xj,yj] = meshgrid(((1:1920)-1920/2)*d,((1:1152)-1152/2)*d);


Delta = pi*z/lambda/f^2.*(xj.^2+yj.^2)+2*pi/f/lambda*(xj.*x+yj.*y);
%% ---- wgs_spots test ----
%%
d = 9.2e-6*2; % SLM pixel size
test_im = zeros(1152,1920);

pts = click_im(test_im,'pts');
%%

z_offset = -61e+1;

x_spt = pts(:,1);
y_spt = pts(:,2);
z_spt = zeros(size(x_spt))+z_offset;
% z_offset = -5.4e-1;


rot = 60;
xyz = [x_spt y_spt z_spt]';
xyz_offset = mean(xyz,2);
xyz = [
     1 0 0 %0 1 1/3 %
    0 cosd(rot)  sind(rot)
     0  sind(rot)  cosd(rot)%0 -1/20 1 % 
    ]*(xyz-xyz_offset) + xyz_offset;

tic
% phase = wgs_spots([reshape(x,[],1),reshape(y,[],1),reshape(z,[],1)],30);
phase = wgs_spots(xyz',30,1);
toc
pause(.5)

hSLMApp.current_mask = phase;
hSLM.project(phase)

%% find optimal x-y plane rotation offset
z_offset = -5.4e-1;

x_spt = pts(:,1)*d;
y_spt = pts(:,2)*d;
z_spt = zeros(size(x_spt))+z_offset;
% z_offset = -5.4e-1;



im_all = [];
nrep =10;
for jj = 1:nrep
rot_all = -70:1:-50;
for ii = 1:length(rot_all)
rot = rot_all(ii);
xyz = [x_spt y_spt z_spt]';
xyz_offset = mean(xyz,2);
xyz = [
     1 0 0 %0 1 1/3 %
    0 cosd(rot)  sind(rot)
     0  sind(rot)  cosd(rot)%0 -1/20 1 % 
    ]*(xyz-xyz_offset) + xyz_offset;

tic
% phase = wgs_spots([reshape(x,[],1),reshape(y,[],1),reshape(z,[],1)],30);
phase = wgs_spots(xyz',30,1);
toc
pause(.5)

hSLMApp.current_mask = phase;
hSLM.project(phase)
pause(.5)
live_im=hDragonflyApp.camApp.snapOneFrame;
im_all(:,:,ii,jj) = live_im; 
end
end
im_all = squeeze(mean(im_all,4));
%% ---- wgs_spots test end ----
%%
figure(99);clf;imshow(abs(fftshift(fft2(exp(1i*phase)))),[]);colormap(jet)
%%

V = sum(exp(1i*2*pi*(2*f+z))/lambda/1i.*d^2/lambda/f.*sum(exp(-1i*Delta),3),3).*slm_mask;
%%

%%
figure;imshow(abs(fftshift(fft2(V))),[]);colormap(jet)

% figure;imshow(imgaussfilt(abs(fftshift(fft2(V))),1),[]);colormap(jet)
%%
phase = angle(V);
%%
% [x,y] = meshgrid(-2:2,-2:2);
% x = x(:)*2e-3;
% y = y(:)*2e-3*0;
% z = zeros(size(x));
x = [-2;0;2]*2e-3;
y = zeros(3,1);
z = zeros(size(x));

slm_mask = false(1152, 1920);
slm_mask(1152/4+1:1152*3/4,1920/4+1:1920*3/4)=true;
%%
tic
phase = wgs_spots([x y z],30,true);
toc
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

phase = reshape(phase(slm_mask),1152/2,1920/2);

phase = imresize(phase,[1152 1920],'nearest');

% figure(100);clf;imshow(mod(pic_fft_out,2*pi),[])
% figure(100);clf;imshow(mod(pic_fft_out,2*pi),[])
% phase = mod(phase+ 2*pi*Gx*X-2*pi*Gy*Y,2*pi); % blaze
phase_out = mod(phase+ 2*pi*Gx*X-2*pi*Gy*Y + k/2/f*(X.^2+Y.^2),2*pi); % blaze + lens
hSLM.project(phase_out)
%%

x = (linspace(-1920/2,1920/2,1920))*9.2e-3;
y = (linspace(-1152/2,1152/2,1152))*9.2e-3;
[X,Y] = meshgrid(x,y);
theta = 45;
Gx = 1/(4/1920*max(X(:)*2));
Gy = 1/(4/1152*max(Y(:)*2));
lambda = .594e-3;
k=2*pi/lambda;


f=3200; %mm

%% search for f with max intens
fAll = 2900:50:4000;
im_all = [];
maxIntens = zeros(10,length(fAll));
for iRep = 1%:10

for ii = 1:length(fAll)
f = fAll(ii);
% figure(101);clf;
% surf(abs(fftshift(fft2(fftshift(abs(input_intensity) .* exp(1i*pic_fft))))),'edgecolor','none');shading interp
% xlim([1200 1600]);ylim([200 400])
pic_fft_out = hSLMApp.current_mask;
pic_fft_out = reshape(pic_fft_out(slm_mask),1152/2,1920/2);

pic_fft_out = imresize(pic_fft_out,[1152 1920],'nearest');
% figure(100);clf;imshow(mod(pic_fft_out,2*pi),[])
% figure(100);clf;imshow(mod(pic_fft_out,2*pi),[])
% pic_fft_out = mod(pic_fft_out+ 2*pi*Gx*X-2*pi*Gy*Y,2*pi); % blaze


% figure(101);clf
% imshow(abs(C),[])
% imshow(abs(fftshift(fft2(fftshift(abs(input_intensity_r) .* exp(1i*angle(A)))))),[]);
% % imshow(input_intensity,[])
% imshow(abs(fftshift(fft2(fftshift(abs(input_intensity) .* exp(1i*pic_fft_out))))),[]);

pic_fft_out = mod(pic_fft_out+ 2*pi*Gx*X-2*pi*Gy*Y + k/2/f*(X.^2+Y.^2),2*pi); % blaze + lens
% phase_c = mod(pic_fft_out+ ...
%                 2*pi*Gx*X-2*pi*Gy*Y + ...
%                 k/2/f*(X.^2+Y.^2) +...
%                 squeeze(zpolys.data(:,:,8))*35 ...
%                 ,2*pi); % blaze + lens + coma
% pic_fft_out = uint8(mat2gray(pic_fft_out)*255);

hSLM.project(pic_fft_out)

pause(.5)
live_im=hDragonflyApp.camApp.snapOneFrame;
im_all(:,:,ii) = live_im;
% maxIntens(iRep,ii) = max(live_im(:));
end
end
figure;moviesc(vm(im_all))
% figure;plot(fAll,mean(maxIntens,1))
%% search for optim aberration correction

slm_mask = false(1152, 1920);
slm_mask(1152/4+1:1152*3/4,1920/4+1:1920*3/4)=true;

x = (linspace(-1920/2,1920/2,1920))*9.2e-3;
y = (linspace(-1152/2,1152/2,1152))*9.2e-3;
[X,Y] = meshgrid(x,y);
theta = 45;
Gx = 1/(4/1920*max(X(:)*2));
Gy = 1/(4/1152*max(Y(:)*2));
lambda = .594e-3;
k=2*pi/lambda;


f=3200; %mm

[xx, yy] = ndgrid(1:1152,1:1920);
xx = (xx-1920/2)/1920;
yy = (yy-1152/2)/1152;
rr = sqrt(xx.^2+yy.^2);
th = atan2(yy,xx);

zpolys = vm(cat(3,...
2*yy,... % tilt 1
2*xx,... % tip 2
sqrt(6)*rr.^2.*sin(2*th),... % oblique astigmatism 3
sqrt(3)*(2*rr.^2-1),...      % defocus 4
sqrt(6)*rr.^2.*cos(2*th),... % vertical astigmatism 5
sqrt(8)*rr.^3.*sin(3*th),...        % vertical trefoil 6
sqrt(8)*(3*rr.^3-2*rr).*sin(th),... % vertical coma 7
sqrt(8)*(3*rr.^3-2*rr).*cos(th),... % horizontal coma 8
sqrt(8)*rr.^3.*cos(3*th),...        % oblique trefoil 9
sqrt(10)*rr.^4.*sin(4*th),...             % Oblique quadrafoil 10
sqrt(10)*(4*rr.^4-3*rr.^2).*sin(2*th),... % Oblique secondary astigmatism 11
sqrt(5)*(6*rr.^4-6*rr.^2+1),...           % Primary spherical 12
sqrt(10)*(4*rr.^4-3*rr.^2).*cos(2*th),... % Vertical secondary astigmatism 13
sqrt(10)*rr.^4.*cos(4*th)));              % Vertical quadrafoil 14

coef_all = 0:1:15;
im_all = [];
for ii=1:length(coef_all)
    coef = coef_all(ii);
pic_fft_out = hSLMApp.current_mask;
% pic_fft_out = reshape(pic_fft_out(slm_mask),1152/2,1920/2);
% pic_fft_out = imresize(pic_fft_out,[1152 1920],'nearest');
phase_c = mod(pic_fft_out+ ...
                squeeze(zpolys.data(:,:,9))*coef ...
                ,2*pi); % blaze + lens + coma
% pic_fft_out = mod(pic_fft_out+ 2*pi*Gx*X-2*pi*Gy*Y + k/2/f*(X.^2+Y.^2),2*pi); % blaze + lens            
hSLM.project(phase_c)
pause(.5)
live_im=hDragonflyApp.camApp.snapOneFrame;
im_all(:,:,ii) = live_im; 
end
%% test spherical aberration
[xx, yy] = ndgrid(1:1152,1:1920);
xx = (xx-1920/2)/1920;
yy = (yy-1152/2)/1152;
rr = sqrt(xx.^2+yy.^2);
th = atan2(yy,xx);

zpolys = vm(cat(3,...
2*yy,... % tilt 1
2*xx,... % tip 2
sqrt(6)*rr.^2.*sin(2*th),... % oblique astigmatism 3
sqrt(3)*(2*rr.^2-1),...      % defocus 4
sqrt(6)*rr.^2.*cos(2*th),... % vertical astigmatism 5
sqrt(8)*rr.^3.*sin(3*th),...        % vertical trefoil 6
sqrt(8)*(3*rr.^3-2*rr).*sin(th),... % vertical coma 7
sqrt(8)*(3*rr.^3-2*rr).*cos(th),... % horizontal coma 8
sqrt(8)*rr.^3.*cos(3*th),...        % oblique trefoil 9
sqrt(10)*rr.^4.*sin(4*th),...             % Oblique quadrafoil 10
sqrt(10)*(4*rr.^4-3*rr.^2).*sin(2*th),... % Oblique secondary astigmatism 11
sqrt(5)*(6*rr.^4-6*rr.^2+1),...           % Primary spherical 12
sqrt(10)*(4*rr.^4-3*rr.^2).*cos(2*th),... % Vertical secondary astigmatism 13
sqrt(10)*rr.^4.*cos(4*th)));     % Vertical quadrafoil 14

pic_fft_out = hSLMApp.current_mask;

phase_c = mod(pic_fft_out+ ...
                squeeze(zpolys.data(:,:,10))*9 ...
                ,2*pi); % blaze + lens + coma
hSLM.project(phase_c)

%%
hSLM.project(hSLMApp.current_mask)