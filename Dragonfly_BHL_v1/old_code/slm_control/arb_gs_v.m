

test_im = zeros(1152,1920);
% [ROIs,~]=clicky_faster(test_im);
ROIs = click_im(test_im,'rect');

[X,Y]=ndgrid(1:1920,1:1152);

in = cellfun(@(pixel_pos) inpolygon(X,Y,pixel_pos(:,1),pixel_pos(:,2)),ROIs,'uniformoutput',false);
x_spt = [];
y_spt = [];    
area =100*100;
for i=1:length(in)
    in_i = find(in{i});
    [x_spt_i,y_spt_i] = ind2sub(size(in{i}),in_i(1:area:end));
    
    x_spt_i = x_spt_i-1920/2;
    y_spt_i = y_spt_i-1152/2;
    
    x_spt = [x_spt; x_spt_i];
    y_spt = [y_spt; y_spt_i];
end
z_spt = zeros(size(x_spt));
%%

theta = linspace(2*pi/10,2*pi,10)';
x_spt = cos(theta)*300;
y_spt = sin(theta)*300;

z_spt = zeros(size(x_spt));

%%

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
27 % Horizontal coma
0 % Oblique trefoil
0 % Oblique quadrafoil
0 % Oblique secondary astigmatism
0 % Primary spherical
0 % Vertical secondary astigmatism
0 % Vertical quadrafoil
];
%%
pic_out = gs(commonz, zpolys, [x_spt y_spt z_spt]', 5);

x = (linspace(-1920/2,1920/2,1920))*9.2e-3;
y = (linspace(-1152/2,1152/2,1152))*9.2e-3;
[X,Y] = meshgrid(x,y);
theta = 45;
Gx = 1/(4/1920*max(X(:)*2));
Gy = 1/(4/1152*max(Y(:)*2));
lambda = .594e-3;
k=2*pi/lambda;
f=3050; %mm

pic_out = mod(angle(pic_out)+ 2*pi*Gx*X-2*pi*Gy*Y + k/2/f*(X.^2+Y.^2),2*pi); % blaze + lens

pic_out = uint8(mat2gray(pic_out)*255);
% figure(100);clf;imshow(pic_fft_uint8,[])
imwrite(pic_out,'C:\Program Files\Meadowlark Optics\Blink 1920 HDMI\Image Files\arb_gs_v.bmp')