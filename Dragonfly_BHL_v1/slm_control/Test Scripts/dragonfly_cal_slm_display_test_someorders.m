% %% display slm calibration mask
% excDotsRowColSLM = [
%     599   127 0
%     617   154 0
%     197   172 0
%    -36    245 0
%    -281   121 0
%    -458   102 0
%    -436   222 0
%     ]';
% excDotsRowColSLM = bsxfun(@plus,excDotsRowColSLM,-[0 50 -0]');
% excDotsRowColSLM = bsxfun(@rdivide,excDotsRowColSLM,[.5 -.5 1]'*2);
% excDotsRowColSLM = bsxfun(@plus,excDotsRowColSLM,+[0 250 0]');
% xyz = excDotsRowColSLM;
% ttd = [
%     0 1 1/3 % 0 sind(60) cosd(60)
%     1 0 0
%     0 -1/20 1 % cosd(60)/100 sind(60)
%     ]*mean(xyz);
% %% zernike even 4n modes
% [xx, yy] = ndgrid(1:size(image_0,1),1:size(image_0,2));
% xx = (xx-1920/2)/1920;
% yy = (yy-1152/2)/1152;
% rr = sqrt(xx.^2+yy.^2);
% th = atan2(yy,xx);
% zeroimg = zeros(size(xx));
% zpolys = vm(cat(3,...
% 2*yy,... % tilt
% 2*xx,... % tip
% sqrt(6)*rr.^2.*sin(2*th),... % oblique astigmatism
% sqrt(3)*(2*rr.^2-1),... % defocus
% sqrt(6)*rr.^2.*cos(2*th),... % vertical astigmatism
% sqrt(8)*rr.^3.*sin(3*th),... % vertical trefoil
% sqrt(8)*(3*rr.^3-2*rr).*sin(th),... % vertical coma
% sqrt(8)*(3*rr.^3-2*rr).*cos(th),... % horizontal coma
% sqrt(8)*rr.^3.*cos(3*th),... % oblique trefoil
% sqrt(10)*rr.^4.*sin(4*th),... % Oblique quadrafoil
% sqrt(10)*(4*rr.^4-3*rr.^2).*sin(2*th),... % Oblique secondary astigmatism
% sqrt(5)*(6*rr.^4-6*rr.^2+1),... % Primary spherical
% sqrt(10)*(4*rr.^4-3*rr.^2).*cos(2*th),... % Vertical secondary astigmatism
% sqrt(10)*rr.^4.*cos(4*th))); % Vertical quadrafoil

%%
commonz = [
0 % Tilt (Y; vertical tilt)
0 % Tip (X; horizontal tilt)
0 % Oblique astigmatism
20*5 % Defocus
5*10 % Vertical astigmatism
0 % Vertical trefoil
0 % Vertical coma
0 % Horizontal coma
0 % Oblique trefoil
0 % Oblique quadrafoil
0 % Oblique secondary astigmatism
10*5 % Primary spherical
0 % Vertical secondary astigmatism
0*10 % Vertical quadrafoil
    ];

cplxit = zeroimg;
for it = 1:size(ttd,2)
    phit = zeroimg;
    phit = phit + sum(zpolys([1 2 4])*ttd(:,it));
    phit = phit + sum(zpolys*commonz);
%     phit = phit + sum(zpolys*updatedcoeffs);
    cplxit = cplxit + exp(1i*(phit+24*pi/4));
end



slmph1 = angle(cplxit);
slmph1 = slmph1*128/pi+128;
slmph1 = max(min(slmph1,255),0);
slmph1 = uint8(slmph1);
imshow(slmph1,[]) 
fullscreen(rot90(slmph1),3)
