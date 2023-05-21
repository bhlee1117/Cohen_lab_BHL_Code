% Example usage of BlinkHdmiSdk.dll
% Meadowlark Optics Spatial Light Modulators
% Sept 21 2017

% Load the DLL BlinkHdmiSdk.dll
% should all be located in the same directory as the program referencing the
% library
% Matlab only supports C-style headers so this is a 'sanitized' version of
% the normal header file
addpath 'C:\Program Files\Meadowlark Optics\Blink 1920 HDMI\SDK'
loadlibrary('BlinkHdmiSdk.dll', 'BlinkHdmiSdk_matlab.h');


% Matlab automatically escapes backslashes (unlike most languages)
% Linear is a default LUT, you should replace this path to your customer
% LUT that was delivered with your SLM. This calibrates the nonlinear
% response of the LC to voltage.
lut_file = 'C:\Program Files\Meadowlark Optics\Blink 1920 HDMI\LUT Files\linear.blt';

% Arrays for image data. Read them in and format the colors so that they
% are correct for loading to the SLM.
image_0 = imread('C:\Program Files\Meadowlark Optics\Blink 1920 HDMI\Image Files\mlo.bmp');
% image_1 = imread('C:\Program Files\Meadowlark Optics\Blink 1920 HDMI\Image Files\half_wedge.bmp');
image_0 = rot90(image_0, 3);
% image_1 = rot90(image_1, 3);

% Array for LUT data
lut = libpointer('uint16Ptr', zeros(256));

sdk = calllib('BlinkHdmiSdk', 'Create_SDK', 0);

disp('Blink SDK was successfully constructed');

% Read the lookup table into memory
calllib('BlinkHdmiSdk', 'Read_lut', sdk, lut_file, lut);

% Load the lookup table to the controller. Make sure that your com port is
% correct. You can open Blink to validate the COM port. 
calllib('BlinkHdmiSdk', 'Load_lut', sdk, 'COM4', lut);

% For RGB operation use the following lines of code to specify the color in
% the image to use for the various SLMs, and load a different LUT to the 
% different SLMs. If only using a single SLM, red bits are the default and 
% the Set Channel functions doesn't need to be used.
% calllib('BlinkHdmiSdk', 'Set_channel', sdk, 'COM18', 0);
% calllib('BlinkHdmiSdk', 'Load_lut', sdk, 'COM14', lut);
% calllib('BlinkHdmiSdk', 'Set_channel', sdk, 'COM14', 1);
% calllib('BlinkHdmiSdk', 'Load_lut', sdk, 'COM15', lut);
% calllib('BlinkHdmiSdk', 'Set_channel', sdk, 'COM15', 2);

% Loop between our images
% for n = 1:3
% 	calllib('BlinkHdmiSdk', 'Write_image', sdk, image_0, 1);
% 	pause(0.5) % This is in seconds
%  	calllib('BlinkHdmiSdk', 'Write_image', sdk, image_1, 1);
%  	pause(0.5) % This is in seconds
% end
%% quadratic
[xx, yy] = ndgrid(1:size(image_0,1),1:size(image_0,2));
xx = (xx-1920/2)/1920;
yy = (yy-1152/2)/1152;
rr = sqrt(xx.^2+yy.^2);
th = atan2(yy,xx);
zeroimg = zeros(size(xx));
zpolys = vm(cat(3,...
2*yy,... % tilt
2*xx,... % tip
sqrt(6)*rr.^2.*sin(2*th),... % oblique astigmatism
sqrt(3)*(2*rr.^2-1),... % defocus
sqrt(6)*rr.^2.*cos(2*th),... % vertical astigmatism
sqrt(8)*rr.^3.*sin(3*th),... % vertical trefoil
sqrt(8)*(3*rr.^3-2*rr).*sin(th),... % vertical coma
sqrt(8)*(3*rr.^3-2*rr).*cos(th),... % horizontal coma
sqrt(8)*rr.^3.*cos(3*th),... % oblique trefoil
sqrt(10)*rr.^4.*sin(4*th),... % Oblique quadrafoil
sqrt(10)*(4*rr.^4-3*rr.^2).*sin(2*th),... % Oblique secondary astigmatism
sqrt(5)*(6*rr.^4-6*rr.^2+1),... % Primary spherical
sqrt(10)*(4*rr.^4-3*rr.^2).*cos(2*th),... % Vertical secondary astigmatism
sqrt(10)*rr.^4.*cos(4*th))); % Vertical quadrafoil

%% display one spot
excDotsRowColSLM = [
   0 50 0
    ]';
% excDotsRowColSLM = bsxfun(@plus,excDotsRowColSLM,-[0 150 0]');
% excDotsRowColSLM = bsxfun(@rdivide,excDotsRowColSLM,[.5 -.5 1]');
% excDotsRowColSLM = bsxfun(@plus,excDotsRowColSLM,+[0 250 0]');
xyz = excDotsRowColSLM;
ttd = [
    0 1 1/3 % 0 sind(60) cosd(60)
    1 0 0
    0 -1/20 1 % cosd(60)/100 sind(60)
    ]*xyz;
commonz = [
0 % Tilt (Y; vertical tilt)
0 % Tip (X; horizontal tilt)
3 % Oblique astigmatism
-21 % Defocus
-3*5 % Vertical astigmatism
0 % Vertical trefoil
0 % Vertical coma
0 % Horizontal coma
0 % Oblique trefoil
0 % Oblique quadrafoil
0 % Oblique secondary astigmatism
0 % Primary spherical
0 % Vertical secondary astigmatism
0 % Vertical quadrafoil
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
slmph1 = slmph1*128/pi+128 + 0*randn(size(slmph1));
slmph1 = max(min(slmph1,255),0);
slmph1 = uint8(slmph1);
clf
imshow(slmph1,[]) 
% colorbar

% prefix = [datestr(now,'HHMMSS') '_slm'];
% save([prefix '.mat'],'dotsparams','cplxit','slmph1')
calllib('BlinkHdmiSdk', 'Write_image', sdk, slmph1, 1);

%% load camera
vid = videoinput('hamamatsu', 1, 'MONO16_2048x2048_SlowMode');
src = getselectedsource(vid);
%% set ROI
% square centered ROI
roi_nr = 256;
roi_nc = roi_nr;
roi_or = (2048-roi_nr)/2;
roi_oc = (2048-roi_nc)/2;
vid.ROIPosition = [roi_or roi_oc roi_nr roi_nc];
%%
preview(vid);
vid.FramesPerTrigger = 600;
src.ExposureTime = 0.05;
%%
tic
start(vid);
for it = 1:14
    calllib('BlinkHdmiSdk', 'Write_image', sdk, slmph1, 1);
    pause(.1)
    calllib('BlinkHdmiSdk', 'Write_image', sdk, 0*slmph1, 1);
    pause(.1)
end
toc
img = vm(squeeze(getdata(vid)));
toc

% imshow(img,[])
% moviesc(img)
plot(img.frameAverage)

colorbar
drawnow
%%
stoppreview(vid)
