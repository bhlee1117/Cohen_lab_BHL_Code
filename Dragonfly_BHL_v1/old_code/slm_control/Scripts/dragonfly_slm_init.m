% % % % % % Example usage of BlinkHdmiSdk.dll
% % % % % % Meadowlark Optics Spatial Light Modulators
% % % % % % Sept 21 2017
% % % % % 
% % % % % % Load the DLL BlinkHdmiSdk.dll
% % % % % % should all be located in the same directory as the program referencing the
% % % % % % library
% % % % % % Matlab only supports C-style headers so this is a 'sanitized' version of
% % % % % % the normal header file
% % % % % addpath 'C:\Program Files\Meadowlark Optics\Blink 1920 HDMI\SDK'
% % % % % loadlibrary('BlinkHdmiSdk.dll', 'BlinkHdmiSdk_matlab.h');
% % % % % 
% % % % % 
% % % % % % Matlab automatically escapes backslashes (unlike most languages)
% % % % % % Linear is a default LUT, you should replace this path to your customer
% % % % % % LUT that was delivered with your SLM. This calibrates the nonlinear
% % % % % % response of the LC to voltage.
% % % % % lut_file = 'C:\Program Files\Meadowlark Optics\Blink 1920 HDMI\LUT Files\linear.blt';
% % % % % 
% % % % % % Arrays for image data. Read them in and format the colors so that they
% % % % % % are correct for loading to the SLM.
% % % % % image_0 = imread('C:\Program Files\Meadowlark Optics\Blink 1920 HDMI\Image Files\mlo.bmp');
% % % % % % image_1 = imread('C:\Program Files\Meadowlark Optics\Blink 1920 HDMI\Image Files\half_wedge.bmp');
% % % % % image_0 = rot90(image_0, 3);
% % % % % % image_1 = rot90(image_1, 3);
% % % % % 
% % % % % % Array for LUT data
% % % % % lut = libpointer('uint16Ptr', zeros(256));
% % % % % 
% % % % % sdk = calllib('BlinkHdmiSdk', 'Create_SDK', 0);
% % % % % isempty(sdk)
% % % % % 
% % % % % disp('Blink SDK was successfully constructed');
% % % % % 
% % % % % % Read the lookup table into memory
% % % % % calllib('BlinkHdmiSdk', 'Read_lut', sdk, lut_file, lut);
% % % % % 
% % % % % % Load the lookup table to the controller. Make sure that your com port is
% % % % % % correct. You can open Blink to validate the COM port. 
% % % % % calllib('BlinkHdmiSdk', 'Load_lut', sdk, 'COM4', lut);
% % % % % 
% % % % % % For RGB operation use the following lines of code to specify the color in
% % % % % % the image to use for the various SLMs, and load a different LUT to the 
% % % % % % different SLMs. If only using a single SLM, red bits are the default and 
% % % % % % the Set Channel functions doesn't need to be used.
% % % % % % calllib('BlinkHdmiSdk', 'Set_channel', sdk, 'COM18', 0);
% % % % % % calllib('BlinkHdmiSdk', 'Load_lut', sdk, 'COM14', lut);
% % % % % % calllib('BlinkHdmiSdk', 'Set_channel', sdk, 'COM14', 1);
% % % % % % calllib('BlinkHdmiSdk', 'Load_lut', sdk, 'COM15', lut);
% % % % % % calllib('BlinkHdmiSdk', 'Set_channel', sdk, 'COM15', 2);
% % % % % 
% % % % % % Loop between our images
% % % % % % for n = 1:3
% % % % % % 	calllib('BlinkHdmiSdk', 'Write_image', sdk, image_0, 1);
% % % % % % 	pause(0.5) % This is in seconds
% % % % % %  	calllib('BlinkHdmiSdk', 'Write_image', sdk, image_1, 1);
% % % % % %  	pause(0.5) % This is in seconds
% % % % % % end
image_0 = zeros(1920,1152);
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

zpolys2 = vm(cat(3,...
sqrt(7)*(-1+12.*rr.^2-30.*rr.^4+20.*rr.^6),...	m = 6, n = 0
sqrt(14)*rr.^4.*(-5+6.*rr.^2).*cos(4*th),...    m = 6, n = 4
sqrt(9)*(1-20.*rr.^2+90.*rr.^4-140.*rr.^6+70.*rr.^8),...	m = 8, n = 0
sqrt(18)*rr.^4.*(15-42.*rr.^2+28.*rr.^4).*cos(4*th),...     m = 8, n = 4
sqrt(18)*rr.^8.*cos(8*th),...                               m = 8, n = 8
sqrt(11)*(-1+30.*rr.^2-210.*rr.^4+560.*rr.^6-630.*rr.^8+252.*rr.^10),...    m = 10, n = 0
sqrt(22)*rr.^4.*(-35+168.*rr.^2-252.*rr.^4+120.*rr.^6).*cos(4*th),...       m = 10, n = 4
sqrt(22)*rr.^8.*(-9+10.*rr.^2).*cos(8*th),...                               m = 10, n = 8
sqrt(13)*(1-42.*rr.^2+420.*rr.^4-1680.*rr.^6+3150.*rr.^8-2772.*rr.^10+924.*rr.^12),...	m = 12, n = 0
sqrt(26)*rr.^4.*(70-504.*rr.^2+1260.*rr.^4-1320.*rr.^6+495.*rr.^8).*cos(4*th),...       m = 12, n = 4
sqrt(26)*rr.^8.*(45-110.*rr.^2+66.*rr.^4).*cos(8*th),...                                m = 12, n = 8
sqrt(26)*rr.^12.*cos(12*th))); %                                                        m = 12, n = 12











