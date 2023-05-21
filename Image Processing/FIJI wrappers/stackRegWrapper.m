function imgsRegistered = stackRegWrapper(imgs, transformationType)

% THIS FUNCTION REQUIRES INSTALLATION OF Fiji
% DOWNLOAD FIJI HERE: http://fiji.sc/Downloads
%
% this function ILLUSTRATES a method to call
% the stackReg function in FIJI (ImageJ) with MIJI 
%
% Thanks to the support team at FIJI,MIJI,StackReg,ImageJ
% for these incredibly easy to use tools!
% 
% useful references
%
%  StackReg :  
%  Fiji: http://fiji.sc/Fiji
%  Miji: http://fiji.sc/Miji
%  ImageJ: http://imagej.nih.gov/ij/
%
% $Revision: 1.0 $ $Date: 2014/02/20 08:00$ $Author: Pangyu Teng $

if nargin < 1
	display('need at least 1 input! (stackRegWrapper.m)');
	return;
end

if nargin < 2
	transformationType = 'Translation';
	disp('Translation used for stackreg registration');
end

% start FIJI without GUI.
Miji(false);

% transfer images from Matlab to FIJI
MIJ.createImage(imgs);
%%alternatively, above line can be replaced with below line and changing the input 'imgs' with the 'filename' of a multi-page tiff
%MIJ.run('Open...', ['path=[' filename ']']);  

% register
MIJ.run('StackReg', ['transformation=' transformationType]);

% transfer image from FIJI to Matlab
imgsRegistered = MIJ.getCurrentImage;

% close image window in FIJI
MIJ.run('Close');
	
% exit FIJI
MIJ.exit