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

%% display slm calibration mask
excDotsRowColSLM = [
    599   127 0
    617   154 0
    197   172 0
   -36    245 0
   -281   121 0
   -458   102 0
   -436   222 0
    ]';
excDotsRowColSLM = bsxfun(@plus,excDotsRowColSLM,-[0 150 0]');
excDotsRowColSLM = bsxfun(@rdivide,excDotsRowColSLM,[.5 -.5 1]');
excDotsRowColSLM = bsxfun(@plus,excDotsRowColSLM,+[0 250 0]');
xyz = excDotsRowColSLM;
ttd = [
    0 1 1/3 % 0 sind(60) cosd(60)
    1 0 0
    0 -1/20 1 % cosd(60)/100 sind(60)
    ]*xyz;
commonz = [
0 % Tilt (Y; vertical tilt)
0 % Tip (X; horizontal tilt)
0 % Oblique astigmatism
-3*5 % Defocus
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
    phit = phit + mean(zpolys([1 2 4])*ttd(:,it));
    phit = phit + mean(zpolys*commonz);
    cplxit = cplxit + exp(1i*phit);
end

slmph1 = angle(cplxit);
slmph1 = slmph1*128/pi+128;
slmph1 = max(min(slmph1,255),0);
slmph1 = uint8(slmph1);
imshow(slmph1',[]),colorbar

% prefix = [datestr(now,'HHMMSS') '_slm'];
% save([prefix '.mat'],'dotsparams','cplxit','slmph1')
calllib('BlinkHdmiSdk', 'Write_image', sdk, slmph1, 1);
%%
calllib('BlinkHdmiSdk', 'Write_image', sdk, 0*slmph1, 1);

%% DMD click calibrate
basepath = 'D:\Code\Dragonfly runtime';
dirlist = dir(fullfile(basepath,'*FFb1*stimcal*.tif'));
[~, ix] = sort(cell2mat({dirlist.datenum}));
plimg = imread(fullfile(basepath,dirlist(ix(end)).name));
figure windowstyle docked
imshow((double(plimg)),[])
% zoom(2)
[xdmd, ydmd] = getpts_o;
%%
stimDotsRowColDMD = [
    599   127
    617   154
    197   172
   -36    245
   -281   121
   -458   102
   -436   222
];
stimDotsRowColDMD = bsxfun(@plus,stimDotsRowColDMD,-[0 170]);
stimDotsRowColDMD = bsxfun(@rdivide,stimDotsRowColDMD,[-2 .5]);
stimDotsRowColDMD = bsxfun(@plus,stimDotsRowColDMD,[1024 768]/2);
stimDotsRowColDMD = round(stimDotsRowColDMD);
plxy = [stimDotsRowColDMD];
[~, idplxy] = sort(plxy(:,1));
clxy = [xdmd ydmd];
[~, idclxy] = sort(clxy(:,1));
plot(plxy(idplxy,2),clxy(idclxy,2),'o')
stimCalMatrix = [clxy(idclxy,:) clxy(idclxy,1)*0+1]\[plxy(idplxy,:) plxy(idplxy,1)*0+1];

%% SLM click calibrate
basepath = 'D:\Code\Dragonfly runtime';
dirlist = dir(fullfile(basepath,'*FFb1*exccal*.tif'));
[~, ix] = sort(cell2mat({dirlist.datenum}));
plimg = imread(fullfile(basepath,dirlist(ix(end)).name));
figure windowstyle docked
imshow((double(plimg)),[])
% zoom(2)
[xslm, yslm] = getpts_o;
%%
plxy = [excDotsRowColSLM(1:2,:)'];
[~, idplxy] = sort(plxy(:,1));
clxy = [-xslm -yslm];
[~, idclxy] = sort(clxy(:,1));
plot(plxy(idplxy,2),clxy(idclxy,2),'o')
excCalMatrix = [clxy(idclxy,:) clxy(idclxy,1)*0+1]\[plxy(idplxy,:) plxy(idplxy,1)*0+1];

%% select clicky spots
dirlist = dir(fullfile(basepath,'*GFP*.tif'));
[~, ix] = sort(cell2mat({dirlist.datenum}));
plimg = imread(fullfile(basepath,dirlist(ix(end)).name));
% clear('lastClickedX','lastClickedY')
if exist('lastClickedX','var')
    if projected
        priorClickedX = [priorClickedX; lastClickedX];
        priorClickedY = [priorClickedY; lastClickedY];
    end
else
    priorClickedX = [];
    priorClickedY = [];
end
projected = false;
figure windowstyle docked
imagesc(plimg)
colormap gray
axis image
zoom(2)
hold on
plot(1024+512*[-1 -1 1 1 -1],1024+72*[-1 1 1 -1 -1])
plot(priorClickedX,priorClickedY,'o')
[lastClickedX, lastClickedY] = getpts_o;
slmlocs = [-lastClickedX -lastClickedY lastClickedX*0+1]*excCalMatrix;
plot(lastClickedX,lastClickedY,'o')
prefix = [datestr(now,'HHMMSS') '_spot_clicky'];
saveas(gcf,fullfile(basepath,[prefix '.fig']))
%%
xyz = bsxfun(@minus,slmlocs,[0 0 0*5]);
ttd = [
    0 1 1/20 % 0 sind(60) cosd(60)
    1 0 0
    0 -1/20 1 % cosd(60)/100 sind(60)
    ]*xyz';
cplxit = zeroimg;
correctedZernikeCoeffs = commonz;

for it = 1:size(ttd,2)
    phit = zeroimg;
    phit = phit + mean(zpolys([1 2 4])*ttd(:,it));
    phit = phit + mean(zpolys*correctedZernikeCoeffs);
    cplxit = cplxit + exp(1i*phit);
end

slmph1 = angle(cplxit);
slmph1 = slmph1*128/pi+128;
slmph1 = max(min(slmph1,255),0);
slmph1 = uint8(slmph1);
clf
imshow(slmph1,[])
colorbar
projected = true;
prefix = [datestr(now,'HHMMSS') '_slm_click_spots'];
save(fullfile(basepath,[prefix '.mat']),'slmlocs','cplxit','slmph1','xyz','ttd','commonz','lastClickedX','lastClickedY','excCalMatrix','stimCalMatrix','x','y','prevdefx','prevdefy')
calllib('BlinkHdmiSdk', 'Write_image', sdk, slmph1, 1);

%% display img and clicked
dirlist = dir(fullfile(basepath,'*.tif'));
[~, ix] = sort(cell2mat({dirlist.datenum}));
plimg = imread(fullfile(basepath,dirlist(ix(end)).name));
figure windowstyle docked
imshow(plimg,[])
hold on
plot(lastClickedX, lastClickedY,'o')
prefix = [datestr(now,'HHMMSS') '_spot'];
saveas(gcf,fullfile(basepath,[prefix '.fig']))


%% analysis
datadir = 'R:\';
pd = dir([datadir '\*']);
pd = pd(cell2mat({pd.isdir}));
pd = pd(arrayfun(@(x)~isequal(x.name,'.'),pd));
pd = pd(arrayfun(@(x)~isequal(x.name,'..'),pd));
pd = pd(arrayfun(@(x)~isequal(x.name,'Temp'),pd));
pd = pd(arrayfun(@(x)~isequal(x.name,'rec'),pd));
[d,ix] = sort(arrayfun(@(x)x.datenum,pd));
sel = ix(end);
pname = fullfile(datadir,pd(sel).name,'Sq_camera.bin');
dcname = fullfile(datadir,pd(sel).name,'Sq_camera.dcimg');
for it = 1:600
    if ~exist(pname,'file')
        pause(.1)
    else
        break
    end
end
for it = 1:600
    if exist(dcname,'file')
        pause(.1)
    else
        break
    end
end

disp(pd(sel).name)
clear mov
mov = vm(fullfile(datadir,pd(sel).name))';
% figure windowstyle docked
% moviesc(mov)
% title(strrep(pd(sel).name,'_','\_'))
% return
%
refimg = imgaussfilt(mov.mean,4,'Padding',0);
locmask = findalllocs(refimg,10,prctile(refimg(:),80));
locIdxs = find(locmask);
[~, ix] = sort(refimg(locmask));
[ixi, ixj] = ind2sub(mov.imsz,locIdxs(ix(end-(1:min(length(lastClickedX),end))+1)));
% clf
% imshow(refimg,[])
% hold on
% plot(ixj,ixi,'o')
%
genroi = 3*[-1,-1;1,-1;1,1;-1,1;-1,-1];
rois = {};
for it = 1:numel(ixi)
    rois = [rois {bsxfun(@plus,genroi,[ixj(it),ixi(it)])}];
end
figure
apply_clicky_faster(mov,rois);
title(strrep(pd(sel).name,'_','\_'))
saveas(gcf,fullfile(datadir,pd(sel).name,[datestr(now,'HHMMSS') '_clicky.fig']))

%%

[rois, intens] = nested_clicky_rects(mov(1000:end));
figure
apply_clicky_faster(rois,mov);
saveas(gcf,fullfile(datadir,pd(sel).name,[datestr(now,'HHMMSS') '_clicky.fig']))
%%
% Select analysis ROIs
% 1. background
% 2. ROI with multiple cells
% 3..N ROIs enclosing each cell
[rois, intens] = clicky_rects(mov);
noiseTrace = imgaussfilt(intens(:,1),eps);
photoTrace = imgaussfilt(intens(:,2),eps);
photoMov = mov.crop_rect([min(rois{2}) max(rois{2})-min(rois{2})]);
photoMovMean = photoMov.mean;
[fx, fy] = gradient(photoMovMean);
oftraces = (tovec(cat(3,photoMovMean,fx,fy))\photoMov.tovec.double)';
dxTrace = imgaussfilt(oftraces(:,2),eps);
dyTrace = imgaussfilt(oftraces(:,3),eps);
figure windowstyle docked
plot([noiseTrace photoTrace dxTrace*100 dyTrace*100])
%% % 3..N ROIs enclosing each cell
it = 4;
cellMov = mov.crop_rect([min(rois{it}) max(rois{it})-min(rois{it})]);
remTraces = [noiseTrace photoTrace dxTrace dyTrace];
mMasks = vm(double(cellMov(1000:end)))/remTraces(1000:end,:)';
newMov = cellMov - mMasks*remTraces';
[u, s, v] = pca_eig(newMov(:,1000:end),5);
[ics, mixmat, sepmat] = sorted_ica(v,5);
a = newMov(:,1000:end)/(ics');
% subplot(121)
% imshow(montage(vm(a,newMov.imsz),[5]),[])
% subplot(122)
stackplot((a\newMov(:,1000:end))')
% title 'brightest spot'
% saveas(gcf,fullfile(datadir,pd(sel).name,[datestr(now,'HHMMSS') '_clicky.fig']))
%%
% [min(rois{1}) max(rois{1})-min(rois{1})]
% image_2 = image_2 + sqrt((rr-1920/2).^2)*-4/100;
% image_2 = image_2 + sqrt((rr-1920/2).^2+0*(cc-1152/2).^2)*-4/100;

% dotsparams = [
% %     0  -20 0 0 0 0
% %     0  50 0 0 0 0 
% %     0  30 0 0 0 0
% %     0  -70 0 0 0 0
% ];
% randparams = [
% 	11 1
% ];
% randparams = bsxfun(@minus,rand(100,2)*diag([1000 10]),[500 0]);
% randparams = [
%     ((-100:2:100)-18)'*[1 0] + 2*rand(101,2) + [0*rand(101,1) 2*rand(101,1)]
%     ((-50:2:50)-18)'*[1 0] + 2*rand(51,2) + [0*rand(51,1) 2*rand(51,1)]
%     ];
% dotsparams = [0*randparams(:,1) randparams 0*randparams(:,[1 1 1 1])];
% dotsparams = [
%      0   -37     3     1     -4    0 0
%      0   -71     7     0     0     0 0
%      0   159     7     0     0     0 0
%      0   177     2     0     0     0 0
%      0  -188     1     0     0     0 0
%      0    -6     5     0     0     0 0
%      0   -33    10     0     0     0 0
%      0    88     3     0     0     0 0
%      0   126     6     0     1     1 2
%      0   154     2     0     1     1 2
% ];
dotsparams = [
%      0   108 1 0 0 0 0
%      0   -38 4 0 0 0 0
%      0   -18 2.5 0 0 0 0
%      0   13 3.5 0 0 0 0
%      0   25 1 0 0 0 0
%      0   40 3 0 0 0 0
%      0   87 4 0 0 0 0
%      0   75 1 0 0 0 0
%      0   -5 4 0 0 0 0
%      0   63 3 0 0 0 0
    0 0 0 0 0 0 0
];

% dotsparams = bsxfun(@plus,dotsparams,[0 0 12 0 40 16 ]);
dotsparams = bsxfun(@plus,dotsparams,[0 0 20 -5 -3 0 -1]);
cplxit = zeroimg;
for it = 1:size(dotsparams,1)
    phit = zeroimg;
    phit = phit + dotsparams(it,1);
    phit = phit + dotsparams(it,2)*10 * xx;
    phit = phit + dotsparams(it,3)*10 * yy;
    phit = phit + dotsparams(it,4)*10 * rr.^2;
    phit = phit + dotsparams(it,5)*10 * xx.^2;
    phit = phit + dotsparams(it,6)*10 * yy.^2;
    phit = phit + dotsparams(it,7)*10 * xx.*yy;
    phit = phit + 2*10 * (6*rr.^4-1*rr.^2+1);
    cplxit = cplxit + exp(1i*phit);
end

slmph1 = angle(cplxit);
slmph1 = slmph1*128/pi+128;
slmph1 = max(min(slmph1,255),0);
slmph1 = uint8(slmph1);
imshow(slmph1,[]),colorbar

% prefix = [datestr(now,'HHMMSS') '_slm'];
% save([prefix '.mat'],'dotsparams','cplxit','slmph1')
calllib('BlinkHdmiSdk', 'Write_image', sdk, slmph1, 1);

%%
% ph1 = zeroimg;
% ph1 = ph1 + 0;
% ph1 = ph1 + xx*-20*10;
% ph1 = ph1 + yy*0*10;
% ph1 = ph1 + rr.^2*10*10;
% cplx1 = exp(1i*ph1);
% 
% ph2 = zeroimg;
% ph2 = ph2 + 21;
% ph2 = ph2 + xx*40*10;
% ph2 = ph2 + yy*-10*10;
% ph2 = ph2 + rr.^2*10*10;
% cplx2 = exp(1i*ph2);

% image_2 = image_2 + cc*3;
% period = 50;
% image_2 = (mod(rr,period)>period/2)*255;
% period 100
% [100 150 166 176 195 211 218 224 236 246 255]
% [max  =  min  =  max  =  min  =  max  =  min]
% image_2 = image_2 - min(image_2(:));
% image_2 = image_2./max(image_2(:));
% image_2 = mod(image_2,256);
% image_2 = uint8(image_2);
% imshow(image_2);
%%
mylut = 'X:\Lab\Labmembers\Vicente Parot\Code\2017-12-02 slm tests\mylut.blt';
dlmwrite(mylut,floor(logspace(1,2048,256)'),'\n');
% calllib('BlinkHdmiSdk', 'Read_lut', sdk, lut_file, lut);
lut = libpointer('uint16Ptr', zeros(256));
calllib('BlinkHdmiSdk', 'Read_lut', sdk, mylut, lut);
calllib('BlinkHdmiSdk', 'Load_lut', sdk, 'COM4', lut);
period = 50;
ph1 = (mod(xx,period)>period/2)*255;
% period 100
% [100 150 166 176 195 211 218 224 236 246 255]
% [max  =  min  =  max  =  min  =  max  =  min]
% image_2 = image_2 - min(image_2(:));
% image_2 = image_2./max(image_2(:));
ph1 = mod(ph1,256);
ph1 = uint8(ph1);
imshow(ph1);
calllib('BlinkHdmiSdk', 'Write_image', sdk, ph1, 1);

%%
% Always call Delete_SDK before exiting
calllib('BlinkHdmiSdk', 'Delete_SDK', sdk);

unloadlibrary('BlinkHdmiSdk');