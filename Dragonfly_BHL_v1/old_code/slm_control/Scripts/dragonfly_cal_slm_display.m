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
excDotsRowColSLM = bsxfun(@plus,excDotsRowColSLM,-[0 50 -0]');
excDotsRowColSLM = bsxfun(@rdivide,excDotsRowColSLM,[.5 -.5 1]'*2);
excDotsRowColSLM = bsxfun(@plus,excDotsRowColSLM,+[0 250 0]');
xyz = excDotsRowColSLM;
ttd = [
    0 1 1/3 % 0 sind(60) cosd(60)
    1 0 0
    0 -1/20 1 % cosd(60)/100 sind(60)
    ]*xyz;
commonz = [
% 0 % Tilt (Y; vertical tilt)
% 0 % Tip (X; horizontal tilt)
% 0 % Oblique astigmatism
% 20*5 % Defocus
% 5*10 % Vertical astigmatism
% 0 % Vertical trefoil
% 0 % Vertical coma
% 0 % Horizontal coma
% 0 % Oblique trefoil
% 0 % Oblique quadrafoil
% 0 % Oblique secondary astigmatism
% 10*5 % Primary spherical
% 0 % Vertical secondary astigmatism
% 0*10 % Vertical quadrafoil

% 0 % Tilt (Y; vertical tilt)
% 40 % Tip (X; horizontal tilt)
% 0 % Oblique astigmatism
% 5*5 % Defocus
% 5*5 % Vertical astigmatism
% 0 % Vertical trefoil
% 0 % Vertical coma
% 0 % Horizontal coma
% 0*10 % Oblique trefoil
% 0*-10 % Oblique quadrafoil
% 0 % Oblique secondary astigmatism
% 1*10 % Primary spherical
% -0 % Vertical secondary astigmatism
% -0*10 % Vertical quadrafoil

% optimized 2018-07-27
  -65 % -68.0259
  140 % 88.5640
   23.7727
   25.6780
    3.0698
   -5.4880
   -3.9627
   54.4698
   -3.5866
  -22.5700
   -5.6261
    9.0160
  -18.5891
   -8.3086
    ];

approxSLM = gs(commonz, zpolys, ttd);

slmph1 = angle(approxSLM);
slmph1 = slmph1*128/pi+128;
slmph1 = max(min(slmph1,255),0);
slmph1 = uint8(slmph1);
imshow(slmph1,[]),colorbar

% prefix = [datestr(now,'HHMMSS') '_slm'];
% save([prefix '.mat'],'dotsparams','cplxit','slmph1')
% calllib('BlinkHdmiSdk', 'Write_image', sdk, slmph1, 1);
    fullscreen(rot90(slmph1),3)