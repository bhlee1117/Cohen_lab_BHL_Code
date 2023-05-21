
cd D:\Code\Dragonfly_yq_v1\Hardware\Kinetix
mex cameraMex.cpp pvcamControl.cpp pvcam_graphicThreads_nomutex.cpp Common.cpp ...
    -ID:\Code\libs\SDL2-vc\include  -LD:\Code\libs\SDL2-vc\lib\x64 -lSDL2 ...
    '-IC:\Program Files\Photometrics\PVCamSDK\Inc' ...
    '-LC:\Program Files\Photometrics\PVCamSDK\Lib\amd64' -lpvcam64

copyfile cameraMex.mexw64 D:\Code\Dragonfly_yq_v1\Hardware\cameraMex.mexw64
cd D:\Code\Dragonfly_yq_v1
%%
cameraMex('startup',0,955,1174,2,24);


inputROI = [2048 2048];
binning = 1;
exposureTime = 15e-3;
[et,rt] = cameraMex('aq_live_restart',inputROI,binning,exposureTime)
% im = cameraMex('aq_snap');
%%
cameraMex('shutdown')
clear cameraMex