%% initialize projector
lab_init_device
%%
device.free;
disp 'device free' 
%% switch back to master mode and Illuminate
            roirows = 18*10+0:55*10+3; % extended central quad 2017-12-18
            roicols = 32*10+0:69*10+5;
% roirows = 1:768; % extended full fov 2017-03-29
% roicols = 14*10+0:87*10+1;
% roirows = 1:768; % extended full fov 2017-12-13
% roicols = 14*10+0:89*10+1;
% roirows = 18*10+9:56*10+9; % extended central half 2017-12-02
% roicols = 13*10+0:87*10+1;
% roirows = 19*10+0:57*10+3; % extended central quad 2017-03-17
% roicols = 32*10+0:69*10+5;
%              roirows = 33*10+0:42*10+3; % rat in vivo roi 2017-11-30 256
%              roicols = 45*10+0:55*10+1;
%              roirows = 33*10+0:42*10+3; % rat in vivo roi 2017-11-30 256
%              roicols = 35*10+0:59*10+1;
% roirows = 29*10+0:57*10+3; % stim roi
% roicols = 39*10+0:69*10+5;
% roirows = 19*10+9:57*10+9; % extended central half 2017-02-13
% roicols = 13*10+0:87*10+1;

pat = zeros(device.height,device.width);

pat(roicols,roirows) = 1;
% pat(:) = 1;
% rng(142857)
% pat = pat.*(rand(size(pat))<.005);
% pat(407*1000+10*100+(1:100)) = 1;
% nnz(pat)
% litidx = find(pat);
% save reg09 litidx

% device.halt;
% device.projcontrol(api.PROJ_MODE,api.MASTER);
device.put(pat*255)

%% switch back to master mode and project selective mask
device.halt;
device.projcontrol(api.PROJ_MODE,api.MASTER);
% new_peaks_mask = circshift(peaks_mask2,[0 0]);
device.put(spots_masks(:,:,3)'*255)
% device.put((checkerboard')*255)

%% switch back to master mode and FLASH
device.halt;
device.projcontrol(api.PROJ_MODE,api.MASTER);
% simple test. set to white, wait 1/4 s, set to black.
device.put(ones(device.height,device.width)*255)
% pause(.25)
device.put(zeros(device.height,device.width))

