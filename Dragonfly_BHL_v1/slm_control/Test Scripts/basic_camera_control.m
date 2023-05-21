%% load camera
imaqreset
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
figure windowstyle docked
%%
preview(vid);
vid.FramesPerTrigger = 600;
src.ExposureTime = 0.105;
%%
tic
start(vid);
for it = 1:14
    fullscreen(rot90(slmph1),3)
    pause(.05)
    fullscreen(0*rot90(slmph1),3)
    pause(.05)
%     calllib('BlinkHdmiSdk', 'Write_image', sdk, slmph1, 1);
%     pause(.1)
%     calllib('BlinkHdmiSdk', 'Write_image', sdk, 0*slmph1, 1);
%     pause(.075)
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
%%
plot(img.frameAverage)
plot(img.frameAverage>threshold)
threshold = 135;
frid = [find(diff(img.frameAverage>threshold)>0)+1 find(diff(img.frameAverage>threshold)<0)];
assert(size(frid,1)==14,'missing frames')
assert(min(frid(:,2)-frid(:,1))>=15,'missing frames')
im = img(frid(:,1)+(6:15)).blnfun(@mean,10);











