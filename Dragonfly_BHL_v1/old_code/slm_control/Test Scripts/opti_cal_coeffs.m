%% NI DAQ init
s = daq.createSession('NI');
addDigitalChannel(s,'Dev1','Port0/Line22','OutputOnly');
addAnalogInputChannel(s,'Dev1','ai0','Voltage');
% outputSingleScan(s,1)
% outputSingleScan(s,0)
% queueOutputData(s,[ones(200,1);zeros(1,1)]);
% s.startForeground;
% clear s

%% camera init
% imaqreset
vid = videoinput('hamamatsu', 1, 'MONO16_2048x2048_SlowMode');
src = getselectedsource(vid);
%% set ROI
% square centered ROI
roi_nr = 256;
roi_nc = 256;
roi_or = (2048-roi_nr)/2;
roi_oc = (2048-roi_nc)/2;
vid.ROIPosition = [roi_oc roi_or roi_nc roi_nr];
hImage = preview(vid);
hImage.CDataMapping = 'scaled';
hImage.Parent.CLim = [0 200];
set(hImage,'cdatamapping','scaled')
hold(hImage.Parent,'on')
plot(hImage.Parent,lastClickedX-cplc/2+plc/2, lastClickedY-cplr/2+plr/2,'o')

%% PSF optimization with zernike coefficients
% refocus and avoid photobleaching
% preview(vid);
src.ExposureTime = .5;
fullscreen(rot90(slmph1),3)

%% reset optim params
% original_xyz_offset = mean(best_xyz_offset);
original_xyz_offset = [0 0 0];
original_commonz_offset = zeros(1,14);
optimParamScale = 10*ones(1,14);
nParams = numel(optimParamScale);
optimSteps = linspace(-1,1,10);
nSteps = numel(optimSteps);
optimImgIdx = sub2ind([roi_nr roi_nc],round(lastClickedY-cplr/2+roi_nr/2),round(lastClickedX-cplc/2+roi_nc/2));
% optimImgIdx = 256*128+128;
nSpots = numel(optimImgIdx);

%% set up timing for auto AQ
fullscreen(0*rot90(slmph1),3)
slmSettlingTime = 0.1; % s
src.ExposureTime = 0.01; % s
nLaserFrames = 10;
laserOnTime = src.ExposureTime*(nLaserFrames + 2)*s.Rate; % ms
vid.FramesPerTrigger = ceil((laserOnTime/1000 + slmSettlingTime)*nSteps/src.ExposureTime*3);

%% reset optim values
best_xyz_offset = ones(nSpots,1)*original_xyz_offset;
best_commonz_offset = ones(nSpots,1)*original_commonz_offset;

%%
figure windowstyle docked
% vid.FramesPerTrigger = ceil((laserOnTime/1000 + slmSettlingTime)*nSteps/src.ExposureTime*2);
outputSingleScan(s,0) % shut off shutter
% outputSingleScan(s,1)
% alldata = zeros(roi_nr,roi_nc,vid.FramesPerTrigger,nParams);
avgFrames = zeros(roi_nr,roi_nc,nSteps,nParams);
deltaOptim = zeros(nSpots,nParams);
allSlmPhases = uint8(zeros(1152,1920,nSteps));
optimStepSortList = ceil(((1:nSteps).*(-1).^(1:nSteps))/2) + ceil(nSteps/2);
% optimStepSortList = optimStepSortList(end:-1:1);
recenter = 1;
if recenter
    itParamList = 1:2;
else
    itParamList = 3:14;
end
% itParamList = [7 11];
for itParam = itParamList % numel(best_commonz_offset) % optimize one parameter each time
    for itStep = 1:nSteps
        test_offset = best_commonz_offset;
        test_offset(:,itParam) = test_offset(:,itParam) + optimParamScale(itParam).*optimSteps(itStep);
%         test_xyz_offset
        xyz = slmlocs - original_xyz_offset;
        ttd = [
            0 1 1/3 % 0 sind(60) cosd(60)
            1 0 0
            0 -1/20 1 % cosd(60)/100 sind(60)
            ]*xyz';
        commonz = test_offset';
        approxSLM = gs(commonz, zpolys, ttd);
        slmph1opt = angle(approxSLM);
        slmph1opt = slmph1opt*128/pi+128;
        slmph1opt = max(min(slmph1opt,255),0);
        slmph1opt = uint8(slmph1opt);
        allSlmPhases(:,:,itStep) = rot90(slmph1opt);
    end
%     moviesc(permute(allframes,[2 1 3]))
    start(vid);
    for itStep = optimStepSortList
        fullscreen(allSlmPhases(:,:,itStep),3)
        pause(slmSettlingTime)
        queueOutputData(s,[ones(round(laserOnTime),1);zeros(1,1)]);
        s.startForeground;
    end
    fprintf('AQ %d%%\n',ceil(itParam/nParams*100));
    bdat = blur(vm(squeeze(getdata(vid))),2);
%     if recenter
%         traces = bdat(optimImgIdx,:)';
%     else
        traces = max(bdat(:,:))';
%     end
    traces = traces - prctile(traces,25);
    aggregateTrace = medfilt2(mean(traces./std(traces),2),[3 1]);
%     agg_diffs = max(tovec(alldata(:,:,:,it_steps)));
    subplot(nSpots+1,nParams,itParam)
    plot(aggregateTrace)
    drawnow
    [~,locs1] = findpeaks(+diff(aggregateTrace),'NPeaks',nSteps,'MinPeakDistance',10,'SortStr','descend');
    locs1 = sort(locs1)+1;
    [~,locs2] = findpeaks(-diff(aggregateTrace),'NPeaks',nSteps,'MinPeakDistance',10,'SortStr','descend');
    locs2 = sort(locs2)-1;
    disp ' '
    disp([locs1 diff([locs1 locs2]')' locs2])
    if any(diff([locs1 locs2]')' < 10)
        continue
    end
    for itSortStep = 1:nSteps
        itStep = optimStepSortList(itSortStep);
        avgFrames(:,:,itStep,itParam) = bdat(locs1(itSortStep):locs2(itSortStep)).mean;
    end
    bstepsmov = vm(avgFrames(:,:,:,itParam));
    plotSteps = linspace(optimSteps(1),optimSteps(end),30);
    for itSpot = 1:nSpots
        if recenter
            objVals = bstepsmov(optimImgIdx(itSpot),:);
        else
            objVals = max(bstepsmov(:,:));
            ihs = 60; % image half side
            mindistance = 128;
            newmov = zeros(2*ihs+1,2*ihs+1,nSteps);
            for itSortStep = 1:nSteps
                [idr, idc] = find(findlocs(bstepsmov(itSortStep).mean,mindistance,1));
                newmov(:,:,itSortStep) = bstepsmov(idr+(-ihs:ihs),idc+(-ihs:ihs),itSortStep).mean;
            end
%             imgbg = mean([mean(mean(newmov(:,[1 end],:)),3) mean(mean(newmov([1 end],:,:),2),3)']);
            imgbg = prctile(newmov(:),95);
            moms = imgmoments(max(newmov-imgbg,0)); 
            objVals = -sqrt(prod(moms(6:7,:)));
        end
        p = polyfit(optimSteps,objVals,2);
        deltaOptim(itSpot,itParam) = -p(2)/p(1)/2;
        if abs(deltaOptim(itSpot,itParam)) > 1 || p(1) > 0
            deltaOptim(itSpot,itParam) = mean(optimSteps(objVals == max(objVals)));
        end
        subplot(nSpots+1,nParams,itParam+(itSpot)*nParams)
        plot(optimSteps,objVals,...
            plotSteps,polyval(p,plotSteps),...
            deltaOptim(itSpot,itParam),polyval(p,deltaOptim(itSpot,itParam)),'o')
    end    
    drawnow
    deltaOptim = min(deltaOptim,+1);
    deltaOptim = max(deltaOptim,-1);
    disp ' '
    disp(deltaOptim)
    best_commonz_offset(:,itParam) = best_commonz_offset(:,itParam) + deltaOptim(:,itParam).*optimParamScale(itParam)
%     itParam = itParam + 1;
end
prefix = [datestr(now,'HHMMSS') '_opti_psf'];
saveas(gcf,fullfile(ramfovpath,[prefix '.fig']))
%%
figure windowstyle docked
        test_offset = best_xyz_offset;
        xyz = slmlocs - test_offset;
        ttd = [
            0 1 1/3 % 0 sind(60) cosd(60)
            1 0 0
            0 -1/20 1 % cosd(60)/100 sind(60)
            ]*xyz';
        approxSLM = gs(commonz, zpolys, ttd, 50, false, true, false, true, true);
%         approxSLM = gs(commonz, zpolys, ttd, 1, false, true, false, true, true, 1./(max(amps)./amps').^1);
        slmph1opt = angle(approxSLM);
%         slmph1opt = (sign(slmph1opt)/2+.5)*128;
        slmph1opt = slmph1opt*128/pi+128;
        slmph1opt = max(min(slmph1opt,255),0);
        slmph1opt = uint8(slmph1opt);
    vid.FramesPerTrigger = ceil((laserOnTime/1000 + slmSettlingTime)*1/src.ExposureTime*1.5);
    start(vid);
        fullscreen(allSlmPhases(:,:,itStep),3)
        pause(slmSettlingTime)
        queueOutputData(s,[ones(round(laserOnTime),1);zeros(1,1)]);
        s.startForeground;
    dat = vm(squeeze(getdata(vid)));
    bdat = blur(dat,2);
    traces = bdat(optimImgIdx,:)';
    traces = traces - prctile(traces,25);
    aggregateTrace = medfilt2(mean(traces./std(traces),2),[3 1]);
    plot(aggregateTrace)
    drawnow
    [~,locs1] = findpeaks(+diff(aggregateTrace),'NPeaks',1,'MinPeakDistance',10,'SortStr','descend');
    locs1 = sort(locs1)+1;
    [~,locs2] = findpeaks(-diff(aggregateTrace),'NPeaks',1,'MinPeakDistance',10,'SortStr','descend');
    locs2 = sort(locs2)-1;
    disp ' '
    disp([locs1 diff([locs1 locs2]')' locs2])
    imshow(dat(locs1(1):locs2(1)).mean,[])
    hold on
    plot(lastClickedX-cplc/2+plc/2, lastClickedY-cplr/2+plr/2,'o')
    prefix = [datestr(now,'HHMMSS') '_spot'];
    saveas(gcf,fullfile(ramfovpath,[prefix '.fig']))
%     amps = ...
mean(bdat(optimImgIdx,locs1(1):locs2(1))')
max(ans)./ans
%     max(amps)./amps
%%
return
% moviefixsc(squeeze(avgFrames(:,:,3,:)))

% best_xyz_offset =[
%     6.0000    4.0719    0.1665
%     6.7978    1.7590    0.2699
%     ]
%%
% plot(linspace(-1,1,100),polyval(p,linspace(-1,1,100)),optimSteps,bstepsmov(optimImgIdx(1),:))

%% release camera and daq
try stoppreview(vid); end
clear src
clear vid
clear s
imaqreset


