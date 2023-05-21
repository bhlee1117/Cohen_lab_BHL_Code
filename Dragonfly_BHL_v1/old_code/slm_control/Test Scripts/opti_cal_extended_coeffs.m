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
hImage.Parent.CLim = [0 100];
set(hImage,'cdatamapping','scaled')
hold(hImage.Parent,'on')
plot(hImage.Parent,lastClickedX-cplc/2+roi_nc/2, lastClickedY-cplr/2+roi_nr/2,'o')

%% PSF optimization with zernike coefficients
% refocus and avoid photobleaching
% preview(vid);
src.ExposureTime = .5;
fullscreen(rot90(slmph1),3)

%% reset optim params
% original_xyz_offset = mean(best_xyz_offset);
original_xyz_offset = [0 0 0];
original_commonz_offset = zeros(1,14);
original_extended_commonz_offset = [commonz' zeros(1,12)];
optimParamScale = [10 10 5*ones(1,12) 1 1 5 1 1 1 1 1 1 5 1 1];
nParams = numel(optimParamScale);
optimSteps = linspace(-1,1,9);
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
% best_commonz_offset = ones(nSpots,1)*original_commonz_offset;
best_extended_commonz_offset = ones(nSpots,1)*original_extended_commonz_offset;

%% state archival 
storagepath = 'R:\state';
if ~exist(storagepath,'file'), mkdir(storagepath), end
%%
figure windowstyle docked
% vid.FramesPerTrigger = ceil((laserOnTime/1000 + slmSettlingTime)*nSteps/src.ExposureTime*2);
outputSingleScan(s,0) % shut off shutter
% outputSingleScan(s,1)
% alldata = zeros(roi_nr,roi_nc,vid.FramesPerTrigger,nParams);
avgBlurredFrames = zeros(roi_nr,roi_nc,nSteps,nParams);
avgOriginalFrames = zeros(roi_nr,roi_nc,nSteps,nParams);
deltaOptim = zeros(nSpots,nParams);
allSlmPhases = uint8(zeros(1152,1920,nSteps));
optimStepSortList = ceil(((1:nSteps).*(-1).^(1:nSteps))/2) + ceil(nSteps/2);
% optimStepSortList = optimStepSortList(end:-1:1);
recenter = 0;
if recenter
    itParamList = [1 2];
else
    itParamList = [3:26];
%     itParamList = [26]*[1 1 1];
%     itParamList = [itParamList; itParamList; itParamList];
%     itParamList = itParamList(:);
end
for itParam = itParamList % numel(best_commonz_offset) % optimize one parameter each time
    temp_commonz = best_extended_commonz_offset;
    temp_zpol = cat(3,zpolys,zpolys2);
    idxMix = setdiff(1:numel(temp_commonz),[1 2 4 itParam]);
    if recenter
        idxStore = 3;
        idxOptim = itParam;
        idxRemove = 5:numel(temp_commonz);
    else
        idxStore = 3;
        idxOptim = 5;
        idxRemove = 6:numel(temp_commonz);
    end
    temp_temp = temp_zpol(itParam);
    temp_zpol(idxStore) = sum(permute(temp_commonz(idxMix),[1 3 2]).*temp_zpol(idxMix));
    temp_zpol(idxOptim) = temp_temp;
    temp_zpol(idxRemove) = [];
    temp_temp = temp_commonz(itParam);
    temp_commonz(idxStore) = 1;
    temp_commonz(idxOptim) = temp_temp;
    temp_commonz(idxRemove) = [];
    for itStep = 1:nSteps
        local_temp_commonz = temp_commonz;
        local_temp_commonz(idxOptim) = local_temp_commonz(idxOptim) + optimParamScale(itParam).*optimSteps(itStep);
        xyz = slmlocs - original_xyz_offset;
        ttd = [
            0 1 1/3 % 0 sind(60) cosd(60)
            1 0 0
            0 -1/20 1 % cosd(60)/100 sind(60)
            ]*xyz';
        approxSLM = gs(local_temp_commonz', temp_zpol, ttd);
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
    nbdat = vm(squeeze(getdata(vid)));
    bdat = blur(nbdat,2);
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
        avgBlurredFrames(:,:,itStep,itParam) = bdat(locs1(itSortStep):locs2(itSortStep)).mean;
        avgOriginalFrames(:,:,itStep,itParam) = nbdat(locs1(itSortStep):locs2(itSortStep)).mean;
    end
    bstepsmov = vm(avgBlurredFrames(:,:,:,itParam));
    stepsmov = vm(avgOriginalFrames(:,:,:,itParam));
    plotSteps = linspace(optimSteps(1),optimSteps(end),30);
    for itSpot = 1:nSpots
        ihs = 60; % image half side
        mindistance = ihs*2;
        newmov = zeros(2*ihs+1,2*ihs+1,nSteps);

        % correct displacement
        [~, idx] = max(bstepsmov.tovec.data);
        [idxr, idxc] = ind2sub(bstepsmov.imsz,idx);
        dxy = [idxc' idxr'];
        bimhsz = bstepsmov.imsz/2;
        dxy = dxy - bimhsz;
        sts = 1:nSteps;
        fitmat = [sts'*[0 1]+1 sts'.^2];
        bbstepsmov = bstepsmov.blur(1);
        bstepsmov = bstepsmov.apply_displacement(fitmat*(fitmat\dxy));
        stepsmov = stepsmov.apply_displacement(fitmat*(fitmat\dxy));
        for itSortStep = 1:nSteps
%             [idr, idc] = find(findlocs(bstepsmov(itSortStep).mean,mindistance,1));
            newmov(:,:,itSortStep) = stepsmov(bimhsz(1)+(-ihs:ihs),bimhsz(2)+(-ihs:ihs),itSortStep).mean;
        end
        
        if recenter
            objVals = bbstepsmov(optimImgIdx(itSpot),:);
        else
            %% max intensity metric
%             objVals = prctile(bstepsmov(:,:),99.9);
%             plot(objVals)
            %% x, y 2nd moment metric
            imgbg = prctile(newmov(:),95);
            moms = imgmoments(max(newmov-imgbg,0)); 
            objVals = -sqrt(prod(moms(6:7,:)));
%             plot(objVals)

            %% x, y, FWHM metric
%             rszf = 10;
%             allyprofs = squeeze(newmov(ceil(end/2),:,:));
%             allyprofs = imresize(allyprofs,size(allyprofs).*[rszf 1]);
%             allyprofs = allyprofs - min(allyprofs);
%             allyprofs = allyprofs./max(allyprofs);
%             
%             allxprofs = squeeze(newmov(:,ceil(end/2),:));
%             allxprofs = imresize(allxprofs,size(allxprofs).*[rszf 1]);
%             allxprofs = allxprofs - min(allxprofs);
%             allxprofs = allxprofs./max(allxprofs);
%             
%             objVals = -(sum(allyprofs>1/2)+sum(allxprofs>1/2))/2/rszf;
%             plot(objVals)
            %% r 2nd moment metric
%             imgbg = prctile(newmov(:),95);
%             objVals = zeros(1,nSteps);
%             for itSortStep = 1:nSteps
%                 projs = radon(max(newmov(:,:,itSortStep)-imgbg,0),0:10:90);
%                 radr = (1:size(projs,1))';
%                 rcoid = mean(radr); % sum(projs.*radr)./sum(projs)
%                 rstds = sqrt(sum(projs.*(radr-rcoid).^2)./sum(projs));
%                 objVals(itSortStep) = -mean(rstds);
%             end
%             plot(objVals)
            %% r FWHM metric
%             imgbg = prctile(newmov(:),80);
%             objVals = zeros(1,nSteps);
%             for itSortStep = 1:nSteps
%                 projs = radon(max(newmov(:,:,itSortStep)-imgbg,0),0:5:90);
%                 projs = projs - min(projs);
%                 projs = projs./max(projs);
%                 objVals(itSortStep) = -mean(sum(projs>1/2));
%             end
%             plot(objVals)
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
    
    for itSortStep = (nSteps-1)/2+1
        local_temp_commonz = best_extended_commonz_offset;
        local_temp_commonz(itParam) = local_temp_commonz(itParam) + optimParamScale(itParam).*optimSteps(itSortStep);
        params = local_temp_commonz;
        image = newmov(:,:,itSortStep);
        fname = fullfile(storagepath,sprintf('state_%s.mat',datestr(now,'YYYYmmDDHHMMSSFFF')));
        save(fname,'params','image')
    end
    for itSortStep = [-(nSteps-1)/2:-1 1:(nSteps-1)/2]+(nSteps-1)/2+1
        local_temp_commonz = best_extended_commonz_offset;
        local_temp_commonz(itParam) = local_temp_commonz(itParam) + optimParamScale(itParam).*optimSteps(itSortStep);
        params = local_temp_commonz;
        image = newmov(:,:,itSortStep);
        fname = fullfile(storagepath,sprintf('test_%s.mat',datestr(now,'YYYYmmDDHHMMSSFFF')));
        save(fname,'params','image')
    end

    deltaOptim = min(deltaOptim,+1);
    deltaOptim = max(deltaOptim,-1);
    disp ' '
    disp(deltaOptim)
    best_extended_commonz_offset(:,itParam) = best_extended_commonz_offset(:,itParam) + deltaOptim(:,itParam).*optimParamScale(itParam)
    scale_learning_rate = .5;
    scale_learning_amp = 2;
%     if abs(deltaOptim(itSpot,itParam)) > 1 || p(1) > 0
%     else
%     if recenter
%     else
        scalefactor = std(objVals-p*optimSteps.^[2;1;0])./std(p(1:2)*optimSteps.^[2;1]);
        if isfinite(scalefactor)
            scalefactor = min(max(scalefactor,.25),4);
            optimParamScale(itParam) = (1-scale_learning_rate).*optimParamScale(itParam) + ...
                scale_learning_rate.*scale_learning_amp.*optimParamScale(itParam).*scalefactor
        else
%             scalefactor = .1;
%             optimParamScale(itParam) = (1-scale_learning_rate).*optimParamScale(itParam) + ...
%                 scale_learning_rate.*scale_learning_amp.*optimParamScale(itParam).*scalefactor
        end
%     end
%     itParam = itParam + 1;

%     params = best_extended_commonz_offset;
%     image = newmov(:,:,itSortStep);
%     fname = fullfile(storagepath,sprintf('state_%s.mat',datestr(now,'YYYYmmDDHHMMSSFFF')));
%     save(fname,'params','image')
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


