%% camera init
% imaqreset
vid = videoinput('hamamatsu', 1, 'MONO16_2048x2048_SlowMode');
src = getselectedsource(vid);
%% set ROI
% square centered ROI
roi_nr = 512;
roi_nc = roi_nr;
roi_or = (2048-roi_nr)/2;
roi_oc = (2048-roi_nc)/2;
vid.ROIPosition = [roi_or roi_oc roi_nr roi_nc];
hImage = preview(vid);
hImage.CDataMapping = 'scaled';
hImage.Parent.CLim = [0 max(reshape(hImage.CData,[],1))];
% set(hImage,'cdatamapping','scaled')
%%
s = daq.createSession('NI');
addDigitalChannel(s,'Dev1','Port0/Line22','OutputOnly');
addAnalogInputChannel(s,'Dev1','ai0','Voltage');
% outputSingleScan(s,1)
% outputSingleScan(s,0)
% queueOutputData(s,[ones(200,1);zeros(1,1)]);
% s.startForeground;
% clear s

%% PSF optimization with zernike coefficients
% refocus and avoid photobleaching
preview(vid);
src.ExposureTime = .5;
fullscreen(rot90(slmph1),3)
%% set up for auto AQ
fullscreen(0*rot90(slmph1),3)
vid.FramesPerTrigger = 150;
src.ExposureTime = 0.03;
%% 
original_xyz_offset = [0 0 5];
ostep_scale = [5 2 5];
osteps = linspace(-1,1,4);
oidx = sub2ind([plr plc],round(lastClickedY-cplr/2+plr/2),round(lastClickedX-cplc/2+plc/2));

%%
% best_xyz_offset = original_xyz_offset;
% best_xyz_offset = best_xyz_offset - delta_ostep.*ostep_scale;
%%
% this code block executes one iteration of optimization of SLM pattern to
% maximize static fluorescence from a group of point targets. x,y,z
% location is independently varied and sequentially measured in multiple
% steps. The resulting fluorescence for each cell along these steps is
% independently fitted with a 2nd order polynomial for finding new
% estimates of optimized target in x, y, z. Thus all cells are optimized in
% parallel.
% camera, slm, are initialized earlier.
% this block is run multiple times in a row until convergence is seen in
% intermediate results. 
clf
outputSingleScan(s,0)
alldata = zeros(roi_nr,roi_nc,vid.FramesPerTrigger,numel(osteps));
avgdata = zeros(roi_nr,roi_nc,numel(ostep_scale),numel(osteps));
allframes = uint8(zeros(1152,1920,numel(ostep_scale)));
no = numel(osteps);
olist = ceil(((1:no).*(-1).^(1:no))/2) + ceil(no/2);
olist = olist(end:-1:1);
for it_sortstep = 1:numel(osteps) % loops over 4-step line search
    it_steps = olist(it_sortstep);
    for it_coeff = 1:numel(ostep_scale) % loops over x,y,z coordinates separately
        test_xyz_offset = best_xyz_offset;
        test_xyz_offset(:,it_coeff) = best_xyz_offset(:,it_coeff) + ostep_scale(it_coeff).*osteps(it_steps);
        
        xyz = bsxfun(@minus,(slmlocs),test_xyz_offset);
        ttd = [
            0 1 1/3 % 0 sind(60) cosd(60)
            1 0 0
            0 -1/20 1 % cosd(60)/100 sind(60)
            ]*xyz';
        approxSLM = gs(commonz, zpolys, ttd);
        slmph1opt = angle(approxSLM);
        slmph1opt = slmph1opt*128/pi+128;
        slmph1opt = max(min(slmph1opt,255),0);
        slmph1opt = uint8(slmph1opt);
        allframes(:,:,it_coeff) = rot90(slmph1opt);
    end
    start(vid);
    for it_coeff = 1:numel(ostep_scale) % loops over x,y,z coordinates separately
        fullscreen(allframes(:,:,it_coeff),3)
        pause(.150)
        queueOutputData(s,[ones(round(src.ExposureTime*24*s.Rate),1);zeros(1,1)]);
        s.startForeground;
    end
    alldata(:,:,:,it_steps) = squeeze(getdata(vid));
    fprintf('AQ %d%%\n',ceil(it_sortstep/numel(osteps)*100));
    bdat = blur(vm(alldata(:,:,:,it_steps)),2);
    traces = bdat(oidx,:)';
    traces = traces - prctile(traces,40);
    agg_diffs = mean(traces./std(traces),2);
%     agg_diffs = max(tovec(alldata(:,:,:,it_steps)));
    plot(agg_diffs)
    drawnow
    [~,locs1] = findpeaks(+diff(agg_diffs),'NPeaks',3,'MinPeakDistance',10,'SortStr','descend');
    locs1 = sort(locs1)+1;
    [~,locs2] = findpeaks(-diff(agg_diffs),'NPeaks',3,'MinPeakDistance',10,'SortStr','descend');
    locs2 = sort(locs2)-1;
    disp ' '
    disp([locs1 diff([locs1 locs2]')' locs2])
    for it_coeff = 1:numel(ostep_scale)
        coeff_frames = locs1(it_coeff):locs2(it_coeff);
        avgdata(:,:,it_coeff,it_steps) = mean(alldata(:,:,coeff_frames,it_steps),3);
    end
%     break
end
%
% fit params
% assume parabolic optim
% find max, optimal scale
% update 
% moviefixsc(squeeze(avgdata(:,:,3,:)))

delta_ostep = zeros(numel(numel(lastClickedY),ostep_scale));
% oidx = oidx(1:2)
for it_coeff = 1:numel(ostep_scale)
    bstepsmov = blur(vm(squeeze(avgdata(:,:,it_coeff,:))),2);
%     bstepsmov(oidx,:)
%     plot(bstepsmov(oidx(1),:)')
    for it_spot = 1:numel(lastClickedY)
%         plot(osteps,bstepsmov(oidx(it_spot),:))
%         title([it_coeff it_spot])
%         drawnow
        pause(.1)
        p = polyfit(osteps,bstepsmov(oidx(it_spot),:),2);
        delta_ostep(it_spot,it_coeff) = -p(2)/p(1)/2;
    end    
end
delta_ostep = min(delta_ostep,+1);
delta_ostep = max(delta_ostep,-1);
disp ' '
disp(delta_ostep)
best_xyz_offset = best_xyz_offset - delta_ostep.*ostep_scale
%%
moviefixsc(squeeze(avgdata(:,:,3,:)))
%%
plot(linspace(-1,1,100),polyval(p,linspace(-1,1,100)),osteps,bstepsmov(oidx(1),:))
%%
%% release camera and daq
try stoppreview(vid); end
clear src
clear vid
clear s
imaqreset
