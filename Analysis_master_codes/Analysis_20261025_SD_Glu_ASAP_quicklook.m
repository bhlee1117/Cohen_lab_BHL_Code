%% Quick look: ASAP7y + iGluSnFR4f culture, spinning disk, blue camera (frames2.bin)
%  Single-camera recording: everything is on the blue cam (Device_Data{4}, 400 x 600 px,
%  12 ms exposure); frames1.bin is only one 2048 x 2048 snap from the orange cam.
%
%  Spinning-disk artifact: the disk scan beats against the rolling shutter, which the Orca
%  reads out from the sensor centre (row 1024 = row 200 of this ROI). This gives horizontal
%  bands that drift through the frame as a chevron about row 200, with a drift speed that
%  changes over the recording. The bands are ADDITIVE (~0.5 counts on the background and on
%  the soma alike, per-pixel coefficient ~1 everywhere) and uniform along each row.
%  Because the drift is not stationary, the period-folded correction of the orange camera
%  (period = 12, Analysis_20260713_V2_iGluSnFR_MC.m) removes nothing here at any period.
%  Instead, per frame and per row, the band is measured on background pixels (cells
%  excluded, so their activity is not subtracted) and subtracted from the whole row.
%  Caveat: a field-wide signal that also lights up the background pixels (e.g. a network
%  burst on dim neurites) would partly be taken as band and removed.
%  Motion correction (NoRMCorre, as for mc2XX.bin) runs after the band correction, because the
%  bands are fixed to sensor rows and would be shifted with the image otherwise.

clear; clc;
fpath = '/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/20261025/223941YD043_SCE_iGluSNFR4f_L23';
camOffset = 100;      % camera dark offset
bgPrctile = 60;       % background = pixels below this percentile of the smoothed mean image
hpWin_s   = 5;        % s, slow part of the row bands kept in the movie (bleaching)
rowSmooth = 5;        % rows, gaussian window smoothing the band profile across rows (noise)
doMC      = true;     % NoRMCorre motion correction after the band correction (~5 min)

%% Load the blue camera movie
load(fullfile(fpath,'output_data.mat'), 'Device_Data');
nCol2 = double(Device_Data{4}.ROI(2));
nRow2 = double(Device_Data{4}.ROI(4));
exposuretime2 = Device_Data{4}.exposuretime;

mov2 = double(readBinMov(fullfile(fpath,'frames2.bin'), nRow2, nCol2)) - camOffset;
meanF = squeeze(mean(mov2,[1 2]));
endFr = find(meanF > 0.5*median(meanF), 1, 'last') - 1;    % blue laser goes off before the camera stops
mov2 = mov2(:,:,2:endFr);
nT = size(mov2,3);
t_glu = (0:nT-1)*exposuretime2;
fprintf('%d x %d px, %d frames (%.1f s at %.1f Hz)\n', nRow2, nCol2, nT, nT*exposuretime2, 1/exposuretime2);

%% Spinning-disk row-band artifact correction
im2 = mean(mov2,3);
imS = imgaussfilt(im2,2);
bgMask = imerode(imS < prctile(imS(:), bgPrctile), strel('disk',2));

rowBG = squeeze(sum(mov2.*bgMask,2)) ./ max(sum(bgMask,2),1);      % nRow2 x nT, background row means
band  = rowBG - movmedian(rowBG, round(hpWin_s/exposuretime2), 2); % fast part = bands
band  = smoothdata(band, 1, 'gaussian', rowSmooth);
mov2_ac = mov2 - reshape(band, nRow2, 1, nT);

%% QC: row kymographs and mean trace, raw vs corrected
% Row means of the fast fluctuations (0.5 s moving mean removed); bands show up as stripes.
fastRow = @(M) squeeze(mean(M,2)) - movmean(squeeze(mean(M,2)), 41, 2);
kRaw = fastRow(mov2);   kAc = fastRow(mov2_ac);
fRaw = meanF(2:endFr) - movmean(meanF(2:endFr), 41);
fAc  = squeeze(mean(mov2_ac,[1 2]));  fAc = fAc - movmean(fAc, 41);
showFr = 1:min(400, nT);
cl = [-1 1]*prctile(abs(kRaw(:)), 99);

figure(1); clf;
tiledlayout(4,3,'TileSpacing','compact','Padding','compact');
nexttile; imshow2(im2, []); title('Mean image');
nexttile; imshow2(bgMask, []); title(sprintf('Background pixels (< %d prctile)', bgPrctile));
nexttile; plot(std(band,0,2), 1:nRow2); axis ij tight; xlabel('Band std (counts)'); ylabel('Row'); title('Band amplitude per row');
nexttile([1 3]); imagesc(showFr, 1:nRow2, kRaw(:,showFr), cl); colormap(gca, turbo); colorbar;
title('Row means of the fast fluctuation, raw'); ylabel('Row');
nexttile([1 3]); imagesc(showFr, 1:nRow2, kAc(:,showFr), cl); colormap(gca, turbo); colorbar;
title('Row means of the fast fluctuation, corrected'); ylabel('Row'); xlabel('Frame');
nexttile([1 3]); plot(t_glu, fRaw, t_glu, fAc); xlim(t_glu([1 end])); legend('raw','corrected');
xlabel('Time (s)'); ylabel('Frame mean - 0.5 s mean (counts)'); title('Whole-frame flicker');

%% Motion correction (NoRMCorre on the band-corrected movie, as in Analysis_20260713_V2_iGluSnFR_MC.m)
if doMC
    mov_ref2 = mean(mov2_ac(:,:,1:500),3);
    % evalc swallows the tic/toc and NoRMCorre progress printing
    [mcLog, mov2_mc, xyField] = evalc('optical_flow_motion_correction_LBH(vm(mov2_ac), mov_ref2, ''normcorre'')');
    if contains(mcLog,'Error during motion correction') || contains(mcLog,'Size not matched')
        fprintf('NoRMCorre fallback triggered:\n%s\n', mcLog);
    end
    mov2_mc  = double(mov2_mc);
    mcTrace2 = xyField.xymean;                    % nT x [x shifts | y shifts], one column per patch
    corrLabel = 'band + motion corrected';
else
    mov2_mc  = mov2_ac;
    mcTrace2 = zeros(nT,2);
    corrLabel = 'band corrected';
end
nPatch = size(mcTrace2,2)/2;
fprintf('Motion: x %.2f to %.2f px, y %.2f to %.2f px (%d patches)\n', ...
    min(mcTrace2(:,1:nPatch),[],'all'), max(mcTrace2(:,1:nPatch),[],'all'), ...
    min(mcTrace2(:,nPatch+1:end),[],'all'), max(mcTrace2(:,nPatch+1:end),[],'all'), nPatch);

figure(2); clf;
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
nexttile; plot(t_glu, mcTrace2(:,1:nPatch)); xlim(t_glu([1 end])); ylabel('x shift (px)');
title(sprintf('NoRMCorre shifts, %d patches', nPatch));
nexttile; plot(t_glu, mcTrace2(:,nPatch+1:end)); xlim(t_glu([1 end])); ylabel('y shift (px)'); xlabel('Time (s)');

%% Quick look 1: raw | band + motion corrected | corrected - movmedian (100 frames), fixed scale
frm = 1:min(1000, nT);                       % frames to show
mov2_sub = mov2_mc(:,:,frm) - movmedian(mov2_mc(:,:,frm), 100, 3);
lim = prctile(im2(:), [1 99.8]);
% mov2_sub sits around 0: its +-3 SD is mapped onto lim so the three panels share one
% colour scale (the colorbar applies to the first two panels only)
subSD = std(mov2_sub(:));
subDisp = mean(lim) + mov2_sub*diff(lim)/(6*subSD);
fprintf('Quick look 1, right panel: full scale = +-%.2f counts (3 SD)\n', 3*subSD);
figure(3); clf;
moviefixsc([mov2(:,:,frm) mov2_mc(:,:,frm) subDisp], lim)

%% Save Quick look 1 as mp4 with time stamp
% raw and corrected share lim; mov2_sub has its own symmetric scale subLim
subLim  = [-3 3]*subSD;
vidFps  = 30;                                % playback rate (recording is 83 Hz)
gap     = 4;                                 % px, white line between panels
toGray  = @(M, L) min(max((M - L(1))/diff(L), 0), 1);
sep     = ones(nRow2, gap);
x0      = [0 1 2]*(nCol2 + gap) + 5;         % label positions of the three panels
vidFile = fullfile(fpath, 'Quicklook1_raw_corrected_movmedsub.mp4');
v = VideoWriter(vidFile, 'MPEG-4'); v.FrameRate = vidFps; v.Quality = 90; open(v);
for j = 1:numel(frm)
    I = [toGray(mov2(:,:,frm(j)), lim) sep toGray(mov2_mc(:,:,frm(j)), lim) sep toGray(mov2_sub(:,:,j), subLim)];
    I = insertText(repmat(I,1,1,3), [x0' [5;5;5]], {sprintf('raw, %.0f-%.0f counts', lim), ...
        sprintf('%s, %.0f-%.0f counts', corrLabel, lim), ...
        sprintf('corrected - movmedian(100), %c%.0f counts', char(177), subLim(2))}, ...
        'FontSize', 16, 'BoxColor', 'black', 'BoxOpacity', 0.6, 'TextColor', 'yellow');
    I = insertText(I, [5 nRow2-32], sprintf('t = %.3f s (frame %d)', t_glu(frm(j)), frm(j)), ...
        'FontSize', 18, 'BoxColor', 'black', 'BoxOpacity', 0.6, 'TextColor', 'white');
    writeVideo(v, I);
end
close(v);
fprintf('Saved %s (%d frames, %.0f fps = %.2fx real time)\n', vidFile, numel(frm), vidFps, vidFps*exposuretime2);

%% Quick look 2: fast fluctuations (0.5 s moving mean removed, gauss 1.5 px), raw | band + motion corrected
% The bands are ~0.5 counts, invisible on the full intensity scale; here they dominate the raw movie.
hp = @(M) imgaussfilt(M - movmean(M, 41, 3), 1.5);
dF_raw = hp(mov2(:,:,frm));
dF_corr = hp(mov2_mc(:,:,frm));
figure(4); clf;
moviefixsc([dF_raw dF_corr], [-3 3]*std(dF_corr(:)))   % +-3 SD of the corrected movie
