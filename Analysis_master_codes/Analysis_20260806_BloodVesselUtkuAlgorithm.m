%% Load file path
clear
clc;
[~, ~, raw] = xlsread(['/Volumes/cohen_lab/Lab/Labmembers' ...
    '/Byung Hun Lee/Data/PrismPCdata_Arrangement.xlsx'], 'Sheet1', 'C5:AA31');

ref_ROI=cellfun(@(x) (str2num(num2str(x))),raw(:,9),'UniformOutput',false);
oblique_ROI=cellfun(@(x) (str2num(num2str(x))),raw(:,10),'UniformOutput',false);
PeriSoma_ROI=cellfun(@(x) (str2num(num2str(x))),raw(:,11),'UniformOutput',false);
basal_ROI=cellfun(@(x) (str2num(num2str(x))),raw(:,12),'UniformOutput',false);
apical_ROI=cellfun(@(x) (str2num(num2str(x))),raw(:,13),'UniformOutput',false);
distal_ROI=cellfun(@(x) (str2num(num2str(x))),raw(:,14),'UniformOutput',false);

fpath=raw(:,1)';
StructureData=raw(:,8);
BadROI=cellfun(@(x) (str2num(num2str(x))),raw(:,17),'UniformOutput',false);
EndFrame=cell2mat(raw(:,15));
ifmotionReject=cell2mat(raw(:,16));
ifdirtRemov=cell2mat(raw(:,18));
Pixelsize=cell2mat(raw(:,6));
NSeesawComponent=cell2mat(raw(:,25));
save_figto='/Volumes/BHL_WD18TB/PP72_PlaceCellResults';
place_bin=150; time_segment=15000; overlap=200;
alignedMovFN = {'STA_Mat_SS','STA_Mat_CS','STA_Mat_dSP'};
bound=6;
title_str={'Basal','Apical','Peri-Soma'};
PlaceFieldList=cellfun(@(x) (str2num(num2str(x))),raw(:,21),'UniformOutput',false);
PlaceFieldBin=cellfun(@(x) (str2num(num2str(x))),raw(:,22),'UniformOutput',false);
set(0,'DefaultFigureWindowStyle','docked')
%foi=[1 4 5 6 8 10 11 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27];
foi=[1 4 5 6 8 10 11 15 16 17 18 19 20 21 22 23 24 25 26 27];

%% ================================================================
%  WHOLE-MOVIE BLOOD-VESSEL CORRECTION  (mc frames 1:358000)
%  Runs bleach-dr on the full motion-corrected movie, chunk by chunk,
%  and writes the blood-corrected movie back out as segmented .bin
%  files (single precision, same chunking as the mc_ShutterReg files).
%
%  Reads mc chunks with loadVoltChunk (photobleach norm + motion
%  residual -> spike-positive residual movie), builds the rank1/spike
%  image from find(Result.spike(1,:)), hands each chunk to
%  run_correction.py, and reads the corrected chunk back.
% ================================================================

% ---- config ----------------------------------------------------
f          = 20;                                  % dataset index (as elsewhere)
venvpy     = '/Users/bhlee1117/Documents/GitHub/bleach-dr/.venv/bin/python';
script     = '/Users/bhlee1117/Documents/GitHub/bleach-dr/run_correction.py';

fpath_mc   = '/Volumes/BHL18TB_D1/20240807/201405BHLm141_N1_VR_LowStim'; % full mc set: mc_ShutterReg01..40.bin + mcTrace01..40.mat
mc_prefix  = 'mc_ShutterReg';                     % movie-file prefix for loadVoltChunk
out_dir    = fpath{f};                            % where corrected chunks + traces are written  (== fpath{20})
scratch_dir= fpath_mc;                            % where the per-chunk in/out HDF5 lives (use a fast/local disk)
corr_prefix= 'mc_BVcorr';                         % output corrected-chunk prefix  -> mc_BVcorr01.bin, ...

Nframes    = 358000;                              % analysed voltage movie length (== mc frames 1:358000)
seg_len    = 15000;                               % frames per mc chunk file

% blood-correction parameters (set these from your sweep_parameters.py run)
bp         = [200 220];    % bandpass [low high] Hz
preset     = 1;            % strategy_preset 0/1/2
r1method   = 'skip';       % 'skip' | 'ols' | 'kalman'
frameRate  = 1000;         % Hz
spikeTargetCount = 1000;   % stop STA accumulation once this many spike frames are pooled

saveFullMovie = true;      % write corrected .bin chunks (~142 GB single). false -> traces only
saveTraces    = true;      % also save corrected traces (corrected projected onto Result.ftprnt)
qc_frames  = 78000:83000;  % QC window: for the chunk(s) covering it, also save the blood component
                           % and a [input | corrected | blood] comparison mp4

% ---- geometry & metadata --------------------------------------
load(fullfile(out_dir,'output_data.mat'),'Device_Data');
sz = double(Device_Data{1,3}.ROI([2 4]));         % sz(1)=nCol(552), sz(2)=nRow(180)
load(fullfile(out_dir,'PC_Result.mat'),'Result'); % Result.spike (roi x T), Result.ftprnt (nRow x nCol x roi)

% ---- segmentation for mc frames 1:Nframes ---------------------
nSeg          = ceil(Nframes/seg_len);            % 24
readFrames    = repmat(seg_len,1,nSeg);
readFrames(end) = Nframes - seg_len*(nSeg-1);     % last chunk partial (8400)
segStart      = [0 cumsum(readFrames)];           % global frame offset per chunk (segStart(j)+1 : segStart(j+1))
qc_segs       = unique(min(max(ceil(qc_frames/seg_len),1),nSeg));  % chunk(s) covering the QC window

spikeIdx = find(Result.spike(1,:));               % spike frames (global 1:Nframes)
spikeIdx = spikeIdx(spikeIdx>=1 & spikeIdx<=Nframes);

in_h5  = fullfile(scratch_dir,'bv_chunk_in.h5');
out_h5 = fullfile(scratch_dir,'bv_chunk_out.h5');

%% ---- PASS 1: build spike (rank1) image = STA of -mov_res -----
fprintf('[BV] Pass 1: building spike image from %d spike frames...\n', numel(spikeIdx));
spikeSum = zeros(sz(2), sz(1));    % (nRow, nCol) = (180, 552)
spikeCnt = 0;
for j = 1:7%nSeg
    movR = loadVoltChunk(fpath_mc, j, sz, readFrames(j), mc_prefix);   % (180,552,nT)
    nT   = size(movR,3);
    loc  = spikeIdx(spikeIdx>segStart(j) & spikeIdx<=segStart(j)+nT) - segStart(j);
    if ~isempty(loc)
        spikeSum = spikeSum + sum(movR(:,:,loc),3);
        spikeCnt = spikeCnt + numel(loc);
    end
    fprintf('  seg %2d/%2d  pooled spikes = %d\n', j, nSeg, spikeCnt);
    clear movR
    if spikeCnt >= spikeTargetCount, break; end
end
assert(spikeCnt>0, 'No spike frames found in Result.spike(1,:) within 1:%d', Nframes);
spikeImg = spikeSum / spikeCnt;                    % (180,552)
spike_out = permute(single(spikeImg), [2 1]);      % -> h5 reads (Ly,Lx)=(180,552)
figure; imagesc(spikeImg); axis image; colormap(gray); title('Rank1 / spike image (STA)');

%% ---- PASS 2: correct each chunk, stream results back ---------
% Split each footprint into odd / even pixels (checkerboard, row+col parity) so
% the corrected movie yields two traces per ROI with independent shot noise but
% the same underlying signal -- the "two channels" used later for a noise-free
% (cross-channel) diagonal in the correlation matrix.
nROIall  = size(Result.ftprnt,3);
[Hh, Ww] = size(Result.ref_im);
[ccg, rrg] = meshgrid(1:Ww, 1:Hh);
oddPix   = mod(rrg+ccg,2)==1;                       % (H x W) checkerboard
ftprnt_vec  = tovec(Result.ftprnt);                % (H*W x nROI) full footprint
ftprnt_odd  = tovec(Result.ftprnt .*  oddPix);     % odd-pixel footprint
ftprnt_even = tovec(Result.ftprnt .* ~oddPix);     % even-pixel footprint
if saveTraces
    tr_corrected     = zeros(Nframes, nROIall);    % combined (all pixels)
    tr_corrected_ch1 = zeros(Nframes, nROIall);    % odd pixels
    tr_corrected_ch2 = zeros(Nframes, nROIall);    % even pixels
end

t_all = tic;
for j = 1:nSeg
    tseg = tic;
    fprintf('[BV] Pass 2: chunk %2d/%2d ...\n', j, nSeg);
    movR = loadVoltChunk(fpath_mc, j, sz, readFrames(j), mc_prefix);   % (180,552,nT), spike-positive residual
    nT   = size(movR,3);

    % --- write chunk + spike image to scratch HDF5 (Python reads (T,Ly,Lx)) ---
    mov_out = permute(single(movR), [2 1 3]);      % (552,180,nT) -> h5 (nT,180,552)
    if exist(in_h5,'file'), delete(in_h5); end
    h5create(in_h5,'/mov',        size(mov_out),  'Datatype','single'); h5write(in_h5,'/mov',        mov_out);
    h5create(in_h5,'/rank1_image',size(spike_out),'Datatype','single'); h5write(in_h5,'/rank1_image',spike_out);
    clear mov_out

    % --- run the Python correction (single pass; chunk already <= seg_len) ---
    % keep the blood component only for the QC chunk(s); --no-blood elsewhere to save time/disk
    if ismember(j, qc_segs), bloodflag = ''; else, bloodflag = '--no-blood'; end
    if exist(out_h5,'file'), delete(out_h5); end
    cmd = sprintf(['"%s" "%s" --input "%s" --output "%s" ' ...
        '--frame-rate %g --bandpass-low %g --bandpass-high %g ' ...
        '--block-size 20 --max-components 30 --strategy-preset %d ' ...
        '--rank1-method %s --device cpu %s'], ...
        venvpy, script, in_h5, out_h5, frameRate, bp(1), bp(2), preset, r1method, bloodflag);
    [status, out] = system(cmd, '-echo');
    assert(status==0, 'Python failed on chunk %d:\n%s', j, out);

    % --- read corrected chunk back, undo permute -> (180,552,nT) ---
    corrected = permute(h5read(out_h5,'/corrected'), [2 1 3]);

    % --- QC: save blood component + [input | corrected | blood] mp4 for the QC chunk ---
    if ismember(j, qc_segs)
        blood = permute(h5read(out_h5,'/blood'), [2 1 3]);
        bfid = fopen(fullfile(out_dir, sprintf('mc_BVblood%02d.bin', j)),'w');
        fwrite(bfid, single(blood), 'single', 'l'); fclose(bfid);   % read back: readBinMovSingle(fn, sz(2), sz(1))

        locqc = qc_frames - segStart(j);
        locqc = locqc(locqc>=1 & locqc<=nT);
        locqc = locqc(1:min(1000,numel(locqc)));      % cap mp4 length (and montage memory)
        if ~isempty(locqc)
            mont   = [movR(:,:,locqc) corrected(:,:,locqc) blood(:,:,locqc)];  % [input | corrected | blood]
            qc_mp4 = fullfile(out_dir, sprintf('BVcorr_QC_seg%02d.mp4', j));
            writeMov(qc_mp4, mont, [], 10, 1, []);   % auto scale
            fprintf('   QC seg %d: saved blood chunk + %s\n', j, qc_mp4);
        end
        clear blood mont
    end

    % --- save corrected chunk as single-precision .bin (read back with readBinMovSingle(fn, sz(2), sz(1))) ---
    if saveFullMovie
        outbin = fullfile(out_dir, sprintf('%s%02d.bin', corr_prefix, j));
        fid = fopen(outbin,'w');
        fwrite(fid, single(corrected), 'single', 'l');   % column-major (nRow,nCol,nT) matches readBinMovSingle
        fclose(fid);
    end

    % --- corrected traces (projection onto footprints; combined + odd/even channels) ---
    if saveTraces
        gidx = (segStart(j)+1):(segStart(j)+nT);
        cvec = tovec(corrected)';                                % (nT x pixels)
        tr_corrected(gidx,:)     = cvec * ftprnt_vec;            % combined
        tr_corrected_ch1(gidx,:) = cvec * ftprnt_odd;           % odd pixels
        tr_corrected_ch2(gidx,:) = cvec * ftprnt_even;          % even pixels
    end

    clear movR corrected
    fprintf('   done in %.1f s\n', toc(tseg));
end
if exist(in_h5,'file'),  delete(in_h5);  end
if exist(out_h5,'file'), delete(out_h5); end

if saveTraces
    save(fullfile(out_dir,'BVcorr_traces.mat'), ...
        'tr_corrected', 'tr_corrected_ch1', 'tr_corrected_ch2', ...
        'bp', 'preset', 'r1method', '-v7.3');
end
fprintf('[BV] All %d chunks done in %.1f min. Output dir: %s\n', nSeg, toc(t_all)/60, out_dir);
% To reload a corrected chunk:  mov = readBinMovSingle(fullfile(out_dir,'mc_BVcorr05.bin'), sz(2), sz(1));

%% ================================================================
%  EXTRACT TRACES FROM ALREADY-SAVED corrected chunks (no re-run)
%  Use this when the blood-vessel correction was already run and the
%  mc_BVcorr##.bin chunks exist. Reads each chunk, projects onto the
%  footprints (combined + odd/even-pixel channels), and saves
%  BVcorr_traces.mat -- identical output to Pass 2, without any Python.
% ================================================================
f          = 20;
out_dir    = fpath{f};                 % folder holding mc_BVcorr##.bin + output_data.mat + PC_Result.mat
corr_prefix= 'mc_BVcorr';

load(fullfile(out_dir,'output_data.mat'),'Device_Data');
sz = double(Device_Data{1,3}.ROI([2 4]));      % sz(1)=nCol(552), sz(2)=nRow(180)
load(fullfile(out_dir,'PC_Result.mat'),'Result');

% odd/even (checkerboard) footprints -- SAME split as Pass 2
nROIall  = size(Result.ftprnt,3);
[Hh, Ww] = size(Result.ref_im);
[ccg, rrg] = meshgrid(1:Ww, 1:Hh);
oddPix   = mod(rrg+ccg,2)==1;
ftprnt_vec  = tovec(Result.ftprnt);            % (H*W x nROI)
ftprnt_odd  = tovec(Result.ftprnt .*  oddPix);
ftprnt_even = tovec(Result.ftprnt .* ~oddPix);

% discover corrected chunks in segment order
D = dir(fullfile(out_dir, [corr_prefix '*.bin']));
assert(~isempty(D), 'No %s*.bin found in %s', corr_prefix, out_dir);
[~, si] = sort({D.name});  D = D(si);

bytesPerFrame = sz(1)*sz(2)*4;                  % single precision (4 bytes)
Ntot = round(sum([D.bytes]) / bytesPerFrame);
tr_corrected     = zeros(Ntot, nROIall);
tr_corrected_ch1 = zeros(Ntot, nROIall);
tr_corrected_ch2 = zeros(Ntot, nROIall);

g0 = 0;  t_all = tic;
for k = 1:numel(D)
    corrected = readBinMovSingle(fullfile(out_dir, D(k).name), sz(2), sz(1));  % (nRow,nCol,nT)
    nT   = size(corrected,3);
    cmat = tovec(corrected);                    % pixels x nT
    idx  = g0 + (1:nT);
    tr_corrected(idx,:)     = (ftprnt_vec'  * cmat)';   % combined
    tr_corrected_ch1(idx,:) = (ftprnt_odd'  * cmat)';   % odd pixels
    tr_corrected_ch2(idx,:) = (ftprnt_even' * cmat)';   % even pixels
    g0 = g0 + nT;
    fprintf('  extracted %-18s %6d frames (total %d)\n', D(k).name, nT, g0);
    clear corrected cmat
end

% trim in case preallocation overshot
tr_corrected     = tr_corrected(1:g0,:);
tr_corrected_ch1 = tr_corrected_ch1(1:g0,:);
tr_corrected_ch2 = tr_corrected_ch2(1:g0,:);

save(fullfile(out_dir,'BVcorr_traces.mat'), ...
    'tr_corrected', 'tr_corrected_ch1', 'tr_corrected_ch2', '-v7.3');
fprintf('[BV] Traces extracted from %d saved chunks (%d frames) in %.1f min -> %s\n', ...
    numel(D), g0, toc(t_all)/60, fullfile(out_dir,'BVcorr_traces.mat'));
