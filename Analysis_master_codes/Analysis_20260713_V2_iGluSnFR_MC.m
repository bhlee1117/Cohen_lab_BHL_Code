%% Load file path
clear
clc;
[~, ~, raw] = xlsread(['/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/Prism_V2+Glu_Data_Arrangement.xlsx'], 'Sheet1', 'C5:AA33');

fpath=raw(:,1)';
V2moviemaxTime=15000;
GlumoviemaxTime=3500;
%% Read Voltage movie & ROI-block motion correction
%
%  Strategy:
%   For each voltage segment, instead of running MC on the full frame
%   (which has large black areas between ROIs), we:
%     a) Extract DMD pattern from Device_Data{6} → bounding boxes per ROI
%     b) Loop over each ROI → crop a padded bounding-box sub-movie
%     c) Run optical_flow_motion_correction_LBH_ROIBlock on the crop
%     d) Paste the corrected crop back into the full frame
%     e) Collect shift traces from all ROIs → concatenate into mcTrace
%        [nFrames x (nROI * nShiftCols)]
%
%  pad_px: pixel padding around each DMD ROI bounding box
%
%  Resume behaviour:
%   A segment is skipped when its motion-corrected movie (mc??.bin / mc2??.bin)
%   AND its shift-trace file (mcTrace??.mat / mc2Trace??.mat) are both already
%   on disk and the .bin holds the expected number of frames. A .bin with the
%   wrong size (e.g. an interrupted write) is redone. Set overwrite_mc = true
%   to force re-running everything.

overwrite_mc = false;   % true → ignore existing files and redo motion correction

% frames stored in a saved .bin (uint16, nRow x nCol per frame)
bin_nframe = @(binfile, nrow, ncol) dir(binfile).bytes / (nrow*ncol*2);

% ── Progress log / phone notification (see mc_log.m) ──────────────────
%  Every message is appended to mc_cfg.logfile — tail that file to see how
%  far the run has got. Start / finish / error messages are also pushed to
%  notify_url so they land on your phone.
%
%  For phone alerts: install the free "ntfy" app, subscribe to a topic name
%  that only you know, and paste it below. The topic name IS the password —
%  anyone who knows it can read the messages on the public ntfy.sh server —
%  so pick a random one and keep the messages free of anything sensitive.
%  A Slack incoming-webhook URL also works and is auto-detected.
%  Leave notify_url = '' to log without notifications.
mc_cfg.logfile      = fullfile(getenv('HOME'), 'Documents', ...
    sprintf('MC_progress_%s.log', char(datetime('now','Format','yyyyMMdd_HHmm'))));
mc_cfg.notify_url   = 'https://ntfy.sh/bhl-matlab-c856d939';   % '' to disable
mc_cfg.notify_title = 'MATLAB motion correction';

file_list    = [20:29];
failed_files = [];
mc_log(mc_cfg, 2, 'Motion correction started: %d files (#%d–#%d)', ...
    numel(file_list), min(file_list), max(file_list));
t_batch = tic;

for f=file_list
    [~, fname_disp] = fileparts(fpath{f});   % short label for progress messages
    mc_log(mc_cfg, 2, 'START file#%d  %s', f, fpath{f});
    t_file = tic;
  try
    Devicedata_filename=fullfile(fpath{f},'output_data.mat');
    Vdata_filename=fullfile(fpath{f},'frames1.bin');
    Gludata_filename=fullfile(fpath{f},'frames2.bin');
    load(Devicedata_filename);

    % set orange cam
    nCol1 = double(Device_Data{3}.ROI([2]));
    nRow1 = double(Device_Data{3}.ROI([4]));
    exposuretime1 = Device_Data{3}.exposuretime;

    o_laser = 200001;
    DAQ_rate=Device_Data{1, 2}.Counter_Inputs(1, 1).rate;

    % voltage camera clock
    cam1_vsyn = Device_Data{1, 2}.Counter_Inputs(1, 1).data;
    start_idx = min(find(cam1_vsyn ==max (cam1_vsyn)));
    end_idx = length(Device_Data{1, 2}.buffered_tasks(1, 2).channels(1, 1).data);
    CamTrigger1_DAQax=find((cam1_vsyn(2:end)-cam1_vsyn(1:end-1))>0);
    segment_size=unique(diff(CamTrigger1_DAQax));
    last_val = cam1_vsyn(start_idx - 1);
    n_to_add = end_idx - start_idx + 1;
    n_segments = ceil(n_to_add / segment_size);
    added_part = repelem((last_val + 1 : last_val + n_segments)', segment_size);
    added_part = added_part(1:n_to_add);
    cam1_vsyn(start_idx:end_idx) = added_part;

    % Prepare
    VmovTimesegments=[(cam1_vsyn(o_laser)+2):V2moviemaxTime:cam1_vsyn(end)];
    VmovTimesegments(end+1)=cam1_vsyn(end);
    nFrame2analyze=VmovTimesegments(end)-VmovTimesegments(1);
    CamTrigger1_DAQaxVec=ind2vec(length(cam1_vsyn),CamTrigger1_DAQax(1):segment_size:cam1_vsyn(end)*segment_size,1);
    CamTrigger1_DAQaxVec=CamTrigger1_DAQaxVec(o_laser:end);
    FirstFrameDAQax=find(CamTrigger1_DAQaxVec>0,1);
    t_vol=[FirstFrameDAQax:segment_size:(FirstFrameDAQax+segment_size*nFrame2analyze)]/DAQ_rate;

    % laser artifacts
    if Device_Data{1, 2}.buffered_tasks(1, 3).channels.name=="Photodiode";
    PD_tr = Device_Data{1, 2}.buffered_tasks(1, 3).channels.data;
    v_cam_vsyn_diff = diff(cam1_vsyn)';
    v_cam_frame_end = find(v_cam_vsyn_diff == 1);
    v_cam_frame_length = mean(diff(v_cam_frame_end));
    v_cam_frame_start = v_cam_frame_end - v_cam_frame_length + 1;
    binned_PD_tr = zeros(length(v_cam_frame_end), 1);
    for ii = 1:length(binned_PD_tr)
        binned_PD_tr(ii) = mean(PD_tr(v_cam_frame_start(ii):v_cam_frame_end(ii)));
    end
    binned_PD_tr = binned_PD_tr(cam1_vsyn(o_laser)+2 :cam1_vsyn(o_laser)+nFrame2analyze+1);
    binned_PD_tr_norm = binned_PD_tr / mean(binned_PD_tr, 1);
    else
        disp('No photodiode channel detected')
        binned_PD_tr_norm=ones(length(t_vol),1);
    end

    % --- Extract DMD pattern → bounding boxes in voltage camera space ---
    % Based on get_blueDMDPatt but targeting the voltage (orange) camera:
    %   - cam index 1  → Device_Data{1,3} is voltage cam (flash), ROI from Device_Data{1,4}
    %   - cam index 2  → Device_Data{1,3} is glu cam (fusion), ROI from Device_Data{1,3}
    % Here we want voltage camera space, so use Device_Data{3}.ROI for the crop.

    inverseTform = invert(Device_Data{1, 6}.tform);
    rawPattern   = double(Device_Data{1, 6}.Target);   % [DMDrows x DMDcols x Npoly]

    % Sensor size for voltage camera
    try
        sensorsize = size(Device_Data{1, 3}.testimage);   % flash cam
    catch
        sensorsize = Device_Data{1, 3}.virtualSensorSize; % fusion cam
    end
    Rfixed_volt = imref2d([sensorsize sensorsize]);

    % Warp DMD pattern into voltage camera space
    dmdpatt= imwarp(rawPattern, inverseTform, 'OutputView', Rfixed_volt);
    bwdmdpatt=bwlabel(dmdpatt);
    revertedImage=zeros(size(dmdpatt,1),size(dmdpatt,2),max(bwdmdpatt(:)));
    for p=1:max(bwdmdpatt(:))
        revertedImage(:,:,p)=bwdmdpatt==p;
    end

    % Crop to voltage camera ROI  (same crop as used when reading frames1.bin)
    ROI_volt = double(Device_Data{3}.ROI([1 3 2 4]));   % [x0 y0 width height]

    Npoly = size(revertedImage, 3);
    pad_px = 20;

    bb = struct('r1',cell(Npoly,1),'r2',cell(Npoly,1),'c1',cell(Npoly,1),'c2',cell(Npoly,1));
    for p = 1:Npoly
        patt_p = imcrop(revertedImage(:,:,p), ROI_volt + [0 0 -1 -1]);
        patt_p = patt_p > 0.5;

        if any(patt_p(:))
            stats  = regionprops(patt_p, 'BoundingBox');
            % BoundingBox: [xmin ymin width height] → convert to row/col
            bbox   = stats(1).BoundingBox;
            c1_raw = bbox(1);
            r1_raw = bbox(2);
            c2_raw = bbox(1) + bbox(3);
            r2_raw = bbox(2) + bbox(4);
        else
            % fallback: use full frame
            c1_raw = 1;  r1_raw = 1;
            c2_raw = nCol1; r2_raw = nRow1;
        end

        bb(p).r1 = max(1,     floor(r1_raw) - pad_px);
        bb(p).r2 = min(nRow1, ceil(r2_raw)  + pad_px);
        bb(p).c1 = max(1,     floor(c1_raw) - pad_px);
        bb(p).c2 = min(nCol1, ceil(c2_raw)  + pad_px);
    end

    mov_ref = [];   % built from first segment

    for vs=1:length(VmovTimesegments)-1
        frm2read=[VmovTimesegments(vs):VmovTimesegments(vs+1)-1];

        % ── Skip segments that are already motion corrected ───────────
        mcbin_file   = [fpath{f} '/mc'      num2str(vs,'%02d') '.bin'];
        mctrace_file = [fpath{f} '/mcTrace' num2str(vs,'%02d') '.mat'];
        if ~overwrite_mc && exist(mcbin_file,'file')==2 && exist(mctrace_file,'file')==2
            nframe_disk = bin_nframe(mcbin_file, nRow1, nCol1);
            if abs(nframe_disk - length(frm2read)) <= 2
                mc_log(mc_cfg, 0, '[file#%d %s] V-seg %d/%d already motion corrected — skipping', ...
                    f, fname_disp, vs, length(VmovTimesegments)-1);
                if isempty(mov_ref)
                    % reference for the remaining segments: mean image of the
                    % already-corrected segment
                    S = load(mctrace_file, 'ave_im');
                    mov_ref = S.ave_im;
                end
                continue
            else
                mc_log(mc_cfg, 0, ['[file#%d %s] V-seg %d/%d: existing mc%02d.bin holds ' ...
                    '%g frames (expected %d) — redoing'], f, fname_disp, vs, ...
                    length(VmovTimesegments)-1, vs, nframe_disk, length(frm2read));
            end
        end

        [mov1, ~] = readBinMov_times(Vdata_filename, nRow1, nCol1, frm2read);
        mov1 = mov1 - 100;

        % spinning disk artifact correction
        period = 12;
        ground_truth   = mean(mov1,3);
        artifact_image = zeros(size(mov1,1), size(mov1,2), period, 'double');
        for jj = 1:period
            artifact_image(:,:,jj) = mean(mov1(:,:,jj:period:end),3);
        end
        artifact_pattern = artifact_image ./ ground_truth;
        corrected_img = zeros(size(mov1), 'double');
        for kk = 1:period
            corrected_img(:,:,kk:period:end) = double(mov1(:,:,kk:period:end)) ./ artifact_pattern(:,:,kk);
        end
        corrected_img(isnan(corrected_img)) = 0;
        corrected_img(isinf(corrected_img)) = 0;

        binned_PD_tr_norm_segment=binned_PD_tr_norm(frm2read-VmovTimesegments(1)+1);
        corrected_img = double(corrected_img) ./ reshape(binned_PD_tr_norm_segment, 1, 1, []);

        % build reference from first segment
        if isempty(mov_ref)
            mov_ref = mean(corrected_img, 3);
        end

        nT      = size(corrected_img, 3);
        mov_out = corrected_img;    % patched in-place per ROI
        mcTrace_all = [];           % [nT x (Npoly * nShiftCols)]

        % ── ROI-block MC ──────────────────────────────────────────────
        for p = 1:Npoly
            mc_log(mc_cfg, 0, '[file#%d %s] V-seg %d/%d  ROI %d/%d  running...', ...
                f, fname_disp, vs, length(VmovTimesegments)-1, p, Npoly);

            r1 = bb(p).r1;  r2 = bb(p).r2;
            c1 = bb(p).c1;  c2 = bb(p).c2;

            sub_mov = corrected_img(r1:r2, c1:c2, :);
            sub_ref = mov_ref(r1:r2, c1:c2);

            sub_mov_vm = vm(sub_mov);
            % evalc swallows the tic/toc and NoRMCorre progress printing
            t_roi = tic;
            [mc_out_log, sub_mc, xyField_p] = evalc([ ...
                'optical_flow_motion_correction_LBH(' ...
                'sub_mov_vm, sub_ref, ''normcorre'')']);
            mc_log(mc_cfg, 0, '[file#%d %s] V-seg %d/%d  ROI %d/%d  done (%.1f s)', ...
                f, fname_disp, vs, length(VmovTimesegments)-1, p, Npoly, toc(t_roi));

            % only surface the captured log when a fallback path was taken
            if contains(mc_out_log,'Error during motion correction') || ...
                    contains(mc_out_log,'Size not matched')
                fprintf('    NoRMCorre fallback triggered:\n%s\n', mc_out_log);
            end

            xyField_p=xyField_p.xymean;

            % paste corrected patch back into full frame
            mov_out(r1:r2, c1:c2, :) = sub_mc;

            % pad shift trace to nT rows if needed
            nT_p = size(xyField_p, 1);
            if nT_p < nT
                xyField_p = [xyField_p; repmat(xyField_p(end,:), nT-nT_p, 1)]; %#ok<AGROW>
            elseif nT_p > nT
                xyField_p = xyField_p(1:nT, :);
            end

            mcTrace_all = [mcTrace_all, xyField_p]; %#ok<AGROW>
        end
        % mcTrace_all: [nT x (Npoly * nShiftCols)]
        % column layout per ROI block: [x-shifts | y-shifts]

        % ── Save ──────────────────────────────────────────────────────
        ave_im     = mean(mov_out, 3);
        mov_out_vm = vm(mov_out);
        mov_out_vm.transpose.savebin([fpath{f} '/mc' num2str(vs,'%02d') '.bin']);

        mcTrace = mcTrace_all;
        save([fpath{f} '/mcTrace' num2str(vs,'%02d') '.mat'], 'mcTrace', 'ave_im');
        mc_log(mc_cfg, 0, '[file#%d %s] saved V-segment %2.0f out of %2.0f', ...
            f, fname_disp, vs, length(VmovTimesegments)-1);

        % ── QC visualisation ──────────────────────────────────────────
        % Figure 10: per-ROI crops — reference | MC mean | difference
        % Figure 11: x/y shift traces per ROI
        nShiftCols = size(mcTrace_all,2) / Npoly;
        t_seg      = (0:nT-1) * exposuretime1;

        figure(10); clf;
        tiledlayout(3, Npoly, 'TileSpacing','compact','Padding','compact');
        for p = 1:Npoly
            r1 = bb(p).r1; r2 = bb(p).r2;
            c1 = bb(p).c1; c2 = bb(p).c2;

            sub_ref_p = mov_ref(r1:r2, c1:c2);
            sub_mc_p  = ave_im(r1:r2, c1:c2);
            diff_p    = sub_ref_p - sub_mc_p;

            nexttile(p);         imshow2(mat2gray(sub_ref_p)); title(sprintf('ROI %d  ref',    p));
            nexttile(Npoly+p);   imshow2(mat2gray(sub_mc_p));  title(sprintf('ROI %d  MC mean',p));
            nexttile(2*Npoly+p); imshow2(mat2gray(diff_p));    title(sprintf('ROI %d  diff',   p));
        end
        sgtitle(sprintf('Segment %d/%d — images', vs, length(VmovTimesegments)-1));

        figure(11); clf;
        tiledlayout(Npoly, 2, 'TileSpacing','compact','Padding','compact');
        for p = 1:Npoly
            col_x = (p-1)*nShiftCols + 1;
            col_y = (p-1)*nShiftCols + 2;

            nexttile;
            plot(t_seg, mcTrace_all(:, col_x));
            ylabel('x-shift (px)'); xlabel('Time (s)');
            title(sprintf('ROI %d  x-shift', p));

            nexttile;
            plot(t_seg, mcTrace_all(:, col_y));
            ylabel('y-shift (px)'); xlabel('Time (s)');
            title(sprintf('ROI %d  y-shift', p));
        end
        sgtitle(sprintf('Segment %d/%d — shift traces', vs, length(VmovTimesegments)-1));
        drawnow;
    end

    clear mov_mc corrected_img mov1

    % Glutamate image motion correction

    cam2_vsyn = Device_Data{1, 2}.Counter_Inputs(1, 2).data;
    cam2_vsyn_trig = find (diff(cam2_vsyn) == 1)+1;
    segment_size = cam2_vsyn_trig (10) - cam2_vsyn_trig(9);
    start_idx = min(find(cam2_vsyn ==max (cam2_vsyn)));
    end_idx = length(Device_Data{1, 2}.buffered_tasks(1, 2).channels(1, 2).data);
    last_val = cam2_vsyn(start_idx - 1);
    n_to_add = end_idx - start_idx + 1;
    n_segments = ceil(n_to_add / segment_size);
    added_part = repelem((last_val + 1 : last_val + n_segments)', segment_size);
    added_part = added_part(1:n_to_add);
    cam2_vsyn(start_idx:end_idx) = added_part;

    mod488 = 200001;
    mod607 = 200001;
    cam2_trig = find (diff(cam2_vsyn) ==1)+1;
    cam2_trig = cam1_vsyn (cam2_trig);
    cam2_trig = cam2_trig (cam2_vsyn(mod488)+2:cam2_vsyn(end))'- (cam1_vsyn(mod607)+2);
    if cam2_vsyn(end)<GlumoviemaxTime
        GluemovTimesegments=[(cam2_vsyn(mod488)+2) cam2_vsyn(end)];
    else
        GluemovTimesegments=[(cam2_vsyn(mod488)+2):GlumoviemaxTime:cam2_vsyn(end)];
        GluemovTimesegments(end+1)=cam2_vsyn(end);
    end

    Blue_Nframe = Device_Data{4}.frames_requested;
    nCol2 = double(Device_Data{4}.ROI([2]));  % ROI on the camera
    nRow2 = double(Device_Data{4}.ROI([4]));  % ROI on the camera
    exposuretime2 = Device_Data{4}.exposuretime;

    mov_ref2 = [];   % built from the first segment of this file

    for vs2=1:length(GluemovTimesegments)-1
        frm2read2=[GluemovTimesegments(vs2):GluemovTimesegments(vs2+1)-1];

        % ── Skip segments that are already motion corrected ───────────
        mc2bin_file   = [fpath{f} '/mc2'      num2str(vs2,'%02d') '.bin'];
        mc2trace_file = [fpath{f} '/mc2Trace' num2str(vs2,'%02d') '.mat'];
        if ~overwrite_mc && exist(mc2bin_file,'file')==2 && exist(mc2trace_file,'file')==2
            nframe_disk2 = bin_nframe(mc2bin_file, nRow2, nCol2);
            if abs(nframe_disk2 - length(frm2read2)) <= 2
                mc_log(mc_cfg, 0, '[file#%d %s] Glu movie %1.0f/%1.0f already motion corrected — skipping', ...
                    f, fname_disp, vs2, length(GluemovTimesegments)-1);
                if isempty(mov_ref2)
                    S2 = load(mc2trace_file, 'ave_im2');
                    mov_ref2 = S2.ave_im2;
                end
                continue
            else
                mc_log(mc_cfg, 0, ['[file#%d %s] Glu movie %1.0f/%1.0f: existing mc2%02d.bin holds ' ...
                    '%g frames (expected %d) — redoing'], f, fname_disp, vs2, ...
                    length(GluemovTimesegments)-1, vs2, nframe_disk2, length(frm2read2));
            end
        end

        [mov2] = double(readBinMov_times(Gludata_filename, nRow2, nCol2,frm2read2));
        mov2 = mov2 - 100;             % baseline
        if isempty(mov_ref2)
            mov_ref2=mean(mov2(:,:,1:500),3);
        end
        mov2=vm(mov2);
        mc_log(mc_cfg, 0, '[file#%d %s] Glu movie %d/%d  running...', ...
            f, fname_disp, vs2, length(GluemovTimesegments)-1);
        % evalc swallows the tic/toc and NoRMCorre progress printing
        t_glu = tic;
        [mc_out_log2, mov2_mc, xyField] = evalc([ ...
            'optical_flow_motion_correction_LBH(' ...
            'mov2, mov_ref2, ''normcorre'')']);
        mc_log(mc_cfg, 0, '[file#%d %s] Glu movie %d/%d  done (%.1f s)', ...
            f, fname_disp, vs2, length(GluemovTimesegments)-1, toc(t_glu));

        if contains(mc_out_log2,'Error during motion correction') || ...
                contains(mc_out_log2,'Size not matched')
            fprintf('    NoRMCorre fallback triggered:\n%s\n', mc_out_log2);
        end

        ave_im2=mean(mov2_mc,3);
        mov2_mc=vm(mov2_mc);
        mov2_mc.transpose.savebin([fpath{f} '/mc2' num2str(vs2,'%02d') '.bin'])

        mcTrace2=xyField; % Normcorre
        save([fpath{f} '/mc2Trace' num2str(vs2,'%02d') '.mat'],'mcTrace2','ave_im2')
        mc_log(mc_cfg, 0, '[file#%d %s] saved Glu movie %1.0f out of %1.0f', ...
            f, fname_disp, vs2, length(GluemovTimesegments)-1);
    end

    mc_log(mc_cfg, 2, 'DONE file#%d (%s) in %.1f min', ...
        f, fname_disp, toc(t_file)/60);

  catch ME
    % keep going with the remaining files, but make the failure loud
    failed_files(end+1) = f; %#ok<SAGROW>
    mc_log(mc_cfg, 3, 'ERROR file#%d (%s): %s', f, fname_disp, ME.message);
    mc_log(mc_cfg, 0, 'stack: %s', ...
        strjoin(arrayfun(@(s) sprintf('%s:%d', s.name, s.line), ...
        ME.stack, 'UniformOutput', false), ' < '));
  end   % try/catch around one file
end

if isempty(failed_files)
    mc_log(mc_cfg, 2, 'ALL DONE: %d files in %.1f min, no errors', ...
        numel(file_list), toc(t_batch)/60);
else
    mc_log(mc_cfg, 3, 'FINISHED with errors: %d/%d files failed (#%s) in %.1f min', ...
        numel(failed_files), numel(file_list), ...
        strjoin(string(failed_files), ', #'), toc(t_batch)/60);
end
fprintf('Progress log: %s\n', mc_cfg.logfile);

%% See Glu movie

f=20
Devicedata_filename=fullfile(fpath{f},'output_data.mat');
Vdata_filename=fullfile(fpath{f},'frames1.bin');
Gludata_filename=fullfile(fpath{f},'frames2.bin');
load(Devicedata_filename);

% set orange cam
nCol1 = double(Device_Data{3}.ROI([2]));  % ROI on the camera
nRow1 = double(Device_Data{3}.ROI([4]));  % ROI on the camera
exposuretime1 = Device_Data{3}.exposuretime;
% ---------- Correct spinning disk artifact ----------
% Discard first some second of frames without 607 ex light, drop last 1s to avoid sync issues
o_laser = 200001; % all the recording started 2s after the acquisition start

% voltage camera clock
cam1_vsyn = Device_Data{1, 2}.Counter_Inputs(1, 1).data;
start_idx = min(find(cam1_vsyn ==max (cam1_vsyn)));
end_idx = length(Device_Data{1, 2}.buffered_tasks(1, 2).channels(1, 1).data);
segment_size = 103*floor(exposuretime1/0.001);
last_val = cam1_vsyn(start_idx - 1);
n_to_add = end_idx - start_idx + 1;
n_segments = ceil(n_to_add / segment_size);
added_part = repelem((last_val + 1 : last_val + n_segments)', segment_size);
added_part = added_part(1:n_to_add);
cam1_vsyn(start_idx:end_idx) = added_part;
VmovTimesegments=[(cam1_vsyn(o_laser)+2):V2moviemaxTime:cam1_vsyn(end)];
nFrame2analyze=VmovTimesegments(end)-VmovTimesegments(1)+1;
t_vol = (0:nFrame2analyze-1)*exposuretime1;

cam2_vsyn = Device_Data{1, 2}.Counter_Inputs(1, 2).data;
cam2_vsyn_trig = find (diff(cam2_vsyn) == 1)+1;
segment_size = cam2_vsyn_trig (10) - cam2_vsyn_trig(9);
start_idx = min(find(cam2_vsyn ==max (cam2_vsyn)));
end_idx = length(Device_Data{1, 2}.buffered_tasks(1, 2).channels(1, 2).data);
last_val = cam2_vsyn(start_idx - 1);
n_to_add = end_idx - start_idx + 1;
n_segments = ceil(n_to_add / segment_size);
added_part = repelem((last_val + 1 : last_val + n_segments)', segment_size);
added_part = added_part(1:n_to_add);
cam2_vsyn(start_idx:end_idx) = added_part;

mod488 = 200001;
mod607 = 200001;
cam2_trig = find (diff(cam2_vsyn) ==1)+1;
cam2_trig = cam1_vsyn (cam2_trig);
cam2_trig = cam2_trig (cam2_vsyn(mod488)+2:cam2_vsyn(end))'- (cam1_vsyn(mod607)+2);
if cam2_vsyn(end)<GlumoviemaxTime
    GluemovTimesegments=[(cam2_vsyn(mod488)+2) cam2_vsyn(end)];
else
    GluemovTimesegments=[(cam2_vsyn(mod488)+2):GlumoviemaxTime:cam2_vsyn(end)];
    GluemovTimesegments(end+1)=cam2_vsyn(end);
end

Blue_Nframe = Device_Data{4}.frames_requested;
nCol2 = double(Device_Data{4}.ROI([2]));  % ROI on the camera
nRow2 = double(Device_Data{4}.ROI([4]));  % ROI on the camera
exposuretime2 = Device_Data{4}.exposuretime;
t_cal = (cam2_trig-1)*exposuretime2;

mov2_mc=[];
for vs2=1:length(GluemovTimesegments)-1
    fprintf('Reading %1.0f out of %1.0f movie...\n', vs2,length(GluemovTimesegments)-1)
    try
        mov_read_tmp=double(readBinMov([fpath{f} '/mc2' num2str(vs2,'%02d') '.bin'], nRow2, nCol2));
    catch
        mov_read_tmp=double(readBinMov_times([fpath{f} '/mc2' num2str(vs2,'%02d') '.bin'], nRow2, nCol2,[1:(GluemovTimesegments(vs2+1)-GluemovTimesegments(vs2)+1)]));
    end
    mov2_mc = cat(3,mov2_mc,mov_read_tmp);
end

mov2_mc_filt=imgaussfilt(mov2_mc,1);
mov2_sub=mov2_mc_filt-movmedian(mov2_mc_filt,100,3);
[V, D, u_sub] = get_eigvector(tovec(mov2_sub)',40);
% figure(2); clf; imshow2_patch(reshape(u_sub,nRow2(1),nCol2(1),[]))
% mov2_sub=SeeResiduals(mov2_sub,V(:,[3 4]));

figure(1);
moviefixsc(mov2_mc)