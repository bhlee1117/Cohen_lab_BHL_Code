
clear
clc;
cd '/Volumes/BHL18TB_D2/Arranged_Data/Prism_OptopatchResult';
[~, ~, raw] = xlsread(['/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/Prism_V2+Glu_Data_Arrangement.xlsx'], 'Sheet1', 'C5:AA31');
load('/Volumes/BHL18TB_D2/20260203_SD_V2+iGluSNFR4/25X_transformationMatrix.mat');
fpath=raw(:,1)';
V2moviemaxTime=15000;
GlumoviemaxTime=5000;
set(0,'DefaultFigureWindowStyle','docked')

%% ROI setting
for f=[13]
    fprintf(['Reading file#%2.0f, %s\n'],f,fpath{f})
    Devicedata_filename=fullfile(fpath{f},'output_data.mat');
    load(Devicedata_filename);

    sz = double(Device_Data{3}.ROI([2 4]));  % ROI on the camera
    Result = struct();
    ref_time=[4500:7000];
    load(fullfile(fpath{f},'mcTrace01.mat'));
    mov_mc=double(readBinMov_times([fpath{f} '/mc' num2str(1,'%02d') '.bin'], sz(2), sz(1),ref_time));
    avgImg=mean(mov_mc,3);

    [u,s,v] = svds(tovec(mov_mc-mean(mov_mc,3)),20);
    reshape_u=reshape(u,sz(2),sz(1),[]);
    Result.bvMask=[];
    [~, Result.bvMask]=get_ROI(max(abs(reshape_u),[],3),Result.bvMask);

    Result.ref_im=mean(mov_mc,3);
    mov_res= mov_mc-mean(mov_mc,3);
    if isstruct(mcTrace)
        mcTrace=movmean(mcTrace.xymean,3,1);
    end
    mov_res = SeeResiduals(mov_res,mcTrace(ref_time,:));
    mov_res = SeeResiduals(mov_res,mcTrace(ref_time,:).^2);
    mov_res = SeeResiduals(mov_res,mcTrace(ref_time,1).*mcTrace(ref_time,end));
    mov_res = mov_res.*double(max(Result.bvMask,[],3)==0);

    [~, ~, icsTrace]=clickyICA(imresize(mov_res,0.5),imresize(mean(mov_mc,3),0.5),10);
    icsImgs=toimg(tovec(mov_res)*icsTrace.intens',size(mov_res,1),size(mov_res,2));

    ROIimg=mean(mov_mc,3);
    %ROIimg=mean(icsImgs(:,:,[1]),3);
    excludeImg=mean(icsImgs(:,:,[end]),3);
    %excludeImg=mean(mov_mc,3);

    figure(3); clf;
    imshow2(ROIimg,[])
    title('set extraction ROI');
    % set extraction ROI
    g=1; ROIpoly=[];
    while g
        h = drawpolygon('Color','r');
        if size(h.Position,1)==1 %no more ROI
            g=0;
        else
            ROIpoly=[ROIpoly; {h.Position}];
            hold all
            plot(h.Position(:,1),h.Position(:,2))
        end
    end
    close(figure(3));
    Result.ROIpoly=ROIpoly;

    % set exclude ROI
    figure; excludeROI=[];
    [~, Result.excludeROI]=get_ROI(excludeImg,excludeROI,'exclude ROI');

    n_comp=6;
    mov_res_reg=mov_res;
    %mov_res_reg=SeeResiduals(mov_res_reg,icsTrace.intens([3 4 5],:));
    mov_filt=imgaussfilt3(mov_res_reg.*double(max(Result.bvMask,[],3)==0).*double(max(Result.excludeROI,[],3)==0),[1.5 1.5 0.1]);
    %mov_filt=mov_filt(:,:,1:500);
    movVec=tovec(mov_filt);
    Npoly=size(Result.ROIpoly,1);
    ftprnt = zeros(size(mov_filt,1)*size(mov_filt,2),Npoly);
    clear mask
    figure(4); clf;
    for p=1:Npoly %each ROIs
        p
        clf; ax2=[];
        tiledlayout(n_comp/2+2,2)
        mask(:,:,p) = poly2mask(Result.ROIpoly{p}(:,1), Result.ROIpoly{p}(:,2), sz(2), sz(1));
        pixelList=find(tovec(squeeze(mask(:,:,p))));
        subMov = movVec(pixelList,:);
        [V D eigTrace]=get_eigvector(subMov);
        nexttile([2 2])
        plot(rescale2(eigTrace(:,1:n_comp),1)+[1:n_comp])

        for n=1:n_comp
            eigImg=NaN(size(mov_filt,1)*size(mov_filt,2),1);
            ax2=[ax2 nexttile([1 1])];
            eigImg(pixelList,1)=V(:,n);
            eigImg=toimg(eigImg,size(mov_filt,1),size(mov_filt,2));
            imshow2(im_merge(cat(3,mat2gray(Result.ref_im),(eigImg)*10),[1 1 1;1 0 0]),[])
            title([num2str(n) ', Fraction: ' num2str(D(n)/sum(D),2)])
        end
        linkaxes(ax2,'xy')
        n_take = input('#components to take: ', 's');
        n_take = str2num(n_take);
        coeff = V(:, n_take) * D(n_take);
        ftprnt(pixelList,p)=coeff;
    end
    close(figure(4));
    Result.ftprnt=toimg(ftprnt,sz(2),sz(1));
    ref_ftprnt=Result.ftprnt;

    % %use ICA/PCA image
    % for p=1:Npoly
    %     mask(:,:,p) = poly2mask(Result.ROIpoly{p}(:,1), Result.ROIpoly{p}(:,2), sz(2), sz(1));
    %     pixelList=find(tovec(squeeze(mask(:,:,p))));
    %     ftprnt(pixelList,p)=-ROIimg(pixelList);
    % end
    % ref_ftprnt=toimg(ftprnt,sz(2),sz(1));
    % Result.ftprnt=ref_ftprnt;

    figure(99); clf;
    show_footprnt(Result.ftprnt,Result.ref_im);
    save(fullfile(fpath{f},'Volt_Result.mat'),'Result','-v7.3')
    disp(['Footprint saved at ' fullfile(fpath{f},'Volt_Result.mat')]);
end

%% Signal extraction from multiple movie files, in streaming mode
bound = 6;
for f = [9 10 12 13]

    time_segment = 15000;
    Result=importdata([fpath{f} '/Volt_Result.mat']);
    Devicedata_filename = fullfile(fpath{f}, 'output_data.mat');
    load(Devicedata_filename);

    sz = double(Device_Data{3}.ROI([2 4]));  % ROI on the camera
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

    %-- Build checkerboard masks for each footprint
    % Checkerboard pattern: odd pixels = (row+col) is even, even pixels = (row+col) is odd
    [rows, cols] = ndgrid(1:sz(2), 1:sz(1));
    checker_mask = mod(rows + cols, 2);  % 0 = "odd" pixels, 1 = "even" pixels

    Npoly = size(Result.ftprnt, 3);
    ftprnt_checker = cell(1, 2);
    for p = 1:Npoly
        fp = Result.ftprnt(:, :, p);
        ftprnt_checker{1}(:, :, p) = fp .* double(checker_mask == 0);  % odd pixels
        ftprnt_checker{2}(:, :, p) = fp .* double(checker_mask == 1);  % even pixels
    end

    Result.traces         = [];
    Result.traces_checker = {[], []};
    Result.mc             = [];
    Result.bvTrace=[];

    for j = 1:length(VmovTimesegments)-1
        load([fpath{f} '/mcTrace' num2str(j, '%02d') '.mat']);
        if isstruct(mcTrace)
            mcTrace=mcTrace.xymean;
        end
        motionTrace = movmean(mcTrace, 5, 1);
        Result.mc   = [Result.mc; motionTrace];

        try
            readFrame = VmovTimesegments(j+1) - VmovTimesegments(j);
            mov_mc = double(readBinMov_times([fpath{f} '/mc' num2str(j, '%02d') '.bin'], sz(2), sz(1), [1:readFrame]));
            %disp('readBinMov_times succeeded')
        catch
            mov_mc = double(readBinMov([fpath{f} '/mc' num2str(j, '%02d') '.bin'], sz(2), sz(1)));
            %disp('Movie loaded via readBinMov')
        end

        %-- Motion correction and artefact regression
        bv_trace = tovec(mov_mc)' * tovec(Result.bvMask);
        V_trace  = tovec(mov_mc)' * tovec(Result.ftprnt);
        bv_trace = squeeze(SeeResiduals(permute(bv_trace, [2 3 1]), V_trace));
        if size(bv_trace,2)~=size(mov_mc,3)
            bv_trace=bv_trace';
        end

        mov_res = mov_mc - mean(mov_mc, 3);
        Result.bvTrace=[Result.bvTrace bv_trace];

        %-- Extract main traces
        Trace2add = -(tovec(mov_mc)' * tovec(Result.ftprnt))';

        %-- Extract checkerboard traces (odd and even pixel subsets)
        for ch = 1:2
            Trace2add_ch{ch} = -(tovec(mov_mc)' * tovec(ftprnt_checker{ch}))';
        end

        %-- Stitch segments with offset correction
        if ~isempty(Result.traces)
            offset = mean(Result.traces(:, end-20:end), 2, 'omitnan') - ...
                mean(Trace2add(:, 1:20), 2, 'omitnan');
            offset_ch{1} = mean(Result.traces_checker{1}(:, end-20:end), 2, 'omitnan') - ...
                mean(Trace2add_ch{1}(:, 1:20), 2, 'omitnan');
            offset_ch{2} = mean(Result.traces_checker{2}(:, end-20:end), 2, 'omitnan') - ...
                mean(Trace2add_ch{2}(:, 1:20), 2, 'omitnan');
        else
            offset     = zeros(Npoly, 1);
            offset_ch  = {zeros(Npoly, 1), zeros(Npoly, 1)};
        end

        Result.traces            = [Result.traces,            Trace2add     + offset];
        Result.traces_checker{1} = [Result.traces_checker{1}, Trace2add_ch{1} + offset_ch{1}];
        Result.traces_checker{2} = [Result.traces_checker{2}, Trace2add_ch{2} + offset_ch{2}];
        Result.t_ax= t_vol;

        fprintf('Trace extracting from movie (%2.0d/%2.0d)...\n',j,length(VmovTimesegments)-1);
    end

    %-- Diagnostic figure
    figure(2*f); clf; tiledlayout(2, 2)
    nexttile([1 1])
    show_footprnt_contour(Result.ftprnt, Result.ref_im);
    nexttile([1 1])
    show_footprnt_contour(Result.bvMask, Result.ref_im)
    nexttile([1 2])
    plot(rescale2(Result.traces, 2)' + [1:size(Result.traces, 1)])
    colormap(turbo(256));
    drawnow;

    save([fpath{f} '/Volt_Result.mat'], 'Result', '-v7.3')
    disp([num2str(f) ' th file is saved.'])
end
%% Signal processing and bleaching/motion correction
exclude_frq  = [241.7 242];  % Monitor frequency to exclude
exclude_frq2 = [55.5  56];   % Motion frequency to exclude
time_bin = 15000;
Fs       = 1000;
nOverlap  = 10;

% Design Butterworth filters
[b,  a ] = butter(4, exclude_frq  / (Fs/2), 'stop');    % Notch: monitor
[b2, a2] = butter(4, exclude_frq2 / (Fs/2), 'stop');    % Notch: motion
[b3, a3] = butter(4, 2             / (Fs/2), 'low');     % Low-pass: sub-threshold
[b4, a4] = butter(4, [5 11]        / (Fs/2), 'bandpass'); % Band-pass: theta

for f = [8 9 10 12 13]

    %-- Load data
    load(fullfile(fpath{f}, 'Volt_Result.mat'), 'Result');
    try
        ref_trace = ref_ROI{f}(1);
    catch
        ref_trace = 1;
    end
    nTime = size(Result.traces, 2);
    nROI  = size(Result.traces, 1);

    %-- Detect spikes on normalized traces (for silent period mask only)
    tr = Result.traces - movmedian(Result.traces, 50, 2);
    tr = tr ./ get_threshold(tr, 1);
    sp = find_spike_bh(tr, 4, 2);

    %-- Build silent period mask (exclude [-100, +250] samples around spikes)
    SilentTau    = -100:250;
    SilentPeriod = ones(1, nTime);
    sp_time_all  = find(max(sp, [], 1))';
    valid_sp     = sum((sp_time_all + SilentTau) < 1 | ...
        (sp_time_all + SilentTau) > nTime, 2) == 0;
    SilentPeriod(sp_time_all(valid_sp) + SilentTau) = NaN;
    t_silent = find(~isnan(SilentPeriod));

    %-- Fit exponential decay on reference ROI (for bleaching normalization only)
    ref_sum = sum(Result.traces(ref_trace, t_silent), 1, 'omitnan');
    [y_fit2, ~, ~] = expfitDM_2(t_silent', -ref_sum', (1:nTime)', [100000 1000]);

    %-- Normalize bleaching: divide all traces by exponential fit before motion regression
    bleach_norm = abs(y_fit2(:)');
    bleach_norm = bleach_norm / mean(bleach_norm);        % rescale so mean = 1

    tr_res            = Result.traces                       ./ bleach_norm;
    tr_res_checker{1} = Result.traces_checker{1}(:,1:nTime) ./ bleach_norm;
    tr_res_checker{2} = Result.traces_checker{2}(:,1:nTime) ./ bleach_norm;

    %-- Plot bleaching fit overview
    figure(f); clf;
    ax3 = nexttile;
    plot(Result.traces(ref_trace, :)); hold on;
    plot(-y_fit2, 'r');
    legend({'Raw', 'Exp fit'});
    yyaxis right;
    plot(tr_res(ref_trace, :), 'k');

    %-- Per-ROI processing
    tN = [1:time_bin:nTime, nTime];
    lwpass_fit2       = NaN(nROI, nTime);
    lwpass_fit_ch2{1} = NaN(nROI, nTime);
    lwpass_fit_ch2{2} = NaN(nROI, nTime);
    mcTrace = squeeze(Result.mc)';

    clear norm_trace norm_trace_check

    for n = 1:nROI
        tr       = tr_res(n, 1:nTime)-median(tr_res(n, 1:nTime),2);
        tr_check = {tr_res_checker{1}(n, 1:nTime)-median(tr_res_checker{1}(n, 1:nTime),2), tr_res_checker{2}(n, 1:nTime)-median(tr_res_checker{2}(n, 1:nTime),2)};

        %-- Regress out motion artefacts in time bins (with overlap + blended average)
        tr_regressed       = zeros(1, nTime);
        tr_weight          = zeros(1, nTime);
        tr_check_regressed = {zeros(1, nTime), zeros(1, nTime)};
        tr_check_weight    = {zeros(1, nTime), zeros(1, nTime)};

        for t = 1:length(tN)-1

            % Core segment
            seg_start = tN(t);
            seg_end   = tN(t+1);

            % Padded segment (clamped to valid range)
            pad_start = max(1,     seg_start - nOverlap);
            pad_end   = min(nTime, seg_end   + nOverlap);
            seg_pad   = pad_start : pad_end;
            nPad      = length(seg_pad);

            mc          = mcTrace(:, seg_pad);
            vesseltrace = Result.bvTrace(:, seg_pad);

            % Regress on padded segment
            tr_pad = squeeze(SeeResiduals(reshape(tr(seg_pad),  1, 1, []), mc));
            tr_pad = squeeze(SeeResiduals(reshape(tr_pad,       1, 1, []), mc.^2));
            tr_pad = squeeze(SeeResiduals(reshape(tr_pad,       1, 1, []), mc(1,:) .* mc(end,:)));
            tr_pad = squeeze(SeeResiduals(reshape(tr_pad,       1, 1, []), vesseltrace));

            for ch = 1:2
                tr_check_pad = squeeze(SeeResiduals(reshape(tr_check{ch}(seg_pad), 1, 1, []), mc));
                tr_check_pad = squeeze(SeeResiduals(reshape(tr_check_pad,          1, 1, []), mc.^2));
                tr_check_pad = squeeze(SeeResiduals(reshape(tr_check_pad,          1, 1, []), mc(1,:) .* mc(end,:)));
                tr_check_pad = squeeze(SeeResiduals(reshape(tr_check_pad,          1, 1, []), vesseltrace));

                % Build linear blend weights for this padded segment
                w = ones(1, nPad);
                left_len  = seg_start - pad_start;   % number of left overlap frames
                right_len = pad_end   - seg_end;      % number of right overlap frames
                if left_len  > 0, w(1:left_len)             = linspace(0, 1, left_len);  end
                if right_len > 0, w(end-right_len+1:end)    = linspace(1, 0, right_len); end

                tr_check_regressed{ch}(seg_pad) = tr_check_regressed{ch}(seg_pad) + tr_check_pad' .* w;
                tr_check_weight{ch}(seg_pad)    = tr_check_weight{ch}(seg_pad)    + w;
            end

            % Build linear blend weights for this padded segment
            w = ones(1, nPad);
            if left_len  > 0, w(1:left_len)          = linspace(0, 1, left_len);  end
            if right_len > 0, w(end-right_len+1:end) = linspace(1, 0, right_len); end

            tr_regressed(seg_pad) = tr_regressed(seg_pad) + tr_pad' .* w;
            tr_weight(seg_pad)    = tr_weight(seg_pad)    + w;
        end

        % Normalize by accumulated weights to get the blended average
        tr = tr_regressed ./ tr_weight;
        for ch = 1:2
            tr_check{ch} = tr_check_regressed{ch} ./ tr_check_weight{ch};
        end
        tr_mc = tr;

        %-- Notch-filter monitor and motion frequencies
        tr             = filtfilt(b,  a,  tr);
        tr             = filtfilt(b2, a2, tr);
        tr_mc_filtered = tr;

        for ch = 1:2
            tr_check{ch} = filtfilt(b,  a,  tr_check{ch});
            tr_check{ch} = filtfilt(b2, a2, tr_check{ch});
        end

        %-- Second-pass slow baseline removal (per ROI)
        norm_trace(n,:)          = tr;
        norm_trace_check{1}(n,:) = tr_check{1};
        norm_trace_check{2}(n,:) = tr_check{2};

        % lwpass_fit2(n, t_silent)       = norm_trace(n, t_silent);
        % lwpass_fit2(n,:)               = movmedian(lwpass_fit2(n,:), 20000, 2, 'omitnan');
        % lwpass_fit_ch2{1}(n, t_silent) = norm_trace_check{1}(n, t_silent);
        % lwpass_fit_ch2{1}(n,:)         = movmedian(lwpass_fit_ch2{1}(n,:), 20000, 2, 'omitnan');
        % lwpass_fit_ch2{2}(n, t_silent) = norm_trace_check{2}(n, t_silent);
        % lwpass_fit_ch2{2}(n,:)         = movmedian(lwpass_fit_ch2{2}(n,:), 20000, 2, 'omitnan');

        % norm_trace(n,:)          = norm_trace(n,:)          - lwpass_fit2(n,:);
        % norm_trace_check{1}(n,:) = norm_trace_check{1}(n,:) - lwpass_fit_ch2{1}(n,:);
        % norm_trace_check{2}(n,:) = norm_trace_check{2}(n,:) - lwpass_fit_ch2{2}(n,:);

        %-- Diagnostic plot for reference ROI
        if n == ref_trace
            nexttile([1 1]);
            plot(rescale2([Result.traces(n,1:nTime); tr_mc; tr_mc_filtered; mcTrace(1,1:nTime)], 2)');
            legend({'Trace', 'MC regressed', 'Monitor filtered', 'Motion trace'});
        end
    end

    %-- Store normalized traces
    Result.normTraces          = norm_trace;
    Result.norm_trace_check{1} = norm_trace_check{1};
    Result.norm_trace_check{2} = norm_trace_check{2};

    figure(2*f); clf;
    plot(rescale2(Result.normTraces,2)'+[1:nROI])

    %-- Save results
    tic;
    save(fullfile(fpath{f}, 'Volt_Result.mat'), 'Result', '-v7.3');
    disp(['Saved,' fullfile(fpath{f}, 'Volt_Result.mat')])
    toc;
end



%% Find Spikes
CS_thres=[5 1];
Sp_thres=[4 2.5];
exclude_frq2=[55.5 56]; %motion
freq_lowhigh2=exclude_frq2/(1000/2);
[b2, a2] = butter(4, freq_lowhigh2, 'stop');

for i=[84]%length(fpath)]
    load([fpath{i} '/OP_Result.mat'])
    Result.normTraces=Result.traces-prctile(Result.traces,25,2);
    Result.normTraces=Result.normTraces./get_threshold(Result.normTraces,1);
    for n=1:size(Result.ftprnt,3)
        Result.normTraces(n,:) = filtfilt(b2, a2, Result.normTraces(n,:));
    end
    tr_ref=Result.normTraces([3],:);
    tr_ref_mean=mean(tr_ref,1,'omitnan');
    tr=Result.normTraces;
    [nROI nTime]=size(Result.normTraces);
    sp=find_spike_bh(tr-movmedian(tr,50,2),Sp_thres(1),Sp_thres(2));
    sp_soma=find_spike_bh(tr_ref-movprc(tr_ref,200,30,2),Sp_thres(1),Sp_thres(2));

    sp_ref_tmp=sum(sp_soma,1)>0;
    sp_ref_1st_ind=find((sp_ref_tmp(2:end)-sp_ref_tmp(1:end-1))==1)+1;
    sp_ref_sumcon=zeros(1,size(sp_soma,2));
    sp_ref_sumcon(1,sp_ref_1st_ind)=sum(sum(reshape(sp_soma(:,sp_ref_1st_ind'+[-1:0]),size(sp_soma,1),[],2),3),1);
    sp_ref=sp_ref_sumcon;

    [~, shift]=max(reshape(tr_ref_mean(find(sp_ref>=1)+[-1:0]'),2,[]),[],1);
    shift=shift-2;
    sp_time_Soma = find(sp_ref>=1)+shift;
    sp_soma=zeros(1,size(tr,2));
    sp_soma(sp_time_Soma)=1;
    sp_soma=[0 (sp_soma(2:end)-sp_soma(1:end-1))==1]; %remove consecutive spikes

    tr_sub=mean(tr_ref,1)-movprc(mean(tr_ref,1),200,20,2);
    tr_sub=get_subthreshold(tr_sub,sp_soma,7,17);

    [trans tr_trace]=detect_transient2(tr_sub,CS_thres,sp_soma,20);
    if isempty(trans.amp)
        CS_trace=zeros(1,nTime);
    else
        transcand=cell2mat(cellfun(@(x) length(x)>1,trans.ISI,'UniformOutput',false));
        meanISI_frnt=cellfun(@(x) mean(x(1:1)),trans.ISI(transcand));
        meanISI_first3=zeros(1,length(trans.length));
        meanISI_first3(transcand)=meanISI_frnt;

        %CS_ind=find(trans.spike_number>2 & trans.mean_ISI<15);
        CS_ind=find(trans.spike_number>2 & meanISI_first3<20);
        CS_trace=ismember(tr_trace,CS_ind);
        CS_spike=sp_soma.*bwlabel(CS_trace);
        [~, CS_spike_time]=unique(CS_spike);
    end

    sp_total=max([sp_soma; sp(2:end,:)],[],1);
    bAP_ind=zeros(1,size(tr,2));
    bAP_ind(unique(find(sp_soma)'+[0:3]))=1;

    SpikeClassMat=zeros(3,size(tr,2));
    SpikeClassMat(1,:)=sp_soma.*(1-CS_trace); %bAPs
    SpikeClassMat(2,CS_spike_time(2:end))=1; %Complex spikes
    SpikeClassMat(3,:)=sp_total.*(1-bAP_ind); %dSpikes

    Result.spike=[sp_soma; sp(2:end,:)];
    Result.SpClass=SpikeClassMat;
    Result.CStrace=CS_trace;
    show_traces_spikes(Result.normTraces,Result.spike,[Result.SpClass(1:2,:); Result.Blue;]);
    save([fpath{i} '/OP_Result.mat'],'Result','-v7.3')
end
%% Label dendrite

f_local=[94 96 97 98 99 102 104:111 123 124 144 145 151:158 164 165 172:175 185 186 190 191 192];

for f=f_local(31:end)
    f
    load([fpath{f} '/OP_Result.mat'])
    Result.roilabel = interactive_ROIlabel(Result.ftprntSplit, Result.ref_im,'region');
    showScaleImage(Result.ftprntSplit>0,Result.roilabel,colormap(turbo));
    save([fpath{f} '/OP_Result.mat'],'Result','-v7.3')
end



%%
load([fpath{89} '/Result.mat']); nROI=size(Result.normTraces,1);
coord_1d=dim_reduce(get_coord(Result.ftprnt));
[~, dist_order]=sort(coord_1d,'ascend');
nTau_bAP=[-20:20];
bAP_ref=find(Result.SpClass(1,:));
prc_normTr=Result.normTraces;
STA_SS=squeeze(mean(reshape(prc_normTr(:,bAP_ref'+nTau_bAP),nROI,[],length(nTau_bAP)),2));
STA_SS=STA_SS-prctile(STA_SS,25,2);
F_ref=mean(STA_SS(:,-nTau_bAP(1)+[10:14]),2);
prc_normTr=prc_normTr./F_ref;
show_traces_spikes(prc_normTr(dist_order,:),Result.spike(dist_order,:),[Result.SpClass; Result.Blue]);
%%
Blueon=find((Result.Blue(2:end)-Result.Blue(1:end-1))>0);
STA_Blue=reshape(prc_normTr(:,Blueon'+[-50:150]),nROI,[],length([-50:150]));
CS_on=find(Result.SpClass(2,:));
STA_CS=reshape(prc_normTr(:,CS_on'+[-100:150]),nROI,[],length([-100:150]));
figure(10); clf;
imagesc(flipud(squeeze(mean(STA_CS(dist_order,:,:,[0 5])
colormap(turbo)

%% STAs
nTau=[-10:20];

for i=[133]%1:length(OP_Result)
    Result.normTrace=Result.traces./get_threshold(Result.traces,1);
    Result.spike=find_spike_bh(Result.normTrace-movmedian(Result.normTrace,300,2),5,3);

    Blue=Result.Blue;
    blueOff = Blue == 0;
    blueOff2 = imerode(blueOff, [ones(1,20), zeros(1, 20)]);
    Blue_di=~blueOff2;
    bwBlue_di=bwlabel(Blue_di);

    ref_ROI=find(sum(Result.spike(1:5,:).*bwBlue_di,2)==max(sum(Result.spike(1:5,:).*bwBlue_di,2)),1);
    nROI=size(Result.normTrace,1);
    tr=Result.normTrace(ref_ROI,:); t=[1:length(tr)];
    spike=Result.spike(ref_ROI,:);


    sp_pulse=[];
    for b=1:max(bwBlue_di)
        t_tmp=find(bwBlue_di==b);
        if max(bwBlue_di)<500
            sp_pulse=[sp_pulse find(spike(t_tmp))+t_tmp(1)-1];
        else
            sp_pulse=[sp_pulse find(spike(t_tmp),1,'first')+t_tmp(1)-1];
        end
    end

    sTau=sp_pulse(1:end-1)'+nTau;
    spikeMat=reshape(Result.normTrace(:,sTau),nROI,length(sp_pulse)-1,[]);
    STAtrace=squeeze(mean(spikeMat,2));
    STAtrace=STAtrace-prctile(STAtrace, 5, 2);
    Result.spikeMat = spikeMat;
    Result.STAtrace = STAtrace;

    if max(bwBlue_di)>15

        load(fullfile(fpath{i},"output_data.mat"))
        load([fpath{i} '/mcTrace' num2str(1,'%02d') '.mat']);
        switch char(CamType(i))
            case 'flash'
                sz=double(Device_Data{1, 4}.ROI([2 4]));
            case 'fusion'
                sz=double(Device_Data{1, 3}.ROI([2 4]));
        end

        CamCounter=Device_Data{1, 2}.Counter_Inputs(1, 1).data;
        CamTrigger=find(CamCounter(2:end)-CamCounter(1:end-1));

        mov_mc=double(readBinMov_times([fpath{i} '/mc_ShutterReg' num2str(1,'%02d') '.bin'],sz(2),sz(1),[1:length(CamTrigger)]));

        mov_res= mov_mc-mean(mov_mc,3);
        bkg = zeros(1, size(mov_mc,3));

        bkg(1,:)=movmedian(get_blueoffTrace(squeeze(mean(mov_mc,[1 2])),Blue,30),3000,'omitnan');
        mov_res = SeeResiduals(mov_res,Result.mc);
        mov_res = SeeResiduals(mov_res,Result.mc.^2);
        mov_res = SeeResiduals(mov_res,Result.mc(:,1).*Result.mc(:,end));
        mov_res= SeeResiduals(mov_res,bkg,1);

        STA_tmp=reshape(mov_res(:,:,sTau),sz(2),sz(1),[],length(nTau));
        %STA_tmp=STA_tmp-mean(STA_tmp(:,:,:,[1:3]),4);
        Result.STAmovie=squeeze(mean(STA_tmp,3));
    end
end

%% Generate SNAPT movie
for i=[133]
    load([fpath{i} '/Result.mat'])
    mask=max(Result.Structure_bin,[],3)>0.01;
    maskSTA=max(-Result.STAmovie,[],3)./Result.ref_im>0.05;
    StrImg=max(Result.Structure,[],3);
    STAmovie=mat2gray(-Result.STAmovie);
    STAmovie=STAmovie-prctile(STAmovie,10,3);
    STAmovie=mat2gray(STAmovie(:,:,6:22));
    tformReg=Result.tform;
    [Result.SNAPT Result.dtimg]=generate_SNAPTmov(mat2gray(STAmovie),mask,StrImg,tformReg);
    [yR xR zR]=size(Result.Structure);
    bluePatt = bwboundaries(imwarp(Result.BlueDMDimg,Result.tform,'OutputView', imref2d([yR xR])));

    figure(20); clf;
    %v = VideoWriter([fpath{i} '/SNAPT_movie'],'MPEG-4');
    v = VideoWriter([fpath{i} '/SNAPT_movie'],'Uncompressed AVI');

    open(v);
    subframeT = 0.025; % ms
    initialT = -2; % ms
    finalT = 2; % ms
    times = initialT:subframeT:finalT;

    for j = 1:length(times)
        clf;
        %set(gca,'units','pixels','position',[200 0 1000 800])
        imshow(Result.SNAPT(:,:,:,j),[])
        pbaspect([size(double(Result.SNAPT(:,:,:,j)),2) size(double(Result.SNAPT(:,:,:,j)),1) 1]),colormap(gray)
        hold all
        plot(bluePatt{1}(:,2),bluePatt{1}(:,1),'color',[0 0.6 1],'linewidth',2)
        axis off
        text(2,20,[num2str(times(j)+0.9) 'ms'], 'FontSize', 20, 'color', [0.99 0.99 0.99])% the value 1. is to adjust timing by eyes
        pause(0.1)
        set(gcf,'color','w')    % Sets background to white
        frame = getframe(gcf);
        writeVideo(v,frame);
        pause(0.1);
    end;
    close(v);
end
%save([save_to 'OP_Result_20240212'],"OP_Result",'fpath','-v7.3')
%%

[a, unqInd] = unique([Mouse NeuronInd] ,'row');

for i=unqInd(20);
    cat_trace=[];
    cat_spike=[];

    SameCellInd=find(Mouse==Mouse(i) & NeuronInd==NeuronInd(i));
    for j=SameCellInd'
        load(fullfile(fpath{j},'OP_Result.mat'),'Result')
        cat_trace=[cat_trace Result.traces];
        cat_spike=[cat_spike Result.spike];
    end

    cat_trace=cat_trace-movmedian(cat_trace,300,2);
    [V D eigTrace]=get_eigvector(cat_trace',10);
    includePC=find(cumsum(D)/sum(D)>0.99,1);
    F0_PCA=sqrt(sum((V(:,[1:includePC]).*sqrt(D([1:includePC]))').^2,2));
    for j=SameCellInd'
        load(fullfile(fpath{j},'OP_Result.mat'),'Result')
        coord_1d=dim_reduce(get_coord(Result.ftprnt));
        [~, Result.dist_order]=sort(coord_1d,'ascend');
        Result.F0_PCA=F0_PCA;
        save([fpath{j} '/OP_Result.mat'],'Result','-v7.3')
    end
    figure(i); clf;
    nexttile([1 1])
    imagesc(cat_trace(Result.dist_order,:)./F0_PCA(Result.dist_order),[-0.02 0.03])
    show_footprnt(Result.ftprnt(:,:,Result.dist_order),Result.ref_im)
end

