clear
clc;
cd '/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/FromBackup/PP72_PlaceCellResults';
[~, ~, raw] = xlsread(['/Volumes/cohen_lab/Lab/Labmembers' ...
    '/Byung Hun Lee/Data/PrismPCdata_Arrangement.xlsx'], 'Sheet1', 'C5:P27');

% [~, ~, NeuronsToUse]=xlsread(['/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/' ...
%     'PlaceCellData_Arrangement.xlsx'], 'Sheet1', 'L8:M46');
%
% NeuronsToUse=cellfun(@(x) (str2num(num2str(x))),NeuronsToUse,'UniformOutput',false);
ref_ROI=cellfun(@(x) (str2num(num2str(x))),raw(:,10),'UniformOutput',false);
fpath=raw(:,1)';
StructureData=raw(:,10);
EndFrame=cell2mat(raw(:,13));
ifmotionReject=cell2mat(raw(:,14));
save_figto='/Volumes/BHL_WD18TB/PP72_PlaceCellResults';
place_bin=150; time_segment=15000; overlap=200;
alignedMovFN = {'STA_Mat_SS','STA_Mat_CS','STA_Mat_dSP'};

%% Motion correction

for f=[23]% 15 4 5 7]%length(fpath)
    f
    load(fullfile(fpath{f},"output_data.mat"))
    sz=double(Device_Data{1, 3}.ROI([2 4]));
    ref_time=[119000:120000];

    %frm_end=max(Device_Data{1, 2}.Counter_Inputs(1, 1).data);
    frm_end=EndFrame(f);
    f_seg=[[1:time_segment:frm_end] frm_end+1]; f_seg(2:end)=f_seg(2:end)-1;
    DAQ_rate=Device_Data{1, 2}.buffered_tasks(1, 1).rate;
    CamDAQ_rate=Device_Data{1, 2}.Counter_Inputs.rate;
    CamTrig=Device_Data{1, 2}.Counter_Inputs.data;
    CamTrig2=find(CamTrig(2:end)-CamTrig(1:end-1)>0);
    Frm_rate=(CamTrig2(2)-CamTrig2(1))/CamDAQ_rate;

    if length(ref_time)>2000
        mov_test=double(readBinMov_times([fpath{f} '/frames1.bin'],sz(2),sz(1),[ref_time(1):ref_time(1)+2000]));
    else
        mov_test=double(readBinMov_times([fpath{f} '/frames1.bin'],sz(2),sz(1),[ref_time(1):ref_time(end)]));
    end
    mov_test=rollingShutter_correction(mov_test,1/Frm_rate,'fusion');
    mov_test=mov_test(:,:,2:end);
    [mov_test_mc,xyField]=optical_flow_motion_correction_LBH(mov_test,mean(mov_test,3),'normcorre');
    mov_test=vm(mov_test);
    mov_test = single(mov_test)./single(max(mov_test.data(:)));
    mov_test = movmean(mov_test,10,3);
    mov_ref = squeeze(median(mov_test,3));

    for j=1:length(f_seg)-1
        try
            mov=double(readBinMov_times([fpath{f} '/frames1.bin'],sz(2),sz(1),[f_seg(j):f_seg(j+1)+overlap]));
        catch % when the image ends
            mov=double(readBinMov_times([fpath{f} '/frames1.bin'],sz(2),sz(1),[f_seg(j):f_seg(j+1)]));
        end

        mov=rollingShutter_correction(mov,1/Frm_rate,'fusion');
        mov=vm(mov(:,:,2:end));
        if j==1
            mov=mov(:,:,[1 1:size(mov,3)]);
        end

        [mov_mc,xyField]=optical_flow_motion_correction_LBH(mov,mov_ref,'normcorre');

        ave_im=mean(mov_mc,3);
        mov_mc=vm(mov_mc);
        mov_mc.transpose.savebin([fpath{f} '/mc_ShutterReg' num2str(j,'%02d') '.bin'])

        %        mcTrace = squeeze(mean(xyField,[1 2])); %optic flow
        mcTrace=xyField; % Normcorre
        save([fpath{f} '/mcTrace' num2str(j,'%02d') '.mat'],'mcTrace','ave_im')

        %  clear mov_mc mov
    end
end

%% Get footprint

for f=[21]%:length(fpath)
    disp(fpath{f}); Result=[];
    DAQ_rate=0.000005;
    load([fpath{f} '/output_data.mat'])
    sz=double(Device_Data{1, 3}.ROI([2 4]));
    ref_time=[2000:3000];
    mov_test=double(readBinMov_times([fpath{f} '/mc_ShutterReg' num2str(10,'%02d') '.bin'],sz(2),sz(1),ref_time));
    avgImg=mean(mov_test,3);
    figure(3); clf;
    imshow2(avgImg,[])
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

    [u,s,v] = svds(tovec(mov_test-mean(mov_test,3)),20);
reshape_u=reshape(u,sz(2),sz(1),[]);
bvMask=[];
[~, Result.bvMask]=get_ROI(max(abs(reshape_u),[],3),bvMask);

    load(fullfile(fpath{f},'mcTrace20.mat'));
    frm_end=max(Device_Data{1, 2}.Counter_Inputs(1, 1).data);
    Result.Blue=Device_Data{1, 2}.buffered_tasks(1, 1).channels(1, 2).data;
    Result.Reward=Device_Data{1, 2}.buffered_tasks(1, 3).channels.data;
    CamCounter=Device_Data{1, 2}.Counter_Inputs(1, 1).data;
    CamTrigger=find(CamCounter(2:end)-CamCounter(1:end-1));

    Result.Blue=Result.Blue(CamTrigger); Result.Reward=Result.Reward(CamTrigger);
    frm_rate=double((CamTrigger(2)-CamTrigger(1))*DAQ_rate);
    mov_mc=double(readBinMov([fpath{f} '/mc_ShutterReg' num2str(20,'%02d') '.bin'],sz(2),sz(1)));
    Result.ref_im=mean(mov_mc,3);
    [clickyROI, original_trace]=clicky(mov_mc);

    mov_res= mov_mc-mean(mov_mc,3);
    mcTrace.xymean=movmean(mcTrace.xymean,3,2);
    mov_res = SeeResiduals(mov_res,mcTrace.xymean);
    mov_res = SeeResiduals(mov_res,mcTrace.xymean.^2);
    mov_res = SeeResiduals(mov_res,mcTrace.xymean(:,1).*mcTrace.xymean(:,3));

    [SeeRes_trace]=apply_clicky(clickyROI,mov_res);
    figure(3); clf;
    plot(rescale(original_trace)); hold all;
    plot(rescale(SeeRes_trace)+1);

    n_comp=6;
    mov_filt=imgaussfilt3(mov_res.*double(max(Result.bvMask,[],3)==0),[2 2 0.1]);
    mov_filt=mov_filt(:,:,1:7000);
    movVec=tovec(mov_filt);
    Npoly=size(Result.ROIpoly,1);
    ftprnt = zeros(size(mov_filt,1)*size(mov_filt,2),Npoly);
    clear mask
    figure(4); 
    for p=1:Npoly %each ROIs
        clf; ax2=[];
    tiledlayout(n_comp/2+2,2)
        mask(:,:,p) = poly2mask(Result.ROIpoly{p}(:,1), Result.ROIpoly{p}(:,2), sz(2), sz(1));
        pixelList=find(tovec(squeeze(mask(:,:,p))));
        subMov = movVec(pixelList,:);
        covMat = subMov*subMov';  % PCA within each region
        [V, D] = eig(covMat);
        D = diag(D);
        D = D(end:-1:1);
        V = V(:,end:-1:1);
        vSign = sign(max(V) - max(-V));  % make the largest value always positive
        V = V.*vSign;
        eigTrace=subMov'*V;
        nexttile([2 2])
        plot(rescale2(eigTrace(:,1:n_comp),1)+[1:n_comp])

        %[icsTrace, ~, sepmat]=sorted_ica(eigTrace(:,1:n_comp),n_comp);
        %plot(rescale2(icsTrace,1)+[1:size(icsTrace,2)])
        %V_ics=V(:,1:n_comp)*sepmat';
        for n=1:n_comp
    eigImg=NaN(size(mov_filt,1)*size(mov_filt,2),1);
    ax2=[ax2 nexttile([1 1])];
    eigImg(pixelList,1)=V(:,n);
    eigImg=toimg(eigImg,size(mov_filt,1),size(mov_filt,2));
    imshow2(im_merge(cat(3,Result.ref_im,eigImg),[1 1 1;1 0 0]),[])
    title([num2str(n) ', Fraction: ' num2str(D(n)/sum(D),2)])
        end
linkaxes(ax2,'xy')
n_take = input('#components to take: ', 's');
n_take = str2num(n_take);
coeff=subMov*mean(eigTrace(:,n_take)*V(:,n_take)',2);
ftprnt(pixelList,p)=coeff;    
    end
close(figure(4));
    Result.ftprnt=toimg(ftprnt,sz(2),sz(1));
    figure; clf;
    show_footprnt(Result.ftprnt,Result.ref_im)
    save(fullfile(fpath{f},'PC_Result.mat'),'Result','fpath','-v7.3')
end
%% Load Virmen data
    
for f=[21]%:length(fpath)
    load([fpath{f} '/output_data.mat'])
    load(fullfile(fpath{f},'PC_Result.mat'),'Result')
    fileList = dir(fullfile(fpath{f}, '*.data'));
    if length(fileList)==1
        fid = fopen(fullfile(fpath{f},fileList.name));
    else
        error('Data file cannot be found');
    end
    VRdata = fread(fid,[12 inf],'double');

    WorldITrack=find(VRdata(2,:)==1); %world 1
    VRdata(5,WorldITrack)=(VRdata(5,WorldITrack)+6)*115/121;

    VRdata=VRdata(:,VRdata(10,:)>0);
    DAQ_rate=Device_Data{1, 2}.Counter_Inputs.rate;
    CamCounter=Device_Data{1, 2}.Counter_Inputs(1, 1).data;
    CamTrigger=find(CamCounter(2:end)-CamCounter(1:end-1));
    t_DAQ=CamTrigger/DAQ_rate;
    t_VR = datetime(datetime(VRdata(1,:),'ConvertFrom','datenum'), 'InputFormat', 'yyyy-MM-dd HH:mm:ss.SSS');
    t_VR= t_VR-t_VR(1);
    t_VR= milliseconds(t_VR)/1000;
    t_VR= t_VR*t_DAQ(end)/t_VR(end);
    VRdata(1,:)=t_VR;
    [Virmen_data_int vel_trace]=virmen_interpolate(VRdata,115,t_DAQ);
    Virmen_data_int(end+1,:)=vel_trace;

    Result.VR=Virmen_data_int;
    %Result.Blue=Result.Blue(:,1:EndFrame(f)).*Result.VR(12,1:EndFrame(f));
    Result.Blue=Result.Blue(:,1:EndFrame(f)).* ...
    ismember(Result.VR(8,1:EndFrame(f)),[13:23]); % From Lap 13 to 23
 
    Result.VR=Result.VR(:,1:EndFrame(f));

    save(fullfile(fpath{f},'PC_Result.mat'),'Result','fpath','-v7.3')
end

%% Signal extraction

for f=[21]%length(fpath)
    load(fullfile(fpath{f},'PC_Result.mat'),'Result');
    load([fpath{f} '/output_data.mat'])
    sz=double(Device_Data{1, 3}.ROI([2 4])); blueDMDcontour=[];
    CamCounter=Device_Data{1, 2}.Counter_Inputs(1, 1).data;
    CamTrigger=find(CamCounter(2:end)-CamCounter(1:end-1));
    Rfixed = imref2d(repmat(Device_Data{1, 3}.virtualSensorSize,1,2));
    inverseTform = invert(Device_Data{1, 6}.tform);
    try
        revertedImage = imwarp(double(Device_Data{6}.pattern_stack(:,:,find(sum(Device_Data{6}.pattern_stack,[1 2])>0))), inverseTform,'OutputView',Rfixed);
    catch
        revertedImage = imwarp(double(Device_Data{6}.Target), inverseTform,'OutputView',Rfixed);
    end
    [blueDMDimg]=imcrop_3d(revertedImage,double(Device_Data{1, 3}.ROI([1 3 2 4]))+[0 0 -1 -1]);
    for d=1:size(blueDMDimg,3)
        blueDMDcontour{d}=bwboundaries(blueDMDimg(:,:,d));
    end
 
    oddmask=checkerboard(1, size(Result.ftprnt,1)/2, size(Result.ftprnt,2)/2) > 0.5;
    oddftprnt=Result.ftprnt.*oddmask;
    evenftprnt=Result.ftprnt.*(1-oddmask);

    Result.blueDMDimg=blueDMDimg;
    Result.blueDMDcontour=blueDMDcontour;
    Result.traces=[]; Result.traces_bvMask=[];
    Result.traces_checker{1}=[]; Result.traces_checker{2}=[];
    Result.tracesPCA=[];
    Result.mcTrace=[];
    Result.im_corr=[];
    bound=5;
    ref_im_vec=tovec(Result.ref_im(bound:end-bound,bound:end-bound));
    ref_im_vec=tovec(Result.ref_im(bound:end-bound,bound:end-bound))/std(ref_im_vec,0,1);    
    %frm_end=max(Device_Data{1, 2}.Counter_Inputs(1, 1).data);
    frm_end=EndFrame(f);
    f_seg=[[1:time_segment:frm_end] frm_end+1]; f_seg(2:end)=f_seg(2:end)-1;
    take_window=repmat([1 time_segment],length(f_seg)-1,1);
    take_window(2:end,1)=take_window(2:end,1)+overlap; take_window(1:end-1,2)=take_window(1:end-1,2)+overlap;
    take_window(end)=mod(f_seg(end),time_segment);
    take_window(take_window==0)=time_segment;
    Blue_on_Seg=unique(ceil(find(Result.Blue)/time_segment));

    for j=1:length(f_seg)-1
        j
        mov_mc=double(readBinMov([fpath{f} '/mc_ShutterReg' num2str(j,'%02d') '.bin'],sz(2),sz(1)));
        load([fpath{f} '/mcTrace' num2str(j,'%02d') '.mat']);

        mov_mc=mov_mc(:,:,[take_window(j,1):take_window(j,2)]);
        %mc= movmean(mcTrace.xymean-movmedian(mcTrace.xymean,500,1),3,1);
        mc=mcTrace.xymean([take_window(j,1):take_window(j,2)],:);
        
        mov_mc_filt=imgaussfilt3(mov_mc,[2 2 0.1]);
        mov_mc_vec=tovec(mov_mc_filt(bound:end-bound,bound:end-bound,:));

        if ismember(j,Blue_on_Seg)
           isfirst_seg=double(j==1);
        [~, blueomitTr]=get_blueoffTrace(squeeze(mean(mov_mc,[1 2])), ...
            Result.Blue(f_seg(j)+take_window(j,1)-isfirst_seg:f_seg(j)+take_window(j,2)-isfirst_seg),30);
        bkg = zeros(1, size(mov_mc,3));    
        bkg(1,:)=movmedian(blueomitTr,4000,'omitnan');
        else
        bkg = zeros(2, size(mov_mc,3));    
        bkg(1,:) = linspace(-1, 1, size(mov_mc,3));  % linear term
        bkg(2,:) = linspace(-1, 1, size(mov_mc,3)).^2;  % quadratic term
        end
        
        mov_res=mov_mc;
        % mov_res= mov_mc-median(mov_mc,3);
        % mov_res = SeeResiduals(mov_res,mc);
        % mov_res = SeeResiduals(mov_res,mc.^2);
        % mov_res = SeeResiduals(mov_res,mc(:,1).*mc(:,end));
        % mov_res= SeeResiduals(mov_res,bkg,1);

        Result.traces=[Result.traces -(tovec(mov_res)'*tovec(Result.ftprnt))'];
        Result.traces_bvMask=[Result.traces_bvMask -(tovec(mov_res.*double(max(Result.bvMask,[],3)==0))'*tovec(Result.ftprnt))'];
        Result.traces_checker{1}=[Result.traces_checker{1} -(tovec(mov_res.*double(max(Result.bvMask,[],3)==0))'*tovec(oddftprnt))'];
        Result.traces_checker{2}=[Result.traces_checker{2} -(tovec(mov_res.*double(max(Result.bvMask,[],3)==0))'*tovec(evenftprnt))'];
%       Result.tracesPCA=[Result.tracesPCA -(tovec(mov_res)'*tovec(Result.pcaMask))'];
        Result.mcTrace=[Result.mcTrace; mc];
        Result.im_corr=[Result.im_corr corr(rescale2(mov_mc_vec,1),ref_im_vec,'type','Spearman')'];  %image correlation

    end

    %Result.traces=Result.traces(:,1:length(CamTrigger)-1);
    %Result.mcTrace=Result.mcTrace(1:length(CamTrigger)-1,:);
    Result.traces=Result.traces(:,1:frm_end);
    Result.mcTrace=Result.mcTrace(1:frm_end,:);
    Result.Blue=Result.Blue(:,1:frm_end);
    Result.Reward=Result.Reward(:,1:frm_end);

    save(fullfile(fpath{f},'PC_Result.mat'),'Result','fpath','-v7.3')
end

%% Clean up and norm
exclude_frq=[241.7 242]; %monitor
%exclude_frq2=[483.5 484]; %monitor
exclude_frq2=[25 65.7]; %motion
time_bin=15000; Fs=1000; ref_trace=[2 4 4 2 3 2 1 1 1 3 1 1 1 1 1 1 1 1 1 3 21 3 23]; %2nd trunk is the reliable trace

for f=[21]%:length(fpath)
    load(fullfile(fpath{f},'PC_Result.mat'),'Result')
    nTime=size(Result.traces,2);
    nROI=size(Result.traces,1);
tr=Result.traces-movmedian(Result.traces,50,2);
tr=tr./get_threshold(tr,1);
sp=find_spike_bh(tr,5,3);

SilentPeriod=ones(1,size(Result.traces,2));
sp_time=find(max(sp,[],1))';
sp_na=sum((find(max(sp,[],1))'+[-10:150])<0 | (find(max(sp,[],1))'+[-10:150])>size(Result.traces,2),2)==0;
SilentPeriod(sp_time(sp_na)+[-10:150])=NaN;

t_fit=find(~isnan(SilentPeriod));
[y_fit2 t_consts coeffY]  = expfitDM_2(t_fit',-Result.traces(1,t_fit)',[1:nTime]',[10000 100000]);
tr_res=squeeze(SeeResiduals(permute(Result.traces_bvMask,[1 3 2]),y_fit2));

tr_res_checker{1}=squeeze(SeeResiduals(permute(Result.traces_checker{1},[1 3 2]),y_fit2));
tr_res_checker{2}=squeeze(SeeResiduals(permute(Result.traces_checker{2},[1 3 2]),y_fit2));

    freq_lowhigh=exclude_frq/(Fs/2);
    [b, a] = butter(4, freq_lowhigh, 'stop');

    freq_lowhigh2=exclude_frq2/(Fs/2);
    [b2, a2] = butter(4, freq_lowhigh2, 'stop');

    sub_pass_frq=[2];
    freq_lowhigh3=sub_pass_frq/(Fs/2);
    [b3, a3] = butter(4, freq_lowhigh3, 'low');

    theta_pass_frq=[5 11];
    freq_lowhigh4=theta_pass_frq/(Fs/2);
    [b4, a4] = butter(4, freq_lowhigh4, 'bandpass');

    % figure(2); clf;
    % [p f]=fft_simple(PC_Result{f}.traces(2,:),1000);
    % ax1=nexttile([1 1]);
    % plot(f(1,:),p')
    % set(gca,'YScale','log')
    %
    % [p f]=fft_simple(PC_Result{f}.mcTrace,1000);
    % ax2=nexttile([1 1]);
    % plot(f(1,:),p')
    % set(gca,'YScale','log')
    % linkaxes([ax1 ax2],'x')

    clear traces_res_filtered noise noise_intp norm_trace sp_height SpHeight_intp sp_time tr_mc_imcorr tr_mc tr norm_trace_check traces_res_filtered_ch
    tN=[1:time_bin:nTime]; tN=[tN nTime];
    sp_time=zeros(nROI,nTime);
    sp_height=zeros(nROI,nTime);
    mcTrace=squeeze(Result.mcTrace)';

    for n=1:size(Result.traces,1)
        tr=tr_res(n,1:nTime);
        tr_check{1}=tr_res_checker{1}(n,1:nTime); tr_check{2}=tr_res_checker{2}(n,1:nTime);

        % regress out motion frequency
        for t=1:length(tN)-1
            tr(tN(t):tN(t+1))=squeeze(SeeResiduals(reshape(tr(tN(t):tN(t+1)),1,1,[]),mcTrace(:,(tN(t):tN(t+1)))));
            tr(tN(t):tN(t+1))=squeeze(SeeResiduals(reshape(tr(tN(t):tN(t+1)),1,1,[]),mcTrace(:,(tN(t):tN(t+1))).^2));
            tr(tN(t):tN(t+1))=squeeze(SeeResiduals(reshape(tr(tN(t):tN(t+1)),1,1,[]),mcTrace(1,(tN(t):tN(t+1))).*mcTrace(end,(tN(t):tN(t+1)))));
            for ch=1:2 % Checkerboard pattern traces
                tr_check{ch}(tN(t):tN(t+1))=squeeze(SeeResiduals(reshape(tr_check{ch}(tN(t):tN(t+1)),1,1,[]),mcTrace(:,(tN(t):tN(t+1)))));
                tr_check{ch}(tN(t):tN(t+1))=squeeze(SeeResiduals(reshape(tr_check{ch}(tN(t):tN(t+1)),1,1,[]),mcTrace(:,(tN(t):tN(t+1))).^2));
                tr_check{ch}(tN(t):tN(t+1))=squeeze(SeeResiduals(reshape(tr_check{ch}(tN(t):tN(t+1)),1,1,[]),mcTrace(1,(tN(t):tN(t+1))).*mcTrace(end,(tN(t):tN(t+1)))));
            end
            tr_mc(tN(t):tN(t+1))=tr(tN(t):tN(t+1));
        end

        % regress out monitor frequency
        traces_res_filtered(n,:) = filtfilt(b, a, tr);
        traces_res_filtered_ch{1}(n,:) = filtfilt(b, a, tr_check{1});
        traces_res_filtered_ch{2}(n,:) = filtfilt(b, a, tr_check{2});
        
        tr_mc_bpMonitorMotion=traces_res_filtered(n,:);

        norm_trace(n,:)=traces_res_filtered(n,:);%-movmedian(traces_res_filtered(n,:),500,2);
        norm_trace_check{1}(n,:)=traces_res_filtered_ch{1}(n,:);%-movmedian(traces_res_filtered(n,:),500,2);
        norm_trace_check{2}(n,:)=traces_res_filtered_ch{2}(n,:);%-movmedian(traces_res_filtered(n,:),500,2);

        if n==ref_trace(f) %compensate spike hight due to bleaching
            for t=1:length(tN)-1
                tr_tmp=norm_trace(n,tN(t):tN(t+1));
                tr_tmp=tr_tmp-movmedian(tr_tmp,300);
                noise(t,n)=get_threshold(tr_tmp,1);
                tr_tmp_norm=tr_tmp./noise(t,n);
                [sp_temp, pks, prom]=find_spike_bh(tr_tmp_norm,12,4);
                sp_time(n,tN(t):tN(t+1))=sp_temp;
            end

            t_fit=find(sp_time(n,:));
            sp_height(n,t_fit)=norm_trace(n,t_fit);
            tx=tN(1:end-1)+time_bin/2;
            %noise_intp(n,:)=movmean(interp1(tx,noise(:,n),[1:size(PC_Result{f}{i}.traces,2)],'linear','extrap'),10000);

            [Ny_fit t_consts coeffY]  = expfitDM_2(tx(~isnan(noise(1:end-1,n)))',noise(~isnan(noise(1:end-1,n)),n),[1:nTime]',10^7);
            noise_intp=Ny_fit;
            %[y_fit t_consts coeffY]  = expfitDM_2(tx(~isnan(sp_height(1:end-1,n)))',sp_height(~isnan(sp_height(1:end-1,n)),n),[1:size(Result{i}.traces,2)]',10^7);
            [y_fit t_consts coeffY]  = expfitDM_2(t_fit',sp_height(n,t_fit)',[1:nTime]',10^7);
            SpHeight_intp=y_fit;


            figure; clf;
            ax1=nexttile([1 1]);
            plot(rescale2([Result.traces(n,1:nTime);tr_mc;tr_mc_bpMonitorMotion;mcTrace(1,1:nTime)],2)')%+[1:4])
            ax2=nexttile([1 1]);
            plot([1:nTime],norm_trace(n,:))
            hold all
            plot(find(sp_time(n,:)),norm_trace(n,find(sp_time(n,:))),'r.')
            plot([1:nTime],y_fit([1:nTime]),'k')
            plot([1:nTime],Ny_fit([1:nTime]),'g')

        end
    end
    Result.SpikeHeight_fit=SpHeight_intp';
    norm_trace=norm_trace./SpHeight_intp';
    norm_trace_check{1}=norm_trace_check{1}./SpHeight_intp';
    norm_trace_check{2}=norm_trace_check{2}./SpHeight_intp';

    %norm_trace=norm_trace;%./(SpHeight_intp./SpHeight_intp(:,1));
    Result.normTraces=norm_trace./get_threshold(norm_trace,1);
    Result.norm_trace_check{1}=norm_trace_check{1}./get_threshold(norm_trace_check{1},1);
    Result.norm_trace_check{2}=norm_trace_check{2}./get_threshold(norm_trace_check{2},1);
    %Result.spike=find_spike_bh(PC_Result{f}.normTraces-movmedian(PC_Result{f}.normTraces,300,2),5,3);

    % for n=1:size(Result.traces,1)
    % Result.subThreshold(n,:) = filtfilt(b3, a3, Result.normTraces(n,:));
    % Result.theta(n,:) = filtfilt(b4, a4, Result.normTraces(n,:));
    % end
    save(fullfile(fpath{f},'PC_Result.mat'),'Result','fpath','-v7.3')
end


%% Classify the spikes
motion_frq=[30 120];
ifmotion_reject=ifmotionReject;
nSpikeThres=1;
for f=[21]
    load(fullfile(fpath{f},'PC_Result.mat'))
    %tr=PC_Result{i}.normTraces(ref_ROI{i},:); %somatic spike
    %sp_ref=max(find_spike_bh(tr-movmedian(tr,100,2),5,3),[],1);
    tr=Result.normTraces;
    tr_ref=Result.normTraces(ref_ROI{f},:);

    nROI=size(tr,1);
    DOI{f}=setdiff([2:nROI],ref_ROI{f});
    sp=find_spike_bh(tr-movmedian(tr,50,2),6,4);
    sp_ref=find_spike_bh(tr_ref-movprc(tr_ref,200,30,2),4,3);
    [wvletTr wvletF] = cwt(Result.im_corr,1000); % Compute the CWT
    motionArtTrace=sum(abs(wvletTr(find(wvletF>(motion_frq(1)) & wvletF<motion_frq(2)),:)),1);
    motionArtTrace=(motionArtTrace-prctile(motionArtTrace,20))./std(motionArtTrace);
    motionReject=zeros(1,size(tr,2));
    if ifmotionReject(f)
        motionReject= motionArtTrace>4;
        motionReject = imdilate(motionReject, strel('square', 200));
        sp(:,motionReject)=0; sp_ref(:,motionReject)=0;
    end
    Result.motionReject=motionReject;
    sp_ref=sum(sp_ref,1);

    [~, shift]=max(reshape(tr(ref_ROI{f}(1),find(sp_ref>nSpikeThres)+[-1:0]'),2,[]),[],1);
    shift=shift-2;
    sp_time_Soma = find(sp_ref>nSpikeThres)+shift;
    sp_soma=zeros(1,size(tr,2));
    sp_soma(sp_time_Soma)=1;
    sp_soma=[0 (sp_soma(2:end)-sp_soma(1:end-1))==1]; %remove consecutive spikes

    tr_sub=mean(tr_ref,1)-movprc(mean(tr_ref,1),200,20,2);
    tr_sub=get_subthreshold(tr_sub,sp_soma,5,10);

    [trans tr_trace]=detect_transient2(tr_sub,[4 1.5],sp_soma,15);
    transcand=cell2mat(cellfun(@(x) length(x)>2,trans.ISI,'UniformOutput',false));
    meanISI_frnt=cellfun(@(x) mean(x(1:2)),trans.ISI(transcand));
    meanISI_first3=NaN(1,length(trans.length));
    meanISI_first3(transcand)=meanISI_frnt;

    %CS_ind=find(trans.spike_number>2 & trans.mean_ISI<15);
    CS_ind=find(trans.spike_number>2 & meanISI_first3<18);
    CS_trace=ismember(tr_trace,CS_ind);
    CS_spike=sp_soma.*bwlabel(CS_trace);
    [~, CS_spike_time]=unique(CS_spike);

    sp_total=max([sp_soma; sp(DOI{f},:)],[],1);
    bAP_ind=zeros(1,size(tr,2));
    bAP_ind(unique(find(sp_soma)'+[-1:3]))=1;

    SpikeClassMat=zeros(3,size(tr,2));
    SpikeClassMat(1,:)=sp_soma.*(1-CS_trace); %bAPs
    SpikeClassMat(2,CS_spike_time(2:end))=1; %Complex spikes
    SpikeClassMat(3,:)=sp_total.*(1-bAP_ind).*(1-CS_trace); %dSpikes
    SpikeClassMat(3,:)=[0 (SpikeClassMat(3,2:end)-SpikeClassMat(3,1:end-1))==1]; %remove consecutive spikes

    Result.spike=[sp_soma; sp(2:end,:)];
    Result.SpClass=SpikeClassMat;
    Result.CStrace=CS_trace;

nROI=size(Result.ftprnt,3);
SkelDend = Skeletonize_dendrite(Result.ref_im,8,0.005,25);
interDendDist=[];
for i=1%:nROI
    for j=1:nROI
        [interDendDist(i,j), ~]=geodesic_distance(SkelDend,get_coord(Result.ftprnt(:,:,i)),get_coord(Result.ftprnt(:,:,j)));
    end
end

coord_1d=dim_reduce(get_coord(Result.ftprnt));
[~, Result.dist_order]=sort(coord_1d,'descend');
Result.geodist=interDendDist(1,:)'.*sign(coord_1d-coord_1d(1));
som_roi=find(Result.dist_order==1);
show_traces_spikes(Result.normTraces(Result.dist_order,:),Result.spike(Result.dist_order,:),[Result.SpClass; double(Result.motionReject)]);

    save(fullfile(fpath{f},'PC_Result.mat'),'Result','fpath','-v7.3')
end

%% Save the STA movies
nTau={[-30:20],[-50:100],[-30:20]}; %SS, CS, dSP
bound=6;
%f_tmp='/Volumes/BHL18TB_D1/20240218/134705BHLm117_FOV2_VR2';
for f=[21]
    load(fullfile(fpath{f},'PC_Result.mat'),'Result')
    StackedSpike=[];
    alignmovlist=dir(fullfile(fpath{f},'STA_Mat*.tiff'));
    alignmov_ind=ones(1,3);
    for l=1:length(alignmovlist)
        delete(fullfile(fpath{f},alignmovlist(l).name));
    end

    for c=1:3 % Get frames of spikes of each classes
        fileInd{c}=ceil(find(Result.SpClass(c,:))/time_segment);
        frameInd{c}=mod(find(Result.SpClass(c,:)),time_segment);
    end
    %line up movie segments at simple spike

    load([fpath{f} '/output_data.mat'])
    sz=double(Device_Data{1, 3}.ROI([2 4]));
    frm_end=max(Device_Data{1, 2}.Counter_Inputs(1, 1).data);
    f_seg=[[1:time_segment:frm_end] frm_end+1]; f_seg(2:end)=f_seg(2:end)-1;
    g=zeros(1,3);
    ss=[0 0 0];
    for j=unique(cell2mat(fileInd))
        j

        mov_mc=double(readBinMov([fpath{f} '/mc_ShutterReg' num2str(j,'%02d') '.bin'],sz(2),sz(1)));
        load([fpath{f} '/mcTrace' num2str(j,'%02d') '.mat']);

        if j==length(f_seg)-1
            mc=mcTrace.xymean;
        else
            mov_mc=mov_mc(:,:,1:time_segment);
            mc=mcTrace.xymean(1:time_segment,:);
        end

        mov_res= mov_mc-mean(mov_mc,3);
        bkg = zeros(2, size(mov_mc,3));
        bkg(1,:) = linspace(-1, 1, size(mov_mc,3));  % linear term
        bkg(2,:) = linspace(-1, 1, size(mov_mc,3)).^2;  % quadratic term
        mov_res = SeeResiduals(mov_res,mc);
        mov_res = SeeResiduals(mov_res,mc.^2);
        mov_res = SeeResiduals(mov_res,mc(:,1).*mc(:,end));
        mov_res= SeeResiduals(mov_res,bkg,1);

        for c=1:3
            STA_movie_align{c}=[];
            for k=find(fileInd{c}==j)
                if frameInd{c}(k)+nTau{c}(1)>1 && frameInd{c}(k)+nTau{c}(end)<size(mov_res,3)
                    mov_seg=mov_res(bound:end-bound,bound:end-bound,frameInd{c}(k)+nTau{c});
                    %mov_seg=mat2gray(mov_seg);
                    STA_movie_align{c} = cat(3,STA_movie_align{c}, mov_seg);
                    g(c)=g(c)+1;
                    StackedSpike{c}(1,g(c))=1;
                    StackedSpike{c}(2,g(c))=time_segment*(j-1)+frameInd{c}(k);
                else
                    STA_movie_align{c} = cat(3,STA_movie_align{c}, zeros(sz(2)-bound*2+1,sz(1)-bound*2+1,length(nTau{c})));
                    g(c)=g(c)+1;
                    StackedSpike{c}(1,g(c))=0;
                    StackedSpike{c}(2,g(c))=NaN;
                end
            end
            if ~isempty(STA_movie_align{c})

                % write_tif_stack(STA_movie_align{c},fullfile(f_tmp,[alignedMovFN{c} '_' num2str(alignmov_ind(c)) '.tiff']))
                % alignmovlist=dir(fullfile(f_tmp,[alignedMovFN{c} '*.tiff']));  
                write_tif_stack(STA_movie_align{c},fullfile(fpath{f},[alignedMovFN{c} '_' num2str(alignmov_ind(c)) '.tiff']))
                alignmovlist=dir(fullfile(fpath{f},[alignedMovFN{c} '*.tiff']));
                if alignmovlist(alignmov_ind(c)).bytes > 2.0*10^9
                    disp('Move on to 2nd tiff')
                    alignmov_ind(c) = alignmov_ind(c) + 1;
                end
            end
        end
        g
    end

    Result.StackedSpike=StackedSpike;
    save(fullfile(fpath{f},'PC_Result.mat'),'Result','fpath','-v7.3')
end

%% Masking blood vessel cross section
f=11;
load(fullfile(fpath{f},'PC_Result.mat')); j=5;
load([fpath{f} '/output_data.mat'])
sz=double(Device_Data{1, 3}.ROI([2 4]));

mov_mc=double(readBinMov([fpath{f} '/mc_ShutterReg' num2str(j,'%02d') '.bin'],sz(2),sz(1)));
load([fpath{f} '/mcTrace' num2str(j,'%02d') '.mat']);

mov_res= mov_mc-mean(mov_mc,3);
bkg = zeros(2, size(mov_mc,3));
bkg(1,:) = linspace(-1, 1, size(mov_mc,3));  % linear term
bkg(2,:) = linspace(-1, 1, size(mov_mc,3)).^2;  % quadratic term
mov_res = SeeResiduals(mov_res,mcTrace.xymean);
mov_res = SeeResiduals(mov_res,mcTrace.xymean.^2);
mov_res = SeeResiduals(mov_res,mcTrace.xymean(:,1).*mcTrace.xymean(:,end));
mov_res= SeeResiduals(mov_res,bkg,1);

[u,s,v] = svds(tovec(mov_res),20);
reshape_u=reshape(u,sz(2),sz(1),[]);
Result.bvMask=[];
[~, Result.bvMask]=get_ROI(max(abs(reshape_u),[],3),Result.bvMask);
%save(fullfile(fpath{f},'PC_Result.mat'),'Result','fpath','-v7.3')

%% Structure segmentation
Struct_valid=find(1-cell2mat(cellfun(@(x) sum(isnan(x)), StructureData, 'UniformOutput', false)));

for i=Struct_valid(5)'
    load(fullfile(fpath{i},'PC_Result.mat'))
    StructureStack=mat2gray(double(tiffreadVolume(StructureData{i})));
    StructureStack(StructureStack==0)=median(StructureStack(:));
    StructureStack=StructureStack(:,:,35:100);
    %StructureStack_med=medfilt2_mov(StructureStack,[15 15]);
    StructureStack_Gauss=imgaussfilt3(StructureStack,[7 7 0.1]);
    %StructureStack_med(StructureStack_med==0)=median(StructureStack_med(:));
    %StructureStack=(StructureStack-StructureStack_med)./StructureStack_med;
    StructureStack_filt=(StructureStack-StructureStack_Gauss);
    StructureStack_filt=mat2gray(StructureStack_filt);
    StructureStack_bin=[]; level=[];
    level = graythresh(StructureStack_filt);
    StructureStack_bin=StructureStack_filt>level*0.97;
    moviefixsc(StructureStack_bin)

    se = strel('sphere', 1);
    StructureStack_bin = imdilate(StructureStack_bin, se);
    bwSeg=bwlabeln(StructureStack_bin);
    segments = regionprops3(bwSeg,'Volume','EquivDiameter');
    segments = table2array(segments);
    %bwlist=find(arrayfun(@(x) x.Volume>4000, segments) & arrayfun(@(x) x.EquivDiameter>10, segments));
    bwlist = find(segments(:,1)>7000 & segments(:,2)>10);

    se = strel('sphere',1);
    dendrite_bin=double(ismember(bwSeg,bwlist));
    dendrite_bin= imdilate(dendrite_bin,se);
    dendrite_bin= imgaussfilt3(dendrite_bin,2);

    figure(3); clf;
    imshow2(max(dendrite_bin,[],3),[])
    g=1; ROIrmv=[];
    while g
        h = drawpolygon('Color','r');
        if size(h.Position,1)==1 %no more ROI
            g=0;
        else
            ROIrmv=[ROIrmv; {h.Position}];
            hold all
            plot(h.Position(:,1),h.Position(:,2))
        end
    end
    ROIrmvmask=roi2mask(ROIrmv,size(dendrite_bin,1),size(dendrite_bin,2));
    close(figure(3));
    dendrite_bin(repmat(ROIrmvmask,1,1,size(dendrite_bin,3)))=0;
    figure(3); clf;
    imshow2(max(dendrite_bin,[],3),[])

    StructureStack_final = double(StructureStack).* dendrite_bin;
    figure(2); clf;
    imshow2(imfuse(mat2gray(max(StructureStack_final,[],3)),mat2gray(max(StructureStack,[],3))),[])


    rot_ang=0;
    Structure_ref=(imrotate(StructureStack_final,rot_ang));
    ref_img=Result.ref_im; ref_img(ref_img<prctile(ref_img(:),20))=median(ref_img(:)); ref_img=ref_img-prctile(ref_img(:),20);
    [RegImg,tformReg]=imReg(ref_img,max(Structure_ref,[],3));
    saveastiff(uint16(mat2gray(Structure_ref)*255), [fpath{i} 'Structure.tiff']);

    Result.Structure=max(Structure_ref,[],3);
    Result.Structure_bin=max(imrotate(dendrite_bin,rot_ang),[],3);
    Result.tform=tformReg;
    save(fullfile(fpath{i},'PC_Result.mat'),'Result','fpath','-v7.3')
end

%% Place field
f=14;
load(fullfile(fpath{f},'PC_Result.mat'))
place_bin=150;
velocity_threshold=-0.002;
cmap=[1 1 1; 1 0 0; 0 0.5 1];
for c=1:3
    PF=PlaceTrigger_average(Result.SpClass(c,:)>0,place_bin,Result.VR,velocity_threshold,115);
    PF=movmean(repmat(PF,1,3),3,2);
    PF=PF(:,place_bin+1:2*place_bin);
    Result.PF_class(:,:,c)=PF;
end

dSP_list=find(Result.SpClass(3,:));
for doi=2:size(Result.normTraces,1)
    dSP_trace=zeros(1,size(Result.normTraces,2));
    dSP_list_sub=ismember(dSP_list,find(Result.spike(doi,:)));
    dSP_trace(1,dSP_list(find(dSP_list_sub)))=1;
    PFdSP=PlaceTrigger_average(dSP_trace,place_bin,Result.VR,velocity_threshold,115);
    PFdSP=movmean(repmat(PFdSP,1,3),3,2);
    PFdSP=PFdSP(:,place_bin+1:2*place_bin);
    Result.Lap_FR(:,:,doi)=PFdSP;
end

figure; clf;
nexttile([[1 1]])
imagesc(im_merge(Result.PF_class,cmap))
nexttile([[1 1]])
imagesc(im_merge(Result.Lap_FR,jet(size(Result.normTraces,1))))

%% pca analysis of spike triggered movie
f=4; c=1; bound_more=6; bound=6;
toi=[16:24];

load(fullfile(fpath{f},'PC_Result.mat'))
alignmovlist=dir(fullfile(fpath{f},[alignedMovFN{c} '*.tiff']));
AlignMov=[];
for l=1:length(alignmovlist)
    AlignMov=cat(3,AlignMov,readtiff(fullfile(fpath{f},alignmovlist(l).name)));
end

sz_align=size(AlignMov);
AlignMov=double(reshape(AlignMov,sz_align(1),sz_align(2),length(nTau{c}),[]));
AlignMov=AlignMov-mean(AlignMov(:,:,1:5,:),3);
% gF = fspecial('gaussian', 7, 2);
% AlignMov = imfilter(AlignMov, gF, 'same', 'replicate');

bvMask=Result.bvMask(bound:end-bound,bound:end-bound,:);
bvMask=bvMask(bound_more:end-bound_more,bound_more:end-bound_more,:);
% cmask=cellmask(bound:end-bound,bound:end-bound,:);
% cmask=cmask(bound_more:end-bound_more,bound_more:end-bound_more,:);

%isolated spike
s_list=[];
for s=1:size(Result.StackedSpike{c},2)
    isnearby=zeros(1,3);
    for cc=1:3
        isnearby(cc)=sum(ismember(Result.StackedSpike{c}(2,s)'+nTau{c},Result.StackedSpike{cc}(2,:)));
    end
    if sum(isnearby)==1
        s_list=[s_list s];
    end
end

%Triggering somatic spike

doi=[21];
dSP_s=find(Result.SpClass(3,:)>0 & sum(Result.spike(doi,:),1)>0);
dSP_s_IndSoma=[];
som_spike=find(Result.spike(1,:)>0);
for s=dSP_s
    isnearby=sum(ismember(s+[0 1 2],som_spike))>0;
    if isnearby
        dSP_s_IndSoma=[dSP_s_IndSoma s];
    end
end
s_list=ismember(Result.StackedSpike{c}(2,:),dSP_s_IndSoma);

%AlignMov_sub=AlignMov(bound_more:end-bound_more,bound_more:end-bound_more,toi,s_list);
AlignMov_sub=mean(AlignMov(bound_more:end-bound_more,bound_more:end-bound_more,toi,s_list),3);
sz_sub=size(AlignMov_sub);
AlignMov_sub=reshape(AlignMov_sub,sz_sub(1)*sz_sub(2)*sz_sub(3),[]);
AlignMov_sub(repmat(tovec(max(bvMask,[],3))>0,sz_sub(3),sz_sub(4)))=0;
% AlignMov_sub(repmat(tovec(cmask==0),sz_sub(3),sz_sub(4)))=0;
% AlignMov_sub=reshape(AlignMov_sub,sz_sub(1),sz_sub(2),sz_sub(3),[]);
%     gF = fspecial('gaussian', 7, 2);
%     AlignMov_sub = imfilter(AlignMov_sub, gF, 'same', 'replicate');
% AlignMov_sub=reshape(AlignMov_sub,sz_sub(1)*sz_sub(2)*sz_sub(3),[]);
[u,s,v] = svds(AlignMov_sub,20);
reshape_u2=reshape(u,sz_sub(1),sz_sub(2),sz_sub(3),[]);

u_movie=[];
for j=1:5
    for i=1:4
        u_movie(sz_sub(1)*(i-1)+1:sz_sub(1)*i,sz_sub(2)*(j-1)+1:sz_sub(2)*j,:)=reshape_u2(:,:,:,j+(i-1)*5);
    end
end
moviefixsc([u_movie; imgaussfilt3(u_movie,[1.5 1.5 0.1])])

[ics, mixmat, sepmat] = sorted_ica(u,size(u,2));
reshape_ics=reshape(ics,sz_sub(1),sz_sub(2),sz_sub(3),[]);
ics_movie=[];
for j=1:5
    for i=1:4
        ics_movie(sz_sub(1)*(i-1)+1:sz_sub(1)*i,sz_sub(2)*(j-1)+1:sz_sub(2)*j,:)=reshape_ics(:,:,:,j+(i-1)*5);
    end
end
moviefixsc([ics_movie; imgaussfilt3(ics_movie,[1.5 1.5 0.1])])

[~, bv_tmp]=get_ROI(max(mean(reshape_u2(:,:,:,[2]),3),[],4));
Result.bvMask(bound+bound_more-2:bound+bound_more-3+sz_sub(1),bound+bound_more-2:bound+bound_more-3+sz_sub(2),end+1:end+size(bv_tmp,3))=bv_tmp;

pca_trace=AlignMov_sub'*u;
ics_trace=AlignMov_sub'*ics; %event
%ics_trace=pca_trace*sepmat';
figure; clf;
tiledlayout(1,2)
nexttile([1 1]);
plot(rescale2(pca_trace,1)+[1:size(pca_trace,2)]-0.5)
axis tight
nexttile([1 1]);
plot(rescale2(ics_trace,1)+[1:size(ics_trace,2)]-0.5)
axis tight

interest_inx=[1 2 3 5 10 11 15]; sign_convert=[1 1 1 1 -1 -1 1];
s_ics=zscore(ics_trace(:,interest_inx).*sign_convert,1)<-2.5;
figure;
show_mov=[];

for inx=interest_inx
    AlignMovVec=reshape(AlignMov(:,:,:,s_list(s_ics(:,interest_inx==inx))),sz_align(1)*sz_align(2)*length(nTau{c}),[]);
    AlignMovVec(repmat(tovec(max(Result.bvMask(bound:end-bound,bound:end-bound,:),[],3))>0,length(nTau{c}),sum(s_ics(:,interest_inx==inx))))=0;
    show_mov=[show_mov; imgaussfilt3(mean(reshape(AlignMovVec,sz_align(1),sz_align(2),length(nTau{c}),[]),4),[2 2 1])];
end
moviefixsc(show_mov,[-1500 1000])

pca_trace_z=zscore(pca_trace,0,1);
ica_trace_z=zscore(ics_trace,0,1);
vr_list=Result.VR(5,Result.StackedSpike{c}(2,s_list));
figure;
for ind=1:6
    s_ind=find(abs(pca_trace_z(:,ind))>1.5);
    plot(repmat(ind,length(s_ind),1),vr_list(s_ind),'.','markersize',12); hold all
end
xlim([0.5 6.5]); %ylim([75 98]);
xlabel('Component #'); ylabel('VR position')

% Correlation with ica image

load([fpath{f} '/output_data.mat'])
sz=double(Device_Data{1, 3}.ROI([2 4])); blueDMDcontour=[];
CamCounter=Device_Data{1, 2}.Counter_Inputs(1, 1).data;
CamTrigger=find(CamCounter(2:end)-CamCounter(1:end-1));

frm_end=EndFrame(f);
f_seg=[[1:time_segment:frm_end] frm_end+1]; f_seg(2:end)=f_seg(2:end)-1;
take_window=repmat([1 time_segment],length(f_seg)-1,1);
take_window(2:end,1)=take_window(2:end,1)+overlap; take_window(1:end-1,2)=take_window(1:end-1,2)+overlap;
take_window(end)=mod(f_seg(end),time_segment);

ics_im_vec=tovec(squeeze(reshape_ics));
ics_im_vec=ics_im_vec./std(ics_im_vec,0,1);
ics_corr=[];

pca_im_vec=tovec(squeeze(reshape_u2));
pca_im_vec=pca_im_vec./std(pca_im_vec,0,1);
pca_corr=[];

for j=1:length(f_seg)-1
    j
    mov_mc=double(readBinMov([fpath{f} '/mc_ShutterReg' num2str(j,'%02d') '.bin'],sz(2),sz(1)));
    load([fpath{f} '/mcTrace' num2str(j,'%02d') '.mat']);

    mov_mc=mov_mc(:,:,[take_window(j,1):take_window(j,2)]);
    mc=mcTrace.xymean([take_window(j,1):take_window(j,2)],:);

    mov_res= mov_mc-mean(mov_mc,3);
    bkg = zeros(2, size(mov_mc,3));
    bkg(1,:) = linspace(-1, 1, size(mov_mc,3));  % linear term
    bkg(2,:) = linspace(-1, 1, size(mov_mc,3)).^2;  % quadratic term
    mov_res = SeeResiduals(mov_res,mc);
    mov_res = SeeResiduals(mov_res,mc.^2);
    mov_res = SeeResiduals(mov_res,mc(:,1).*mc(:,end));
    mov_res= SeeResiduals(mov_res,bkg,1);

    mov_res=mov_res(bound:end-bound,bound:end-bound,:);
    mov_res=mov_res(bound_more:end-bound_more,bound_more:end-bound_more,:);

    mov_mc_vec=tovec(mov_res);
    mov_mc_vec=(mov_mc_vec-mean(mov_mc_vec,1))./std(mov_mc_vec,0,1);
    ics_corr=[ics_corr; mov_mc_vec'*ics_im_vec/(size(mov_mc_vec,1)-1)];
    pca_corr=[pca_corr; mov_mc_vec'*pca_im_vec/(size(mov_mc_vec,1)-1)];

end
Result.ics_corr=ics_corr;
Result.pca_corr=pca_corr;

%% dSpike vs bAP amplitude versus distance
f=11;
geodist=[]; geopath=[];
load(fullfile(fpath{f},'PC_Result.mat'))
nROI=size(Result.ftprnt,3);
SkelDend = Skeletonize_dendrite(Result.ref_im,4,0.02,25);
for roi=1:nROI
    [geodist(roi), geopath(:,:,roi)]=geodesic_distance(SkelDend,get_coord(Result.ftprnt(:,:,1)),get_coord(Result.ftprnt(:,:,roi)));
end
show_footprnt_contour(geopath,Result.ref_im)
Result.geoDist=geodist;
Result.geoPath=geopath;

som_spike=find(Result.SpClass(1,:)>0);
bAP_s=[];
for s=som_spike
    isnearby=sum(ismember(s+nTau{1},som_spike))>1;
    if ~isnearby
        bAP_s=[bAP_s s];
    end
end
%bAP_s=som_spike;

STA_SSmat=reshape(Result.normTraces(:,bAP_s'+nTau{1}),nROI,[],length(nTau{1}));
STA_SSmat=STA_SSmat-mean(mink(STA_SSmat,5,3),3);
STA_SS=squeeze(mean(STA_SSmat,2));

ds=find(Result.SpClass(3,:)>0);
dSP_s=[];
for s=ds
    isnearby=sum(ismember(s+nTau{1},som_spike))>1;
    if ~isnearby
        dSP_s=[dSP_s s];
    end
end

STA_dSPmat=reshape(Result.normTraces(:,dSP_s'+nTau{3}),nROI,[],length(nTau{3}));
STA_dSPmat=STA_dSPmat-mean(STA_dSPmat(:,:,[1:5]),3);
STA_dSP=squeeze(mean(STA_dSPmat,2));

F_ref=mean(STA_SS(:,21+[10:15]),2);
Result.F_ref=F_ref;

Amp_SSmat=max(STA_SSmat(:,:,[20:23]),[],3)./F_ref;
%Amp_SSmat=Amp_SSmat./Amp_SSmat(1,:);

Amp_dSPmat=max(STA_dSPmat(:,:,[20:23]),[],3)./F_ref;
%Amp_dSPmat=Amp_dSPmat./Amp_dSPmat(1,:);

attenuation_SS=[]; attenuation_dSP=[];
ft = fittype('poly1'); % 'poly1' specifies a first-degree polynomial (linear fit)
for s=1:size(Amp_SSmat,2)
    [fitresult, gof] = fit(Result.geoDist', Amp_SSmat(:,s), ft);
    attenuation_SS(s)=fitresult.p1;
end

for s=1:size(Amp_dSPmat,2)
    [fitresult, gof] = fit(Result.geoDist', Amp_dSPmat(:,s), ft);
    attenuation_dSP(s)=fitresult.p1;
end

[gD, Result.distOrder]=sort(Result.geoDist,'ascend');
figure(3); clf;
plot(gD,Amp_SSmat(Result.distOrder,:)','.','color',[0.8 0.5 0.5]); hold all
plot(gD,Amp_dSPmat(Result.distOrder,:)','.','color',[0.3 0.5 0.8]);

errorbar(gD,mean(Amp_SSmat(Result.distOrder,:),2),std(Amp_SSmat(Result.distOrder,:),0,2),'color','r','LineStyle','-','marker','+','MarkerSize',15); hold all
errorbar(gD,mean(Amp_dSPmat(Result.distOrder,:),2),std(Amp_dSPmat(Result.distOrder,:),0,2),'color',[0 0.5 1],'LineStyle','-','marker','+','MarkerSize',15);
axis tight
xlabel('Distance (\mum)');
ylabel('Spike amplitude')

figure(4); clf;
plot(Result.geoDist,max(STA_SS(:,[20:23])./F_ref,[],2),'r')

figure;
imshow2(max(reshape(max(STA_SS(:,[20:23])./F_ref,[],2),1,1,20).*double(Result.ftprnt>0),[],3),[])
colormap('turbo')


%% Generate SNAPT movie

f=14;
load(fullfile(fpath{f},'PC_Result.mat'))
mask=max(Result.Structure_bin,[],3)>0.01;
StrImg=max(Result.Structure,[],3);
STAmovie=zeros(size(Result.ref_im,1),size(Result.ref_im,2),size(AlignMov,3));
STAmovie(bound:end-bound,bound:end-bound,:)=mat2gray(-mean(AlignMov(:,:,:,s_list),4));
STAmovie=STAmovie-prctile(STAmovie,10,3);
STAmovie=mat2gray(STAmovie(:,:,16:25));
tformReg=Result.tform;
[Result.SNAPT Result.dtimg]=generate_SNAPTmov(mat2gray(STAmovie),mask,StrImg,tformReg);

dtimg_Reg=imwarp(Result.dtimg, tformReg, 'OutputView', imref2d(size(StrImg)));
dtimg_Reg=dtimg_Reg.*mask; dtimg_Reg(dtimg_Reg==0)=NaN;
imshow2(dtimg_Reg-prctile(dtimg_Reg(:),5), [0 3]); title('Timing')
hold all
Maskboundary = cell2mat(bwboundaries(mask));
plot(Maskboundary(:,2),Maskboundary(:,1),'r.','markersize',2);
colormap('turbo')


[yR xR zR]=size(Result.Structure);
figure(20); clf;
v = VideoWriter([fpath{f} '/dSP_indSS_SNAPT_movie'],'MPEG-4');

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
    axis off
    text(2,20,[num2str(times(j)+0.9) 'ms'], 'FontSize', 20, 'color', [0.99 0.99 0.99])% the value 1. is to adjust timing by eyes
    pause(0.1)
    set(gcf,'color','w')    % Sets background to white
    frame = getframe(gcf);
    writeVideo(v,frame);
    pause(0.1);
end;
close(v);


%% pca analysis of spike triggered movie
f=14; c=3; bound_more=6; bound=6; suffix=['dSpike_indSomSpike'];
toi=[16:25];

 load(fullfile(fpath{f},'PC_Result.mat'))
alignmovlist=dir(fullfile(fpath{f},[alignedMovFN{c} '*.tiff']));
AlignMov=[];
for l=1:length(alignmovlist)
    l
    AlignMov=cat(3,AlignMov,readtiff(fullfile(fpath{f},alignmovlist(l).name)));
end

sz_align=size(AlignMov);
AlignMov=double(reshape(AlignMov,sz_align(1),sz_align(2),length(nTau{c}),[]));
%AlignMov=AlignMov-median(AlignMov,3);
AlignMov=AlignMov-mean(maxk(AlignMov,7,3),3);
%AlignMov=AlignMov-mean(AlignMov(:,:,1:5,:),3);
% % gF = fspecial('gaussian', 7, 2);
% % AlignMov=imfilter(AlignMov_sub, gF, 'same', 'replicate');

%s_list=[1:sz_align(3)/length(nTau{c})];
s_list=[];
for s=1:size(Result.StackedSpike{c},2)
    isnearby=zeros(1,3);
    for cc=1:3
        isnearby(cc)=sum(ismember(Result.StackedSpike{c}(2,s)'+nTau{c},Result.StackedSpike{cc}(2,:)));
    end
    if sum(isnearby)==1
        s_list=[s_list s];
    end
end

%Triggering somatic spike
% 
doi=[2:size(Result.normTraces,1)];
dSP_s=find(Result.SpClass(3,:)>0 & sum(Result.spike(doi,:),1)>0);
dSP_s_IndSoma=[];
som_spike=find(Result.spike(1,:)>0);
for s=dSP_s
    isnearby=sum(ismember(s+[0 1 2],som_spike))>0;
if isnearby
dSP_s_IndSoma=[dSP_s_IndSoma s];
end
end
s_list=ismember(Result.StackedSpike{c}(2,:),dSP_s_IndSoma);

AlignMov_Vec=reshape(AlignMov(:,:,toi,s_list),sz_align(1)*sz_align(2),[]);
bvMask=tovec(max(Result.bvMask(bound:end-bound,bound:end-bound,:),[],3));
ROImask=tovec(max(Result.ftprnt(bound:end-bound,bound:end-bound,:)>0,[],3));
AlignMov_Vec2=AlignMov_Vec;
%AlignMov_Vec2((bvMask==1 | ROImask==0),:)=0;
AlignMov_Vec2((bvMask==1),:)=0;
AlignMov_Vec((bvMask==1),:)=0;
AlignMov
%AlignMov_Vec=AlignMov_Vec(~bvMask & ROImask,:);

covMat = AlignMov_Vec*AlignMov_Vec';

[V, D] = eig(covMat);
D = diag(D);
D = D(end:-1:1);
V = V(:,end:-1:1);
vSign = sign(max(V) - max(-V));  % make the largest value always positive
V = V.*vSign;
eigTraces = (V'*AlignMov_Vec);

figure(8); clf
nexttile([1 1])
nKeep = 30;
for i=1:nKeep
    plot(rescale(eigTraces(i,:)-movmedian(eigTraces(i,:),30))+i-0.5)
    hold all
end
grid_x=tovec(repmat([length(toi):length(toi):size(eigTraces,2)],3,1));
grid_y=tovec(repmat([0 nKeep+1 NaN]',length(find(s_list>0)),1));
plot(grid_x,grid_y,'linestyle',':','color',[0.7 0.7 0.7])
set(gca, 'YDir','reverse')
nexttile([1 1])
plot(cumsum(D)./sum(D))

saveas(gca,fullfile(fpath{f} ,['EigTrace_' num2str(c) '_' suffix '_' num2str(toi([1 end])) '.fig']))

figure(9); clf;
eigImgs=toimg(AlignMov_Vec2*eigTraces(1:nKeep,:)',sz_align(1),sz_align(2));
ax1=[];
tiledlayout(5,8);
for j=1:10
    ax1=[ax1 nexttile([1 3])];
    imshow2(imgaussfilt(eigImgs(:,:,j),1),[])
    title(num2str(j))
nexttile([1 1])
shadedErrorBar(toi-21,mean(reshape(eigTraces(j,:),length(toi),[]),2),std(reshape(eigTraces(j,:),length(toi),[]),0,2))
end
linkaxes(ax1,'xy');
saveas(gca,fullfile(fpath{f} ,['EigImg_' num2str(c) '_' suffix '_' num2str(toi([1 end])) '.fig']))

% [~, bv_tmp]=get_ROI(max(eigImgs(:,:,[4 5 6 7 8 9 10]),[],3));
% Result.bvMask(bound:end-bound,bound:end-bound,end+1:end+size(bv_tmp,3))=bv_tmp;
% 
keep_ind=[1:10];
eigImgs_vec=tovec(eigImgs);
eigImgs_vec=eigImgs_vec(~bvMask & ROImask,:);

[ics, mixmat, sepmat] = sorted_ica(rescale2(eigImgs_vec(:,keep_ind),2),length(keep_ind));
%[ics, mixmat, sepmat] = sorted_ica(V(:,keep_ind),length(keep_ind));

%icsImgs = tovec(eigImgs(:,:,keep_ind))*sepmat';
icsImgs=zeros(sz_align(1)*sz_align(2),size(ics,2));
icsImgs(~bvMask & ROImask,:)=ics;    
icsImgs = toimg(icsImgs,sz_align(1),sz_align(2));
icsTrace = (AlignMov_Vec'*ics)';

figure(10); clf;
ax2=[];
tiledlayout(ceil(size(ics,2)/2),8)
for j=1:size(ics,2)
    ax2=[ax2 nexttile([1 3])];
    imshow2(imgaussfilt(icsImgs(:,:,j),1),[])
    title(num2str(j))
    nexttile([1 1])
    shadedErrorBar(toi-21,mean(reshape(icsTrace(j,:),length(toi),[]),2),std(reshape(icsTrace(j,:),length(toi),[]),0,2))
end
linkaxes(ax2,'xy');
saveas(gca,fullfile(fpath{f} ,['ICAImg_' num2str(c) '_' suffix '_' num2str(toi([1 end])) '.fig']))

%%
f=14; loadload(fullfile(fpath{f},'PC_Result.mat'))
ROItrace=Result.normTraces(:,~Result.motionReject);

Result.tracesPCA=Result.tracesPCA-median(Result.tracesPCA,2);
Result.tracesPCA=Result.tracesPCA./Result.SpikeHeight_fit;
bAP_frame=zeros(1,size(Result.traces,2));
bAP_frame(tovec(find(Result.spike(1,:))'+[0 1 2]))=1;
Result.tracesPCA(:,bAP_frame==1)=0;
Result.tracesPCA(:,Result.motionReject==1)=0;
Lap_PCA=[];
for n=1:size(Result.tracesPCA,1)
Lap_PCA(:,:,n)=PlaceTrigger_average(Result.tracesPCA(n,:),150,Result.VR,-0.02,115);
end

%% Branch-Branch correlation
f=14; load(fullfile(fpath{f},'PC_Result.mat'))
nROI=size(Result.normTraces,1);
nTau={[-20:20],[-30:300],[-20:20]}; %SS, CS, dSP
noi=setdiff([1:nROI],[20 8 3 6]);

%load aligned movie somatic spike
alignmovlist=dir(fullfile(fpath{f},[alignedMovFN{1} '*.tiff']));
AlignMov=[];
for l=1:length(alignmovlist)
    l
    AlignMov=cat(3,AlignMov,readtiff(fullfile(fpath{f},alignmovlist(l).name)));
end
sz_align=size(AlignMov);
AlignMov=double(reshape(AlignMov,sz_align(1),sz_align(2),length(nTau{1}),[]));
AlignMov=AlignMov-median(AlignMov(:,:,1:20,:),3);
%AlignMov=AlignMov-mean(AlignMov(:,:,1:5,:),3);
%AlignMov=AlignMov-mean(maxk(AlignMov,7,3),3);
%AlignMov=AlignMov-prctile(AlignMov,80,3);

SkelDend = Skeletonize_dendrite(Result.ref_im,4,0.02,25);
interDendDist=[];
for i=1:size(Result.normTraces,1)
    for j=1:size(Result.normTraces,1)
        [interDendDist(i,j), ~]=geodesic_distance(SkelDend,get_coord(Result.ftprnt(:,:,i)),get_coord(Result.ftprnt(:,:,j)));
    end
end
[~, dist_order_noi]=sort(interDendDist(1,noi),'ascend');
[~, dist_order]=sort(interDendDist(1,:),'ascend');
dist_order_inv_noi=mod(find([1:length(noi)]==dist_order_noi'),length(noi));
dist_order_inv_noi(dist_order_inv_noi==0)=length(noi);
dist_order_inv=mod(find([1:nROI]==dist_order'),nROI);
dist_order_inv(dist_order_inv==0)=nROI;

noi=setdiff([1:nROI],[20 8 3 6]);
rois={[6],[21]}; toi_heatmap=[18:19]; toi_pca=[15:19]; toi_movie=[10:25]; 

% Isolated Somatic spike
som_spike=Result.StackedSpike{1}(2,:); bAP_s=[]; 
for s=som_spike
    isnearby=sum(ismember(s+nTau{1},som_spike))>1;
    if ~isnearby & ~isnan(s)
        bAP_s=[bAP_s s];
    end
end
som_list=ismember(Result.StackedSpike{1}(2,:),bAP_s);

% % dSpike inducing Soma spike
% dSP_s=find(Result.SpClass(3,:)>0 & sum(Result.spike(2:end,:),1)>0);
% dSP_s_IndSoma=[];
% som_spike=find(Result.spike(1,:)>0);
% for s=dSP_s
%     isnearby=sum(ismember(s+[0 1 2],som_spike))>0;
%     if isnearby
%         [~, shift]=max(Result.normTraces(1,s+[0 1 2]));
%         dSP_s_IndSoma=[dSP_s_IndSoma s+shift-1];
%     end
% end

%PCA on Soma spike trace
nKeep = 50;
STA_SSmat=reshape(Result.normTraces(:,bAP_s'+nTau{1}),nROI,[],length(nTau{1}));
STA_SSmat=STA_SSmat-mean(mink(STA_SSmat,10,3),3);
%STA_SSmat=STA_SSmat-median(STA_SSmat,3);
STA_SSmat=STA_SSmat./max(tovec(permute(STA_SSmat,[1 3 2])),[],1);
STA_SS=squeeze(mean(STA_SSmat,2));
F_ref=mean(STA_SS(:,21+[10:14]),2);
STA_SSmat=STA_SSmat./F_ref;
STA_SSmatVec=tovec(permute(STA_SSmat(noi,:,toi_pca),[1 3 2]));
covMat=STA_SSmatVec*STA_SSmatVec';
[V, D] = eig(covMat);
D = diag(D);
D = D(end:-1:1);
V = V(:,end:-1:1);
vSign = sign(max(V) - max(-V));  % make the largest value always positive
V = V.*vSign;
% STA_SSeigImg = (V'*STA_SSmatVec);
% STA_SSeigImg=toimg(STA_SSeigImg',length(noi),length(toi_pca));
% STA_SSeigTr=tovec(STA_SSeigImg)'*STA_SSmatVec';
% 
STA_SSeigTr = (V'*STA_SSmatVec);
STA_SSeigImg=STA_SSeigTr(1:nKeep,:)*STA_SSmatVec';
[STA_SScoV, STA_SScoV_P]=corr(STA_SSeigImg',STA_SSmatVec);
STA_SSeigImg=toimg(STA_SSeigImg',length(noi),length(toi_pca));

figure(13); clf;
pos_bin=150; pos_range=[115 130]; 
for i=2:9
PCs=[1 i];
nexttile([1 1])
cmap_pos=turbo(pos_range(2)-pos_range(1)+1);
cmap_pos=[repmat(cmap_pos(1,:),pos_range(1)-1,1); cmap_pos; repmat(cmap_pos(end,:),pos_bin-pos_range(2),1)];
Sp_PosBin=ceil(Result.VR(5,bAP_s)/115*pos_bin);
sp_PosRange = Sp_PosBin>pos_range(1) & Sp_PosBin<pos_range(2);
scatter(STA_SSeigTr(PCs(1),sp_PosRange),STA_SSeigTr(PCs(2),sp_PosRange),[],'r','filled'); hold all
scatter(STA_SSeigTr(PCs(1),~sp_PosRange),STA_SSeigTr(PCs(2),~sp_PosRange),[],'k','filled');
xlabel(['PC #' num2str(PCs(1))])
ylabel(['PC #' num2str(PCs(2))])
end





% imagesc(abs(STA_SScoV.*(STA_SScoV_P<0.01)));
% colormap('turbo')

STA_SScoVSig=STA_SScoV.*(STA_SScoV_P<0.01);
CorrEvents=abs(STA_SScoVSig)>0.5;

Corr_movie=[]; Eigmov_corr=[];
AlignMov_sub=reshape(AlignMov(:,:,toi_movie,som_list),sz_align(1)*sz_align(2)*length(toi_movie),[]);
AlignMov_sub=AlignMov_sub./min(AlignMov_sub,[],1);
%AlignMov_sub=AlignMov_sub./prctile(AlignMov_sub,0.01,1);
for i=1:size(CorrEvents,1)
Eigmov_corr(:,:,:,i)=reshape(mean(AlignMov_sub(:,find(CorrEvents(i,:))),2),sz_align(1),sz_align(2),[]);
end
mean_STA=reshape(mean(AlignMov_sub,2),sz_align(1),sz_align(2),[]);
for j=1:5
    for i=1:2   
        Corr_movie(sz_align(1)*(i-1)+1:sz_align(1)*i,sz_align(2)*(j-1)+1:sz_align(2)*j,:)=Eigmov_corr(:,:,:,j+(i-1)*5);%-mean_STA;
    end
end
figure(7); clf;
nexttile([1 1])
%moviefixsc([Corr_movie; imgaussfilt3(Corr_movie,[2 2 0.1])],[-0.01 0.02])
moviefixsc([Corr_movie; imgaussfilt3(Corr_movie,[2 2 0.1])],[-0.01 0.05])

% for i=1:nKeep
%     %plot(rescale(STA_SSeigTr(i,:)-movmedian(STA_SSeigTr(i,:),30))+i-0.5)
%     hold all
% end
figure(8); clf;
nexttile([1 1])
plot(cumsum(D./sum(D))); hold all
bar([1:length(D)],D./sum(D)); xlim([0.5 nKeep+0.5]); ylim([0 1])

figure(9); clf;
for i=1:12
    nexttile([1 1])
    %imagesc(STA_SSeigImg(dist_order,:,i))
    STA_eigCorrImg(:,:,i)=squeeze(mean(STA_SSmat(:,find(CorrEvents(i,:)),:),2));
    imagesc(STA_eigCorrImg(dist_order,[7:25],i),[0.5 1.5]);
    title(['Pattern#' num2str(i) ', N=' num2str(sum(CorrEvents(i,:)))])
    %STA_eigImgSub=mean(STA_SSeigImg(:,toi_heatmap-toi_pca(1)+1,i),2,'omitnan');
    ftmap=double(Result.ftprnt(:,:,dist_order)>0.07).*reshape(mean(STA_eigCorrImg(dist_order,toi_heatmap,i),2),1,1,[]);
    ftprnt_XY=get_coord(Result.ftprnt(:,:,dist_order));
    ftmap(ftmap==0)=NaN;
    STA_SSeigImg_ftprnt(:,:,i)=max(ftmap,[],3);
    nexttile([1 1])
    imagesc(STA_SSeigImg_ftprnt(:,:,i),[0.5 1.5])
    text(ftprnt_XY(:,1),ftprnt_XY(:,2),num2str([1:nROI]'),'color','w')
    colormap('turbo')
end

pcoi={[1],[2],[3]};
pcoi_vennmat=[];
for p=1:length(pcoi)
    pcoi_vennmat(:,p)=max(double(CorrEvents(pcoi{p},:)),[],1)+1;
end
V=venn_data(pcoi_vennmat);
figure(10); clf;
[~, S]=venn(V,'FaceColor',{'r','y','b'},'FaceAlpha',{0.5,0.6,0.7},'EdgeColor','black');
text(S.ZoneCentroid(:,1),S.ZoneCentroid(:,2),num2str(V'),'HorizontalAlignment','center');
text(S.Position(:,1),S.Position(:,2)+S.Radius'/2,cellfun(@num2str,pcoi,'UniformOutput',false),'color','r','HorizontalAlignment','center')
axis equal off

sp_corr_trace=[];
for i=1:10
sp_corr_trace(i,:)=zeros(1,size(Result.normTraces,2));
sp_corr_trace(i,bAP_s(find(CorrEvents(i,:))))=1;
Lap_sp_Corr(:,:,i)=PlaceTrigger_average(sp_corr_trace(i,:),150,Result.VR,0,115);  
end

figure(11); clf;
for i=1:8
    nexttile([1 1])
imagesc(ringmovMean(Lap_sp_Corr(:,:,i),5))
colormap('turbo')
colorbar
title(num2str(i))
end


figure(12); clf;
rois={[13 8 9],[16 19 17]}; g=1;
cmap=distinguishable_colors(4);
cmap_light=distinguishable_colors(4)+0.3; cmap_light(cmap_light>1)=1;
for i=[1 2 3 7]
CorrContour=permute(STA_SSmat(dist_order,find(CorrEvents(i,:)),:),[1 3 2]);
MContour=mean(CorrContour,3);
CorrContour=reshape(CorrContour,nROI,[]);
%plot(mean(CorrContour(rois{1},:),1),mean(CorrContour(rois{2},:),1),'color',cmap_light(g,:)); hold all
plot(mean(MContour(rois{1},:),1),mean(MContour(rois{2},:),1),'marker','.','markersize',8,'color',cmap(g,:)); hold all
xlabel(['ROI ' num2str(rois{1})])
ylabel(['ROI ' num2str(rois{2})])
g=g+1;
end

%%
figure; toi_corr=[15:25];
imshow2(min(mean(AlignMov,4),[],3),[])
[x y]=ginput;
pointftprnt=zeros(size(AlignMov,1),size(AlignMov,2),length(x));
for i=1:length(x)
pointftprnt(round(y(i)),round(x(i)),i)=1;
SE = strel("disk",5);
pointftprnt(:,:,i)=imdilate(pointftprnt(:,:,i),SE);
end
pointTr=reshape(AlignMov(:,:,toi_corr,:),sz_align(1)*sz_align(2),[])'*tovec(pointftprnt);
AlignMov_time=reshape(AlignMov(:,:,toi_corr,:),sz_align(1)*sz_align(2),[]);
pointCorr=corr(AlignMov_time',pointTr);
pointCorr=reshape(pointCorr,sz_align(1),sz_align(2),[]);

figure(12); clf;
for i=1:size(pointCorr,3)
    nexttile([1 1])
    imshow2(pointCorr(:,:,i),[])
    [x y]=ginput(2);
    lineProfile = improfile(pointCorr(:,:,i), x, y);
    plot(lineProfile);
    nexttile([1 1])
    imshow2(im_merge(cat(3,pointCorr(:,:,i),pointftprnt(:,:,i)),[1 0.7 0.7; 0 0.5 1]),[])
    hold all
    plot(x,y,'r')
end
%%

subthAmp=[];
subthAmp(1,:)=mean(STA_SSmat(rois{1},:,toi),3);
subthAmp(2,:)=mean(mean(STA_SSmat(rois{2},:,toi),1),3);
figure;
cmap=jet(20);
subthAmp2=mean(STA_SSmat(:,:,toi),3);
figure(4); clf; 
Bias=[];
for i=1:size(Result.normTraces,1)
    %doi=setdiff([1:size(Result.normTraces,1)],i);
    for j=1:size(Result.normTraces,1)
%[interDendDist(j), ~]=geodesic_distance(SkelDend,get_coord(Result.ftprnt(:,:,i)),get_coord(Result.ftprnt(:,:,doi(j))));
[fitresult, gof] = fit(subthAmp2(i,:)', subthAmp2(j,:)', ft);
%Bias(i,doi(j))=fitresult.p1;
Bias(i,j)=corr(subthAmp2(i,:)', subthAmp2(j,:)');
    end
nexttile([1 1])
imagesc(max((Result.ftprnt>0.07).*reshape([Bias(i,:)],1,1,[]),[],3))
colormap(turbo)
axis equal tight off
hold all
c_ref=get_coord(Result.ftprnt(:,:,i));
text(c_ref(1),c_ref(2),num2str(i),'color','w')
end

STA_dSPmat=reshape(Result.normTraces(:,dSP_s_IndSoma'+nTau{1}),nROI,[],length(nTau{1}));
STA_dSPmat=STA_dSPmat-mean(STA_dSPmat(:,:,1:5),3);
STA_dSPmat=STA_dSPmat./Result.F_ref;
dSPAmp=[];
dSPAmp(1,:)=mean(STA_dSPmat(rois{1},:,toi),3);
dSPAmp(2,:)=mean(mean(STA_dSPmat(rois{2},:,toi),1),3);
dSPAmp2=mean(STA_dSPmat(:,:,toi),3);

figure(5); clf; 
Bias_dSP=[];
for i=1:size(Result.normTraces,1)
    %doi=setdiff([1:size(Result.normTraces,1)],i);
    for j=1:size(Result.normTraces,1)
[interDendDist(i,j), ~]=geodesic_distance(SkelDend,get_coord(Result.ftprnt(:,:,i)),get_coord(Result.ftprnt(:,:,j)));
[fitresult, gof] = fit(dSPAmp2(i,:)', dSPAmp2(j,:)', ft);
%Bias(i,doi(j))=fitresult.p1;
Bias_dSP(i,j)=corr(dSPAmp2(i,:)', dSPAmp2(j,:)');
    end
nexttile([1 1])
imagesc(max((Result.ftprnt>0.07).*reshape([Bias_dSP(i,:)],1,1,[]),[],3))
colormap(turbo)
axis equal tight off
hold all
c_ref=get_coord(Result.ftprnt(:,:,i));
text(c_ref(1),c_ref(2),num2str(i),'color','w')
end


figure;
nexttile([1 1]);
plot(subthAmp(1,:),subthAmp(2,:),'.'); hold all
plot(dSPAmp(1,:),dSPAmp(2,:),'o');
axis equal tight 
xlabel('Branch 1 amplitude')
ylabel('Branch 2 amplitude');
xlim([0 6])
ylim([0 6])
x_edges = 0:0.25:10;
y_edges = 0:0.25:10;
nexttile([1 1]);
h = histogram2(subthAmp(1,:),subthAmp(2,:),x_edges,y_edges,'DisplayStyle','tile');
colormap(jet)
cb = colorbar();
cb.Label.String = 'Bin Count';
xlabel('Branch 1 amplitude')
ylabel('Branch 2 amplitude');
axis equal
xlim([0 9])
ylim([0 9])
% 
% nexttile([1 1]);
% x_edges = 0:0.5:10;
% y_edges = 0:0.5:10;
% h = histogram2(dSPAmp(1,:),dSPAmp(2,:),x_edges,y_edges,'DisplayStyle','tile');
% colormap(jet)
% cb = colorbar();
% cb.Label.String = 'Bin Count';
% xlabel('Branch 1 amplitude')
% ylabel('Branch 2 amplitude');
% axis equal
% xlim([0 9])
% ylim([0 9])
%%
f=14; loadload(fullfile(fpath{f},'PC_Result.mat'))
dSP_s=find(Result.SpClass(3,:)>0 & sum(Result.spike(2:end,:),1)>0);
dSP_s_IndSoma=[];
som_spike=find(Result.spike(1,:)>0);
for s=dSP_s
    isnearby=sum(ismember(s+[0 1 2],som_spike))>0;
    if isnearby
        [~, shift]=max(Result.normTraces(1,s+[0 1 2]));
        dSP_s_IndSoma=[dSP_s_IndSoma s+shift-1];
    end
end

Result.normTraces
