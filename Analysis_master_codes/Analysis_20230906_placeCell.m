clear
clc;
sourcePath='/Volumes/BHL_WD18TB/20230815_20230903_VRs';
cd(sourcePath)
[fpath] = uigetfile_n_dir(); %only Treadmill data
%% Parameter setting and get cell coordinate
block_size=25;

for i=1:2%:length(fpath)
    clear bs ROI
    load(fullfile(fpath{i},"output_data.mat"))
    sz=double(Device_Data{1, 3}.ROI([2 4]));
    mov_test=double(readBinMov_times([fpath{i} '/frames1.bin'],sz(2),sz(1),[3000:4000]));
    avgImg=mean(mov_test,3);
    [CellCoordinate, ~]=Cell_segment_circle_25x_VU(avgImg,0.85);
    CellCoordinate=cell_detection_manual(avgImg,CellCoordinate,[]);
    figure(1); clf;
    imshow2(avgImg,[]); hold all
    for n=1:size(CellCoordinate,1)
        ROI(n,:)=round([CellCoordinate(n,1)-block_size CellCoordinate(n,1)+block_size ...
            CellCoordinate(n,2)-block_size CellCoordinate(n,2)+block_size]);
        bs(n)=block_size;

        % if ROI go further than outbound, reduce window size
        while sum([ROI(n,1) > 1, ROI(n,2)<sz(1), ROI(n,3)> 1, ROI(n,4)<sz(2)])<4
            bs(n)=bs(n)-1;
            ROI(n,:)=round([CellCoordinate(n,1)-bs(n) CellCoordinate(n,1)+bs(n) ...
                CellCoordinate(n,2)-bs(n) CellCoordinate(n,2)+bs(n)]);
        end
        ROI_box=[ROI(n,[1 3]); ROI(n,[1 4]); ROI(n,[2 4]); ROI(n,[2 3]); ROI(n,[1 3])];
        plot(ROI_box(:,1),ROI_box(:,2),'r')
    end

    plot(CellCoordinate(:,1),CellCoordinate(:,2),'ro','markersize',15)
    saveas(gca,[fpath{i} '/ROIs.png'])
    save([fpath{i} '/Analysis_parameter.mat'],'CellCoordinate','bs','ROI','-v7.3')
end
%% motion correction
time_size=150000; %segment size

for i=2:length(fpath)
    clear mcTrace_block ave_im
    %load the center positions
    load(fullfile(fpath{i},"output_data.mat"))
    load(fullfile(fpath{i},'Analysis_parameter.mat'))
    sz=double(Device_Data{1, 3}.ROI([2 4]));
    CamCounter=Device_Data{1, 2}.Counter_Inputs(1, 1).data;
    CamTrigger=find(CamCounter(2:end)-CamCounter(1:end-1));

    ref_time=[10000:11000];

    % make time segment
    t_seg=[[1:time_size:length(CamTrigger)] length(CamTrigger)+1];
    t_seg=[t_seg(1:end-1)' t_seg(2:end)'+100];
    t_seg(end)=length(CamTrigger);

    for n=1:size(CellCoordinate,1)

        % load ROI reference
        mov_seg=double(readBinMov_times_ROI([fpath{i} '/frames1.bin'],sz(1),sz(2),ref_time,ROI(n,:)));

        % motion correction of reference
        options_rigid = NoRMCorreSetParms('d1',size(mov_seg,1),'d2',size(mov_seg,2),'bin_width',200,'max_shift',bs(n),'us_fac',50,'init_batch',200);
        tic; [mov_seg,shifts1,template1,options_rigid] = normcorre(mov_seg,options_rigid); toc

        % make reference image
        mov_seg=  vm(mov_seg);
        mov_seg = single(mov_seg)./single(max(mov_seg.data(:)));
        mov_seg = movmean(mov_seg,10,3);
        mov_ref = squeeze(median(mov_seg,3));

        % load time segment of ROI
        for t=1:size(t_seg,1)
            mcTrace_block{n,t}=[];
            mov=vm(readBinMov_times_ROI([fpath{i} '/frames1.bin'],sz(1),sz(2),[t_seg(t,1):t_seg(t,2)],ROI(n,:)));
            [mov_mc,xyField]=optical_flow_motion_correction_LBH_ROIBlock(mov,double(mov_ref),'normcorre');

            ave_im{n,t}=mean(mov_mc,3);
            mov_mc=vm(mov_mc);
            mov_mc.transpose.savebin([fpath{i} '/mc' num2str(n,'%02d') '_' num2str(t,'%02d') '.bin'])

            %mcTrace_block{n,t} = squeeze(mean (xyField,[1 2]))'; %optic flow
            mcTrace_block{n,t}=xyField'; % Normcorre
        end


    end
    save([fpath{i} '/Analysis_parameter.mat'],'CellCoordinate','mcTrace_block','ave_im','ROI','bs','-v7.3')
    disp(['MC done:' fpath{i}])

end

%% Signal Extraction

time_size=150000; %segment size

for i=1:length(fpath)
    %load files
    load(fullfile(fpath{i},'Analysis_parameter.mat'))
    load([fpath{i} '/output_data.mat'])
    fileList = dir(fullfile(fpath{i}, '*.data'));
    if length(fileList)==1
        fid = fopen(fullfile(fpath{i},fileList.name));
        VRdata = fread(fid,[12 inf],'double');
    else
        error('Data file cannot be found');
    end

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

    Result{i}.VR=Virmen_data_int;
    Result{i}.centers=CellCoordinate;
    Result{i}.frm_rate=double((CamTrigger(2)-CamTrigger(1))*DAQ_rate);
    Result{i}.fpath=fpath{i};
    disp(['loading  ' fpath{i}])

    sz=double(Device_Data{1, 3}.ROI([2 4]));
    mov_test=double(readBinMov_times([fpath{i} '/frames1.bin'],sz(2),sz(1),[3000:4000]));
    Result{i}.FOV=mean(mov_test,3);

    try
        Result{i}.Blue=Device_Data{1, 2}.buffered_tasks(1, 1).channels(1, 2).data(CamTrigger).*Virmen_data_int(12,:);
    end
    try
        Result{i}.Reward=Device_Data{1, 2}.buffered_tasks(1, 3).channels.data(CamTrigger);
    end

    t_seg=[[1:time_size:length(CamTrigger)] length(CamTrigger)+1];
    t_seg=[t_seg(1:end-1)' t_seg(2:end)'+100];
    t_seg(end)=length(CamTrigger);

    Result{i}.traces=[];
    Result{i}.traces_res=[];
    Result{i}.im_corr=[];
    Result{i}.mcTrace=[];

    for n=1:size(Result{i}.centers,1)

        mov=double(readBinMov([fpath{i} '/mc' num2str(n,'%02d') '_' num2str(1,'%02d') '.bin'], ...
            bs(n)*2+1,bs(n)*2+1));
        mov_res= mov-mean(mov,3);
        bkg = zeros(2, size(mov,3));
        bkg(1,:) = linspace(-1, 1, size(mov,3));  % linear term
        bkg(2,:) = linspace(-1, 1, size(mov,3)).^2;  % quadratic term
        mov_res=SeeResiduals(mov_res,mcTrace_block{n,1});
        mov_res=SeeResiduals(mov_res,mcTrace_block{n,1}.^2);
        mov_res=SeeResiduals(mov_res,mcTrace_block{n,1}(1,:).*mcTrace_block{n,1}(2,:));
        mov_res= SeeResiduals(mov_res,bkg,1);


        Result{i}.ref_im{n}=mean(mov,3);
        ref_im_vec=tovec(Result{i}.ref_im{n});
        ref_im_vec=(ref_im_vec-mean(ref_im_vec,1))./std(ref_im_vec,0,1);
        mov_mc_vec=tovec(mov);
        mov_mc_vec=(mov_mc_vec-mean(mov_mc_vec,1))./std(mov_mc_vec,0,1);
        Result{i}.im_corr{n}=[];
        Result{i}.im_corr{n}=[Result{i}.im_corr{n} sum(mov_mc_vec.*ref_im_vec,1)/(size(mov_mc_vec,1)-1)];

        Result{i}.c_ftprnt(1:2*bs(n)+1,1:2*bs(n)+1,n)=mask_footprint([bs(n)+0.5 bs(n)+0.5],mov_res(:,:,1000:end-1000),[],20);
        Cellpsf = fspecial('gaussian', 2*bs(n)+1, 8);
        Result{i}.c_ftprnt(1:2*bs(n)+1,1:2*bs(n)+1,n)=Result{i}.c_ftprnt(1:2*bs(n)+1,1:2*bs(n)+1,n).*Cellpsf;
        Result{i}.c_ftprnt(1:2*bs(n)+1,1:2*bs(n)+1,n)=imgaussfilt(Result{i}.c_ftprnt(1:2*bs(n)+1,1:2*bs(n)+1,n),1);

        for t=2:size(t_seg,1)
            Result{i}.traces(n,t_seg(t-1,1):t_seg(t-1,2)-101)=-(tovec(mov_res(:,:,1:time_size))'*tovec(Result{i}.c_ftprnt(1:2*bs(n)+1,1:2*bs(n)+1,n)))';
            Result{i}.mcTrace(:,t_seg(t-1,1):t_seg(t-1,2)-101,n)=mcTrace_block{n,t-1}(:,1:time_size);

            mov=double(readBinMov([fpath{i} '/mc' num2str(n,'%02d') '_' num2str(t,'%02d') '.bin'], ...
                bs(n)*2+1,bs(n)*2+1));

            mov_res= mov-mean(mov,3);
            bkg = zeros(2, size(mov,3));
            bkg(1,:) = linspace(-1, 1, size(mov,3));  % linear term
            bkg(2,:) = linspace(-1, 1, size(mov,3)).^2;  % quadratic term
            mov_res=SeeResiduals(mov_res,mcTrace_block{n,t});
            mov_res=SeeResiduals(mov_res,mcTrace_block{n,t}.^2);
            mov_res= SeeResiduals(mov_res,bkg,1);

            mov_mc_vec=tovec(mov);
            mov_mc_vec=(mov_mc_vec-mean(mov_mc_vec,1))./std(mov_mc_vec,0,1);
            Result{i}.im_corr{n}=[Result{i}.im_corr{n} sum(mov_mc_vec.*ref_im_vec,1)/(size(mov_mc_vec,1)-1)];

        end
        Result{i}.traces(n,t_seg(end,1):t_seg(end,2))=-(tovec(mov_res)'*tovec(Result{i}.c_ftprnt(1:2*bs(n)+1,1:2*bs(n)+1,n)))';
        Result{i}.mcTrace(:,t_seg(end,1):t_seg(end,2),n)=mcTrace_block{n,end};


        mcT=squeeze(Result{i}.mcTrace(:,:,n));
        Result{i}.traces_res(n,:)=squeeze(SeeResiduals(reshape(Result{i}.traces(n,:),1,1,[]),mcT));
        Result{i}.traces_res(n,:)=squeeze(SeeResiduals(reshape(Result{i}.traces_res(n,:),1,1,[]),mcT.^2'));
        Result{i}.traces_res(n,:)=squeeze(SeeResiduals(reshape(Result{i}.traces_res(n,:),1,1,[]),mcT(1,:).*mcT(2,:)));
        Result{i}.AvgImg=mean(mov,3);
    end

end

for i=1:length(fpath)
    Result{i}.spike=zeros(size(Result{i}.traces));
    tmp=squeeze(Result{i}.traces_res) - movmedian(squeeze(Result{i}.traces_res),250,2); tmp=tmp./get_threshold(tmp,1);
    Result{i}.spike=(Result{i}.spike+find_spike_bh(tmp,5.5,2))>0;
end

save(fullfile(sourcePath,'20230908_Result.mat'),'Result','fpath','-v7.3')
%% Clean up and bleach correction
exclude_frq=[20 60];
time_bin=10000; Fs=1000;

freq_lowhigh=exclude_frq/(Fs/2);
[b, a] = butter(4, freq_lowhigh, 'stop');

for i=1%:length(fpath)
    tN=[1:time_bin:size(Result{i}.traces,2)]; tN=[tN size(Result{i}.traces,2)];
    for n=1:size(Result{i}.traces,1)

        mcTrace=squeeze(Result{i}.mcTrace(:,:,n));
            tr=Result{i}.traces_res(n,:);
            
            % regress out motion frequency
            for t=1:length(tN)-1
                tr(tN(t):tN(t+1))=squeeze(SeeResiduals(reshape(tr(tN(t):tN(t+1)),1,1,[]),mcTrace(:,(tN(t):tN(t+1)))));
                tr(tN(t):tN(t+1))=squeeze(SeeResiduals(reshape(tr(tN(t):tN(t+1)),1,1,[]),mcTrace(:,(tN(t):tN(t+1))).^2));
                tr(tN(t):tN(t+1))=squeeze(SeeResiduals(reshape(tr(tN(t):tN(t+1)),1,1,[]),mcTrace(1,(tN(t):tN(t+1))).*mcTrace(2,(tN(t):tN(t+1)))));
            end

        % regress out motion frequency
        traces_res_filtered(n,:) = filtfilt(b, a, tr);
        norm_trace(n,:)=traces_res_filtered(n,:);%-movmedian(traces_res_filtered(n,:),500,2);

        for t=1:length(tN)-1
            noise(t,n)=get_threshold(norm_trace(n,tN(t):tN(t+1)),1);
        end

        noise_intp(n,:)=interp1(tN(1:end-1),noise(:,n),[1:size(Result{i}.traces,2)]);
    end
    
    norm_trace=norm_trace./noise_intp;        
    Result{i}.normTraces=norm_trace;

end


%% Get PC















