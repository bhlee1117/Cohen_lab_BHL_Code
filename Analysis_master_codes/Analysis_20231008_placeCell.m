clear
clc;
sourcePath='/Volumes/BHL_WD18TB/20231007_BHLm081_82_PlaceCell';
cd(sourcePath)
[fpath] = uigetfile_n_dir(); %only Treadmill data
%% Parameter setting and get cell coordinate
block_size=25;

for i=1:length(fpath)
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

for i=1:length(fpath)
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

        Result{i}.c_ftprnt(1:2*bs(n)+1,1:2*bs(n)+1,n)=mask_footprint([bs(n)+0.5 bs(n)+0.5],mov_res(:,:,1000:end-1000),[],18);
        Cellpsf = fspecial('gaussian', 2*bs(n)+1, 8);
        Result{i}.c_ftprnt(1:2*bs(n)+1,1:2*bs(n)+1,n)=Result{i}.c_ftprnt(1:2*bs(n)+1,1:2*bs(n)+1,n).*Cellpsf;
        Result{i}.c_ftprnt(1:2*bs(n)+1,1:2*bs(n)+1,n)=imgaussfilt(Result{i}.c_ftprnt(1:2*bs(n)+1,1:2*bs(n)+1,n),2);

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

save(fullfile(sourcePath,'20231008_Result.mat'),'Result','fpath','-v7.3')
%%
load(fullfile(sourcePath,'20231008_Result.mat'))
%% Clean up and bleach correction

exclude_frq=[241.7 242]; %monitor
%exclude_frq2=[483.5 484]; %monitor
exclude_frq2=[20 65]; %motion
time_bin=10000; Fs=1000;

freq_lowhigh=exclude_frq/(Fs/2);
[b, a] = butter(4, freq_lowhigh, 'stop');

freq_lowhigh2=exclude_frq2/(Fs/2);
[b2, a2] = butter(4, freq_lowhigh2, 'stop');

  pass_frq=[5];
freq_lowhigh3=pass_frq/(Fs/2);
[b3, a3] = butter(4, freq_lowhigh3, 'low');

for i=1:length(fpath)
    clear traces_res_filtered noise noise_intp norm_trace sp_height
    tN=[1:time_bin:size(Result{i}.traces,2)]; tN=[tN size(Result{i}.traces,2)];
    for n=1:size(Result{i}.traces,1)

        mcTrace=squeeze(Result{i}.mcTrace(:,:,n));
            tr=Result{i}.traces(n,:);
            
            % regress out motion frequency
            for t=1:length(tN)-1
                tr(tN(t):tN(t+1))=squeeze(SeeResiduals(reshape(tr(tN(t):tN(t+1)),1,1,[]),mcTrace(:,(tN(t):tN(t+1)))));
                tr(tN(t):tN(t+1))=squeeze(SeeResiduals(reshape(tr(tN(t):tN(t+1)),1,1,[]),mcTrace(:,(tN(t):tN(t+1))).^2));
                tr(tN(t):tN(t+1))=squeeze(SeeResiduals(reshape(tr(tN(t):tN(t+1)),1,1,[]),mcTrace(1,(tN(t):tN(t+1))).*mcTrace(2,(tN(t):tN(t+1)))));
            end

        % regress out motion frequency
        traces_res_filtered(n,:) = filtfilt(b, a, tr);
        traces_res_filtered(n,:) = filtfilt(b2, a2, traces_res_filtered(n,:));
        norm_trace(n,:)=traces_res_filtered(n,:);%-movmedian(traces_res_filtered(n,:),500,2);

       for t=1:length(tN)-1
        tr_tmp=norm_trace(n,tN(t):tN(t+1));
        tr_tmp=tr_tmp-movmedian(tr_tmp,300);
        noise(t,n)=get_threshold(tr_tmp,1);
        tr_tmp=tr_tmp./noise(t,n);
        [sp_temp, pks, prom]=find_spike_bh(tr_tmp,4,3);
        
        dec_fac=0.5;
        while isempty(pks)    
        [sp_temp, pks, prom]=find_spike_bh(tr_tmp,4-dec_fac,3-dec_fac);
        dec_fac=dec_fac+0.5;
        end

        %sp_height(t,n)=median(pks);%/noise(t,n);
        sp_height(t,n)=prctile(pks,80);%/noise(t,n);
    end
    tx=tN(1:end-1)+time_bin/2;
    noise_intp(n,:)=movmean(interp1(tx,noise(:,n),[1:size(Result{i}.traces,2)],'linear','extrap'),10000);

%     [Ny_fit t_consts coeffY]  = expfitDM_2(tx(~isnan(noise(:,n)))',noise(~isnan(noise(:,n)),n),[1:size(Result{i}.traces,2)]',10^7);
%     noise_intp(n,:)=Ny_fit;
    [y_fit t_consts coeffY]  = expfitDM_2(tx(~isnan(sp_height(:,n)))',sp_height(~isnan(sp_height(:,n)),n),[1:size(Result{i}.traces,2)]',10^7);
    SpHeight_intp(n,:)=y_fit;

    end
    
    norm_trace=norm_trace./noise_intp;        
    %norm_trace=norm_trace;%./(SpHeight_intp./SpHeight_intp(:,1));        
    Result{i}.normTraces=norm_trace./get_threshold(norm_trace,1);
    Result{i}.spike=find_spike_bh(Result{i}.normTraces-movmedian(Result{i}.normTraces,300,2),5,3);

    
for n=1:size(Result{i}.traces,1)
Result{i}.subThreshold(n,:) = filtfilt(b3, a3, Result{i}.normTraces(n,:));
%Result{i}.subThreshold(n,:) = movmean(Result{i}.normTraces(n,:),1000);
end

end

%%
    goi=[1 2]; place_bin=50;

    tr_cat=[]; sp_cat=[]; vr_cat=[];
    for g=goi
    [sz_fprntY, sz_fprntX, nNeurons]=size(Result{g}.c_ftprnt);
    sz=size(Result{g}.FOV);
    
    Result{g}.c_fprnt_merge=zeros(sz(1),sz(2),nNeurons);
    [~, n_sort]=sort(Result{g}.centers(:,1),'ascend');
    Result{g}.n_sort = n_sort;
    tr_cat=[tr_cat Result{g}.normTraces(n_sort,:)];
    sp_cat=[sp_cat Result{g}.spike(n_sort,:)];
    if ~isempty(vr_cat)
    addVR=Result{g}.VR;
    addVR(8,:)=vr_cat(8,end)+addVR(8,:);
    else
        addVR=Result{g}.VR;
    end
    vr_cat=[vr_cat addVR];
    
    for n=1:length(n_sort)
    offset=round(Result{g}.centers(n_sort(n),:)+[-floor(sz_fprntY/2):1:floor(sz_fprntY/2); -floor(sz_fprntX/2):1:floor(sz_fprntX/2)]');
    Result{g}.c_fprnt_merge(offset(:,2),offset(:,1),n)=Result{g}.c_ftprnt(:,:,n_sort(n))';
    [Result{g}.Lap_FR(:,:,n) Result{g}.Lap_V]=PlaceTrigger_average(Result{g}.spike(n_sort(n),:),place_bin,Result{g}.VR,0.003,115);
    [Result{g}.Lap_F(:,:,n) Result{g}.Lap_V]=PlaceTrigger_average(Result{g}.normTraces(n_sort(n),:),place_bin,Result{g}.VR,0.003,115);
    [Result{g}.Lap_sub(:,:,n) Result{g}.Lap_V]=PlaceTrigger_average(Result{g}.subThreshold(n_sort(n),:),place_bin,Result{g}.VR,0.003,115);
    end
    figure;
    show_footprnt(Result{g}.c_fprnt_merge,mat2gray(Result{g}.FOV))
    end 
        


    %%

        %lap_seg=[5 9;10 12; 13 17;26 28;30 32];
    lap_seg=[5 9;13 17;30 32];
    cmap=distinguishable_colors(size(lap_seg,1));
    show_shift=10; stim_bin=place_bin*0.3;
    legend_string={'Lap 5-9','Lap 13-17','Lap 30-32'};

        figure(2); clf
        tiledlayout(2,3)
        Lap_FR=[];
            for g=goi
            [~, Lap_start]=unique(Result{g}.VR(8,:));
            VRworld=Result{g}.VR(2,Lap_start);
            Lap_FR=[Lap_FR; Result{g}.Lap_FR(VRworld==2,:,:)];
            end
            Lap_FR_filt=movmean(repmat(Lap_FR,1,3),5,2,'omitnan');
            for n=1:nNeurons
                nexttile(n)
                imagesc(Lap_FR_filt(:,place_bin+1:place_bin*2,n))
                xlabel('position (cm)')
                set(gca,'xtick',[0:place_bin/2:place_bin],'xticklabel',[0:200/2:200])
                ylabel('Laps')
                colormap('turbo')
        
                nexttile(n+nNeurons)
                for i=1:size(lap_seg,1)
                plot(mean(Lap_FR_filt(lap_seg(i,1):lap_seg(i,2),place_bin-show_shift:place_bin*2+show_shift,n),1,'omitnan'),'color',cmap(i,:)); hold all
                end
                line([stim_bin stim_bin]+show_shift,[0 10],'color','r','linewidth',1)
                xlabel('position (cm)')
                set(gca,'xtick',show_shift+[0:place_bin/2:place_bin],'xticklabel',[0:200/2:200])
                ylabel('Firing Rate (Hz)')
                legend(legend_string)
                axis tight
            end
        
            
        figure(3); clf
        tiledlayout(2,3)
        Lap_F=[];
            for g=goi
            [~, Lap_start]=unique(Result{g}.VR(8,:));
            VRworld=Result{g}.VR(2,Lap_start);
            Lap_F=[Lap_F; Result{g}.Lap_F(VRworld==2,:,:)];
            end
            Lap_F_filt=movmean(repmat(Lap_F,1,3),5,2,'omitnan');
            for n=1:nNeurons
                nexttile(n)
                imagesc(Lap_F_filt(:,place_bin+1:place_bin*2,n),[-900 1500])
                xlabel('position (cm)')
                set(gca,'xtick',[0:place_bin/2:place_bin],'xticklabel',[0:200/2:200])
                ylabel('Laps')
                colormap('turbo')
        
                nexttile(n+nNeurons)
                for i=1:size(lap_seg,1)
                plot(mean(Lap_F_filt(lap_seg(i,1):lap_seg(i,2),place_bin-show_shift:place_bin*2+show_shift,n),1,'omitnan'),'color',cmap(i,:)); hold all
                end
                line([stim_bin stim_bin]+show_shift,[-700 700],'color','r','linewidth',1)
                xlabel('position (cm)')
                set(gca,'xtick',show_shift+[0:place_bin/2:place_bin],'xticklabel',[0:200/2:200])
                ylabel('\DeltaF')
                legend(legend_string)
            end
        
            figure(4); clf
        tiledlayout(2,3)
        Lap_sub=[];
            for g=goi
            [~, Lap_start]=unique(Result{g}.VR(8,:));
            VRworld=Result{g}.VR(2,Lap_start);
            Lap_sub=[Lap_sub; Result{g}.Lap_sub(VRworld==2,:,:)];
            end
            Lap_sub_filt=movmean(repmat(Lap_sub,1,3),5,2,'omitnan');
            for n=1:nNeurons
                nexttile(n)
                imagesc(Lap_sub_filt(:,place_bin+1:place_bin*2,n),[-900 1500])
                colormap('turbo')
                xlabel('position (cm)')
                set(gca,'xtick',[0:place_bin/2:place_bin],'xticklabel',[0:200/2:200])
                ylabel('Laps')
                nexttile(n+nNeurons)
        for i=1:size(lap_seg,1)
        plot(mean(Lap_sub_filt(lap_seg(i,1):lap_seg(i,2),place_bin-show_shift:place_bin*2+show_shift,n),1,'omitnan'),'color',cmap(i,:)); hold all
        end
        line([stim_bin stim_bin]+show_shift,[-700 700],'color','r','linewidth',1)
        xlabel('position (cm)')
                set(gca,'xtick',show_shift+[0:place_bin/2:place_bin],'xticklabel',[0:200/2:200])
                ylabel('\DeltaF')
                legend(legend_string)
    end

    show_traces_align_Position(tr_cat,sp_cat,34,[-3000:1500],vr_cat,[2])
%% find CS

    goi=[1 2]; place_bin=50;

spike_time=find(spike); spike_height=tr(spike_time);

pass_frq=[15];
Fs=1000;

freq_lowhigh=pass_frq/(Fs/2);
[b, a] = butter(4, freq_lowhigh, 'low');
tr_pass = filtfilt(b, a, tr);
figure(2); clf; show_t=[750000:850000];
plot(t(show_t),tr(show_t)); hold all
plot(t(show_t),tr_pass(show_t));
plot(t(show_t),zeros(length(show_t),1));

[trans tr_trace]=detect_transient(tr_pass,[3 0.5],spike);
CS_ind=find(trans.spike_number>2 & trans.mean_ISI<40);
CS_trace=ismember(tr_trace,CS_ind);
figure(3); clf;
plot(t,tr)
hold all
tr_nan=tr; tr_nan(CS_trace==0)=NaN;
plot(t,tr_nan,'r')

CS_trace_spike=spike.*bwlabel(CS_trace);




    
    for g=goi
    [sz_fprntY, sz_fprntX, nNeurons]=size(Result{g}.c_ftprnt);
    sz=size(Result{g}.FOV);
    
    Result{g}.c_fprnt_merge=zeros(sz(1),sz(2),nNeurons);
    [~, n_sort]=sort(Result{g}.centers(:,1),'ascend');
    Result{g}.n_sort = n_sort;
    tr_cat=[tr_cat Result{g}.normTraces(n_sort,:)];
    sp_cat=[sp_cat Result{g}.spike(n_sort,:)];
    if ~isempty(vr_cat)
    addVR=Result{g}.VR;
    addVR(8,:)=vr_cat(8,end)+addVR(8,:);
    else
        addVR=Result{g}.VR;
    end
    vr_cat=[vr_cat addVR];
    
    for n=1:length(n_sort)
    offset=round(Result{g}.centers(n_sort(n),:)+[-floor(sz_fprntY/2):1:floor(sz_fprntY/2); -floor(sz_fprntX/2):1:floor(sz_fprntX/2)]');
    Result{g}.c_fprnt_merge(offset(:,2),offset(:,1),n)=Result{g}.c_ftprnt(:,:,n_sort(n))';
    [Result{g}.Lap_FR(:,:,n) Result{g}.Lap_V]=PlaceTrigger_average(Result{g}.spike(n_sort(n),:),place_bin,Result{g}.VR,0.003,115);
    [Result{g}.Lap_F(:,:,n) Result{g}.Lap_V]=PlaceTrigger_average(Result{g}.normTraces(n_sort(n),:),place_bin,Result{g}.VR,0.003,115);
    [Result{g}.Lap_sub(:,:,n) Result{g}.Lap_V]=PlaceTrigger_average(Result{g}.subThreshold(n_sort(n),:),place_bin,Result{g}.VR,0.003,115);
    end
    figure;
    show_footprnt(Result{g}.c_fprnt_merge,mat2gray(Result{g}.FOV))
    end 

%%
i=1;
show_traces_spikes(Result{i}.normTraces,Result{i}.spike,[Result{i}.VR(5,:)/30; Result{i}.Blue]);
figure;
for n=1:size(Result{i}.Lap_FR,3)
nexttile([1 1])
fr=repmat(Result{i}.Lap_FR(:,:,n),1,3);
fr=movmean(fr,5,2);
s=size(Result{i}.Lap_FR(:,:,n));
imagesc(fr(:,s(2)+1:s(2)*2))
colormap('turbo')
title(num2str(n))
end










