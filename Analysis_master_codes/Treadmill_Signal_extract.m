clear
[fpath] = uigetfile_n_dir;

%%
for i=5%:length(fpath)
    load([fpath{i} '/output_data.mat'])
    sz=double(Device_Data{1, 3}.ROI([2 4]));
    ref_time=[2000:3000];
    mov_test=double(readBinMov_times([fpath{i} '/mc' num2str(1,'%02d') '.bin'],sz(2),sz(1),ref_time));
    avgImg=mean(mov_test,3);
    [centers, radii]=Cell_segment_circle_10x_VU(avgImg,0.85);
    %Result{i}.centers=cell_detection_manual(mean(mov_test,3),centers,[0 9000]);
    Result{i}.centers=cell_detection_manual(mean(mov_test,3),Result{1}.centers+[3 0],[0 9000]);
end


%% Signal extraction


DAQ_rate=0.000005;

for i=1:length(fpath)


    %load device data

    disp(['loading  ' fpath{i}])
    load([fpath{i} '/output_data.mat'])
    load(fullfile(fpath{i},'mcTrace01.mat'));

    frm_end=max(Device_Data{1, 2}.Counter_Inputs(1, 1).data);
    f_seg=[[1:10000:frm_end] frm_end+1];
    try
    Result{i}.Blue=Device_Data{1, 2}.buffered_tasks(1, 1).channels(1, 2).data;
    end
    Result{i}.Reward=Device_Data{1, 2}.buffered_tasks(1, 3).channels.data;
    CamCounter=Device_Data{1, 2}.Counter_Inputs(1, 1).data;
    CamTrigger=find(CamCounter(2:end)-CamCounter(1:end-1));
    frm_rate=double((CamTrigger(2)-CamTrigger(1))*DAQ_rate);
    sz=double(Device_Data{1, 3}.ROI([2 4]));
    mov_mc=double(readBinMov([fpath{i} '/mc' num2str(1,'%02d') '.bin'],sz(2),sz(1)));

    % initialize variables
    Result{i}.ref_im=mean(mov_mc,3);
    Result{i}.traces=[];
    Result{i}.traces_res=[];
    Result{i}.mcTrace=[];
    Result{i}.im_corr=[];

    shifts_r = squeeze(cat(3,mcTrace(:).shifts));
    shifts_nr = cat(ndims(mcTrace(1).shifts)+1,mcTrace(:).shifts);
    shifts_nr = reshape(shifts_nr,[],ndims(mov_mc)-1,size(mov_mc,3));
    shifts_x = squeeze(shifts_nr(:,2,:))';
    shifts_y = squeeze(shifts_nr(:,1,:))';
    mcTrace=[shifts_x shifts_y];

    mov_res= mov_mc-mean(mov_mc,3);
    mov_res = SeeResiduals(mov_res,mcTrace);
    mov_res = SeeResiduals(mov_res,mcTrace.^2);
    mov_res = SeeResiduals(mov_res,mcTrace(:,3).*mcTrace(:,9));

    % calculate footprint
    if frm_end<10000; end_frame=size(mov_mc,3); else end_frame=10000; end
    Result{i}.c_ftprnt=mask_footprint(Result{i}.centers,movmean(mov_res(:,:,1000:end),10,3),[],6);
    N=size(Result{i}.c_ftprnt,3);
    Result{i}.coord=get_coord(Result{i}.c_ftprnt);
    Result{i}.traces=[Result{i}.traces -(tovec(mov_res(:,:,1:end_frame))'*tovec(Result{i}.c_ftprnt))'];
    Result{i}.mcTrace=[Result{i}.mcTrace; mcTrace(1:end_frame,:)];

    % calculate correlation
    mov_mc_vec=tovec(mov_mc(:,:,1:end_frame)); mov_mc_vec=(mov_mc_vec-mean(mov_mc_vec,1))./std(mov_mc_vec,0,1);
    ref_im_vec=tovec(Result{i}.ref_im); ref_im_vec=(ref_im_vec-mean(ref_im_vec))./std(ref_im_vec);
    Result{i}.im_corr=sum(mov_mc_vec.*ref_im_vec,1)/(size(mov_mc_vec,1)-1);

    for j=2:length(f_seg)-1
        mov_mc=double(readBinMov([fpath{i} '/mc' num2str(j,'%02d') '.bin'],sz(2),sz(1)));
        load([fpath{i} '/mcTrace' num2str(j,'%02d') '.mat']);

        shifts_r = squeeze(cat(3,mcTrace(:).shifts));
        shifts_nr = cat(ndims(mcTrace(1).shifts)+1,mcTrace(:).shifts);
        shifts_nr = reshape(shifts_nr,[],ndims(mov_mc)-1,size(mov_mc,3));
        shifts_x = squeeze(shifts_nr(:,2,:))';
        shifts_y = squeeze(shifts_nr(:,1,:))';
        mcTrace=[shifts_x shifts_y];

        try mov_mc=mov_mc(:,:,1:end_frame); catch mov_mc=mov_mc; end
        try mcTrace=mcTrace(1:end_frame,:); catch mcTrace=mcTrace; end
        mov_res= mov_mc-mean(mov_mc,3);
        mov_res = SeeResiduals(mov_res,mcTrace);
        mov_res = SeeResiduals(mov_res,mcTrace.^2);
        mov_res = SeeResiduals(mov_res,mcTrace(:,3).*mcTrace(:,9));

        Result{i}.traces=[Result{i}.traces -(tovec(mov_res)'*tovec(Result{i}.c_ftprnt))'];
        Result{i}.mcTrace=[Result{i}.mcTrace; mcTrace];

    end

    Result{i}.traces=Result{i}.traces(:,1:length(CamTrigger));
    Result{i}.mcTrace=Result{i}.mcTrace(1:length(CamTrigger),:);

    Result{i}.traces_res=squeeze(SeeResiduals(reshape(Result{i}.traces,N,1,[]),Result{i}.mcTrace'));
    Result{i}.traces_res=squeeze(SeeResiduals(reshape(Result{i}.traces_res,N,1,[]),Result{i}.mcTrace.^2'));
    Result{i}.traces_res=squeeze(SeeResiduals(reshape(Result{i}.traces_res,N,1,[]),(Result{i}.mcTrace(:,3).*Result{i}.mcTrace(:,6))'));

    

end

%%

figure;
tiledlayout((length(Result)-1)*2,1)

for i=2:length(Result)
    fileList = dir(fullfile(fpath{i}, '*.data'));
    if length(fileList)==1
        fid = fopen(fullfile(fpath{i},fileList.name));
        VRdata = fread(fid,[12 inf],'double');
    else
        error('There is more than one .data file');
    end

    VRdata(:,1:4078)=[]; VRdata(:,end)=[];
    VRdata=VRdata([6 3 2 5 10 4 7 1 9 11 12],:);
    load([fpath{i} '/output_data.mat'])

    l=find(VRdata(8,2:end)-VRdata(8,1:end-1)>0);
    lap=[[1; (l+1)'] [(l)'; size(VRdata,2)]];
    cumTrack=[];
    VRdata(5,:)=VRdata(5,:)-min(VRdata(5,:));
    cumTrack=[cumTrack VRdata(5,lap(1,1):lap(1,2))];
    for l2=2:size(lap,1)
    cumTrack=[cumTrack VRdata(5,lap(l2,1):lap(l2,2))+cumTrack(lap(l2,1)-1)];    
    end
 
    DAQ_rew=rescale(Device_Data{1, 2}.buffered_tasks(1, 3).channels.data)>0.7;
    RewardOn_DAQ=find([0 (DAQ_rew(2:end)-DAQ_rew(1:end-1))==1]);
    RewardOn_VR=find([0 (VRdata(6,2:end)-VRdata(6,1:end-1))==1]);
   
    NewTimebin=[];
    for tt=1:length(RewardOn_VR)-1
        N=length(RewardOn_DAQ(tt):RewardOn_DAQ(tt+1));
        Gap=(RewardOn_VR(tt+1)-RewardOn_VR(tt))/(N-1);
        NewTimebin=[NewTimebin [RewardOn_VR(tt):Gap:RewardOn_VR(tt+1)]];
    end

        Pos=interp1([1:size(VRdata,2)],VRdata(5,:),NewTimebin,'linear','extrap');
        Lap=round(interp1([1:size(VRdata,2)],VRdata(8,:),NewTimebin,'linear','extrap'));
        CumPos=round(interp1([1:size(VRdata,2)],cumTrack,NewTimebin,'linear','extrap'));
        
        truncate=length(NewTimebin);
        DAQ_Vir=[t_DAQ(1:truncate); DAQ_rew(1:truncate); Pos; Lap; CumPos];
    

    CamCounter=Device_Data{1, 2}.Counter_Inputs(1, 1).data;
    CamTrigger=find(CamCounter(2:end)-CamCounter(1:end-1));


    Result{i}.Virmen=DAQ_Vir(:,CamTrigger);

    tmp=squeeze(Result{i}.traces_res) - movmedian(squeeze(Result{i}.traces_res),10,2); tmp=tmp./get_threshold(tmp,1);
    Result{i}.spike=find_spike_bh(tmp,4,3);
    tmp=squeeze(Result{i}.traces_res) - movmedian(squeeze(Result{i}.traces_res),150,2); tmp=tmp./get_threshold(tmp,1);
    Result{i}.spike=(Result{i}.spike+find_spike_bh(tmp,3.2,1.5))>0;

end






%%
VirDAQ_ratio=1.027250016585005;
figure;
tiledlayout((length(Result)-1)*2,1)


for i=2:length(Result)
    fileList = dir(fullfile(fpath{i}, '*.data'));
    if length(fileList)==1
        fid = fopen(fullfile(fpath{i},fileList.name));
        VRdata = fread(fid,[12 inf],'double');
    else
        error('There is more than one .data file');
    end

    VRdata(:,1:4078)=[]; VRdata(:,end)=[];
    VRdata=VRdata([6 3 2 5 10 4 7 1 9 11 12],:);
    load([fpath{i} '/output_data.mat'])

    l=find(VRdata(8,2:end)-VRdata(8,1:end-1)>0);
    lap=[[1; (l+1)'] [(l)'; size(VRdata,2)]];
    cumTrack=[];
    VRdata(5,:)=VRdata(5,:)-min(VRdata(5,:));
    cumTrack=[cumTrack VRdata(5,lap(1,1):lap(1,2))];
    for l2=2:size(lap,1)
    cumTrack=[cumTrack VRdata(5,lap(l2,1):lap(l2,2))+cumTrack(lap(l2,1)-1)];    
    end
 
    DAQ_rew=rescale(Device_Data{1, 2}.buffered_tasks(1, 3).channels.data)>0.7;
    DAQ_rate=Device_Data{1, 2}.buffered_tasks(1, 3).rate;
    t_DAQ=[1:length(DAQ_rew)]/DAQ_rate;

    t_VR = datetime(datetime(VRdata(1,:),'ConvertFrom','datenum'), 'InputFormat', 'yyyy-MM-dd HH:mm:ss.SSS');
    t_VR= t_VR-t_VR(1);
    t_VR= milliseconds(t_VR)/1000;
    VR_pos_interp=interp1(t_VR*VirDAQ_ratio,VRdata(5,:),t_DAQ,'linear');
    %VR_pos_interp=interp1(t_VR,VRdata(5,:),t_DAQ,'linear');
    VR_rew=[0 (VRdata(6,2:end)-VRdata(6,1:end-1))>0];
    VR_rew_interp=interp1(t_VR,VR_rew,t_DAQ,'linear')==1;
    VR_rew_interp=[0 (VR_rew_interp(2:end)-VR_rew_interp(1:end-1))==1];
    rewind=find(VR_rew_interp==1)'+[0:9000]; 
    VR_rew_interp(rewind(:))=1;

    VR_reward=zeros(1,length(VR_pos_interp));

    bw=bwlabel(VR_pos_interp>110*0.8);
    for b=1:max(bw)
        tmp=find(bw==b);
        VR_reward(tmp(1)+[0:500])=1;
    end

    [DAQ_VirCorr lag]=xcorr(VR_rew_interp,DAQ_rew,1e7);
    [~, maxLag]=max(DAQ_VirCorr);
    DAQ_Virlag=lag(maxLag)/DAQ_rate;
    Pos=interp1(t_VR*VirDAQ_ratio-DAQ_Virlag,VRdata(5,:),t_DAQ,'linear');
    Lap=round(interp1(t_VR*VirDAQ_ratio-DAQ_Virlag,VRdata(8,:),t_DAQ,'linear'));
    CumPos=round(interp1(t_VR*VirDAQ_ratio-DAQ_Virlag,cumTrack,t_DAQ,'linear'));
    
    DAQ_Vir=[t_DAQ; DAQ_rew; Pos; Lap; CumPos];
    

    CamCounter=Device_Data{1, 2}.Counter_Inputs(1, 1).data;
    CamTrigger=find(CamCounter(2:end)-CamCounter(1:end-1));


    Result{i}.Virmen=DAQ_Vir(:,CamTrigger);

    tmp=squeeze(Result{i}.traces_res) - movmedian(squeeze(Result{i}.traces_res),10,2); tmp=tmp./get_threshold(tmp,1);
    Result{i}.spike=find_spike_bh(tmp,4,3);
    tmp=squeeze(Result{i}.traces_res) - movmedian(squeeze(Result{i}.traces_res),150,2); tmp=tmp./get_threshold(tmp,1);
    Result{i}.spike=(Result{i}.spike+find_spike_bh(tmp,3.2,1.5))>0;

end


