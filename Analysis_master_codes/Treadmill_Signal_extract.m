clear
[fpath] = uigetfile_n_dir;

%%
for i=5%:length(fpath)
    load([fpath{i} '/output_data.mat'])
    sz=double(Device_Data{1, 3}.ROI([2 4]));
    ref_time=[2000:3000];
    mov_test=double(readBinMov_times([fpath{i} '/frames1.bin'],sz(2),sz(1),ref_time));
    [centers, radii]=Cell_segment_circle_10x_VU(avgImg,0.85);
    %Result{i}.centers=cell_detection_manual(mean(mov_test,3),centers,[0 9000]);
    Result{i}.centers=cell_detection_manual(mean(mov_test,3),Result{1}.centers+[6 -0.5],[0 9000]);
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
    Blue=Device_Data{1, 2}.buffered_tasks(1, 1).channels(1, 2).data;
    Reward=Device_Data{1, 2}.buffered_tasks(1, 3).channels.data;  
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
    
    mov_res= mov_mc-mean(mov_mc,3); 
    mov_res = SeeResiduals(mov_res,mcTrace);
    mov_res = SeeResiduals(mov_res,mcTrace.^2);
    mov_res = SeeResiduals(mov_res,mcTrace(:,1).*mcTrace(:,2));
    
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

    try mov_mc=mov_mc(:,:,1:end_frame); catch mov_mc=mov_mc; end
    try mcTrace=mcTrace(1:end_frame,:); catch mcTrace=mcTrace; end
    mov_res= mov_mc-mean(mov_mc,3); 
    mov_res = SeeResiduals(mov_res,mcTrace);
    mov_res = SeeResiduals(mov_res,mcTrace.^2);
    mov_res = SeeResiduals(mov_res,mcTrace(:,1).*mcTrace(:,2));

    Result{i}.traces=[Result{i}.traces -(tovec(mov_res)'*tovec(Result{i}.c_ftprnt))'];
    Result{i}.mcTrace=[Result{i}.mcTrace; mcTrace];

    end

    Result{i}.traces_res=squeeze(SeeResiduals(reshape(Result{i}.traces,N,1,[]),Result{i}.mcTrace'));
    Result{i}.traces_res=squeeze(SeeResiduals(reshape(Result{i}.traces_res,N,1,[]),Result{i}.mcTrace.^2'));
    Result{i}.traces_res=squeeze(SeeResiduals(reshape(Result{i}.traces_res,N,1,[]),(Result{i}.mcTrace(:,1).*Result{i}.mcTrace(:,2))'));

end

%%
for i=3:length(Result)
    fileList = dir(fullfile(fpath{i}, '*.data'));
    if length(fileList)==1
    fid = fopen(fullfile(fpath{i},fileList.name));
    VRdata = fread(fid,[10 inf],'double');
    else
        error('There is more than one .data file');
    end
    Result.trace_res
    tmp=squeeze(Result{i}.traces_res_hi) - movmedian(squeeze(Result{i}.traces_res_hi),4,2); tmp=tmp./get_threshold(tmp,1);
    Result{i}.spike=find_spike_bh(tmp,4,3);
    tmp=squeeze(Result{i}.traces_res_hi_filtered) - movmedian(squeeze(Result{i}.traces_res_hi_filtered),30,2); tmp=tmp./get_threshold(tmp,1);
    Result{i}.spike=(Result{i}.spike+find_spike_bh(tmp,3.2,1.5))>0;

end


