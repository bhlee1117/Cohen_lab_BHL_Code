clear
[fpath] = uigetfile_n_dir;
save_name='pcResult_20230714.mat';
%%
for i=14%:length(fpath)
    disp(fpath{i})
    load([fpath{i} '/output_data.mat'])
    sz=double(Device_Data{1, 3}.ROI([2 4]));
    ref_time=[2000:3000];
    mov_test=double(readBinMov_times([fpath{i} '/mc' num2str(1,'%02d') '.bin'],sz(2),sz(1),ref_time));
    avgImg=mean(mov_test,3);
    [centers, radii]=Cell_segment_circle_10x_VU(avgImg,0.85);
    %Result{i}.centers=cell_detection_manual(mean(mov_test,3),centers,[0 6000]);
    shift=[1 0]; %[x y]
    Result{i}.centers=cell_detection_manual(mean(mov_test,3),Result{i-1}.centers+shift,[0 9000]); 
end
save(save_name,'Result','fpath','-v7.3')

%% Signal extraction


DAQ_rate=0.000005;

for i=1%:length(fpath)


    %load device data

    disp(['loading  ' fpath{i}])
    load([fpath{i} '/output_data.mat'])
    load(fullfile(fpath{i},'mcTrace01.mat'));

    frm_end=max(Device_Data{1, 2}.Counter_Inputs(1, 1).data);
    f_seg=[[1:10000:frm_end] frm_end+1];
    try
        Result{i}.Blue=Device_Data{1, 2}.buffered_tasks(1, 1).channels(1, 2).data;
    end
    try
        Result{i}.Reward=Device_Data{1, 2}.buffered_tasks(1, 3).channels.data;
    end
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

    %     shifts_r = squeeze(cat(3,mcTrace(:).shifts));
    %     shifts_nr = cat(ndims(mcTrace(1).shifts)+1,mcTrace(:).shifts);
    %     shifts_nr = reshape(shifts_nr,[],ndims(mov_mc)-1,size(mov_mc,3));
    %     shifts_x = squeeze(shifts_nr(:,2,:))';
    %     shifts_y = squeeze(shifts_nr(:,1,:))';
    %     mcTrace=[shifts_x shifts_y];

    mov_res= mov_mc-mean(mov_mc,3);
    mov_res = SeeResiduals(mov_res,mcTrace);
    mov_res = SeeResiduals(mov_res,mcTrace.^2);
    mov_res = SeeResiduals(mov_res,mcTrace(:,1).*mcTrace(:,2));

    % calculate footprint
    if frm_end<10000; end_frame=size(mov_mc,3); else end_frame=10000; end
    Result{i}.c_ftprnt=mask_footprint(Result{i}.centers,movmean(mov_res(:,:,1000:end),10,3),[],6);
    for n=1:size(Result{i}.c_ftprnt,3)
        Result{i}.c_ftprnt(:,:,n)=imgaussfilt(Result{i}.c_ftprnt(:,:,n),0.6);
    end
    N=size(Result{i}.c_ftprnt,3);
    Result{i}.coord=get_coord(Result{i}.c_ftprnt);
    Result{i}.traces=[Result{i}.traces -(tovec(mov_res(:,:,1:end_frame))'*tovec(Result{i}.c_ftprnt))'];
    Result{i}.mcTrace=[Result{i}.mcTrace; mcTrace(1:end_frame,:)];

    % calculate correlation
    mov_mc_vec=tovec(mov_mc(:,:,1:end_frame)); mov_mc_vec=(mov_mc_vec-mean(mov_mc_vec,1))./std(mov_mc_vec,0,1);
    ref_im_vec=tovec(Result{i}.ref_im); ref_im_vec=(ref_im_vec-mean(ref_im_vec))./std(ref_im_vec);
    Result{i}.im_corr=[Result{i}.im_corr sum(mov_mc_vec.*ref_im_vec,1)/(size(mov_mc_vec,1)-1)];

    for j=2:length(f_seg)-1
        mov_mc=double(readBinMov([fpath{i} '/mc' num2str(j,'%02d') '.bin'],sz(2),sz(1)));
        load([fpath{i} '/mcTrace' num2str(j,'%02d') '.mat']);

        %         shifts_r = squeeze(cat(3,mcTrace(:).shifts));
        %         shifts_nr = cat(ndims(mcTrace(1).shifts)+1,mcTrace(:).shifts);
        %         shifts_nr = reshape(shifts_nr,[],ndims(mov_mc)-1,size(mov_mc,3));
        %         shifts_x = squeeze(shifts_nr(:,2,:))';
        %         shifts_y = squeeze(shifts_nr(:,1,:))';
        %         mcTrace=[shifts_x shifts_y];

        try mov_mc=mov_mc(:,:,1:end_frame); catch mov_mc=mov_mc; end
        try mcTrace=mcTrace(1:end_frame,:); catch mcTrace=mcTrace; end
        mov_mc_vec=tovec(mov_mc);
        mov_mc_vec=(mov_mc_vec-mean(mov_mc_vec,1))./std(mov_mc_vec,0,1);

        mov_res= mov_mc-mean(mov_mc,3);
        mov_res = SeeResiduals(mov_res,mcTrace);
        mov_res = SeeResiduals(mov_res,mcTrace.^2);
        mov_res = SeeResiduals(mov_res,mcTrace(:,1).*mcTrace(:,2));

        Result{i}.traces=[Result{i}.traces -(tovec(mov_res)'*tovec(Result{i}.c_ftprnt))'];
        Result{i}.mcTrace=[Result{i}.mcTrace; mcTrace];
        Result{i}.im_corr=[Result{i}.im_corr sum(mov_mc_vec.*ref_im_vec,1)/(size(mov_mc_vec,1)-1)];

    end

    Result{i}.traces=Result{i}.traces(:,1:length(CamTrigger));
    Result{i}.mcTrace=Result{i}.mcTrace(1:length(CamTrigger),:);

    Result{i}.traces_res=squeeze(SeeResiduals(reshape(Result{i}.traces,N,1,[]),Result{i}.mcTrace'));
    Result{i}.traces_res=squeeze(SeeResiduals(reshape(Result{i}.traces_res,N,1,[]),Result{i}.mcTrace.^2'));
    Result{i}.traces_res=squeeze(SeeResiduals(reshape(Result{i}.traces_res,N,1,[]),(Result{i}.mcTrace(:,1).*Result{i}.mcTrace(:,2))'));



end

save(save_name,'Result','fpath','-v7.3')


%%
lap_dist=115;

for i=1:length(Result)

    fileList = dir(fullfile(fpath{i}, '*.data'));
    if ~isempty(fileList)
        figure;
        ax1=[];
        tiledlayout(2,1);
        if length(fileList)==1
            fid = fopen(fullfile(fpath{i},fileList.name));
            VRdata = fread(fid,[12 inf],'double');
        else
            error('There is more than one .data file');
        end

        load([fpath{i} '/output_data.mat'])

        % make cumTrack and segment laps
        l=find(VRdata(8,2:end)-VRdata(8,1:end-1)>0);
        lap=[[1; (l+1)'] [(l)'; size(VRdata,2)]];
        cumTrack=[];
        VRdata(5,:)=VRdata(5,:)-min(VRdata(5,:));
        cumTrack=[cumTrack VRdata(5,lap(1,1):lap(1,2))];
        for l2=2:size(lap,1)
            %cumTrack=[cumTrack VRdata(5,lap(l2,1):lap(l2,2))+cumTrack(lap(l2,1)-1)];
            cumTrack=[cumTrack VRdata(5,lap(l2,1):lap(l2,2))+lap_dist*(l2-1)];
        end

        % load from DAQ
        DAQ_rew=rescale(Device_Data{1, 2}.buffered_tasks(1, 3).channels.data)>0.7;
        DAQ_rate=Device_Data{1, 2}.buffered_tasks(1, 3).rate;
        t_DAQ=[1:length(DAQ_rew)]/DAQ_rate;
        CamCounter=Device_Data{1, 2}.Counter_Inputs(1, 1).data;
        CamTrigger=find(CamCounter(2:end)-CamCounter(1:end-1));

        VRdata_IsAc=VRdata(:,find(VRdata(10,:)==1));
        cumTrack=cumTrack(find(VRdata(10,:)==1));
        t_VR = datetime(datetime(VRdata_IsAc(1,:),'ConvertFrom','datenum'), 'InputFormat', 'yyyy-MM-dd HH:mm:ss.SSS');
        t_VR= t_VR-t_VR(1);
        t_VR= milliseconds(t_VR)/1000;
        t_VR= t_VR*t_DAQ(end)/max(t_VR);

        clear VRdata_IsAc_itp
        VRdata_IsAc_itp(1,:)=t_DAQ;

        for v=2:size(VRdata_IsAc,1)
            VRdata_IsAc_itp(v,:)=interp1(t_VR,VRdata_IsAc(v,:),t_DAQ,'linear','extrap');
        end
        VRdata_IsAc_itp(end+1,:)=interp1(t_VR,cumTrack,t_DAQ,'linear','extrap');

        ax1=[ax1 nexttile([1 1])];
        plot(t_DAQ,DAQ_rew,t_VR,VRdata_IsAc(11,:),'--')
        legend('Rig','Virmen')
        title(fpath{i},'Interpreter','none')
        ax1=[ax1 nexttile([1 1])];
        for v=[5 11 12]
            plot(VRdata_IsAc_itp(1,:),rescale(VRdata_IsAc_itp(v,:)));
            hold all
        end
        legend('Pos','Reward','blue')
        linkaxes(ax1,'x');
        Result{i}.Virmen=VRdata_IsAc_itp(:,CamTrigger)';
        
    end

end

save(save_name,'Result','fpath','-v7.3')

