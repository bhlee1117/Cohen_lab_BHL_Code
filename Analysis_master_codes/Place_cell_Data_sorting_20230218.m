
%addpath(genpath('/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Code'));
clear
[fpath] = uigetfile_n_dir();
[fpath_const] = uigetfile_n_dir();

%%
for i=5%:length(fpath_const)

    Sz = importdata([fpath_const{i} '/experimental_parameters.txt']);
    load(fullfile(fpath_const{i},'mcTrace.mat'));
    sz1=Sz.data(1); sz2=Sz.data(2);

    mov_mc=double(readBinMov_times([fpath_const{i} '/mc.bin'],sz2,sz1,[2000:5000]));
    im_G=imgaussfilt(mean(mov_mc,3),2);
    [centers radii]=Cell_segment_circle_2x(im_G,0.93);
    %Result_opto{i}.centers=cell_detection_manual(mean(mov_mc,3),centers,[0 9000]);
    Result_opto{i}.centers=cell_detection_manual(mean(mov_mc,3),Result{1}.centers+[-4 -18],[0 9000]);
end
%%
for i=3%:length(fpath)

    Sz = importdata([fpath{i} '/experimental_parameters.txt']);
    load(fullfile(fpath{i},'settings.mat'));
    sz1=Sz.data(1); sz2=Sz.data(2);
    [a,b]=system(sprintf('GetFileInfo "%s"',fullfile(fpath{i},'settings.mat'))); s=strfind(b,'modified:')+10; crdat=b(s:s+18);
    image_start=datestr(datenum(crdat));
    dt=datetime(image_start(end-8:end),'InputFormat','HH:mm:ss'); dt.Format = 'HH:mm:ss';
    fnm=dir(fullfile(fpath{i}, '*.csv'));

    a=find(DAQ_waves.amplitude(1,[1:500])); frm_rate=1/((a(2)-a(1)-2)*1e-5);
    frm_end=sum(DAQ_waves.amplitude(1,:)); f_seg=[[1:10000:frm_end] frm_end+1];
    [Arduino_data reward_pos lap_dist]=match_treadmill_DAQ(fullfile(fpath{i},fnm(2).name),1e-5,DAQ_data,(1/frm_rate+0.00002),2); % time, treadmill, Reward, Run
    Arduino_data(:,1)=Arduino_data(:,1)*(1/frm_rate*1000)/(1/frm_rate*1000+0.02);
    aa=find(bwlabel(Arduino_data(:,4))==1);

    frm_end=sum(DAQ_waves.amplitude(1,:)); f_seg=[[1:10000:frm_end] frm_end+1];

    mov_mc=double(readBinMov_times([fpath{i} '/mc01.bin'],sz2,sz1,[8000:8600]));
    im_G=imgaussfilt(mean(mov_mc,3),2);
    [centers radii]=Cell_segment_circle_2x(im_G,0.93);
    %Result{i}.centers=cell_detection_manual(mean(mov_mc,3),centers,[0 9000]);
    %Result{i}.centers=cell_detection_manual(mean(mov_mc,3),Result{1}.centers+[-2 -2],[0 9000]);
    Result{i}.centers=cell_detection_manual(mean(mov_mc,3),Result{1}.centers+[-1.5 -3],[0 9000]);
end


%%
for i=2:3%3:length(fpath_const)
    load(fullfile(fpath_const{i},'settings.mat'));
    load(fullfile(fpath_const{i},'mcTrace.mat'));
    Sz = importdata([fpath_const{i} '/experimental_parameters.txt']);
    sz1=Sz.data(1); sz2=Sz.data(2);
    a=find(DAQ_waves.amplitude(1,[1:500])); frm_rate=1/((a(2)-a(1)-2)*1e-5);

    mov_mc=double(readBinMov([fpath_const{i} '/mc.bin'],sz2,sz1));
    mov_res= mov_mc-mean(mov_mc,3);

    mov_res = SeeResiduals(mov_res,mcTrace);
    mov_res = SeeResiduals(mov_res,mcTrace.^2);
    mov_res = SeeResiduals(mov_res,mcTrace(:,1).*mcTrace(:,2));

    Result_opto{i}.c_ftprnt=mask_footprint(Result_opto{i}.centers,movmean(mov_res(:,:,1000:end),10,3),[],6);
    Result_opto{i}.ref_im=mean(mov_mc,3);
    Result_opto{i}.traces=[]; Result_opto{i}.traces_hi=[]; Result_opto{i}.traces_bin=[]; Result_opto{i}.traces_bin_hi=[];
    Result_opto{i}.traces_res=[]; Result_opto{i}.traces_res_hi=[]; Result_opto{i}.mcTrace=[];
    Result_opto{i}.im_corr=[];

    Result_opto{i}.frm_rate=frm_rate;
    Result_opto{i}.DMDPatt=dmd_mask_sequence_rois;
    Result_opto{i}.coord=get_coord(Result_opto{i}.c_ftprnt);
    Result_opto{i}.traces=[Result_opto{i}.traces -(tovec(mov_mc)'*tovec(Result_opto{i}.c_ftprnt))'];
    Result_opto{i}.traces_bin=[Result_opto{i}.traces_bin -(tovec(mov_mc)'*tovec(Result_opto{i}.c_ftprnt>0))'];
    mov_mc_vec=tovec(mov_mc); mov_mc_vec=(mov_mc_vec-mean(mov_mc_vec,1))./std(mov_mc_vec,0,1);
    ref_im_vec=tovec(Result_opto{i}.ref_im); ref_im_vec=(ref_im_vec-mean(ref_im_vec))./std(ref_im_vec);
    Result_opto{i}.im_corr=sum(mov_mc_vec.*ref_im_vec,1)/(size(mov_mc_vec,1)-1);
    Result_opto{i}.mcTrace=[Result_opto{i}.mcTrace; mcTrace];

    [V, D, W] = eig(corrcoef(zscore(Result_opto{i}.traces,0,2)'));
    D = diag(D); D = D(end:-1:1); V = V(:,end:-1:1);
    PCA_trace=zscore(V(:,1)'*Result_opto{i}.traces,0,2);
    ns=size(Result_opto{i}.traces,1);

    Result_opto{i}.traces_res=squeeze(SeeResiduals(reshape(Result_opto{i}.traces,ns,1,[]),Result_opto{i}.mcTrace'));
    Result_opto{i}.traces_res=squeeze(SeeResiduals(reshape(Result_opto{i}.traces_res,ns,1,[]),Result_opto{i}.mcTrace.^2'));
    Result_opto{i}.traces_res=squeeze(SeeResiduals(reshape(Result_opto{i}.traces_res,ns,1,[]),(Result_opto{i}.mcTrace(:,1).*Result_opto{i}.mcTrace(:,2))'));
    %Result{i}.traces_res=squeeze(SeeResiduals(reshape(Result{i}.traces_res,ns,1,[]),PCA_trace'));
    Result_opto{i}.traces_res=moving_residual(Result_opto{i}.traces_res,Result_opto{i}.mcTrace,PCA_trace,300);
    Result_opto{i}.traces_res_hi=Result_opto{i}.traces_res - movmedian(Result_opto{i}.traces_res,150,2);

    Result_opto{i}.traces_hi=squeeze(Result_opto{i}.traces) - movmedian(squeeze(Result_opto{i}.traces),150,2);
    Result_opto{i}.traces_bin_hi=squeeze(Result_opto{i}.traces_bin) - movmedian(squeeze(Result_opto{i}.traces_bin),200,2);

    Result_opto{i}.noise=get_threshold(Result_opto{i}.traces,1); Result_opto{i}.noise_res=get_threshold(Result_opto{i}.traces_res,1);
    Result_opto{i}.traces_hi=Result_opto{i}.traces_hi./Result_opto{i}.noise;
    Result_opto{i}.traces_bin_hi=Result_opto{i}.traces_bin_hi./get_threshold(Result_opto{i}.traces_bin_hi,1);
    Result_opto{i}.traces_res_hi=Result_opto{i}.traces_res_hi./Result_opto{i}.noise_res;

    Result_opto{i}.spike=[]; WaveT=[];
    %tmp=squeeze(Result_opto{i}.traces_res) - movmedian(squeeze(Result_opto{i}.traces_res),4,2); tmp=tmp./get_threshold(tmp,1);
    %Result_opto{i}.spike=find_spike_bh(tmp,4,3);
    tmp=squeeze(Result_opto{i}.traces_res) - movmedian(squeeze(Result_opto{i}.traces_res),100,2); tmp=tmp./get_threshold(tmp,1);
    %Result_opto{i}.spike=(Result_opto{i}.spike+find_spike_bh(tmp,3,1.5))>0;
    Result_opto{i}.spike=find_spike_bh(tmp,3.2,1.5);
 
end

save(['Treadmill_BHLm008_20230223_result.mat'],'Result','Result_opto','fpath','fpath_const','-v7.3')
%%
for i=%1:length(fpath)
    load(fullfile(fpath{i},'settings.mat'));
    load(fullfile(fpath{i},'mcTrace01.mat'));
    Sz = importdata([fpath{i} '/experimental_parameters.txt']);
    sz1=Sz.data(1); sz2=Sz.data(2);

    [a,b]=system(sprintf('GetFileInfo "%s"',fullfile(fpath{i},'settings.mat'))); s=strfind(b,'modified:')+10; crdat=b(s:s+18);
    image_start=datestr(datenum(crdat));
    dt=datetime(image_start(end-8:end),'InputFormat','HH:mm:ss'); dt.Format = 'HH:mm:ss';
    fnm=dir(fullfile(fpath{i}, '*.csv'));

    a=find(DAQ_waves.amplitude(1,[1:500])); frm_rate=1/((a(2)-a(1)-2)*1e-5);
    frm_end=sum(DAQ_waves.amplitude(1,:)); f_seg=[[1:10000:frm_end] frm_end+1];
    [Arduino_data reward_pos lap_dist]=match_treadmill_DAQ(fullfile(fpath{i},fnm(2).name),1e-5,DAQ_data,(1/frm_rate+0.00002),2); % time, treadmill, Reward, Run
    Arduino_data(:,1)=Arduino_data(:,1)*(1/frm_rate*1000)/(1/frm_rate*1000+0.02);
    Result{i}.Arduino=Arduino_data;

    mov_mc=double(readBinMov([fpath{i} '/mc' num2str(11,'%02d') '.bin'],sz2,sz1));
    mov_res= mov_mc(:,:,1:10000)-mean(mov_mc,3);

    mov_res = SeeResiduals(mov_res,mcTrace(1:10000,:));
    mov_res = SeeResiduals(mov_res,mcTrace(1:10000,:).^2);
    mov_res = SeeResiduals(mov_res,mcTrace(1:10000,1).*mcTrace(1:10000,2));

    Result{i}.c_ftprnt=mask_footprint(Result{i}.centers,movmean(mov_res(:,:,1000:end),10,3),[],6);
    Result{i}.ref_im=mean(mov_mc,3);
    Result{i}.traces=[]; Result{i}.traces_hi=[]; Result{i}.traces_bin=[]; Result{i}.traces_bin_hi=[];
    Result{i}.traces_res=[]; Result{i}.traces_res_hi=[]; Result{i}.mcTrace=[];
    Result{i}.im_corr=[];

    Result{i}.frm_rate=frm_rate;
    Result{i}.DMDPatt=dmd_mask_sequence_rois;
    Result{i}.coord=get_coord(Result{i}.c_ftprnt);
    Result{i}.traces=[Result{i}.traces -(tovec(mov_mc(:,:,1:10000))'*tovec(Result{i}.c_ftprnt))'];
    Result{i}.traces_bin=[Result{i}.traces_bin -(tovec(mov_mc(:,:,1:10000))'*tovec(Result{i}.c_ftprnt>0))'];
    mov_mc_vec=tovec(mov_mc(:,:,1:10000)); mov_mc_vec=(mov_mc_vec-mean(mov_mc_vec,1))./std(mov_mc_vec,0,1);
    ref_im_vec=tovec(Result{i}.ref_im); ref_im_vec=(ref_im_vec-mean(ref_im_vec))./std(ref_im_vec);
    Result{i}.im_corr=sum(mov_mc_vec.*ref_im_vec,1)/(size(mov_mc_vec,1)-1);

    for j=2:length(f_seg)-1
        mov_mc=double(readBinMov([fpath{i} '/mc' num2str(j,'%02d') '.bin'],sz2,sz1));
        try mov_mc=mov_mc(:,:,1:10000); catch mov_mc=mov_mc; end

        load([fpath{i} '/mcTrace' num2str(j,'%02d') '.mat']);
        mov_mc_vec=tovec(mov_mc); mov_mc_vec=(mov_mc_vec-mean(mov_mc_vec,1))./std(mov_mc_vec,0,1);
        Result{i}.im_corr=[Result{i}.im_corr sum(mov_mc_vec.*ref_im_vec,1)/(size(mov_mc_vec,1)-1)];

        Result{i}.traces=[Result{i}.traces -(tovec(mov_mc)'*tovec(Result{i}.c_ftprnt))'];
        Result{i}.traces_bin=[Result{i}.traces_bin -(tovec(mov_mc)'*tovec(Result{i}.c_ftprnt>0))'];
        j
    end
    for j=1:length(f_seg)-1
        load([fpath{i} '/mcTrace' num2str(j,'%02d') '.mat']);
        try mcTrace=mcTrace(1:10000,:); catch mcTrace=mcTrace; end
        Result{i}.mcTrace=[Result{i}.mcTrace; mcTrace];
    end

    [V, D, W] = eig(corrcoef(zscore(Result{i}.traces,0,2)'));
    D = diag(D); D = D(end:-1:1); V = V(:,end:-1:1);
    PCA_trace=zscore(V(:,1)'*Result{i}.traces,0,2);
    ns=size(Result{i}.traces,1);

    Result{i}.traces_res=squeeze(SeeResiduals(reshape(Result{i}.traces,ns,1,[]),Result{i}.mcTrace'));
    Result{i}.traces_res=squeeze(SeeResiduals(reshape(Result{i}.traces_res,ns,1,[]),Result{i}.mcTrace.^2'));
    Result{i}.traces_res=squeeze(SeeResiduals(reshape(Result{i}.traces_res,ns,1,[]),(Result{i}.mcTrace(:,1).*Result{i}.mcTrace(:,2))'));
    %Result{i}.traces_res=squeeze(SeeResiduals(reshape(Result{i}.traces_res,ns,1,[]),PCA_trace'));
    Result{i}.traces_res=moving_residual(Result{i}.traces_res,Result{i}.mcTrace,PCA_trace,300);
    Result{i}.traces_res_hi=Result{i}.traces_res - movmedian(Result{i}.traces_res,150,2);

    Result{i}.traces_hi=squeeze(Result{i}.traces) - movmedian(squeeze(Result{i}.traces),150,2);
    Result{i}.traces_bin_hi=squeeze(Result{i}.traces_bin) - movmedian(squeeze(Result{i}.traces_bin),200,2);

    Result{i}.noise=get_threshold(Result{i}.traces,1); Result{i}.noise_res=get_threshold(Result{i}.traces_res,1);
    Result{i}.traces_hi=Result{i}.traces_hi./Result{i}.noise;
    Result{i}.traces_bin_hi=Result{i}.traces_bin_hi./get_threshold(Result{i}.traces_bin_hi,1);
    Result{i}.traces_res_hi=Result{i}.traces_res_hi./Result{i}.noise_res;

    Result{i}.spike=[]; WaveT=[];
    fpass=[15 55]; [b, a] = butter(4, fpass/(Result{i}.frm_rate/2), 'stop');
    for n=1:size(traces,1)
    Result{i}.traces_res_hi_filtered(n,:) = filtfilt(b, a, double(Result{1}.traces_res_hi(n,:)));
    end

    tmp=squeeze(Result{i}.traces_res_hi) - movmedian(squeeze(Result{i}.traces_res_hi),4,2); tmp=tmp./get_threshold(tmp,1);
    Result{i}.spike=find_spike_bh(tmp,4,3);
    tmp=squeeze(Result{i}.traces_res_hi_filtered) - movmedian(squeeze(Result{i}.traces_res_hi_filtered),30,2); tmp=tmp./get_threshold(tmp,1);
    Result{i}.spike=(Result{i}.spike+find_spike_bh(tmp,3.2,1.5))>0;
    
    traces=Result{i}.traces_res-movmedian(Result{i}.traces_res,200,2);
    traces=traces./get_moving_threshold(traces,1,2000);
   
Result{i}.CS_list=[]; Result{i}.CS_spike=[];
for n=1:size(traces,1)
[l Result{i}.CS_spike(n,:)]=find_CS(traces(n,:),Result{i}.spike(n,:),12,0.05);
Result{i}.CS_list=[Result{i}.CS_list; [repmat(n,size(l,1),1) l]];
end

end
%save(['Treadmill_BHLm008_20230222_result.mat'],'Result','Result_opto','fpath','fpath_const','-v7.3')
%%