clear
clc;

[fpath] = uigetfile_n_dir(); %only Treadmill data
save_name='pcResult_20230830.mat';
%% Parameter setting and get cell coordinate
block_size=20;
DAQ_rate=0.000005;

for i=1:length(fpath)

load(fullfile(fpath{i},"output_data.mat"))
sz=double(Device_Data{1, 3}.ROI([2 4]));
mov_test=double(readBinMov_times([fpath{i} '/frames1.bin'],sz(2),sz(1),[3000:4000]));
    avgImg=mean(mov_test,3);
    [CellCoordinate, ~]=Cell_segment_circle_25x_VU(avgImg,0.85);
    CellCoordinate=cell_detection_manual(avgImg,CellCoordinate,[0 9000]); 
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
end
%%
time_size=150000; %segment size

for i=1:length(fpath)
    clear mcTrace_block ave_im
    %load the center positions
    load(fullfile(fpath{i},"output_data.mat"))
    sz=double(Device_Data{1, 3}.ROI([2 4]));
    Blue=Device_Data{1, 2}.buffered_tasks(1, 1).channels(1, 2).data;
    CamCounter=Device_Data{1, 2}.Counter_Inputs(1, 1).data;
    CamTrigger=find(CamCounter(2:end)-CamCounter(1:end-1));
    Blue=Blue(CamTrigger);
    ref_time=[10000:11000];

    % make time segment
    t_seg=[[1:time_size:length(CamTrigger)] length(CamTrigger)+1];
    t_seg=[t_seg(1:end-1)' t_seg(2:end)'+100];
    t_seg(end)=length(CamTrigger);

    for n=1:size(CellCoordinate,1)
        
        % load ROI reference
        mov_seg=double(readBinMov_times_ROI([fpath{i} '/frames1.bin'],sz(1),sz(2),ref_time,ROI(n,:)));

        % motion correction of reference
        options_rigid = NoRMCorreSetParms('d1',size(mov_seg,1),'d2',size(mov_seg,2),'bin_width',200,'max_shift',block_size,'us_fac',50,'init_batch',200);
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
            [mov_mc,xyField]=optical_flow_motion_correction_LBH(mov,double(mov_ref),'optic_flow');

            ave_im{n,t}=mean(mov_mc,3);
            mov_mc=vm(mov_mc);
            mov_mc.transpose.savebin([fpath{i} '/mc' num2str(n,'%02d') '_' num2str(t,'%02d') '.bin'])

            mcTrace_block{n,t} = squeeze(mean (xyField,[1 2]))'; %optic flow
            %mcTrace_block{n,t}=xyField; % Normcorre
        end


    end
    save([fpath{i} '/mcTrace.mat'],'CellCoordinate','mcTrace_block','ave_im','bs','-v7.3')
    disp(['MC done:' fpath{i}])

end


%% Signal Extraction

for i=1:length(fpath)
    load(fullfile(fpath{i},'mcTrace.mat'))
    Result{i}.centers=CellCoordinate;
    Result{i}.fpath=fpath{i};
    disp(['loading  ' fpath{i}])
    load([fpath{i} '/output_data.mat'])
    sz=double(Device_Data{1, 3}.ROI([2 4]));
mov_test=double(readBinMov_times([fpath{i} '/frames1.bin'],sz(2),sz(1),[3000:4000]));
Result{i}.FOV=mean(mov_test,3);

    try
        Result{i}.Blue=Device_Data{1, 2}.buffered_tasks(1, 1).channels(1, 2).data;
    end
    try
        Result{i}.Reward=Device_Data{1, 2}.buffered_tasks(1, 3).channels.data;
    end

    CamCounter=Device_Data{1, 2}.Counter_Inputs(1, 1).data;
    CamTrigger=find(CamCounter(2:end)-CamCounter(1:end-1));
    Result{i}.frm_rate=double((CamTrigger(2)-CamTrigger(1))*DAQ_rate);
    sz=double(Device_Data{1, 3}.ROI([2 4]));
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
        mov_res= SeeResiduals(mov_res,bkg,1);


        Result{i}.ref_im(1:2*bs(n)+1,1:2*bs(n)+1,n)=mean(mov,3);
        ref_im_vec=tovec(Result{i}.ref_im(1:2*bs(n)+1,1:2*bs(n)+1,n));
        ref_im_vec=(ref_im_vec-mean(ref_im_vec,1))./std(ref_im_vec,0,1);
        mov_mc_vec=tovec(mov);
        mov_mc_vec=(mov_mc_vec-mean(mov_mc_vec,1))./std(mov_mc_vec,0,1);
        Result{i}.im_corr{n}=[];
        Result{i}.im_corr{n}=[Result{i}.im_corr{n} sum(mov_mc_vec.*ref_im_vec,1)/(size(mov_mc_vec,1)-1)];


        Result{i}.c_ftprnt(1:2*bs(n)+1,1:2*bs(n)+1,n)=mask_footprint([bs(n)+0.5 bs(n)+0.5],movmean(mov_res(:,:,1000:end-1000),10,3),[],8);
        Result{i}.c_ftprnt(1:2*bs(n)+1,1:2*bs(n)+1,n)=imgaussfilt(Result{i}.c_ftprnt(1:2*bs(n)+1,1:2*bs(n)+1,n),0.6);

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
        Result{i}.AvgImg=mean(mov_mc,3);
    end

end

for i=1:length(fpath)
    Result{i}.spike=zeros(size(Result{i}.traces));
    tmp=squeeze(Result{i}.traces_res) - movmedian(squeeze(Result{i}.traces_res),250,2); tmp=tmp./get_threshold(tmp,1);
    Result{i}.spike=(Result{i}.spike+find_spike_bh(tmp,5.5,2))>0;
end

save(save_name,'Result','fpath','-v7.3')

%%

mperVR=6.44*1e-3;

Lap_FR=[]; Lap_V=[];  
Reward_lap=[];
rewardPos=0.8; LapDist=115;
DAQ_rate=0.000005;


for i=1:length(VirFile)

        fileList = dir(fullfile(fpath{i}, '*.data'));
        if length(fileList)==1
            fid = fopen(fullfile(fpath{i},fileList.name));
            VRdata = fread(fid,[12 inf],'double');
        else
            error('There is more than one .data file');
        end

        load([fpath{i} '/output_data.mat'])
    CamCounter=Device_Data{1, 2}.Counter_Inputs(1, 1).data;
    CamTrigger=find(CamCounter(2:end)-CamCounter(1:end-1));
    


    VRdata=Dat{i}(:,find(Dat{i}(10,:)));
    %VRdata(7,:)=zscore(abs(VRdata(7,:)-movmedian(VRdata(7,:),20)));
    %HF=abs(VRdata(7,:)-movmedian(VRdata(7,:),20));
    HF=abs(VRdata(9,:)-movmedian(VRdata(9,:),500));
    VRdata(9,:)=(HF-mean(HF))/std(HF);
    [pks lick]=findpeaks(VRdata(9,:));
   
    Lick_trace=zeros(1,size(VRdata,2));
    Lick_trace(lick(pks>2))=1;
   
    t_VR = datetime(datetime(VRdata(1,:),'ConvertFrom','datenum'), 'InputFormat', 'yyyy-MM-dd HH:mm:ss.SSS');
    t_VR= t_VR-t_VR(1);
    t_VR= milliseconds(t_VR)/1000;
    t_DAQ=CamTrigger*DAQ_rate;*t_VR(end);

    %[Lap_FR Lap_V]=get_LapFR_VU(double(VRdata(7,:)>0.1),100,VRdata([1 1 2 3 4 5 6 8],:),0.005,1);
    mov_trace_hi=-(mov_trace-movmedian(mov_trace,300,2));
    mov_trace_hi=mov_trace_hi./get_threshold(mov_trace_hi,1);
    sp=find_spike_bh(mov_trace_hi,4,3);
    %mov_trace_hi=-movmean((mov_trace-movmedian(mov_trace,500,2)),20,2);
    


    clear FR
    for n=1:size(mov_trace_hi)
    %[FR(:,:,n) V]=get_LapFR_VU(double(mov_trace_hi(n,:)),t_DAQ,30,VRdata,0.005,1,115);
    [FR(:,:,n) V]=get_LapFR_VU(double(movsum(sp(n,:),1000)),t_DAQ,50,VRdata,0.001,1,115);
    end
    [LickFR V]=get_LapFR_VU(double(Lick_trace),t_VR,50,VRdata,0,1,115);
Lap_FR=[Lap_FR; FR];
Lap_V=[Lap_V; V*mperVR;];

end
%%
for noi=1:size(ROI,1)
figure(noi); clf;
tiledlayout(2,6)
ax1=[];
ax1=[ax1 nexttile([1 1])];
imagesc(Lap_FR(:,:,noi)); hold all;
title('Firing rate')
ax1=[ax1 nexttile([1 1])];
plot(mean(Lap_FR(:,:,noi),1,'omitnan'))

ax1=[ax1 nexttile([1 1])];
imagesc(LickFR); hold all;
title('Lick rate')
ax1=[ax1 nexttile([1 1])];
plot(mean(LickFR,1,'omitnan'))

ax1=[ax1 nexttile([1 1])];
imagesc(Lap_V,[0 0.01]); hold all;
colormap('turbo')
title('Speed')
colorbar

ax1=[ax1 nexttile([1 1])];
plot(mean(Lap_V,1,'omitnan'))

linkaxes(ax1,'x')
axis tight

nexttile([1 6])
t_DAQ_scaled=[t_VR(1): (t_VR(end)-t_VR(1))/(size(mov_trace,2)-1):t_VR(end)];
plot(t_DAQ_scaled,mov_trace_hi(noi,:),'k',t_DAQ_scaled(find(sp(noi,:))),mov_trace_hi(noi,find(sp(noi,:))),'r.')
hold all
plot(t_VR,VRdata(12,:)*15)
plot(t_VR,rescale(VRdata(5,:))*15)

filename=char(VirFile);
saveas(gca,[filename(1:end-5) '_' num2str(noi) '.fig'])

end


