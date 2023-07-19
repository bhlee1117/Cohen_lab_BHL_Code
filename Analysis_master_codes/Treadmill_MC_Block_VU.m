clear
clc;

[fpath] = uigetfile_n_dir(); %only Treadmill data
%%
[fpath_optopatch] = uigetfile_n_dir(); %load optopatch data

%%
fop=[1 2 2 2];
spike_threshold=4.2;

for i=1:length(fpath_optopatch)
    load(fullfile(fpath_optopatch{i},"output_data.mat"))
    load(fullfile(fpath_optopatch{i},"mcTrace01.mat"))
    sz=double(Device_Data{1, 3}.ROI([2 4]));
    Blue=Device_Data{1, 2}.buffered_tasks(1, 1).channels(1, 2).data;
    CamCounter=Device_Data{1, 2}.Counter_Inputs(1, 1).data;
    CamTrigger=find(CamCounter(2:end)-CamCounter(1:end-1));
    Blue=Blue(CamTrigger);
    mov_opto=double(readBinMov([fpath_optopatch{i} '/mc01.bin'],sz(2),sz(1)));
    AvgImg{i}=mean(mov_opto,3);

    [c, radii]=Cell_segment_circle_10x_VU(mean(mov_opto,3),0.85);
    c=cell_detection_manual(mean(mov_opto,3),c,[0 6000]);
    mov_opto=SeeResiduals(mov_opto,mcTrace);
    c_ftprnt=mask_footprint(c,movmean(mov_opto(:,:,1000:end),10,3),[],7);
    traces=-(tovec(mov_opto)'*tovec(c_ftprnt))';
    traces_hi=traces(:,1:length(CamTrigger))-movmedian(traces(:,1:length(CamTrigger)),150,2);
    traces_hi=traces_hi./get_threshold(traces_hi,1);
    [sp]=find_spike_bh(traces_hi,3,2);
    sp(:,Blue==0)=0;
    spike_height=traces(:,1:length(CamTrigger))-movmedian(traces(:,1:length(CamTrigger)),15,2);
    spike_height=spike_height./get_threshold(spike_height,1);
    spike_height(~sp)=NaN;
    spike_height=mean(spike_height,2,'omitnan');
    goodCell=find(spike_height>spike_threshold & sum(sp,2,'omitnan')>15);
    Centers_optopatch{i}=get_coord(c_ftprnt(:,:,goodCell)>0);

    CellCoordinate=Centers_optopatch{i};
    ref_im=AvgImg{i};
    save(fullfile(fpath_optopatch{i},'Centers.mat'),'CellCoordinate','ref_im')

    show_footprnt(c_ftprnt(:,:,goodCell),AvgImg{i})
    title(fpath_optopatch{i},'Interpreter','none')
    show_traces_spikes(traces_hi(goodCell,:),sp(goodCell,:),Blue)
    title(fpath_optopatch{i},'Interpreter','none')
end

%%
block_size=15;
time_size=150000;
DAQ_rate=0.000005;

for i=1:length(fpath)
    clear mcTrace_block ave_im
    %load the center positions
    load(fullfile(fpath_optopatch{fop(i)},'Centers.mat'))
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

    % load the segment from VR movie
    mov_test=double(readBinMov_times([fpath{i} '/frames1.bin'],sz(2),sz(1),[3000:4000]));
    % calculate xcorrelation btw optopatch and VR
    xcorrRefIm=normxcorr2(mean(mov_test,3),ref_im);
    [~, ind]=max(xcorrRefIm(:));
    [shifty shiftx]=ind2sub(size(xcorrRefIm,1),ind);
    shifty=sz(2)-shifty; shiftx=sz(1)-shiftx;
    % shift
    CellCoordinate=CellCoordinate+[shiftx shifty];
    % remove cells in outside of boundary
    rmv=find([CellCoordinate(:,1) < 1 | CellCoordinate(:,1) > sz(1) | CellCoordinate(:,2) < 1 | CellCoordinate(:,2) > sz(2)]);
    CellCoordinate(rmv,:)=[];

    for n=1:size(CellCoordinate,1)
        % set ROI
        ROI=round([CellCoordinate(n,1)-block_size CellCoordinate(n,1)+block_size ...
            CellCoordinate(n,2)-block_size CellCoordinate(n,2)+block_size]);
        bs(n)=block_size;

        % if ROI go further than outbound, reduce window size
        while sum([ROI(1) > 1, ROI(2)<sz(1), ROI(3)> 1, ROI(4)<sz(2)])<4
            bs(n)=bs(n)-1;
            ROI=round([CellCoordinate(n,1)-bs(n) CellCoordinate(n,1)+bs(n) ...
                CellCoordinate(n,2)-bs(n) CellCoordinate(n,2)+bs(n)]);
        end
        
        % load ROI reference
        mov_seg=readBinMov_times_ROI([fpath{i} '/frames1.bin'],sz(1),sz(2),ref_time,ROI);

        % motion correction of reference
        options_rigid = NoRMCorreSetParms('d1',size(mov_seg,1),'d2',size(mov_seg,2),'bin_width',200,'max_shift',block_size,'us_fac',50,'init_batch',200);
        tic; [mov_seg,shifts1,template1,options_rigid] = normcorre(mov_seg,options_rigid); toc
        
        % make reference image
        mov_seg=vm(mov_seg);
        mov_seg = single(mov_seg)./single(max(mov_seg.data(:)));
        mov_seg = movmean(mov_seg,10,3);
        mov_ref = squeeze(median(mov_seg,3));

        % load time segment of ROI
        for t=1:size(t_seg,1)
            mcTrace_block{n,t}=[];
            mov=vm(readBinMov_times_ROI([fpath{i} '/frames1.bin'],sz(1),sz(2),[t_seg(t,1):t_seg(t,2)],ROI));
            [mov_mc,xyField]=optical_flow_motion_correction_LBH(mov,double(mov_ref),'normcorre');

            ave_im{n,t}=mean(mov_mc,3);
            mov_mc=vm(mov_mc);
            mov_mc.transpose.savebin([fpath{i} '/mc' num2str(n,'%02d') '_' num2str(t,'%02d') '.bin'])

            %mcTrace = squeeze(mean (xyField,[1 2])); %optic flow
            mcTrace_block{n,t}=xyField; % Normcorre
        end


    end
    save([fpath{i} '/mcTrace.mat'],'CellCoordinate','mcTrace_block','ave_im','bs','-v7.3')
    disp(['MC done:' fpath{i}])



end

