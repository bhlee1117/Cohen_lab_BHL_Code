clear
clc;

[fpath] = uigetfile_n_dir(); %only Treadmill data
%% Parameter setting and get cell coordinate
block_size=20;

for i=6:length(fpath)

load(fullfile(fpath{i},"output_data.mat"))
sz=double(Device_Data{1, 3}.ROI([2 4]));
mov_test=double(readBinMov_times([fpath{i} '/frames1.bin'],sz(2),sz(1),[3000:4000]));
    avgImg=mean(mov_test,3);
    [CellCoordinate{i}, ~]=Cell_segment_circle_25x_VU(avgImg,0.85);
    CellCoordinate{i}=cell_detection_manual(avgImg,CellCoordinate{i},[0 10000]); 
figure(1); clf;
    imshow2(avgImg,[]); hold all
for n=1:size(CellCoordinate{i},1)
    ROI{i}(n,:)=round([CellCoordinate{i}(n,1)-block_size CellCoordinate{i}(n,1)+block_size ...
            CellCoordinate{i}(n,2)-block_size CellCoordinate{i}(n,2)+block_size]);
        bs{i}(n)=block_size;

        % if ROI go further than outbound, reduce window size
        while sum([ROI{i}(n,1) > 1, ROI{i}(n,2)<sz(1), ROI{i}(n,3)> 1, ROI{i}(n,4)<sz(2)])<4
            bs{i}(n)=bs{i}(n)-1;
            ROI{i}(n,:)=round([CellCoordinate{i}(n,1)-bs{i}(n) CellCoordinate{i}(n,1)+bs{i}(n) ...
                CellCoordinate{i}(n,2)-bs{i}(n) CellCoordinate{i}(n,2)+bs{i}(n)]);
        end
        ROI_box=[ROI{i}(n,[1 3]); ROI{i}(n,[1 4]); ROI{i}(n,[2 4]); ROI{i}(n,[2 3]); ROI{i}(n,[1 3])];
    plot(ROI_box(:,1),ROI_box(:,2),'r')
end
    
    plot(CellCoordinate{i}(:,1),CellCoordinate{i}(:,2),'ro','markersize',15)
    saveas(gca,[fpath{i} '/ROIs.png'])
end
%%
time_size=150000; %segment size
DAQ_rate=0.000005;

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

    for n=1:size(CellCoordinate{i},1)
        
        % load ROI reference
        mov_seg=double(readBinMov_times_ROI([fpath{i} '/frames1.bin'],sz(1),sz(2),ref_time,ROI{i}(n,:)));

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
            mov=vm(readBinMov_times_ROI([fpath{i} '/frames1.bin'],sz(1),sz(2),[t_seg(t,1):t_seg(t,2)],ROI{i}(n,:)));
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

