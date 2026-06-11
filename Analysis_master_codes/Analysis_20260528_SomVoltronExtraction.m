%% Set the path
clear
clc;
cd '/Users/bhlee1117/Documents/GitHub/Cohen_lab_BHL_Code/Analysis_master_codes';
[~, ~, raw] = xlsread(['/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/' ...
    'NaVInactivationData_Arrangement.xlsx'], 'Sheet1', 'C5:T462');

save_to='/Volumes/BHL18TB_D2/Arranged_Data/Prism_OptopatchResult';
fpath=raw(:,1);
set(0,'DefaultFigureWindowStyle','docked')
time_segment=25000;
%% MC

for i=[439:447]
    i
    load(fullfile(fpath{i},"output_data.mat"))
    ref_time=[6000:7000];

    sz=double(Device_Data{1, 3}.ROI([2 4]));
    info = dir([fpath{i} '/frames1.bin']);
    bytes_per_pixel = 2;  % uint16 = 2 bytes
    frm_end= info.bytes / (sz(1) * sz(2) * bytes_per_pixel);

    %frm_end =max(Device_Data{1, 2}.Counter_Inputs(1, 1).data);
    f_seg=[[1:time_segment:frm_end] frm_end];
    CamDAQ_rate=Device_Data{1, 2}.Counter_Inputs.rate;
    CamTrig=Device_Data{1, 2}.Counter_Inputs.data;
    CamTrig2=find(CamTrig(2:end)-CamTrig(1:end-1)>0);
    Frm_rate=(CamTrig2(2)-CamTrig2(1))/CamDAQ_rate;

    disp(['N frames: ' num2str(frm_end)]);

    mov_test=double(readBinMov_times([fpath{i} '/frames1.bin'],sz(2),sz(1),[ref_time(1):ref_time(end)]));
    mov_test=mov_test(:,:,2:end);
    [mov_test_mc,xyField]=optical_flow_motion_correction_LBH(mov_test,mean(mov_test,3),'normcorre');
    mov_test=vm(mov_test);
    mov_test = single(mov_test)./single(max(mov_test.data(:)));
    mov_test = movmean(mov_test,10,3);
    mov_ref = squeeze(median(mov_test,3));

    for j=1:length(f_seg)-1

        mov=double(readBinMov_times([fpath{i} '/frames1.bin'],sz(2),sz(1),[f_seg(j):f_seg(j+1)]));
        mov=vm(mov);
        [mov_mc,xyField]=optical_flow_motion_correction_LBH(mov,mov_ref,'normcorre');

        ave_im=mean(mov_mc,3);
        mov_mc=vm(mov_mc);
        mov_mc.transpose.savebin([fpath{i} '/mc_ShutterReg' num2str(j,'%02d') '.bin'])

        % mcTrace = squeeze(mean(xyField,[1 2])); %optic flow
        mcTrace=xyField; % Normcorre
        save([fpath{i} '/mcTrace' num2str(j,'%02d') '.mat'],'mcTrace','ave_im')

        %  clear mov_mc mov
    end
    sd_mov=std(double(mov),0,3); sd_mov_mc=std(double(mov_mc),0,3);
    figure(i); clf;
    nexttile([1 1]); imshow2(sd_mov,[]); title('before mc')
    nexttile([1 1]); imshow2(sd_mov_mc,[]); title('after mc')
    nexttile([1 1]); imshow2(imfuse(sd_mov,sd_mov_mc),[]);
    title(fpath{i},'Interpreter',  'none')
    saveas(gca,[char(fpath{i}) '/' 'MC_result.fig'])
end

%% Extract voltage signal
bound = 6;
for f=[439:447]
    fprintf('Processing %2.0f ...',f);
    Result=[];
    load(fullfile(fpath{f},"output_data.mat"))
    sz=double(Device_Data{1, 3}.ROI([2 4]));
    CamCounter=Device_Data{1, 2}.Counter_Inputs(1, 1).data;
    CamTrigger=find(CamCounter(2:end)-CamCounter(1:end-1));

    info = dir([fpath{f} '/frames1.bin']);
    bytes_per_pixel = 2;  % uint16 = 2 bytes
    frm_end= info.bytes / (sz(1) * sz(2) * bytes_per_pixel);
    CamTrigger=interp1([1:length(CamTrigger)],CamTrigger,[1:frm_end],'linear','extrap');

    %frm_end=max(Device_Data{1, 2}.Counter_Inputs(1, 1).data);
    f_seg=[[1:time_segment:frm_end] frm_end+1];
    readFrame=diff(f_seg);
    Result.Blue=Device_Data{1, 2}.buffered_tasks(1, 1).channels(1, 2).data;
    Result.Blue=Result.Blue(CamTrigger);
    [~, Result.bluePatt]=get_blueDMDPatt(Device_Data,'normal');
    totalNbin=ceil(frm_end/time_segment);

    Result.traces=[];
    Result.mc=[];
    MINIALIresult=[];

    for j=1:totalNbin
        mov_path=[fpath{f} '/mc_ShutterReg' num2str(j,'%02d') '.bin'];
        mc_path=[fpath{f} '/mcTrace' num2str(j,'%02d') '.mat'];
        load(mc_path);
        motionTrace=movmean(mcTrace.xymean([1:readFrame(j)],:),5,1);
        Result.mc=[Result.mc; motionTrace];

        try
            mov_mc=double(readBinMov_times(mov_path,sz(2),sz(1),[1:readFrame(j)]));
            disp('readBinMov_times failed, move on to readBinMov')
        catch
            mov_mc=double(readBinMov(mov_path,sz(2),sz(1)));
            disp('Movie loaded!')
        end

        if ~isfield(Result,'ref_im')
            Result.ref_im=mean(mov_mc,3);
        end

        bkg = zeros(1, size(mov_mc,3));
        mean_F=squeeze(mean(mov_mc,[1 2]));
        [~, blueOff]=get_blueoffTrace(mean_F,[Result.Blue(f_seg(j):f_seg(j+1)-1)],70);
        [y_fit]=expfitDM_2(find(blueOff)',mean_F(find(blueOff)),[1:size(mov_mc,3)]',1000);
        bkg(1,:)=y_fit;

        mov_mc2=mov_mc./reshape(bkg,1,1,[]);
        mov_res= mov_mc2-mean(mov_mc2,3);
        mov_res = SeeResiduals(mov_res,motionTrace);
        mov_res = SeeResiduals(mov_res,motionTrace.^2);
        mov_res = SeeResiduals(mov_res,motionTrace(:,1).*motionTrace(:,end));

        opts=[];
        opts.zSeedMin=3.8;
        opts.r_max = 8.4;
        opts.doPlot = false;
        opts.uiSeedThresh  = false;
        opts.uiClusterSize = false;
        opts.uiSynFilter   = false;
        opts.synFilter = struct( ...
            'ampRange',   [-inf inf],  ...
            'countRange', [3  inf],  ...
            'areaRange',  [-inf inf]);

        if isempty(MINIALIresult)
            MINIALIresult=extractVoltronST(-mov_res,mov_mc,opts);
            Result.ftprnt=MINIALIresult.S_glu;
            Result.traces=MINIALIresult.T_glu;
            Result.extractionResult=MINIALIresult;
        else
            opts= struct;
            opts.mode = 'project';
            MINIALIresulttmp= extractVoltronST(-mov_res, mov_mc, 0.001, opts, MINIALIresult.S_glu);
            Result.traces=[Result.traces MINIALIresulttmp.T_glu];
        end
    end

    figure(100); clf; tiledlayout(2,2)
    ax1=nexttile([1 1]);
    show_footprnt_contour(Result.ftprnt,Result.ref_im);
    colormap(ax1,gray(256));
    nexttile([1 1])
    imshow2(im_merge(Result.ftprnt,lines(size(Result.ftprnt,3)))*10,[])
    nexttile(3,[1 2])
    plot(rescale2(Result.traces,2)'+[1:size(Result.traces,1)]);
    drawnow;
    sgtitle(sprintf(fpath{f}), 'Interpreter', 'none')
    saveas(gcf,fullfile(fpath{f},'extraction_result.fig'));
    saveas(gcf,fullfile(fpath{f},'extraction_result.png'),'png');
    save([fpath{f} '/OP_Result.mat'],'Result','-v7.3')
    disp([num2str(f) ' th file is saved.'])
end
%% find spikes

for f=[447]
    fprintf(['Processing ' fpath{f} '...\n']);
    VoltResult=importdata([fpath{f} '/OP_Result.mat'],'-v7.3');
    [nROI nTime]=size(VoltResult.traces);
    zscoredTr=(VoltResult.extractionResult.dFF_glu);
    zscoredTr=(zscoredTr-median(zscoredTr(:,VoltResult.Blue==0),2))./mad(zscoredTr(:,VoltResult.Blue==0), 0, 2);
    zscoredTr_hi=zscoredTr-movmedian(zscoredTr,300,2);
    isvalideCell=zeros(nROI,1);
    spikeMat=zeros(nROI,nTime);
    CStrace=zeros(nROI,nTime);
    SpClassMat=zeros(nROI,nTime);
    for n=1:nROI 
        figure(10); clf; tiledlayout(1,4);
        cmap_ftprnt=repmat(0.7,nROI,3);
        cmap_ftprnt(n,:)=[1 0 0];
        nexttile([1 1]);
        imshow2(im_merge(VoltResult.ftprnt,cmap_ftprnt)*10,[])
        nexttile([1 3])
        plot(zscoredTr(n,:));
        kk=input(sprintf('Cell# %2.0f / %2.0f, Is valid?\n',n, nROI));
        isvalideCell(n)=(kk==1);
        if isvalideCell(n)
            spikeMat(n,:)=find_spike_bh(zscoredTr_hi(n,:),3,1);
            subtr=movmean(zscoredTr_hi(n,:),20,2);
            subtr=subtr-movmedian(subtr,300,2);
            [trans, tr_trace] = detect_transient2(subtr, [3 0.5], spikeMat(n,:), 20);
            if isempty(trans.amp)
                CS_trace_tmp = zeros(1, nTime);
            else
                transcand = cell2mat(cellfun(@(x) length(x) > 1, trans.ISI, 'UniformOutput', false));
                meanISI_frnt = cellfun(@(x) mean(x(1:min(3-1,end))), trans.ISI(transcand));
                meanISI_first3 = zeros(1,length(trans.length));
                meanISI_first3(transcand) = meanISI_frnt;

                CS_ind = find(trans.spike_number >= 3 & meanISI_first3 <= 20 & trans.int > 30);
                CS_trace_tmp = ismember(tr_trace, CS_ind);
            end
            CStrace(n,:)=CS_trace_tmp;
            figure(10); clf;
            plot(zscoredTr(n,:),'k'); hold all
            scatter(find(spikeMat(n,:)),zscoredTr(n,find(spikeMat(n,:))),'ro'); 
            scatter(find(CS_trace_tmp),zscoredTr(n,find(CS_trace_tmp)),'b.'); 
           
            zz=input(sprintf('Need spike correction?\n',n));
            if zz>0
            spframe=interactive_spike_finder(zscoredTr_hi(n,:));
            [spikeMat(n,:)]=ind2vec(nTime,spframe,1);
            [CStrace(n,:)]=interactive_cs_finder(spikeMat(n,:),zscoredTr(n,:));
            SpClassMattmp=[];
            SpClassMattmp(1,:)=spikeMat(n,:).*(CStrace(n,:)==0);
            SpClassMattmp(2,:)=spikeMat(n,:).*(CStrace(n,:)>0);
            SpClassMattmp=interactiveFrameCheck(zscoredTr(n,:),SpClassMattmp,200,'trace');
            addedCStrace=find(SpClassMattmp(2,:))'+[0:50];
            addedCStrace(addedCStrace>nTime)=[];
            CStrace(n,addedCStrace)=1;
            disp('Correction done.');
            else
                SpClassMattmp=[];
                SpClassMattmp(1,:)=spikeMat(n,:).*(CStrace(n,:)==0);
                SpClassMattmp(2,:)=spikeMat(n,:).*(CStrace(n,:)>0);
            end
            SpClassMat(n,:)=SpClassMattmp(1,:)+SpClassMattmp(2,:)*2;
        end
    end
    VoltResult.isvalideCell=isvalideCell;
    VoltResult.spike=spikeMat;
    VoltResult.CStrace=CStrace;
    VoltResult.SpClass=SpClassMat;
    Result=VoltResult;
    figure(11); clf;
    set(gcf, 'Color', [1 1 1]);
    show_traces_spikes(zscoredTr(isvalideCell>0,:),spikeMat(isvalideCell>0,:),sum(SpClassMat==2,1));
    drawnow;
    saveas(gcf,fullfile(fpath{f},'SpikeDetectionResult.png'),'png')
    save([fpath{f} '/OP_Result.mat'],'Result','-v7.3')
    disp('Spike information added to result.');
end