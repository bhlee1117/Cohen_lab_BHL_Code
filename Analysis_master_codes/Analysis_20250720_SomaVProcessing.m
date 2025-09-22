% Analysis on AAV expression sample and plot, YQ601
% 2025/07/20, Byung Hun Lee
clear
clc;
[~, ~, raw] = xlsread(['/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/NaVInactivationData_Arrangement.xlsx'], ...
    'Sheet1', 'C5:Q415');
fpath=raw(:,1);
Mouse=cell2mat(raw(:,2));
NeuronInd=cell2mat(raw(:,5));
CamType=raw(:,3);
StructureData=raw(:,10);
StimROI=raw(:,6);
StimWfn=raw(:,7);
place_bin=150; time_segment=35000; overlap=200;
set(0,'DefaultFigureWindowStyle','docked')
%% Motion correction
figure; clf;
for i=366:length(fpath)
    load([fpath{i} '/output_data.mat'])
    try
        sz=double(Device_Data{1, 3}.ROI([2 4]));
    catch
        sz=double(Device_Data{1, 4}.ROI([2 4]));
    end
    mov=double(readBinMov([fpath{i} '/frames1.bin'],sz(2),sz(1)));
    %mov=double(readBinMov_times([fpath{i} '/frames1.bin'],sz(2),sz(1),[1:30000]));

    mov_test=mov(:,:,150:350);
    try mov_test = single(mov_test)./single(max(mov_test.data(:)));
    catch disp('change to vm')

        mov_test=vm(mov_test); mov_test = single(mov_test)./single(max(mov_test.data(:))); end
    mov_test = movmean(mov_test,10,3);
    mov_ref = squeeze(median(mov_test,3));
    try
        [mov_mc,xyField]=optical_flow_motion_correction_LBH(vm(mov),mov_ref,'normcorre');
    catch
        [mov_mc,xyField]=optical_flow_motion_correction_LBH_ROIBlock(vm(mov),mov_ref,'normcorre');
    end
    mov_mc=vm(mov_mc);
    mov_mc.transpose.savebin([fpath{i} '/mc.bin'])
    %mcTrace = squeeze(mean(xyField,[1 2]));
    save([fpath{i} '/mcTrace.mat'],'xyField')

    nexttile([1 1])
    imagesc(im_merge(cat(3,mat2gray(std(mov,0,3)),mat2gray(std(double(mov_mc),0,3))),[1 0 0; 1 1 1]));
    axis tight off;
    drawnow;
end

%moviefixsc(mov)
%hold all
%plot(dmd_current_roi{1}(:,1),dmd_current_roi{1}(:,2))
%Blue=DAQ_waves.amplitude(4,round([1:size(mov,3)]*1.27*1e-3/1e-5));

%% ROI setting
for i=143:163%length(fpath)
    Result=[];
    load(fullfile(fpath{i},"output_data.mat"))
    disp([fpath{i} ' , file:' num2str(i)])
    cd(fpath{i})

    try
        sz=double(Device_Data{1, 3}.ROI([2 4]));
        Result.CamType='fusion';
        camind=3;
    catch
        sz=double(Device_Data{1, 4}.ROI([2 4]));
        Result.CamType='flash';
        camind=4;
    end
    ref_time=[2500:5000];
    load(fullfile(fpath{i},'mcTrace.mat'));
    mov_mc=double(readBinMov_times([fpath{i} '/mc.bin'],sz(2),sz(1),ref_time));
    avgImg=mean(mov_mc,3);

    [u,s,v] = svds(tovec(mov_mc-mean(mov_mc,3)),20);
    reshape_u=reshape(u,sz(2),sz(1),[]);
    Result.bvMask=[];
    [~, Result.bvMask]=get_ROI(max(abs(reshape_u),[],3),Result.bvMask);

    Result.ref_im=mean(mov_mc,3);
    mov_res= mov_mc-mean(mov_mc,3);
    mov_res = SeeResiduals(mov_res,mcTrace(ref_time,:));
    mov_res = SeeResiduals(mov_res,mcTrace(ref_time,:).^2);
    mov_res = SeeResiduals(mov_res,mcTrace(ref_time,1).*mcTrace(ref_time,end));
    mov_res = mov_res.*double(max(Result.bvMask,[],3)==0);

    ROIimg=mean(mov_mc,3);

    figure(3); clf;
    imshow2(ROIimg,[])
    title('set extraction ROI');
    % set extraction ROI
    g=1; ROIpoly=[];
    while g
        h = drawpolygon('Color','r');
        if size(h.Position,1)==1 %no more ROI
            g=0;
        else
            ROIpoly=[ROIpoly; {h.Position}];
            hold all
            plot(h.Position(:,1),h.Position(:,2))
        end
    end
    close(figure(3));
    Result.ROIpoly=ROIpoly;

    n_comp=6;
    mov_res_reg=mov_res;
    mov_filt=imgaussfilt3(mov_res_reg.*double(max(Result.bvMask,[],3)==0),[1.5 1.5 0.1]);
    movVec=tovec(mov_filt);
    Npoly=size(Result.ROIpoly,1);
    ftprnt = zeros(size(mov_filt,1)*size(mov_filt,2),Npoly);
    clear mask
    figure(4); clf;
    for p=1:Npoly %each ROIs
        p
        clf; ax2=[];
        tiledlayout(n_comp/2+2,2)
        mask(:,:,p) = poly2mask(Result.ROIpoly{p}(:,1), Result.ROIpoly{p}(:,2), sz(2), sz(1));
        pixelList=find(tovec(squeeze(mask(:,:,p))));
        subMov = movVec(pixelList,:);
        [V D eigTrace]=get_eigvector(subMov);
        nexttile([2 2])
        plot(rescale2(eigTrace(:,1:n_comp),1)+[1:n_comp])

        for n=1:n_comp
            eigImg=NaN(size(mov_filt,1)*size(mov_filt,2),1);
            ax2=[ax2 nexttile([1 1])];
            eigImg(pixelList,1)=V(:,n);
            eigImg=toimg(eigImg,size(mov_filt,1),size(mov_filt,2));
            imshow2(im_merge(cat(3,mat2gray(Result.ref_im),(eigImg)*10),[1 1 1;1 0 0]),[])
            title([num2str(n) ', Fraction: ' num2str(D(n)/sum(D),2)])
        end
        linkaxes(ax2,'xy')
        n_take = input('#components to take: ', 's');
        n_take = str2num(n_take);
        coeff=subMov*mean(eigTrace(:,n_take)*V(:,n_take)',2);
        ftprnt(pixelList,p)=coeff;
    end
    close(figure(4));
    Result.ftprnt=toimg(ftprnt,sz(2),sz(1));
    ref_ftprnt=Result.ftprnt;

    figure(99); clf;
    show_footprnt(Result.ftprnt,Result.ref_im)
    save(fullfile(fpath{i},'SomOP_Result.mat'),'Result','-v7.3')
end

%% Signal extraction from multiple movie files, in streaming mode
bound = 5;
for f=164:365%366:length(fpath)

    disp([fpath{f} ' , file:' num2str(f)])
    load([fpath{f} '/SomOP_Result.mat'])
    load(fullfile(fpath{f},"output_data.mat"))
    if ~isfield(Result,'CamType')
        Result.CamType='flash';
    end
    switch char(Result.CamType)
        case 'flash'
            sz=double(Device_Data{1, 4}.ROI([2 4]));
        case 'fusion'
            sz=double(Device_Data{1, 3}.ROI([2 4]));
    end

    CamCounter=Device_Data{1, 2}.Counter_Inputs(1, 1).data;
    CamTrigger=find(CamCounter(2:end)-CamCounter(1:end-1));

    Result.Blue=Device_Data{1, 2}.buffered_tasks(1, 1).channels(1, 2).data;
    Result.Blue=Result.Blue(CamTrigger);

    Result.traces=[];
    Result.traces_bvMask=[];
    Result.mc=[];

    ref_im_vec=tovec(Result.ref_im(bound:end-bound,bound:end-bound));
    load([fpath{f} '/mcTrace.mat']);
    if ~isfield(xyField,'xymean')
        Shift2=xyField; clear xyField;
        shifts_r = squeeze(cat(3,Shift2(:).shifts));
        shifts_nr = cat(ndims(Shift2(1).shifts)+1,Shift2(:).shifts);
        shifts_nr = reshape(shifts_nr,[],2,length(Shift2));
        shifts_x = squeeze(shifts_nr(:,2,:))';
        shifts_y = squeeze(shifts_nr(:,1,:))';
        xyField.xymean=[shifts_x; shifts_y];
        xyField.xymean=xyField.xymean';
    end
    motionTrace=movmean(xyField.xymean,5,1);
    %motionTrace=movmean(mcTrace,5,1);
    Result.mc=[Result.mc; motionTrace];

    mov_mc=double(readBinMov([fpath{f} '/mc.bin'],sz(2),sz(1)));
    disp('Movie loaded!')
    if length(Result.Blue)<size(mov_mc,3)
        mov_mc=mov_mc(:,:,1:length(Result.Blue));
        Result.mc=Result.mc(1:length(Result.Blue),:);
    end

    mov_res= mov_mc-mean(mov_mc,3);
    bkg = zeros(1, size(mov_mc,3));
    mean_F=squeeze(mean(mov_mc(bound:end-bound,bound:end-bound,:),[1 2]));
    [~, blueOff]=get_blueoffTrace(mean_F,[Result.Blue],70);
    [y_fit]=expfitDM_2(find(blueOff)',mean_F(find(blueOff)),[1:size(mov_mc,3)]',1000);
    bkg(1,:)=y_fit;
    mov_res = SeeResiduals(mov_res,Result.mc);
    mov_res = SeeResiduals(mov_res,Result.mc.^2);
    mov_res = SeeResiduals(mov_res,Result.mc(:,1).*Result.mc(:,end));
    mov_res= SeeResiduals(mov_res,bkg,1);
    %mov_res=mov_res.*double(max(Result.bvMask,[],3)==0);

    if ~isfield(Result,'ftprnt')
        Result.ftprnt=Result.c_ftprnt;
    end
    Result.traces=[Result.traces, -(tovec(mov_res)'*tovec(Result.ftprnt))'];
    Result.traces=(Result.traces-prctile(Result.traces,30,2))./(y_fit/max(y_fit))';
    if isfield(Result,'bvMask')
        Result.traces_bvMask=[Result.traces_bvMask, -(tovec(mov_res)'*tovec(Result.ftprnt))'];
        Result.traces_bvMask=Result.traces_bvMask-prctile(Result.traces_bvMask,30,2)./(y_fit/max(y_fit))';
    else
        Result.bvMask=[];
        Result.traces_bvMask=Result.traces;
    end


    clf; tiledlayout(2,2)
    nexttile([1 1])
    show_footprnt_contour(Result.ftprnt,Result.ref_im); colormap(gray(256));
    nexttile([1 1])
    show_footprnt_contour(Result.bvMask,Result.ref_im); colormap(gray(256));
    nexttile([1 2])
    plot(rescale2(Result.traces_bvMask,2)'+[1:size(Result.traces_bvMask,1)]);
    saveas(gca,[fpath{f},'/Traces.png']);
    drawnow;
    save([fpath{f} '/SomOP_Result.mat'],'Result','-v7.3')
    disp([num2str(f) ' th file is saved.'])
end

%%
isgoodNeuron=[];
%%
for f=371:length(fpath)
    f;
    load([fpath{f} '/SomOP_Result.mat'])
    [nROI nTime]=size(Result.traces);
    Result.normtrace=Result.traces./get_threshold(Result.traces,1);
    if 1%~isfield(Result,'spike')
        tr=Result.normtrace-movmedian(Result.normtrace,50,2);
        Result.spike=find_spike_bh(tr,12,3);
    end
    %Result.spike=Result.spike(:,1:nTime);
    for n=1:nROI
        if length(find(Result.spike(n,:)))>2
            y_fit=expfitDM_2(find(Result.spike(n,:))',Result.normtrace(n,find(Result.spike(n,:)))',[1:nTime]',10000);
        else
            y_fit=ones(nTime,1);
        end
        lowpass=NaN(1,nTime); lowpass2=NaN(1,nTime);
        [~, blueOff]=get_blueoffTrace(Result.normtrace(n,:),[Result.Blue],200,200);
        lowpass(blueOff)=Result.normtrace(n,blueOff>0);
        lowpass2(blueOff)=movmedian(lowpass(blueOff),1000);
        lowpass2=interpolateNaN(lowpass2);
        lowpass(isnan(lowpass))=lowpass2(isnan(lowpass));
        lowpass=movmean(lowpass,500,'omitnan');

        Result.normtrace(n,:)=Result.normtrace(n,:)-lowpass;
        clf; tiledlayout(2,1);
        nexttile([1 1]);
        plot(Result.normtrace(n,:)); hold all
        plot(find(Result.spike(n,:)),Result.normtrace(n,find(Result.spike(n,:))),'ro');
        plot(y_fit);
        plot(lowpass);
        %Result.normtrace(n,:)=Result.normtrace(n,:)./(y_fit/max(y_fit))';
        nexttile([1 1]);
        plot(Result.normtrace(n,:)); hold all
        plot(find(Result.spike(n,:)),Result.normtrace(n,find(Result.spike(n,:))),'ro');
        isgood=input(['File: ' num2str(f) ' is Good? \n']);
        if ~ismember(isgood,[0 1])
            break;
        else
            isgoodNeuron=[isgoodNeuron; [f n isgood]];
        end
        Result.isGood(n)=isgood;
    end
    save([fpath{f} '/SomOP_Result.mat'],'Result','-v7.3')
end
%%
CS_thres=[5 1.5];
for f=411:length(fpath)
    load([fpath{f} '/SomOP_Result.mat'])
    f
    [nROI nTime]=size(Result.normtrace);
    for n=find(Result.isGood)
        figure(201); clf;
        n
        tr_ref=Result.normtrace(n,:);
        %sp_soma=Result.spike(1,:);
        sp_soma=find_spike_bh(tr_ref-movmedian(tr_ref,50,2),2.5,10);

        tr_hi=mean(tr_ref,1)-movprc(mean(tr_ref,1),200,20,2);
        tr_sub=get_subthreshold(tr_hi,sp_soma,7,17);

        %find CS
        [trans tr_trace]=detect_transient2(tr_sub,CS_thres,sp_soma,20);
        if isempty(trans.amp)
            CS_trace=zeros(1,nTime);
        else
            transcand=cell2mat(cellfun(@(x) length(x)>1,trans.ISI,'UniformOutput',false));
            meanISI_frnt=cellfun(@(x) mean(x(1:2)),trans.ISI(transcand));
            meanISI_first3=zeros(1,length(trans.length));
            meanISI_first3(transcand)=meanISI_frnt;

            %CS_ind=find(trans.spike_number>2 & trans.mean_ISI<15);
            CS_ind=find(trans.spike_number>2 & meanISI_first3<20);
            CS_trace=ismember(tr_trace,CS_ind);
            bwCS=bwlabel(CS_trace);
            CS_spike=sp_soma.*bwCS;
            [~, CS_spike_time]=unique(CS_spike);
            for b=1:max(bwCS)
                bfrm=find(bwCS==b);
                CS_trace(bfrm(1):CS_spike_time(bwCS(CS_spike_time)==b)-5)=0;
            end
        end

        SpikeClassMat=zeros(2,size(tr_ref,2));
        SpikeClassMat(1,:)=sp_soma.*(1-CS_trace); %bAPs
        SpikeClassMat(2,CS_spike_time(2:end))=1; %Complex spikes

        Result.spike(n,:)=[sp_soma];
        bwCS=bwlabel(CS_trace);
        CStrace=ismember(bwCS,bwCS(find(SpikeClassMat(2,:)>0)+5));
        Result.CStrace(n,:)=CS_trace;

        SpikeClassMat(1,:)=sp_soma.*(1-CS_trace); %bAPs
        Result.SpClass(:,:,n)=SpikeClassMat;

        TraceonlyCS=NaN(1,nTime);
        TraceonlyCS(1,Result.CStrace(n,:)>0)=Result.normtrace(n,Result.CStrace(n,:)>0);

        plot(Result.normtrace(n,:)); hold all;
        plot(TraceonlyCS,'r');
        plot(find(Result.spike(n,:)>0),Result.normtrace(n,Result.spike(n,:)>0),'ro');

        isGood=input('Is good to go? \n');
        while isGood==0
            sp_soma = interactiveFrameCheck(tr_ref, sp_soma, 300,'trace');
            tr_hi=mean(tr_ref,1)-movprc(mean(tr_ref,1),200,20,2);
            tr_sub=get_subthreshold(tr_hi,sp_soma,7,17);

            %find CS
            [trans tr_trace]=detect_transient2(tr_sub,CS_thres,sp_soma,20); 
            if isempty(trans.amp)
                CS_trace=zeros(1,nTime);
            else
                transcand=cell2mat(cellfun(@(x) length(x)>1,trans.ISI,'UniformOutput',false));
                meanISI_frnt=cellfun(@(x) mean(x(1:2)),trans.ISI(transcand));
                meanISI_first3=zeros(1,length(trans.length));
                meanISI_first3(transcand)=meanISI_frnt;

                %CS_ind=find(trans.spike_number>2 & trans.mean_ISI<15);
                CS_ind=find(trans.spike_number>2 & meanISI_first3<20);
                CS_trace=ismember(tr_trace,CS_ind);
                bwCS=bwlabel(CS_trace);
                CS_spike=sp_soma.*bwCS;
                [~, CS_spike_time]=unique(CS_spike);
                for b=1:max(bwCS)
                    bfrm=find(bwCS==b);
                    CS_trace(bfrm(1):CS_spike_time(bwCS(CS_spike_time)==b)-5)=0;
                end
            end

            SpikeClassMat=zeros(2,size(tr_ref,2));
            SpikeClassMat(1,:)=sp_soma.*(1-CS_trace); %bAPs
            SpikeClassMat(2,CS_spike_time(2:end))=1; %Complex spikes

            SpikeClassMat = interactiveFrameCheck(tr_ref, SpikeClassMat, 300,'trace');
            CS_trace(unique(find(SpikeClassMat(2,:))'+[-3:50]))=1;
            CS_trace(unique(find(SpikeClassMat(1,:))'+[-3:5]))=0;
            CS_trace=CS_trace(1:nTime);
            Result.spike(n,:)=[sp_soma];
            bwCS=bwlabel(CS_trace);
            bIndwithCS=bwCS(find(SpikeClassMat(2,:)>0)+5); bIndwithCS(bIndwithCS==0)=[];
            CS_trace=ismember(bwCS,bIndwithCS);
            Result.CStrace(n,:)=CS_trace;

            SpikeClassMat(1,:)=sp_soma.*(1-CS_trace); %bAPs
            Result.SpClass(:,:,n)=SpikeClassMat;

            TraceonlyCS=NaN(1,nTime);
            TraceonlyCS(1,Result.CStrace(n,:)>0)=Result.normtrace(n,Result.CStrace(n,:)>0);

            figure(f); clf;

            plot(Result.normtrace(n,:)); hold all;
            plot(TraceonlyCS,'r');
            plot(find(Result.spike(n,:)>0),Result.normtrace(n,Result.spike(n,:)>0),'ro');
            plot(find(Result.SpClass(2,:,n)>0),Result.normtrace(n,find(Result.SpClass(2,:,n)>0)),'co')
            plot(find(Result.SpClass(1,:,n)>0),Result.normtrace(n,find(Result.SpClass(1,:,n)>0)),'yo')
            plot(Result.CStrace(n,:));
            isGood=input('Is good to go? \n');
        end
        disp('done');
    end
    save([fpath{f} '/SomOP_Result.mat'],'Result','-v7.3')
end