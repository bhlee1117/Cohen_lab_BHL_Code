
clear
clc;
cd '/Volumes/BHL18TB_D2/Arranged_Data/Prism_OptopatchResult';
[~, ~, raw] = xlsread(['/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/' ...
    'Prism_OptopatchData_Arrangement.xlsx'], 'Sheet1', 'B5:K165');

save_to='/Volumes/BHL18TB_D2/Arranged_Data/Prism_OptopatchResult';
fpath=raw(:,1);
Mouse=cell2mat(raw(:,2));
NeuronInd=cell2mat(raw(:,5));
CamType=raw(:,3);
StructureData=raw(:,10);
%%
figure;
[~, unqInd] = unique([Mouse NeuronInd] ,'row');
%tiledlayout(4,4);
for i=unqInd'
    nexttile([1 1])
    load(fullfile(fpath{i},"output_data.mat"))
    switch char(CamType(i))
        case 'flash'
            sz=double(Device_Data{1, 4}.ROI([2 4]));
        case 'fusion'
            sz=double(Device_Data{1, 3}.ROI([2 4]));
    end
    ref_time=[2000:3000];
    mov_test=double(readBinMov_times([fpath{i} '/frames1.bin'],sz(2),sz(1),[ref_time(1):ref_time(end)]));
    imshow2(mean(mov_test,3),[])
    title(['Mouse #' num2str(Mouse(i)) '-Neuron#' num2str(NeuronInd(i))])
end

%% MC

for i=[160]%length(fpath)

    load(fullfile(fpath{i},"output_data.mat"))

    ref_time=[6000:7000]; overlap=200;
    time_segment=25000;

    frm_end=max(Device_Data{1, 2}.Counter_Inputs(1, 1).data);
    f_seg=[[1:time_segment:frm_end] frm_end+1]; f_seg(2:end)=f_seg(2:end)-1;
    CamDAQ_rate=Device_Data{1, 2}.Counter_Inputs.rate;
    CamTrig=Device_Data{1, 2}.Counter_Inputs.data;
    CamTrig2=find(CamTrig(2:end)-CamTrig(1:end-1)>0);
    Frm_rate=(CamTrig2(2)-CamTrig2(1))/CamDAQ_rate;

    switch char(CamType(i))
        case 'flash'
            sz=double(Device_Data{1, 4}.ROI([2 4]));
        case 'fusion'
            sz=double(Device_Data{1, 3}.ROI([2 4]));
    end

    mov_test=double(readBinMov_times([fpath{i} '/frames1.bin'],sz(2),sz(1),[ref_time(1):ref_time(end)]));
    mov_test=mov_test(:,:,2:end);
    [mov_test_mc,xyField]=optical_flow_motion_correction_LBH(mov_test,mean(mov_test,3),'normcorre');
    mov_test=vm(mov_test);
    mov_test = single(mov_test)./single(max(mov_test.data(:)));
    mov_test = movmean(mov_test,10,3);
    mov_ref = squeeze(median(mov_test,3));

    for j=1:length(f_seg)-1

        mov=double(readBinMov_times([fpath{i} '/frames1.bin'],sz(2),sz(1),[f_seg(j):f_seg(j+1)]));


        switch char(CamType(i))
            case 'flash'
                mov=rollingShutter_correction(mov,1/Frm_rate,'flash');
            case 'fusion'
                mov=rollingShutter_correction(mov,1/Frm_rate,'fusion');
        end

        mov=vm(mov(:,:,2:end));
        if j==1
            mov=mov(:,:,[1 1:size(mov,3)]);
        end

        [mov_mc,xyField]=optical_flow_motion_correction_LBH(mov,mov_ref,'normcorre');

        ave_im=mean(mov_mc,3);
        mov_mc=vm(mov_mc);
        mov_mc.transpose.savebin([fpath{i} '/mc_ShutterReg' num2str(j,'%02d') '.bin'])

        %        mcTrace = squeeze(mean(xyField,[1 2])); %optic flow
        mcTrace=xyField; % Normcorre
        save([fpath{i} '/mcTrace' num2str(j,'%02d') '.mat'],'mcTrace','ave_im')

        %  clear mov_mc mov
    end
    sd_mov=std(double(mov),0,3); sd_mov_mc=std(double(mov_mc),0,3);
    figure; clf;
    nexttile([1 1]); imshow2(sd_mov,[]); title('before mc')
    nexttile([1 1]); imshow2(sd_mov_mc,[]); title('after mc')
    nexttile([1 1]); imshow2(imfuse(sd_mov,sd_mov_mc),[]);
    title(fpath{i},'Interpreter',  'none')
    saveas(gca,[char(fpath{i}) '/' 'MC_result.fig'])
end

%% ROI setting

[a, unqInd] = unique([Mouse NeuronInd] ,'row');
figure;
for i=unqInd([9])'
    nexttile([1 1])
    load(fullfile(fpath{i},"output_data.mat"))
    disp(fpath{i})

    switch char(CamType(i))
        case 'flash'
            sz=double(Device_Data{1, 4}.ROI([2 4]));
        case 'fusion'
            sz=double(Device_Data{1, 3}.ROI([2 4]));
    end

    ref_time=[2000:5000];
    mov_test=double(readBinMov_times([fpath{i} '/mc_ShutterReg' num2str(1,'%02d') '.bin'],sz(2),sz(1),ref_time));
    avgImg=mean(mov_test,3);
    figure(3); clf;
    imshow2(avgImg,[])
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
    Result.ref_im=avgImg;
        % blood vessel masking
[u,s,v] = svds(tovec(mov_test-mean(mov_test,3)),20);
reshape_u=reshape(u,sz(2),sz(1),[]);
bvMask=[];
[~, bvMask]=get_ROI(max(abs(reshape_u),[],3),bvMask);

    save(fullfile(fpath{i},'Result.mat'),'Result','-v7.3')

    SameCellInd=find(Mouse==Mouse(i) & NeuronInd==NeuronInd(i));
    
    for j=SameCellInd'
        Result.ROIpoly=ROIpoly;
        mov_mc=double(readBinMov_times([fpath{j} '/mc_ShutterReg' num2str(1,'%02d') '.bin'],sz(2),sz(1),[3000:4000]));
        Result.ref_im=mean(mov_mc,3);
        [offsetY offsetX]=calculate_shift(avgImg,mean(mov_mc,3))
        ROIpoly_shift=cellfun(@(x) x+[offsetX offsetY],Result.ROIpoly,'UniformOutput',false);
        Result.ROIpoly=ROIpoly_shift;        
        
        [x, y] = meshgrid(1:sz(1), 1:sz(2));
        Result.bvMask=[];
        for n=1:size(bvMask,3)
        Result.bvMask(:,:,n) = interp2(x, y, bvMask(:,:,n), x+offsetX, y+offsetY, 'linear', 0);
        end
        save([fpath{j} '/Result.mat'],'Result','-v7.3')
    end
end

%% Set Footprint
bound=6;
for i=unqInd([9])'
    load([fpath{i} '/Result.mat'])
    load(fullfile(fpath{i},"output_data.mat"))
    switch char(CamType(i))
        case 'flash'
            sz=double(Device_Data{1, 4}.ROI([2 4]));
        case 'fusion'
            sz=double(Device_Data{1, 3}.ROI([2 4]));
    end

    ref_time=[1000:9000];
    load(fullfile(fpath{i},['/mcTrace' num2str(1,'%02d') '.mat']));

    mov_mc=double(readBinMov_times([fpath{i} '/mc_ShutterReg' num2str(1,'%02d') '.bin'],sz(2),sz(1),ref_time));
    averageImg=mean(mov_mc,3);
    Result.ref_im=mean(mov_mc,3);
    Result.mc=mcTrace.xymean;
    mov_res= mov_mc-mean(mov_mc,3);
    mcTrace.xymean=movmean(mcTrace.xymean,3,2);
    mov_res = SeeResiduals(mov_res,mcTrace.xymean(ref_time,:));
    mov_res = SeeResiduals(mov_res,mcTrace.xymean(ref_time,:).^2);
    mov_res = SeeResiduals(mov_res,mcTrace.xymean(ref_time,1).*mcTrace.xymean(ref_time,2));

    n_comp=6;
    mov_filt=imgaussfilt3(mov_res.*double(max(Result.bvMask,[],3)==0),[2 2 0.1]);
    %mov_filt=imgaussfilt3(mov_res,[2 2 0.1]);
    movVec=tovec(mov_filt);
    Npoly=size(Result.ROIpoly,1);
    ftprnt = zeros(size(mov_filt,1)*size(mov_filt,2),Npoly);
    clear mask
    figure(4);
    for p=1:Npoly %each ROIs
        clf; ax2=[];
        tiledlayout(n_comp/2+2,2)
        mask(:,:,p) = poly2mask(Result.ROIpoly{p}(:,1), Result.ROIpoly{p}(:,2), sz(2), sz(1));
        pixelList=find(tovec(squeeze(mask(:,:,p))));
        subMov = movVec(pixelList,:);
        covMat = subMov*subMov';  % PCA within each region
        [V, D] = eig(covMat);
        D = diag(D);
        D = D(end:-1:1);
        V = V(:,end:-1:1);
        vSign = sign(max(V) - max(-V));  % make the largest value always positive
        V = V.*vSign;
        eigTrace=subMov'*V;
        nexttile([2 2])
        plot(rescale2(eigTrace(:,1:n_comp),1)+[1:n_comp])
        for n=1:n_comp
            eigImg=NaN(size(mov_filt,1)*size(mov_filt,2),1);
            ax2=[ax2 nexttile([1 1])];
            eigImg(pixelList,1)=V(:,n);
            eigImg=toimg(eigImg,size(mov_filt,1),size(mov_filt,2));
            imshow2(im_merge(cat(3,Result.ref_im,eigImg),[1 1 1;1 0 0]),[])
            title([num2str(n) ', Fraction: ' num2str(D(n)/sum(D),2)])
        end
        linkaxes(ax2,'xy')
        n_take = input('#components to take: ', 's');
        n_take = str2num(n_take);
        coeff=subMov*mean(eigTrace(:,n_take)*V(:,n_take)',2);
        ftprnt(pixelList,p)=coeff;
    end

    Result.ftprnt=toimg(ftprnt,sz(2),sz(1));
    Result.ftprnt(1:bound,:,:)=0;
    Result.ftprnt(:,1:bound,:)=0;
    Result.ftprnt(end-bound:end,:,:)=0;
    Result.ftprnt(:,end-bound:end,:)=0;
    ref_ftprnt=Result.ftprnt;

    figure; clf;
    SameCellInd=find(Mouse==Mouse(i) & NeuronInd==NeuronInd(i));
    for j=SameCellInd'
        load([fpath{j} '/Result.mat'])
        mov_mc=double(readBinMov_times([fpath{j} '/mc_ShutterReg' num2str(1,'%02d') '.bin'],sz(2),sz(1),[3000:4000]));
        [offsetY offsetX]=calculate_shift(averageImg,mean(mov_mc,3));
        [x, y] = meshgrid(1:sz(1), 1:sz(2));
        for n=1:size(ref_ftprnt,3)
        Result.ftprnt(:,:,n) = interp2(x, y, ref_ftprnt(:,:,n), x+offsetX, y+offsetY, 'linear', 0);
        end
        Result.ref_im=mean(mov_mc,3);
        coord_1d=dim_reduce(get_coord(Result.ftprnt));
[~, Result.dist_order]=sort(coord_1d,'descend');

        save([fpath{j} '/Result.mat'],'Result','-v7.3')
        nexttile([1 1])
        show_im=mean(mov_mc,3);
        show_im(show_im<15)=median(show_im(:));
        show_im=show_im-imgaussfilt(show_im,15);
        show_im(show_im<prctile(show_im(:),20))=prctile(show_im(:),20);
        imshow2(imfuse(max(Result.ftprnt,[],3),show_im),[])
    end
end


%% Signal extraction
bound=5;

for i=[72:75]%length(fpath)]

    load([fpath{i} '/Result.mat'])
    load(fullfile(fpath{i},"output_data.mat"))
    load([fpath{i} '/mcTrace' num2str(1,'%02d') '.mat']);
    switch char(CamType(i))
        case 'flash'
            sz=double(Device_Data{1, 4}.ROI([2 4]));
        case 'fusion'
            sz=double(Device_Data{1, 3}.ROI([2 4]));
    end

    Result.traces=[];
    Result.traces_bvMask=[];
    Result.mc=movmean(mcTrace.xymean,5,1);
    Result.im_corr=[];

    ref_im_vec=tovec(Result.ref_im(bound:end-bound,bound:end-bound));

    CamCounter=Device_Data{1, 2}.Counter_Inputs(1, 1).data;
    CamTrigger=find(CamCounter(2:end)-CamCounter(1:end-1));

    Result.Blue=Device_Data{1, 2}.buffered_tasks(1, 1).channels(1, 2).data;
    Result.Blue=Result.Blue(CamTrigger);
    Rfixed = imref2d(repmat(Device_Data{1, 3}.virtualSensorSize,1,2));
    inverseTform = invert(Device_Data{1, 6}.tform);
    revertedImage = imwarp(double(Device_Data{1, 6}.Target), inverseTform,'OutputView',Rfixed);
    [blueDMDimg blueDMDcontour]=imcrop(revertedImage,double(Device_Data{1, 3}.ROI([1 3 2 4]))+[0 0 -1 -1]);
    Result.BlueDMD=blueDMDcontour;
    Result.BlueDMDimg=blueDMDimg;

    mov_mc=double(readBinMov_times([fpath{i} '/mc_ShutterReg' num2str(1,'%02d') '.bin'],sz(2),sz(1),[1:length(CamTrigger)]));

    mov_mc_vec=tovec(mov_mc(bound:end-bound,bound:end-bound,:));
    mov_mc_vec=(mov_mc_vec-mean(mov_mc_vec,1))./std(mov_mc_vec,0,1);

    mov_res= mov_mc-mean(mov_mc,3);
    bkg = zeros(1, size(mov_mc,3));

    %     bkg(1,:) = linspace(-1, 1, size(mov_mc,3));  % linear term
    %     bkg(2,:) = linspace(-1, 1, size(mov_mc,3)).^2;  % quadratic term
    mean_F=squeeze(mean(mov_mc(bound:end-bound,bound:end-bound,:),[1 2]));
    [~, blueOff]=get_blueoffTrace(mean_F,[Result.Blue],70);
    [y_fit]=expfitDM_2(find(blueOff)',mean_F(find(blueOff)),[1:size(mov_mc,3)]',1000);
    bkg(1,:)=y_fit;
    mov_res = SeeResiduals(mov_res,Result.mc);
    mov_res = SeeResiduals(mov_res,Result.mc.^2);
    mov_res = SeeResiduals(mov_res,Result.mc(:,1).*Result.mc(:,end));
    mov_res= SeeResiduals(mov_res,bkg,1);

    Result.traces=[-(tovec(mov_res)'*tovec(Result.ftprnt))'];
    if isfield(Result,'bvMask')
    Result.traces_bvMask=[-(tovec(mov_res.*double(max(Result.bvMask,[],3)==0))'*tovec(Result.ftprnt))'];
    end
    Result.im_corr=[Result.im_corr corr(rescale2(mov_mc_vec,1),ref_im_vec,'type','Spearman')'];  %image correlation
    %Result.im_corr=[sum(mov_mc_vec.*ref_im_vec,1)/(size(mov_mc_vec,1)-1)];  %image correlation
    save([fpath{i} '/Result.mat'],'Result','-v7.3')
end

%% Clean up and norm
exclude_frq=[241.7 242]; %monitor
%exclude_frq2=[483.5 484]; %monitor
exclude_frq2=[55.5 56]; %motion
time_bin=15000; Fs=1000;

for f=[160]%:length(fpath)
    load(fullfile(fpath{f},'Result.mat'),'Result')
    ref_trace=1;
    nTime=size(Result.traces,2);
    nROI=size(Result.traces,1);
tr=Result.traces-movmedian(Result.traces,50,2);
tr=tr./get_threshold(tr,1);
sp=find_spike_bh(tr,5,3);

SilentPeriod=ones(1,size(Result.traces,2));
sp_time=find(max(sp,[],1))';
sp_na=sum((find(max(sp,[],1))'+[-10:150])<0 | (find(max(sp,[],1))'+[-10:150])>size(Result.traces,2),2)==0;
SilentPeriod(sp_time(sp_na)+[-10:150])=NaN;

lwpass_fit=NaN(nROI,nTime); lwpass_fit_ch{1}=NaN(nROI,nTime); lwpass_fit_ch{2}=NaN(nROI,nTime);
t_silent=find(~isnan(SilentPeriod) & Result.Blue==0);
%[y_fit2 t_consts coeffY]  = expfitDM_2(t_fit',-Result.traces(1,t_fit)',[1:nTime]',[100000 1000]);
lwpass_fit=movmedian(get_blueoffTrace(Result.traces(1,:),ind2vec([size(Result.traces,2)],t_silent,0,1),30),1000,'omitnan');

tr_res=squeeze(SeeResiduals(permute(Result.traces_bvMask,[1 3 2]),lwpass_fit));
% tr_res=Result.traces_bvMask-lwpass_fit;
% tr_res_checker{1}=Result.traces_checker{1}-lwpass_fit_ch{1};
% tr_res_checker{2}=Result.traces_checker{2}-lwpass_fit_ch{2};

    freq_lowhigh=exclude_frq/(Fs/2);
    [b, a] = butter(4, freq_lowhigh, 'stop');

    freq_lowhigh2=exclude_frq2/(Fs/2);
    [b2, a2] = butter(4, freq_lowhigh2, 'stop');

    sub_pass_frq=[2];
    freq_lowhigh3=sub_pass_frq/(Fs/2);
    [b3, a3] = butter(4, freq_lowhigh3, 'low');

    theta_pass_frq=[5 11];
    freq_lowhigh4=theta_pass_frq/(Fs/2);
    [b4, a4] = butter(4, freq_lowhigh4, 'bandpass');

    % figure(2); clf;
    % [p f]=fft_simple(PC_Result{f}.traces(2,:),1000);
    % ax1=nexttile([1 1]);
    % plot(f(1,:),p')
    % set(gca,'YScale','log')
    %
    % [p f]=fft_simple(PC_Result{f}.mcTrace,1000);
    % ax2=nexttile([1 1]);
    % plot(f(1,:),p')
    % set(gca,'YScale','log')
    % linkaxes([ax1 ax2],'x')

    clear traces_res_filtered noise noise_intp norm_trace sp_height SpHeight_intp sp_time tr_mc_imcorr tr_mc tr norm_trace_check traces_res_filtered_ch
    tN=[1:time_bin:nTime]; tN=[tN nTime];
    lwpass_fit=NaN(nROI,nTime); lwpass_fit_ch{1}=NaN(nROI,nTime); lwpass_fit_ch{2}=NaN(nROI,nTime);
    sp_time=zeros(nROI,nTime);
    sp_height=zeros(nROI,nTime);
    mcTrace=squeeze(Result.mcTrace)';

    for n=1:size(Result.traces,1)
        tr=tr_res(n,1:nTime);
        tr_check{1}=tr_res_checker{1}(n,1:nTime); tr_check{2}=tr_res_checker{2}(n,1:nTime);

        % regress out motion frequency
        for t=1:length(tN)-1
            tr(tN(t):tN(t+1))=squeeze(SeeResiduals(reshape(tr(tN(t):tN(t+1)),1,1,[]),mcTrace(:,(tN(t):tN(t+1)))));
            tr(tN(t):tN(t+1))=squeeze(SeeResiduals(reshape(tr(tN(t):tN(t+1)),1,1,[]),mcTrace(:,(tN(t):tN(t+1))).^2));
            tr(tN(t):tN(t+1))=squeeze(SeeResiduals(reshape(tr(tN(t):tN(t+1)),1,1,[]),mcTrace(1,(tN(t):tN(t+1))).*mcTrace(end,(tN(t):tN(t+1)))));
            %tr(tN(t):tN(t+1))=squeeze(SeeResiduals(reshape(tr(tN(t):tN(t+1)),1,1,[]),Result.im_corr(:,(tN(t):tN(t+1)))));
            for ch=1:2 % Checkerboard pattern traces
                tr_check{ch}(tN(t):tN(t+1))=squeeze(SeeResiduals(reshape(tr_check{ch}(tN(t):tN(t+1)),1,1,[]),mcTrace(:,(tN(t):tN(t+1)))));
                tr_check{ch}(tN(t):tN(t+1))=squeeze(SeeResiduals(reshape(tr_check{ch}(tN(t):tN(t+1)),1,1,[]),mcTrace(:,(tN(t):tN(t+1))).^2));
                tr_check{ch}(tN(t):tN(t+1))=squeeze(SeeResiduals(reshape(tr_check{ch}(tN(t):tN(t+1)),1,1,[]),mcTrace(1,(tN(t):tN(t+1))).*mcTrace(end,(tN(t):tN(t+1)))));
            end
            tr_mc(tN(t):tN(t+1))=tr(tN(t):tN(t+1));
        end

        % regress out monitor frequency
        traces_res_filtered(n,:) = filtfilt(b, a, tr);
        traces_res_filtered_ch{1}(n,:) = filtfilt(b, a, tr_check{1});
        traces_res_filtered_ch{2}(n,:) = filtfilt(b, a, tr_check{2});

        traces_res_filtered(n,:) = filtfilt(b2, a2, traces_res_filtered(n,:));
        traces_res_filtered_ch{1}(n,:) = filtfilt(b2, a2, traces_res_filtered_ch{1}(n,:));
        traces_res_filtered_ch{2}(n,:) = filtfilt(b2, a2, traces_res_filtered_ch{2}(n,:));
        
        tr_mc_bpMonitorMotion=traces_res_filtered(n,:);

        norm_trace(n,:)=traces_res_filtered(n,:);%-movmedian(traces_res_filtered(n,:),500,2);
        norm_trace_check{1}(n,:)=traces_res_filtered_ch{1}(n,:);%-movmedian(traces_res_filtered(n,:),500,2);
        norm_trace_check{2}(n,:)=traces_res_filtered_ch{2}(n,:);%-movmedian(traces_res_filtered(n,:),500,2);

        lwpass_fit(n,t_silent)=norm_trace(n,t_silent); lwpass_fit(n,:)=movmedian(lwpass_fit(n,:),20000,2,'omitnan');
        lwpass_fit_ch{1}(n,t_silent)=norm_trace_check{1}(n,t_silent); lwpass_fit_ch{1}(n,:)=movmedian(lwpass_fit_ch{1}(n,:),20000,2,'omitnan');
        lwpass_fit_ch{2}(n,t_silent)=norm_trace_check{2}(n,t_silent); lwpass_fit_ch{2}(n,:)=movmedian(lwpass_fit_ch{2}(n,:),20000,2,'omitnan');

        norm_trace(n,:)=norm_trace(n,:)-lwpass_fit(n,:);
        norm_trace_check{1}(n,:)=norm_trace_check{1}(n,:)-lwpass_fit_ch{1}(n,:);
        norm_trace_check{2}(n,:)=norm_trace_check{2}(n,:)-lwpass_fit_ch{2}(n,:);

        if n==ref_trace %compensate spike hight due to bleaching
            for t=1:length(tN)-1
                tr_tmp=norm_trace(n,tN(t):tN(t+1));
                tr_tmp=tr_tmp-movmedian(tr_tmp,300);
                noise(t,n)=get_threshold(tr_tmp,1);
                tr_tmp_norm=tr_tmp./noise(t,n);
                [sp_temp, pks, prom]=find_spike_bh(tr_tmp_norm,6,3);
                sp_time(n,tN(t):tN(t+1))=sp_temp;
            end

            t_fit=find(sp_time(n,:));
            sp_height(n,t_fit)=norm_trace(n,t_fit);
            tx=tN(1:end-1)+time_bin/2;
            %noise_intp(n,:)=movmean(interp1(tx,noise(:,n),[1:size(PC_Result{f}{i}.traces,2)],'linear','extrap'),10000);

            % [Ny_fit t_consts coeffY]  = expfitDM_2(tx(~isnan(noise(1:end-1,n)))',noise(~isnan(noise(1:end-1,n)),n),[1:nTime]',[10^6 10^4]);
            % noise_intp=Ny_fit;
            %[y_fit t_consts coeffY]  = expfitDM_2(tx(~isnan(sp_height(1:end-1,n)))',sp_height(~isnan(sp_height(1:end-1,n)),n),[1:size(Result{i}.traces,2)]',10^7);
            [y_fit t_consts coeffY]  = expfitDM_2(t_fit',sp_height(n,t_fit)',[1:nTime]',[10^7]);
            SpHeight_intp=y_fit;

            figure(99); clf;
            ax1=nexttile([1 1]);
            plot(rescale2([Result.traces(n,1:nTime);tr_mc;tr_mc_bpMonitorMotion;mcTrace(1,1:nTime)],2)')%+[1:4])
            ax2=nexttile([1 1]);
            plot([1:nTime],norm_trace(n,:))
            hold all
            plot(find(sp_time(n,:)),norm_trace(n,find(sp_time(n,:))),'r.')
            plot([1:nTime],y_fit([1:nTime]),'k')
            plot([1 nTime],[0 0],'g')

        end
    end
    Result.SpikeHeight_fit=SpHeight_intp';
    norm_trace=norm_trace./SpHeight_intp';
    norm_trace_check{1}=norm_trace_check{1}./SpHeight_intp';
    norm_trace_check{2}=norm_trace_check{2}./SpHeight_intp';

    %norm_trace=norm_trace;%./(SpHeight_intp./SpHeight_intp(:,1));
    Result.normTraces=norm_trace;%./get_threshold(norm_trace,1);
    Result.norm_trace_check{1}=norm_trace_check{1};%./get_threshold(norm_trace_check{1},1);
    Result.norm_trace_check{2}=norm_trace_check{2};%./get_threshold(norm_trace_check{2},1);
    %Result.spike=find_spike_bh(PC_Result{f}.normTraces-movmedian(PC_Result{f}.normTraces,300,2),5,3);

    % for n=1:size(Result.traces,1)
    % Result.subThreshold(n,:) = filtfilt(b3, a3, Result.normTraces(n,:));
    % Result.theta(n,:) = filtfilt(b4, a4, Result.normTraces(n,:));
    % end
    save(fullfile(fpath{f},'PC_Result.mat'),'Result','fpath','-v7.3')
end


%% Structure registration

Struct_valid=find(1-cell2mat(cellfun(@(x) sum(isnan(x)), StructureData, 'UniformOutput', false)));

for i=unqInd([48])'
    load(fullfile(fpath{i},'OP_Result.mat'))
    StructureStack=mat2gray(double(tiffreadVolume(StructureData{i})));
    StructureStack(StructureStack==0)=median(StructureStack(:));
    StructureStack=StructureStack(:,:,1:180);
    %StructureStack_med=medfilt2_mov(StructureStack,[15 15]);
    illumination_field=imgaussfilt(max(StructureStack,[],3),50);
    StructureStack=StructureStack./illumination_field;
    StructureStack_Gauss=imgaussfilt3(StructureStack,[6 6 0.1]);
    %StructureStack_med(StructureStack_med==0)=median(StructureStack_med(:));
    %StructureStack=(StructureStack-StructureStack_med)./StructureStack_med;
    StructureStack_filt=(StructureStack-StructureStack_Gauss);
    StructureStack_filt=mat2gray(StructureStack_filt);
    StructureStack_bin=[]; level=[];
    level = graythresh(StructureStack_filt);
    StructureStack_bin=StructureStack_filt>level*1.19;
    moviefixsc(StructureStack_bin)

    se = strel('sphere', 1);
    StructureStack_bin = imdilate(StructureStack_bin, se);
    bwSeg=bwlabeln(StructureStack_bin);
    segments = regionprops3(bwSeg,'Volume','EquivDiameter');
    segments = table2array(segments);
    %bwlist=find(arrayfun(@(x) x.Volume>4000, segments) & arrayfun(@(x) x.EquivDiameter>10, segments));
    bwlist = find(segments(:,1)>7000 & segments(:,2)>10);

    se = strel('sphere',1);
    dendrite_bin=double(ismember(bwSeg,bwlist));
    dendrite_bin= imdilate(dendrite_bin,se);
    dendrite_bin= imgaussfilt3(dendrite_bin,2);

    figure(3); clf;
    imshow2(max(dendrite_bin,[],3),[])
    g=1; ROIrmv=[];
    while g
        h = drawpolygon('Color','r');
        if size(h.Position,1)==1 %no more ROI
            g=0;
        else
            ROIrmv=[ROIrmv; {h.Position}];
            hold all
            plot(h.Position(:,1),h.Position(:,2))
        end
    end
    ROIrmvmask=roi2mask(ROIrmv,size(dendrite_bin,1),size(dendrite_bin,2));
    close(figure(3));
    dendrite_bin(repmat(ROIrmvmask,1,1,size(dendrite_bin,3)))=0;
    figure(3); clf;
    imshow2(max(dendrite_bin,[],3),[])

    StructureStack_final = double(StructureStack).* dendrite_bin;
    figure(2); clf;
    imshow2(imfuse(mat2gray(max(StructureStack_final,[],3)),mat2gray(max(StructureStack,[],3))),[])

    rot_ang=90;
    Structure_ref=(imrotate(StructureStack_final,rot_ang));
    ref_img=Result.ref_im; ref_img(ref_img<prctile(ref_img(:),20))=median(ref_img(:)); ref_img=ref_img-prctile(ref_img(:),20);
    [RegImg,tformReg]=imReg(max(Structure_ref,[],3),ref_img);

    Result.Structure=RegImg;
    Result.Structure_bin=imwarp(max(imrotate(dendrite_bin,rot_ang),[],3), tformReg, 'OutputView', imref2d(size(ref_img)));
    Result.tform=tformReg;
    save(fullfile(fpath{i},'OP_Result.mat'),'Result','-v7.3')

    SameCellInd=find(Mouse==Mouse(i) & NeuronInd==NeuronInd(i));
    
    for j=SameCellInd'
        load(fullfile(fpath{j},'OP_Result.mat'))
        Result.tform=tformReg;
        [offsetY offsetX]=calculate_shift(ref_img,Result.ref_im)
        Result.tform.T(3,1:2)=Result.tform.T(3,1:2)-[offsetX offsetY];
        Result.Structure=imwarp(max(Structure_ref,[],3), Result.tform, 'OutputView', imref2d(size(ref_img)));
        Result.Structure_bin=imwarp(max(imrotate(dendrite_bin,rot_ang),[],3), Result.tform, 'OutputView', imref2d(size(ref_img)));
        figure; imshow2(imfuse(Result.ref_im,Result.Structure),[]);
        save([fpath{j} '/OP_Result.mat'],'Result','-v7.3')
    end
end
%%
for i=[160]%length(fpath)]
    load([fpath{i} '/Result.mat'])
    Result.normTraces=Result.traces-prctile(Result.traces,25,2);
    Result.normTraces=Result.normTraces./get_threshold(Result.normTraces,1);
    tr_ref=Result.normTraces(1,:);
    tr=Result.normTraces;
    nROI=size(Result.normTraces,1);
    sp=find_spike_bh(tr-movmedian(tr,50,2),5,4);
    sp_soma=find_spike_bh(tr_ref-movprc(tr_ref,200,30,2),6,3);

    tr_sub=mean(tr_ref,1)-movprc(mean(tr_ref,1),200,20,2);
    tr_sub=get_subthreshold(tr_sub,sp_soma,5,10);

    [trans tr_trace]=detect_transient2(tr_sub,[6 1.2],sp_soma,20);
    transcand=cell2mat(cellfun(@(x) length(x)>2,trans.ISI,'UniformOutput',false));
    meanISI_frnt=cellfun(@(x) mean(x(1:2)),trans.ISI(transcand));
    meanISI_first3=zeros(1,length(trans.length));
    meanISI_first3(transcand)=meanISI_frnt;

    %CS_ind=find(trans.spike_number>2 & trans.mean_ISI<15);
    CS_ind=find(trans.spike_number>2 & meanISI_first3<20);
    CS_trace=ismember(tr_trace,CS_ind);
    CS_spike=sp_soma.*bwlabel(CS_trace);
    [~, CS_spike_time]=unique(CS_spike);

    sp_total=max([sp_soma; sp(2:end,:)],[],1);
    bAP_ind=zeros(1,size(tr,2));
    bAP_ind(unique(find(sp_soma)'+[0:3]))=1;

    SpikeClassMat=zeros(3,size(tr,2));
    SpikeClassMat(1,:)=sp_soma.*(1-CS_trace); %bAPs
    SpikeClassMat(2,CS_spike_time(2:end))=1; %Complex spikes
    SpikeClassMat(3,:)=sp_total.*(1-bAP_ind); %dSpikes

    Result.spike=[sp_soma; sp(2:end,:)];
    Result.SpClass=SpikeClassMat;
    Result.CStrace=CS_trace;
    show_traces_spikes(Result.normTraces,Result.spike,[Result.SpClass; Result.Blue]);
    save([fpath{i} '/Result.mat'],'Result','-v7.3')
end

%%
load([fpath{89} '/Result.mat']); nROI=size(Result.normTraces,1);
coord_1d=dim_reduce(get_coord(Result.ftprnt));
[~, dist_order]=sort(coord_1d,'ascend');
nTau_bAP=[-20:20];
bAP_ref=find(Result.SpClass(1,:));
prc_normTr=Result.normTraces;
STA_SS=squeeze(mean(reshape(prc_normTr(:,bAP_ref'+nTau_bAP),nROI,[],length(nTau_bAP)),2));
STA_SS=STA_SS-prctile(STA_SS,25,2);
F_ref=mean(STA_SS(:,-nTau_bAP(1)+[10:14]),2);
prc_normTr=prc_normTr./F_ref;
show_traces_spikes(prc_normTr(dist_order,:),Result.spike(dist_order,:),[Result.SpClass; Result.Blue]);
%%
Blueon=find((Result.Blue(2:end)-Result.Blue(1:end-1))>0);
STA_Blue=reshape(prc_normTr(:,Blueon'+[-50:150]),nROI,[],length([-50:150]));
CS_on=find(Result.SpClass(2,:));
STA_CS=reshape(prc_normTr(:,CS_on'+[-100:150]),nROI,[],length([-100:150]));
figure(10); clf;
imagesc(flipud(squeeze(mean(STA_CS(dist_order,:,:,[0 5])
colormap(turbo)

%% STAs
nTau=[-10:20];

for i=[133]%1:length(OP_Result)
Result.normTrace=Result.traces./get_threshold(Result.traces,1);
Result.spike=find_spike_bh(Result.normTrace-movmedian(Result.normTrace,300,2),5,3);

Blue=Result.Blue;
blueOff = Blue == 0;
blueOff2 = imerode(blueOff, [ones(1,20), zeros(1, 20)]);
Blue_di=~blueOff2;
bwBlue_di=bwlabel(Blue_di);

ref_ROI=find(sum(Result.spike(1:5,:).*bwBlue_di,2)==max(sum(Result.spike(1:5,:).*bwBlue_di,2)),1);
nROI=size(Result.normTrace,1);
tr=Result.normTrace(ref_ROI,:); t=[1:length(tr)];
spike=Result.spike(ref_ROI,:);   


sp_pulse=[];
for b=1:max(bwBlue_di)
t_tmp=find(bwBlue_di==b);
if max(bwBlue_di)<500
sp_pulse=[sp_pulse find(spike(t_tmp))+t_tmp(1)-1];
else
sp_pulse=[sp_pulse find(spike(t_tmp),1,'first')+t_tmp(1)-1];
end
end

sTau=sp_pulse(1:end-1)'+nTau;
spikeMat=reshape(Result.normTrace(:,sTau),nROI,length(sp_pulse)-1,[]);
STAtrace=squeeze(mean(spikeMat,2));
STAtrace=STAtrace-prctile(STAtrace, 5, 2);
Result.spikeMat = spikeMat;
Result.STAtrace = STAtrace;

if max(bwBlue_di)>15

load(fullfile(fpath{i},"output_data.mat"))
load([fpath{i} '/mcTrace' num2str(1,'%02d') '.mat']);
switch char(CamType(i))
    case 'flash'
sz=double(Device_Data{1, 4}.ROI([2 4]));
    case 'fusion'
sz=double(Device_Data{1, 3}.ROI([2 4]));        
end

CamCounter=Device_Data{1, 2}.Counter_Inputs(1, 1).data;
CamTrigger=find(CamCounter(2:end)-CamCounter(1:end-1));

mov_mc=double(readBinMov_times([fpath{i} '/mc_ShutterReg' num2str(1,'%02d') '.bin'],sz(2),sz(1),[1:length(CamTrigger)]));

mov_res= mov_mc-mean(mov_mc,3);
bkg = zeros(1, size(mov_mc,3));

bkg(1,:)=movmedian(get_blueoffTrace(squeeze(mean(mov_mc,[1 2])),Blue,30),3000,'omitnan');
mov_res = SeeResiduals(mov_res,Result.mc);
mov_res = SeeResiduals(mov_res,Result.mc.^2);
mov_res = SeeResiduals(mov_res,Result.mc(:,1).*Result.mc(:,end));
mov_res= SeeResiduals(mov_res,bkg,1);

STA_tmp=reshape(mov_res(:,:,sTau),sz(2),sz(1),[],length(nTau));
%STA_tmp=STA_tmp-mean(STA_tmp(:,:,:,[1:3]),4);
Result.STAmovie=squeeze(mean(STA_tmp,3));
end
end

%% Generate SNAPT movie
for i=[133]
    load([fpath{i} '/Result.mat'])
mask=max(Result.Structure_bin,[],3)>0.01;
maskSTA=max(-Result.STAmovie,[],3)./Result.ref_im>0.05;
StrImg=max(Result.Structure,[],3);
STAmovie=mat2gray(-Result.STAmovie);
STAmovie=STAmovie-prctile(STAmovie,10,3);
STAmovie=mat2gray(STAmovie(:,:,6:22));
tformReg=Result.tform;
[Result.SNAPT Result.dtimg]=generate_SNAPTmov(mat2gray(STAmovie),mask,StrImg,tformReg);
[yR xR zR]=size(Result.Structure);
bluePatt = bwboundaries(imwarp(Result.BlueDMDimg,Result.tform,'OutputView', imref2d([yR xR])));

figure(20); clf;
%v = VideoWriter([fpath{i} '/SNAPT_movie'],'MPEG-4');
v = VideoWriter([fpath{i} '/SNAPT_movie'],'Uncompressed AVI');

open(v);
subframeT = 0.025; % ms
initialT = -2; % ms
finalT = 2; % ms
times = initialT:subframeT:finalT;

for j = 1:length(times)
    clf;
    %set(gca,'units','pixels','position',[200 0 1000 800])
    imshow(Result.SNAPT(:,:,:,j),[])
    pbaspect([size(double(Result.SNAPT(:,:,:,j)),2) size(double(Result.SNAPT(:,:,:,j)),1) 1]),colormap(gray)
    hold all
    plot(bluePatt{1}(:,2),bluePatt{1}(:,1),'color',[0 0.6 1],'linewidth',2)
    axis off
    text(2,20,[num2str(times(j)+0.9) 'ms'], 'FontSize', 20, 'color', [0.99 0.99 0.99])% the value 1. is to adjust timing by eyes       
    pause(0.1)
    set(gcf,'color','w')    % Sets background to white
    frame = getframe(gcf);
    writeVideo(v,frame);
    pause(0.1);
end;
close(v);
end
%save([save_to 'OP_Result_20240212'],"OP_Result",'fpath','-v7.3')
