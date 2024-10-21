clear
clc;
cd '/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/FromBackup/PP72_PlaceCellResults';
[~, ~, raw] = xlsread(['/Volumes/cohen_lab/Lab/Labmembers' ...
    '/Byung Hun Lee/Data/PrismPCdata_Arrangement.xlsx'], 'Sheet1', 'C5:R31');

% [~, ~, NeuronsToUse]=xlsread(['/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/' ...
%     'PlaceCellData_Arrangement.xlsx'], 'Sheet1', 'L8:M46');
%
% NeuronsToUse=cellfun(@(x) (str2num(num2str(x))),NeuronsToUse,'UniformOutput',false);
ref_ROI=cellfun(@(x) (str2num(num2str(x))),raw(:,10),'UniformOutput',false);
fpath=raw(:,1)';
StructureData=raw(:,9);
EndFrame=cell2mat(raw(:,13));
ifmotionReject=cell2mat(raw(:,14));
ifdirtRemoval=cell2mat(raw(:,16));
save_figto='/Volumes/BHL_WD18TB/PP72_PlaceCellResults';
place_bin=150; time_segment=15000; overlap=200;
alignedMovFN = {'STA_Mat_SS','STA_Mat_CS','STA_Mat_dSP'};

%% Motion correction

for f=[27]% 15 4 5 7]%length(fpath)
    f
    load(fullfile(fpath{f},"output_data.mat"))
    sz=double(Device_Data{1, 3}.ROI([2 4]));
    ref_time=[119000:120000];

    %frm_end=max(Device_Data{1, 2}.Counter_Inputs(1, 1).data);
    frm_end=EndFrame(f);
    f_seg=[[1:time_segment:frm_end] frm_end+1]; f_seg(2:end)=f_seg(2:end)-1;
    DAQ_rate=Device_Data{1, 2}.buffered_tasks(1, 1).rate;
    CamDAQ_rate=Device_Data{1, 2}.Counter_Inputs.rate;
    CamTrig=Device_Data{1, 2}.Counter_Inputs.data;
    CamTrig2=find(CamTrig(2:end)-CamTrig(1:end-1)>0);
    Frm_rate=(CamTrig2(2)-CamTrig2(1))/CamDAQ_rate;

    if length(ref_time)>2000
        mov_test=double(readBinMov_times([fpath{f} '/frames1.bin'],sz(2),sz(1),[ref_time(1):ref_time(1)+2000]));
    else
        mov_test=double(readBinMov_times([fpath{f} '/frames1.bin'],sz(2),sz(1),[ref_time(1):ref_time(end)]));
    end
    mov_test=rollingShutter_correction(mov_test,1/Frm_rate,'fusion');
    mov_test=mov_test(:,:,2:end);
    [mov_test_mc,xyField]=optical_flow_motion_correction_LBH(mov_test,mean(mov_test,3),'normcorre');
    mov_test=vm(mov_test);
    mov_test = single(mov_test)./single(max(mov_test.data(:)));
    mov_test = movmean(mov_test,10,3);
    mov_ref = squeeze(median(mov_test,3));

    for j=1:length(f_seg)-1
        try
            mov=double(readBinMov_times([fpath{f} '/frames1.bin'],sz(2),sz(1),[f_seg(j):f_seg(j+1)+overlap]));
        catch % when the image ends
            mov=double(readBinMov_times([fpath{f} '/frames1.bin'],sz(2),sz(1),[f_seg(j):f_seg(j+1)]));
        end

        mov=rollingShutter_correction(mov,1/Frm_rate,'fusion');
        mov=vm(mov(:,:,2:end));
        if j==1
            mov=mov(:,:,[1 1:size(mov,3)]);
        end

        [mov_mc,xyField]=optical_flow_motion_correction_LBH(mov,mov_ref,'normcorre');

        ave_im=mean(mov_mc,3);
        mov_mc=vm(mov_mc);
        mov_mc.transpose.savebin([fpath{f} '/mc_ShutterReg' num2str(j,'%02d') '.bin'])

        %        mcTrace = squeeze(mean(xyField,[1 2])); %optic flow
        mcTrace=xyField; % Normcorre
        save([fpath{f} '/mcTrace' num2str(j,'%02d') '.mat'],'mcTrace','ave_im')

        %  clear mov_mc mov
    end
end

%% Get footprint

for f=[18]%:length(fpath)
    disp(fpath{f}); Result=[];
    DAQ_rate=0.000005;
    load([fpath{f} '/output_data.mat'])
    sz=double(Device_Data{1, 3}.ROI([2 4]));
    ref_time=[2000:3000];
    mov_test=double(readBinMov_times([fpath{f} '/mc_ShutterReg' num2str(10,'%02d') '.bin'],sz(2),sz(1),ref_time));
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

    [u,s,v] = svds(tovec(mov_test-mean(mov_test,3)),20);
reshape_u=reshape(u,sz(2),sz(1),[]);
bvMask=[];
[~, Result.bvMask]=get_ROI(max(abs(reshape_u),[],3),bvMask);

    load(fullfile(fpath{f},'mcTrace15.mat'));
    frm_end=max(Device_Data{1, 2}.Counter_Inputs(1, 1).data);
    Result.Blue=Device_Data{1, 2}.buffered_tasks(1, 1).channels(1, 2).data;
    Result.Reward=Device_Data{1, 2}.buffered_tasks(1, 3).channels.data;
    CamCounter=Device_Data{1, 2}.Counter_Inputs(1, 1).data;
    CamTrigger=find(CamCounter(2:end)-CamCounter(1:end-1));

    Result.Blue=Result.Blue(CamTrigger); Result.Reward=Result.Reward(CamTrigger);
    frm_rate=double((CamTrigger(2)-CamTrigger(1))*DAQ_rate);
    mov_mc=double(readBinMov([fpath{f} '/mc_ShutterReg' num2str(15,'%02d') '.bin'],sz(2),sz(1)));
    Result.ref_im=mean(mov_mc,3);
    [clickyROI, original_trace]=clicky(mov_mc);

    mov_res= mov_mc-mean(mov_mc,3);
    mcTrace.xymean=movmean(mcTrace.xymean,3,2);
    mov_res = SeeResiduals(mov_res,mcTrace.xymean);
    mov_res = SeeResiduals(mov_res,mcTrace.xymean.^2);
    mov_res = SeeResiduals(mov_res,mcTrace.xymean(:,1).*mcTrace.xymean(:,3));

    [SeeRes_trace]=apply_clicky(clickyROI,mov_res);
    figure(3); clf;
    plot(rescale(original_trace)); hold all;
    plot(rescale(SeeRes_trace)+1);

    n_comp=6;
    mov_filt=imgaussfilt3(mov_res.*double(max(Result.bvMask,[],3)==0),[2 2 0.1]);
    mov_filt=mov_filt(:,:,4000:8000);
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

        %[icsTrace, ~, sepmat]=sorted_ica(eigTrace(:,1:n_comp),n_comp);
        %plot(rescale2(icsTrace,1)+[1:size(icsTrace,2)])
        %V_ics=V(:,1:n_comp)*sepmat';
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
close(figure(4));
    Result.ftprnt=toimg(ftprnt,sz(2),sz(1));
    figure; clf;
    show_footprnt(Result.ftprnt,Result.ref_im)
    save(fullfile(fpath{f},'PC_Result.mat'),'Result','fpath','-v7.3')
end

%% Dirt Removal
if ifdirtRemoval(f)

load(fullfile(fpath{20},'disk_ref'),'z2');

for f=[27]%length(fpath)
    load(fullfile(fpath{f},'PC_Result.mat'),'Result');
    load([fpath{f} '/output_data.mat'])
    sz=double(Device_Data{1, 3}.ROI([2 4])); blueDMDcontour=[];
    CamCounter=Device_Data{1, 2}.Counter_Inputs(1, 1).data;
    CamTrigger=find(CamCounter(2:end)-CamCounter(1:end-1));

    bound=10;
    frm_end=EndFrame(f);
    f_seg=[[1:time_segment:frm_end] frm_end+1]; f_seg(2:end)=f_seg(2:end)-1;
    take_window=repmat([1 time_segment],length(f_seg)-1,1);
    take_window(2:end,1)=take_window(2:end,1)+overlap; take_window(1:end-1,2)=take_window(1:end-1,2)+overlap;
    take_window(end)=mod(f_seg(end),time_segment);
    take_window(take_window==0)=time_segment;
    Result.dirtTrace=[];

    for j=1:length(f_seg)-1
        j
        mov_mc=double(readBinMov([fpath{f} '/mc_ShutterReg' num2str(j,'%02d') '.bin'],sz(2),sz(1)));
        load([fpath{f} '/mcTrace' num2str(j,'%02d') '.mat']);

        mov_mc=mov_mc(:,:,[take_window(j,1):take_window(j,2)]);
        mc=mcTrace.xymean([take_window(j,1):take_window(j,2)],:);
        Blue_on_Seg=unique(ceil(find(Result.Blue)/time_segment));

        if ismember(j,Blue_on_Seg)
            isfirst_seg=double(j==1);
            [~, blueomitTr]=get_blueoffTrace(squeeze(mean(mov_mc,[1 2])), ...
                Result.Blue(f_seg(j)+take_window(j,1)-isfirst_seg:f_seg(j)+take_window(j,2)-isfirst_seg),30);
            bkg = zeros(1, size(mov_mc,3));
            bkg(1,:)=movmedian(blueomitTr,4000,'omitnan');
        else
            bkg = zeros(2, size(mov_mc,3));
            bkg(1,:) = linspace(-1, 1, size(mov_mc,3));  % linear term
            bkg(2,:) = linspace(-1, 1, size(mov_mc,3)).^2;  % quadratic term
        end

        mov_res=mov_mc;
        mov_res= mov_mc-median(mov_mc,3);
        mov_res = SeeResiduals(mov_res,mc);
        mov_res = SeeResiduals(mov_res,mc.^2);
        mov_res = SeeResiduals(mov_res,mc(:,1).*mc(:,end));
        mov_res= SeeResiduals(mov_res,bkg,1);

        mov_int=imresize(mov_res(bound:end-bound,bound:end-bound,:),1/2);
        mov_int_filt=imgaussfilt3(mov_int,[8 8 5])-imgaussfilt3(mov_int,[3 3 5]);

        match_img=[]; match_data=[];
        for k=1:size(mov_int_filt,3)
            [ match_data{k}, match_img(:,:,k)] = matchPattern(mov_int_filt(:,:,k), z2, 0.3, 2);
        end
        match_data=cellfun(@(x) 2*(x+[21 21 0]),match_data,'UniformOutput',false);
        detectedPtl=trackParticle(match_data,5,20);
        traveldist=cell2mat(cellfun(@(x) sum(sqrt(sum((x(end,1:2)-x(1,1:2)).^2,2)),1),detectedPtl,'UniformOutput',false));
        traveltime=cell2mat(cellfun(@(x) size(x,1),detectedPtl,'UniformOutput',false));
        valid_ptl=traveldist>50 & traveltime>70;
        filteredPtl=detectedPtl(valid_ptl); 
        dirtMov=zeros(size(mov_res));

        figure(2); clf; filteredPtl_exp=[];
        imshow2(Result.ref_im(bound:end-bound,bound:end-bound),[]); hold all
        for p=1:length(filteredPtl)
            t_exp=[round(filteredPtl{p}(1,3) - range(filteredPtl{p}(:,3))*0.4):1:round(filteredPtl{p}(end,3) + range(filteredPtl{p}(:,3))*0.4)];

            p1 = polyfit(filteredPtl{p}(:,3), filteredPtl{p}(:,1), 1);
            p2 = polyfit(filteredPtl{p}(:,3), filteredPtl{p}(:,2), 1);

            filteredPtl_exp{p}(:,1)=round(p1(1)*t_exp+p1(2));
            filteredPtl_exp{p}(:,2)=round(p2(1)*t_exp+p2(2));
            filteredPtl_exp{p}(:,3)=t_exp;
            outofrngInd=find(filteredPtl_exp{p}(:,1)+bound-1<1 | filteredPtl_exp{p}(:,1)+bound-1>sz(2) | filteredPtl_exp{p}(:,2)+bound-1<1 | filteredPtl_exp{p}(:,2)+bound-1>sz(1) | t_exp' < 1 | t_exp' > size(mov_res,3));
            filteredPtl_exp{p}(outofrngInd,:)=[];
            plot(filteredPtl_exp{p}(:,2)+bound-1,filteredPtl_exp{p}(:,1)+bound-1,'r'); 
            plot(filteredPtl{p}(:,2)+bound-1,filteredPtl{p}(:,1)+bound-1,'wo'); 
            for dd=1:size(filteredPtl_exp{p},1)
            dirtMov(filteredPtl_exp{p}(dd,1)+bound-1,filteredPtl_exp{p}(dd,2)+bound-1,filteredPtl_exp{p}(dd,3))=1;
            end
        end
        drawnow
        se = strel('disk',30);
        dirtMov_dilate = imdilate(dirtMov, se);
        Result.dirtTrace=[Result.dirtTrace (tovec(dirtMov_dilate)'*tovec(Result.ftprnt))'];
    end
    save(fullfile(fpath{f},'PC_Result.mat'),'Result','fpath','-v7.3')
end
end
%% Load Virmen data
    
for f=[18]%:length(fpath)
    load([fpath{f} '/output_data.mat'])
    load(fullfile(fpath{f},'PC_Result.mat'),'Result')
    fileList = dir(fullfile(fpath{f}, '*.data'));
    if length(fileList)==1
        fid = fopen(fullfile(fpath{f},fileList.name));
    else
        error('Data file cannot be found');
    end
    VRdata = fread(fid,[12 inf],'double');

    WorldITrack=find(VRdata(2,:)==1); %world 1
    VRdata(5,WorldITrack)=(VRdata(5,WorldITrack)+6)*115/121;

    VRdata=VRdata(:,VRdata(10,:)>0);
    DAQ_rate=Device_Data{1, 2}.Counter_Inputs.rate;
    CamCounter=Device_Data{1, 2}.Counter_Inputs(1, 1).data;
    CamTrigger=find(CamCounter(2:end)-CamCounter(1:end-1));

    Result.Blue=Device_Data{1, 2}.buffered_tasks(1, 1).channels(1, 2).data;
    Result.Reward=Device_Data{1, 2}.buffered_tasks(1, 3).channels.data;
    Result.Blue=Result.Blue(CamTrigger); Result.Reward=Result.Reward(CamTrigger);


    t_DAQ=CamTrigger/DAQ_rate;
    t_VR = datetime(datetime(VRdata(1,:),'ConvertFrom','datenum'), 'InputFormat', 'yyyy-MM-dd HH:mm:ss.SSS');
    t_VR= t_VR-t_VR(1);
    t_VR= milliseconds(t_VR)/1000;
    t_VR= t_VR*t_DAQ(end)/t_VR(end);
    VRdata(1,:)=t_VR;
    [Virmen_data_int vel_trace]=virmen_interpolate(VRdata,115,t_DAQ);
    Virmen_data_int(end+1,:)=vel_trace;

    Result.VR=Virmen_data_int;
    %Result.Blue=Result.Blue(:,1:EndFrame(f)).*Result.VR(12,1:EndFrame(f));
    Result.Blue=Result.Blue(:,1:EndFrame(f)).* ...
    ismember(Result.VR(8,1:EndFrame(f)),[12:21]); % From Lap 12 to 21
 
    Result.VR=Result.VR(:,1:EndFrame(f));

    save(fullfile(fpath{f},'PC_Result.mat'),'Result','fpath','-v7.3')
end

%% Signal extraction

for f=[18]%length(fpath)
    load(fullfile(fpath{f},'PC_Result.mat'),'Result');
    load([fpath{f} '/output_data.mat'])
    sz=double(Device_Data{1, 3}.ROI([2 4])); blueDMDcontour=[];
    CamCounter=Device_Data{1, 2}.Counter_Inputs(1, 1).data;
    CamTrigger=find(CamCounter(2:end)-CamCounter(1:end-1));
    Rfixed = imref2d(repmat(Device_Data{1, 3}.virtualSensorSize,1,2));
    inverseTform = invert(Device_Data{1, 6}.tform);
    try
        revertedImage = imwarp(double(Device_Data{6}.pattern_stack(:,:,find(sum(Device_Data{6}.pattern_stack,[1 2])>0))), inverseTform,'OutputView',Rfixed);
    catch
        revertedImage = imwarp(double(Device_Data{6}.Target), inverseTform,'OutputView',Rfixed);
    end
    [blueDMDimg]=imcrop_3d(revertedImage,double(Device_Data{1, 3}.ROI([1 3 2 4]))+[0 0 -1 -1]);
    for d=1:size(blueDMDimg,3)
        blueDMDcontour{d}=bwboundaries(blueDMDimg(:,:,d));
    end
 
    oddmask=checkerboard(1, size(Result.ftprnt,1)/2, size(Result.ftprnt,2)/2) > 0.5;
    oddftprnt=Result.ftprnt.*oddmask;
    evenftprnt=Result.ftprnt.*(1-oddmask);

    Result.blueDMDimg=blueDMDimg;
    Result.blueDMDcontour=blueDMDcontour;
    Result.traces=[]; Result.traces_bvMask=[];
    Result.traces_checker{1}=[]; Result.traces_checker{2}=[];
    Result.tracesPCA=[];
    Result.mcTrace=[];
    Result.im_corr=[];
    bound=5;
    ref_im_vec=tovec(Result.ref_im(bound:end-bound,bound:end-bound));
    ref_im_vec=tovec(Result.ref_im(bound:end-bound,bound:end-bound))/std(ref_im_vec,0,1);    
    %frm_end=max(Device_Data{1, 2}.Counter_Inputs(1, 1).data);
    frm_end=EndFrame(f);
    f_seg=[[1:time_segment:frm_end] frm_end+1]; f_seg(2:end)=f_seg(2:end)-1;
    take_window=repmat([1 time_segment],length(f_seg)-1,1);
    take_window(2:end,1)=take_window(2:end,1)+overlap; take_window(1:end-1,2)=take_window(1:end-1,2)+overlap;
    take_window(end)=mod(f_seg(end),time_segment);
    take_window(take_window==0)=time_segment;
    Blue_on_Seg=unique(ceil(find(Result.Blue)/time_segment));

    for j=1:length(f_seg)-1
        j
        mov_mc=double(readBinMov([fpath{f} '/mc_ShutterReg' num2str(j,'%02d') '.bin'],sz(2),sz(1)));
        load([fpath{f} '/mcTrace' num2str(j,'%02d') '.mat']);

        mov_mc=mov_mc(:,:,[take_window(j,1):take_window(j,2)]);
        %mc= movmean(mcTrace.xymean-movmedian(mcTrace.xymean,500,1),3,1);
        mc=mcTrace.xymean([take_window(j,1):take_window(j,2)],:);
        
        mov_mc_filt=imgaussfilt3(mov_mc,[2 2 0.1]);
        mov_mc_vec=tovec(mov_mc_filt(bound:end-bound,bound:end-bound,:));

        if ismember(j,Blue_on_Seg)
           isfirst_seg=double(j==1);
        [~, blueomitTr]=get_blueoffTrace(squeeze(mean(mov_mc,[1 2])), ...
            Result.Blue(f_seg(j)+take_window(j,1)-isfirst_seg:f_seg(j)+take_window(j,2)-isfirst_seg),30,30);
        bkg = zeros(1, size(mov_mc,3));    
        bkg(1,:)=movmedian(blueomitTr,4000,'omitnan');
        else
        bkg = zeros(2, size(mov_mc,3));    
        bkg(1,:) = linspace(-1, 1, size(mov_mc,3));  % linear term
        bkg(2,:) = linspace(-1, 1, size(mov_mc,3)).^2;  % quadratic term
        end
        
        mov_res=mov_mc;
        % mov_res= mov_mc-median(mov_mc,3);
        % mov_res = SeeResiduals(mov_res,mc);
        % mov_res = SeeResiduals(mov_res,mc.^2);
        % mov_res = SeeResiduals(mov_res,mc(:,1).*mc(:,end));
        % mov_res= SeeResiduals(mov_res,bkg,1);

        Result.traces=[Result.traces -(tovec(mov_res)'*tovec(Result.ftprnt))'];
        Result.traces_bvMask=[Result.traces_bvMask -(tovec(mov_res.*double(max(Result.bvMask,[],3)==0))'*tovec(Result.ftprnt))'];
        Result.traces_checker{1}=[Result.traces_checker{1} -(tovec(mov_res.*double(max(Result.bvMask,[],3)==0))'*tovec(oddftprnt))'];
        Result.traces_checker{2}=[Result.traces_checker{2} -(tovec(mov_res.*double(max(Result.bvMask,[],3)==0))'*tovec(evenftprnt))'];
%       Result.tracesPCA=[Result.tracesPCA -(tovec(mov_res)'*tovec(Result.pcaMask))'];
        Result.mcTrace=[Result.mcTrace; mc];
        Result.im_corr=[Result.im_corr corr(rescale2(mov_mc_vec,1),ref_im_vec,'type','Spearman')'];  %image correlation

    end

    %Result.traces=Result.traces(:,1:length(CamTrigger)-1);
    %Result.mcTrace=Result.mcTrace(1:length(CamTrigger)-1,:);
    Result.traces=Result.traces(:,1:frm_end);
    Result.mcTrace=Result.mcTrace(1:frm_end,:);
    Result.Blue=Result.Blue(:,1:frm_end);
    Result.Reward=Result.Reward(:,1:frm_end);

    save(fullfile(fpath{f},'PC_Result.mat'),'Result','fpath','-v7.3')
end

%% Clean up and norm
exclude_frq=[241.7 242]; %monitor
%exclude_frq2=[483.5 484]; %monitor
exclude_frq2=[55.5 56]; %motion
time_bin=15000; Fs=1000; ref_trace=[2 4 4 2 3 2 1 1 1 3 1 1 1 1 1 1 1 1 1 3 23 3 23]; %2nd trunk is the reliable trace

for f=[18]%:length(fpath)
    load(fullfile(fpath{f},'PC_Result.mat'),'Result')
    ref_trace=ref_ROI{f}(1);
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
lwpass_fit(:,t_silent)=Result.traces(:,t_silent); lwpass_fit=movmedian(lwpass_fit,20000,2,'omitnan');
lwpass_fit_ch{1}(:,t_silent)=Result.traces_checker{1}(:,t_silent); lwpass_fit_ch{1}=movmedian(lwpass_fit_ch{1},20000,2,'omitnan');
lwpass_fit_ch{2}(:,t_silent)=Result.traces_checker{2}(:,t_silent); lwpass_fit_ch{2}=movmedian(lwpass_fit_ch{2},20000,2,'omitnan');

tr_res=squeeze(SeeResiduals(permute(Result.traces_bvMask,[1 3 2]),lwpass_fit));
tr_res_checker{1}=squeeze(SeeResiduals(permute(Result.traces_checker{1},[1 3 2]),lwpass_fit));
tr_res_checker{2}=squeeze(SeeResiduals(permute(Result.traces_checker{2},[1 3 2]),lwpass_fit));
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

%% Classify the spikes
motion_frq=[30 120];
ifmotion_reject=ifmotionReject;
nSpikeThres=1; % set 
for f=[26]
    load(fullfile(fpath{f},'PC_Result.mat'))
    %tr=PC_Result{i}.normTraces(ref_ROI{i},:); %somatic spike
    %sp_ref=max(find_spike_bh(tr-movmedian(tr,100,2),5,3),[],1);
    tr=Result.normTraces./get_threshold(Result.normTraces,1);
    tr_ref=Result.normTraces(ref_ROI{f},:)./get_threshold(Result.normTraces(ref_ROI{f},:),1);
    tr_ref_mean=mean(tr_ref,1);

    nROI=size(tr,1);
    DOI{f}=setdiff([2:nROI],ref_ROI{f});
    sp=find_spike_bh(tr-movmedian(tr,50,2),6,4);
    sp_ref=find_spike_bh(tr_ref-movprc(tr_ref,200,30,2),5,3); % Set somatic spike threshold
    [wvletTr wvletF] = cwt(Result.im_corr,1000); % Compute the CWT
    motionArtTrace=sum(abs(wvletTr(find(wvletF>(motion_frq(1)) & wvletF<motion_frq(2)),:)),1);
    motionArtTrace=(motionArtTrace-prctile(motionArtTrace,20))./std(motionArtTrace);
    motionReject=zeros(1,size(tr,2));
    if ifmotionReject(f)
        motionReject= motionArtTrace>6;
        motionReject = imdilate(motionReject, strel('square', 200));
        sp(:,motionReject)=0; sp_ref(:,motionReject)=0;
    end
    if ifdirtRemoval(f)
        sp=sp.*(Result.dirtTrace==0);
    end
    Result.motionReject=motionReject;
    sp_ref=sum([sp_ref; sp(ref_ROI{f},:)],1);

    [~, shift]=max(reshape(tr_ref_mean(find(sp_ref>nSpikeThres)+[-1:0]'),2,[]),[],1);
    shift=shift-2;
    sp_time_Soma = find(sp_ref>nSpikeThres)+shift;
    sp_soma=zeros(1,size(tr,2));
    sp_soma(sp_time_Soma)=1;
    sp_soma=[0 (sp_soma(2:end)-sp_soma(1:end-1))==1]; %remove consecutive spikes

    tr_sub=tr_ref_mean-movprc(tr_ref_mean,200,20,2);
    tr_sub=get_subthreshold(tr_sub,sp_soma,5,10);

    [trans tr_trace]=detect_transient2(tr_sub,[5 1.5],sp_soma,15); % Set Complex spike threshold
    transcand=cell2mat(cellfun(@(x) length(x)>1,trans.ISI,'UniformOutput',false));
    meanISI_frnt=cellfun(@(x) mean(x(1:2)),trans.ISI(transcand));
    meanISI_first3=NaN(1,length(trans.length));
    meanISI_first3(transcand)=meanISI_frnt;

    %CS_ind=find(trans.spike_number>2 & trans.mean_ISI<15);
    CS_ind=find(trans.spike_number>2 & meanISI_first3<18);
    CS_trace=ismember(tr_trace,CS_ind);
    CS_spike=sp_soma.*bwlabel(CS_trace);
    [~, CS_spike_time]=unique(CS_spike);

    sp_total=sum([sp_soma; sp(DOI{f},:)],1);
    bAP_ind=zeros(1,size(tr,2));
    bAP_ind(unique(find(sp_soma)'+[-1:3]))=1;

    SpikeClassMat=zeros(3,size(tr,2));
    SpikeClassMat(1,:)=sp_soma.*(1-CS_trace); %bAPs
    SpikeClassMat(2,CS_spike_time(2:end))=1; %Complex spikes
    SpikeClassMat(3,:)=sp_total.*(1-bAP_ind).*(1-CS_trace); %dSpikes
    SpikeClassMat(3,:)=SpikeClassMat(3,:)>1;
    SpikeClassMat(3,:)=[0 (SpikeClassMat(3,2:end)-SpikeClassMat(3,1:end-1))==1]; %remove consecutive spikes

    Result.spike=[sp_soma; sp(2:end,:)];
    Result.spike(ref_ROI{f},:)=repmat(sp_soma,length(ref_ROI{f}),1);
    Result.SpClass=SpikeClassMat;
    Result.CStrace=CS_trace;

nROI=size(Result.ftprnt,3);
SkelDend = Skeletonize_dendrite(Result.ref_im,8,0.005,25);
interDendDist=[];
for i=1%:nROI
    for j=1:nROI
        [interDendDist(i,j), ~]=geodesic_distance(SkelDend,get_coord(Result.ftprnt(:,:,i)),get_coord(Result.ftprnt(:,:,j)));
    end
end

coord_1d=dim_reduce(get_coord(Result.ftprnt));
[~, Result.dist_order]=sort(coord_1d,'descend');
Result.geodist=interDendDist(1,:)'.*sign(coord_1d-coord_1d(1));
som_roi=find(Result.dist_order==1);
%show_traces_spikes(Result.normTraces(Result.dist_order,:),Result.spike(Result.dist_order,:),[Result.SpClass; double(Result.motionReject)]);
show_traces_spikes(Result.normTraces./get_threshold(Result.normTraces,1),Result.spike,[Result.SpClass; double(Result.motionReject)]);
save(fullfile(fpath{f},'PC_Result.mat'),'Result','fpath','-v7.3')
end

%% Save the STA movies
nTau={[-40:20],[-50:50],[-30:20]}; %SS, CS, dSP
bound=6;
%f_tmp='/Volumes/BHL18TB_D1/20240218/134705BHLm117_FOV2_VR2';
for f=[26]
    load(fullfile(fpath{f},'PC_Result.mat'),'Result')
    StackedSpike=[];
    alignmovlist=dir(fullfile(fpath{f},'STA_Mat*.tiff'));
    alignmov_ind=ones(1,3);
    for l=1:length(alignmovlist)
        delete(fullfile(fpath{f},alignmovlist(l).name));
    end

    for c=1:3 % Get frames of spikes of each classes
        fileInd{c}=ceil(find(Result.SpClass(c,:))/time_segment);
        frameInd{c}=mod(find(Result.SpClass(c,:)),time_segment);
    end
    %line up movie segments at simple spike

    load([fpath{f} '/output_data.mat'])
    sz=double(Device_Data{1, 3}.ROI([2 4]));
    frm_end=max(Device_Data{1, 2}.Counter_Inputs(1, 1).data);
    f_seg=[[1:time_segment:frm_end] frm_end+1]; f_seg(2:end)=f_seg(2:end)-1;
    g=zeros(1,3);
    ss=[0 0 0];
    for j=unique(cell2mat(fileInd))
        j

        mov_mc=double(readBinMov([fpath{f} '/mc_ShutterReg' num2str(j,'%02d') '.bin'],sz(2),sz(1)));
        load([fpath{f} '/mcTrace' num2str(j,'%02d') '.mat']);

        if j==length(f_seg)-1
            mc=mcTrace.xymean;
        else
            mov_mc=mov_mc(:,:,1:time_segment);
            mc=mcTrace.xymean(1:time_segment,:);
        end

        mov_res= mov_mc-mean(mov_mc,3);
        bkg = zeros(2, size(mov_mc,3));
        bkg(1,:) = linspace(-1, 1, size(mov_mc,3));  % linear term
        bkg(2,:) = linspace(-1, 1, size(mov_mc,3)).^2;  % quadratic term
        mov_res = SeeResiduals(mov_res,mc);
        mov_res = SeeResiduals(mov_res,mc.^2);
        mov_res = SeeResiduals(mov_res,mc(:,1).*mc(:,end));
        mov_res= SeeResiduals(mov_res,bkg,1);

        for c=1:3
            STA_movie_align{c}=[];
            for k=find(fileInd{c}==j)
                if frameInd{c}(k)+nTau{c}(1)>1 && frameInd{c}(k)+nTau{c}(end)<size(mov_res,3)
                    mov_seg=mov_res(bound:end-bound,bound:end-bound,frameInd{c}(k)+nTau{c});
                    %mov_seg=mat2gray(mov_seg);
                    STA_movie_align{c} = cat(3,STA_movie_align{c}, mov_seg);
                    g(c)=g(c)+1;
                    StackedSpike{c}(1,g(c))=1;
                    StackedSpike{c}(2,g(c))=time_segment*(j-1)+frameInd{c}(k);
                else
                    STA_movie_align{c} = cat(3,STA_movie_align{c}, zeros(sz(2)-bound*2+1,sz(1)-bound*2+1,length(nTau{c})));
                    g(c)=g(c)+1;
                    StackedSpike{c}(1,g(c))=0;
                    StackedSpike{c}(2,g(c))=NaN;
                end
            end
            if ~isempty(STA_movie_align{c})

                % write_tif_stack(STA_movie_align{c},fullfile(f_tmp,[alignedMovFN{c} '_' num2str(alignmov_ind(c)) '.tiff']))
                % alignmovlist=dir(fullfile(f_tmp,[alignedMovFN{c} '*.tiff']));  
                write_tif_stack(STA_movie_align{c},fullfile(fpath{f},[alignedMovFN{c} '_' num2str(alignmov_ind(c)) '.tiff']))
                alignmovlist=dir(fullfile(fpath{f},[alignedMovFN{c} '*.tiff']));
                if alignmovlist(alignmov_ind(c)).bytes > 2.0*10^9
                    disp('Move on to 2nd tiff')
                    alignmov_ind(c) = alignmov_ind(c) + 1;
                end
            end
        end
        g
    end

    Result.StackedSpike=StackedSpike;
    save(fullfile(fpath{f},'PC_Result.mat'),'Result','fpath','-v7.3')
end
