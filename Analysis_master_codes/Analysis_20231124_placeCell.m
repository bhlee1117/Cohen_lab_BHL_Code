clear
clc;
cd '/Volumes/BHL18TB_D2/YQ601_PlaceCellResults';
[~, fpath, raw] = xlsread(['/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/' ...
    'PlaceCellData_Arrangement.xlsx'], 'Sheet1', 'C8:C42');

[~, ~, NeuronsToUse]=xlsread(['/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/' ...
    'PlaceCellData_Arrangement.xlsx'], 'Sheet1', 'L8:M42');

NeuronsToUse=cellfun(@(x) (str2num(num2str(x))),NeuronsToUse,'UniformOutput',false);

save_figto='/Volumes/BHL_WD18TB/YQ601_PlaceCellResults';
place_bin=150;
velocity_threshold=0.002;
%% Parameter setting and get cell coordinate
block_size=30;

for i=1:length(fpath)
    clear bs ROI
    load(fullfile(fpath{i},"output_data.mat"))
    sz=double(Device_Data{1, 3}.ROI([2 4]));
    mov_test=double(readBinMov_times([fpath{i} '/frames1.bin'],sz(2),sz(1),[3000:4000]));
    avgImg=mean(mov_test,3);
    [CellCoordinate, ~]=Cell_segment_circle_25x_VU(avgImg,0.85);
    CellCoordinate=cell_detection_manual(avgImg,CellCoordinate,[]);
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
    saveas(gca,[fpath{i} '/ROIs.png'])
    save([fpath{i} '/Analysis_parameter.mat'],'CellCoordinate','bs','ROI','-v7.3')
end
%% motion correction
time_size=150000; %segment size
overlap=500;

for i=1:length(fpath)
    clear mcTrace_block ave_im
    %load the center positions
    load(fullfile(fpath{i},"output_data.mat"))
    load(fullfile(fpath{i},'Analysis_parameter.mat'))
    sz=double(Device_Data{1, 3}.ROI([2 4]));
    CamCounter=Device_Data{1, 2}.Counter_Inputs(1, 1).data;
    CamTrigger=find(CamCounter(2:end)-CamCounter(1:end-1));

    ref_time=[10000:11000];

    % make time segment
    t_seg=[[1:time_size:length(CamTrigger)] length(CamTrigger)+1];
    t_seg=[t_seg(1:end-1)' t_seg(2:end)'+overlap-1];
    t_seg(2:end,1)=t_seg(2:end,1)-overlap;
    t_seg(end)=length(CamTrigger);
    % t_seg=t_seg+119;


    for n=1:size(CellCoordinate,1)

        % load ROI reference
        mov_seg=double(readBinMov_times_ROI([fpath{i} '/frames1.bin'],sz(1),sz(2),ref_time,ROI(n,:)));

        % motion correction of reference
        options_rigid = NoRMCorreSetParms('d1',size(mov_seg,1),'d2',size(mov_seg,2),'bin_width',200,'max_shift',bs(n),'us_fac',50,'init_batch',200);
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
            [mov_mc,xyField]=optical_flow_motion_correction_LBH_ROIBlock(mov,double(mov_ref),'normcorre');

            ave_im{n,t}=mean(mov_mc,3);
            mov_mc=vm(mov_mc);
            mov_mc.transpose.savebin([fpath{i} '/mc' num2str(n,'%02d') '_' num2str(t,'%02d') '.bin'])

            %mcTrace_block{n,t} = squeeze(mean (xyField,[1 2]))'; %optic flow
            mcTrace_block{n,t}=xyField'; % Normcorre
        end


    end
    save([fpath{i} '/Analysis_parameter.mat'],'CellCoordinate','mcTrace_block','ave_im','ROI','bs','-v7.3')
    disp(['MC done:' fpath{i}])

end
%% Get footprint

time_size=150000; %segment size
overlap=500;

for i=9%1:length(fpath)
    %load files
    cd(fpath{i});
    load(fullfile(fpath{i},'Analysis_parameter.mat'))
    load([fpath{i} '/output_data.mat'])
    load(fullfile(fpath{i},'PC_Result.mat'),'Result');


    fileList = dir(fullfile(fpath{i}, '*.data'));
    if length(fileList)==1
        fid = fopen(fullfile(fpath{i},fileList.name));
        VRdata = fread(fid,[12 inf],'double');
    else
        error('Data file cannot be found');
    end

    WorldITrack=find(VRdata(2,:)==1); %world 1
    VRdata(5,WorldITrack)=(VRdata(5,WorldITrack)+6)*115/121;

    Result.centers=CellCoordinate;
    Result.fpath=fpath{i};
    disp(['loading  ' fpath{i}])

    sz=double(Device_Data{1, 3}.ROI([2 4]));
    mov_test=double(readBinMov_times([fpath{i} '/frames1.bin'],sz(2),sz(1),[10000:12000]));
    Result.FOV=mean(mov_test,3);

    CamCounter=Device_Data{1, 2}.Counter_Inputs(1, 1).data;
    CamTrigger=find(CamCounter(2:end)-CamCounter(1:end-1));
    t_seg=[[1:time_size:length(CamTrigger)] length(CamTrigger)+1];
    t_seg=[t_seg(1:end-1)' t_seg(2:end)'+overlap-1];
    t_seg(2:end,1)=t_seg(2:end,1)-overlap;
    t_seg(end)=length(CamTrigger);

    t_seg2=[[1:time_size:length(CamTrigger)] length(CamTrigger)+1];
    t_seg2=[t_seg2(1:end-1)' t_seg2(2:end)'];
    t_seg2(1:end-1,2)=t_seg2(1:end-1,2)-1;
    t_seg2(end)=length(CamTrigger);

    extract_seg = repmat([1 time_size],size(t_seg,1),1);
    extract_seg(2:end,:)=extract_seg(2:end,:)+overlap;
    extract_seg(end,2)=mod(length(CamTrigger),time_size)+overlap;

    Result.traces=[];
    Result.traces_res=[];
    Result.im_corr=[];
    Result.mcTrace=[];
    Result.c_ftprnt=[];
    Result.ref_im=[];
    for n=1:size(Result.centers,1)

        mov=double(readBinMov([fpath{i} '/mc' num2str(n,'%02d') '_' num2str(1,'%02d') '.bin'], ...
            bs(n)*2+1,bs(n)*2+1));
        mov=mov(:,:,extract_seg(1,1):extract_seg(1,2));
        mcTrace_block{n,1}=mcTrace_block{n,1}(:,extract_seg(1,1):extract_seg(1,2));

        mov_res= mov-mean(mov,3);
        bkg = zeros(2, size(mov,3));
        bkg(1,:) = linspace(-1, 1, size(mov,3));  % linear term
        bkg(2,:) = linspace(-1, 1, size(mov,3)).^2;  % quadratic term
        mcTrace_block{n,1}=movmean(mcTrace_block{n,1},10,2);
        mov_res=SeeResiduals(mov_res,mcTrace_block{n,1});
        mov_res=SeeResiduals(mov_res,mcTrace_block{n,1}.^2);
        mov_res=SeeResiduals(mov_res,mcTrace_block{n,1}(1,:).*mcTrace_block{n,1}(2,:));
        mov_res= SeeResiduals(mov_res,bkg,1);

        [u,s,v] = svds(tovec(mov_res(:,:,100000:110000)),20);
        reshape_u=reshape(u,size(mov_res,2),size(mov_res,1),[]);
        bvMask=[]; Result.bvMask{n}=[];
        [~, Result.bvMask{n}]=get_ROI(max(abs(reshape_u),[],3),bvMask);
        plot(rescale2(tovec(Result.bvMask{n})'*tovec(mov_res),2)'+[1:size(Result.bvMask{n},3)])

        [~, pcaTrace, icsTrace]=clickyICA(imresize(mov_res,0.5),imresize(mean(mov,3),0.5),5);
        icsImgs=toimg(tovec(mov_res)*icsTrace.intens',size(mov_res,1),size(mov_res,2));
        pcaImgs=toimg(tovec(mov_res)*pcaTrace.intens',size(mov_res,1),size(mov_res,2));

        %ROIimg=icsImgs(:,:,1);
        ROIimg=sqrt(mean(mov_res(:,:,2:end).*mov_res(:,:,1:end-1),3));
        excludeImg=mean(icsImgs(:,:,[3]),3);
        %excludeImg=mean(mov_mc,3);

        figure(13); clf;
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
        close(figure(13));
        Result.ROIpoly=ROIpoly;

        % set exclude ROI
        figure; Result.excludeROI{n}=[];
        [~, Result.excludeROI{n}]=get_ROI(excludeImg,[],'exclude ROI');
        % if ifregressROI(f)
        % regressTrace=tovec(mov_res)'*tovec(Result.excludeROI);
        % mov_res = SeeResiduals(mov_res,regressTrace);
        % end

        n_comp=6;
        mov_res_reg=mov_res;
        %mov_res_reg=SeeResiduals(mov_res_reg,icsTrace.intens(5,:));
        mov_filt=imgaussfilt3(mov_res_reg.*double(max(Result.bvMask{n},[],3)==0).*double(max(Result.excludeROI{n},[],3)==0),[1.5 1.5 0.1]);
        mov_filt=mov_filt(:,:,1:110000);
        movVec=tovec(mov_filt);
        Npoly=size(Result.ROIpoly,1);
        ftprnt = zeros(size(mov_filt,1)*size(mov_filt,2),Npoly);
        clear mask
        figure(4);
        for p=1:Npoly %each ROIs
            clf; ax2=[];
            tiledlayout(n_comp/2+2,2)
            mask(:,:,p) = poly2mask(Result.ROIpoly{p}(:,1), Result.ROIpoly{p}(:,2), size(mov_res,2), size(mov_res,1));
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
            for np=1:n_comp
                eigImg=NaN(size(mov_filt,1)*size(mov_filt,2),1);
                ax2=[ax2 nexttile([1 1])];
                eigImg(pixelList,1)=V(:,np);
                eigImg=toimg(eigImg,size(mov_filt,1),size(mov_filt,2));
                %imshow2(im_merge(cat(3,mean(mov,3),eigImg),[1 1 1;1 0 0]),[])
                imshow2(eigImg,[])
                title([num2str(np) ', Fraction: ' num2str(D(np)/sum(D),2)])
            end
            linkaxes(ax2,'xy')
            n_take = input('#components to take: ', 's');
            n_take = str2num(n_take);
            coeff=subMov*mean(eigTrace(:,n_take)*V(:,n_take)',2);
            ftprnt(pixelList,p)=coeff;
        end
        close(figure(4));

        Result.c_ftprnt(1:2*bs(n)+1,1:2*bs(n)+1,n)=toimg(ftprnt,size(mov_res,2),size(mov_res,1));
        figure(6+n); clf; tiledlayout(2,2);
        show_footprnt(Result.c_ftprnt,mean(mov,3))
        tr_prnt=tovec(Result.c_ftprnt)'*tovec(mov_res(:,:,50000:100000));
        tr_poly=tovec(double(Result.c_ftprnt>0))'*tovec(mov_res(:,:,50000:100000));

        nexttile([1 2])
        plot(zscore(-tr_prnt)); hold all
        plot(zscore(-tr_poly));

        plot(rescale2(mcTrace_block{n,1},2)'+2);
        drawnow;

        % Cellpsf = fspecial('gaussian', 2*bs(n)+1, 10);
        % Result.c_ftprnt(1:2*bs(n)+1,1:2*bs(n)+1,n)=Result.c_ftprnt(1:2*bs(n)+1,1:2*bs(n)+1,n).*Cellpsf;
    end

    save(fullfile(fpath{i},'PC_Result.mat'),'Result','-v7.3')
    load(fullfile(fpath{i},'PC_Result.mat'),'Result');
backupServer(fpath{i},'BHL18TB_D2','cohen_lab/Lab/Labmembers/Byung Hun Lee/Data','PC_Result.mat')
end
%% Signal Extraction

time_size=150000; %segment size
overlap=500;

for i=9%1:length(fpath)
    %load files
    Result=[];
    load(fullfile(fpath{i},'Analysis_parameter.mat'))
    load([fpath{i} '/output_data.mat'])
    fileList = dir(fullfile(fpath{i}, '*.data'));
    if length(fileList)==1
        fid = fopen(fullfile(fpath{i},fileList.name));
        VRdata = fread(fid,[12 inf],'double');
    else
        error('Data file cannot be found');
    end

    WorldITrack=find(VRdata(2,:)==1); %world 1
    VRdata(5,WorldITrack)=(VRdata(5,WorldITrack)+6)*115/121;

    VRdata=VRdata(:,VRdata(10,:)>0);
    DAQ_rate=Device_Data{1, 2}.Counter_Inputs.rate;
    CamCounter=Device_Data{1, 2}.Counter_Inputs(1, 1).data;
    CamTrigger=find(CamCounter(2:end)-CamCounter(1:end-1));
    t_DAQ=CamTrigger/DAQ_rate;
    t_VR = datetime(datetime(VRdata(1,:),'ConvertFrom','datenum'), 'InputFormat', 'yyyy-MM-dd HH:mm:ss.SSS');
    t_VR= t_VR-t_VR(1);
    t_VR= milliseconds(t_VR)/1000;
    t_VR= t_VR*t_DAQ(end)/t_VR(end);
    VRdata(1,:)=t_VR;
    [Virmen_data_int vel_trace]=virmen_interpolate(VRdata,115,t_DAQ);
    Virmen_data_int(end+1,:)=vel_trace;

    Result.VR=Virmen_data_int;
    Result.centers=CellCoordinate;
    Result.frm_rate=double((CamTrigger(2)-CamTrigger(1))*DAQ_rate);
    Result.fpath=fpath{i};
    disp(['loading  ' fpath{i}])

    sz=double(Device_Data{1, 3}.ROI([2 4]));
    mov_test=double(readBinMov_times([fpath{i} '/frames1.bin'],sz(2),sz(1),[10000:11000]));
    Result.FOV=mean(mov_test,3);

    Result.Blue=Device_Data{1, 2}.buffered_tasks(1, 1).channels(1, 2).data(CamTrigger).*Virmen_data_int(12,:);
    Result.Reward=Device_Data{1, 2}.buffered_tasks(1, 3).channels.data(CamTrigger);

    t_seg=[[1:time_size:length(CamTrigger)] length(CamTrigger)+1];
    t_seg=[t_seg(1:end-1)' t_seg(2:end)'+overlap-1];
    t_seg(2:end,1)=t_seg(2:end,1)-overlap;
    t_seg(end)=length(CamTrigger);

    t_seg2=[[1:time_size:length(CamTrigger)] length(CamTrigger)+1];
    t_seg2=[t_seg2(1:end-1)' t_seg2(2:end)'];
    t_seg2(1:end-1,2)=t_seg2(1:end-1,2)-1;
    t_seg2(end)=length(CamTrigger);

    extract_seg = repmat([1 time_size],size(t_seg,1),1);
    extract_seg(2:end,:)=extract_seg(2:end,:)+overlap;
    extract_seg(end,2)=mod(length(CamTrigger),time_size)+overlap;

    Result.traces=[];
    Result.traces_res=[];
    Result.im_corr=[];
    Result.mcTrace=[];
    Result.c_ftprnt=[];
    Result.ref_im=[];
    for n=1:size(Result.centers,1)

        mov=double(readBinMov([fpath{i} '/mc' num2str(n,'%02d') '_' num2str(1,'%02d') '.bin'], ...
            bs(n)*2+1,bs(n)*2+1));
        mov=mov(:,:,extract_seg(1,1):extract_seg(1,2));
        mcTrace_block{n,1}=mcTrace_block{n,1}(:,extract_seg(1,1):extract_seg(1,2));

        %mov_res= mov-mean(mov,3);
        %bkg = zeros(2, size(mov,3));
        % bkg(1,:) = linspace(-1, 1, size(mov,3));  % linear term
        % bkg(2,:) = linspace(-1, 1, size(mov,3)).^2;  % quadratic term
        %         mov_res=SeeResiduals(mov_res,mcTrace_block{n,1});
        %         mov_res=SeeResiduals(mov_res,mcTrace_block{n,1}.^2);
        %         mov_res=SeeResiduals(mov_res,mcTrace_block{n,1}(1,:).*mcTrace_block{n,1}(2,:));
        %mov_res= SeeResiduals(mov_res,bkg,1);


        Result.ref_im{n}=mean(mov,3);
        ref_im_vec=tovec(Result.ref_im{n});
        ref_im_vec=(ref_im_vec-mean(ref_im_vec,1))./std(ref_im_vec,0,1);
        Result.im_corr{n}=[];

        Result.c_ftprnt(1:2*bs(n)+1,1:2*bs(n)+1,n)=mask_footprint([bs(n)+0.5 bs(n)+0.5],mov_res(:,:,1000:end-1000),[],18);
        Cellpsf = fspecial('gaussian', 2*bs(n)+1, 10);
        Result.c_ftprnt(1:2*bs(n)+1,1:2*bs(n)+1,n)=Result.c_ftprnt(1:2*bs(n)+1,1:2*bs(n)+1,n).*Cellpsf;
        Result.c_ftprnt(1:2*bs(n)+1,1:2*bs(n)+1,n)=imgaussfilt(Result.c_ftprnt(1:2*bs(n)+1,1:2*bs(n)+1,n),2);

        for t=2:size(t_seg,1)
            Result.traces(n,t_seg2(t-1,1):t_seg2(t-1,2))=-(tovec(mov)'*tovec(Result.c_ftprnt(1:2*bs(n)+1,1:2*bs(n)+1,n)))';
            Result.mcTrace(:,t_seg2(t-1,1):t_seg2(t-1,2),n)=mcTrace_block{n,t-1};
            mov_mc_vec=tovec(mov);
            mov_mc_vec=(mov_mc_vec-mean(mov_mc_vec,1))./std(mov_mc_vec,0,1);
            Result.im_corr{n}=[Result.im_corr{n} sum(mov_mc_vec.*ref_im_vec,1)/(size(mov_mc_vec,1)-1)];

            mov=double(readBinMov([fpath{i} '/mc' num2str(n,'%02d') '_' num2str(t,'%02d') '.bin'], ...
                bs(n)*2+1,bs(n)*2+1));

            mov=mov(:,:,extract_seg(t,1):extract_seg(t,2));
            mcTrace_block{n,t}=mcTrace_block{n,t}(:,extract_seg(t,1):extract_seg(t,2));

            % mov_res= mov-mean(mov,3);
            % bkg = zeros(2, size(mov,3));
            % bkg(1,:) = linspace(-1, 1, size(mov,3));  % linear term
            % bkg(2,:) = linspace(-1, 1, size(mov,3)).^2;  % quadratic term
            % mov_res=SeeResiduals(mov_res,mcTrace_block{n,t});
            % mov_res=SeeResiduals(mov_res,mcTrace_block{n,t}.^2);
            % mov_res= SeeResiduals(mov_res,bkg,1);

        end
        Result.traces(n,t_seg2(end,1):t_seg2(end,2))=-(tovec(mov)'*tovec(Result.c_ftprnt(1:2*bs(n)+1,1:2*bs(n)+1,n)))';
        Result.mcTrace(:,t_seg2(end,1):t_seg2(end,2),n)=mcTrace_block{n,end};
        mov_mc_vec=tovec(mov);
        mov_mc_vec=(mov_mc_vec-mean(mov_mc_vec,1))./std(mov_mc_vec,0,1);
        Result.im_corr{n}=[Result.im_corr{n} sum(mov_mc_vec.*ref_im_vec,1)/(size(mov_mc_vec,1)-1)];


        % mcT=squeeze(Result.mcTrace(:,:,n));
        % Result.traces_res(n,:)=squeeze(SeeResiduals(reshape(Result.traces(n,:),1,1,[]),mcT));
        % Result.traces_res(n,:)=squeeze(SeeResiduals(reshape(Result.traces_res(n,:),1,1,[]),mcT.^2'));
        % Result.traces_res(n,:)=squeeze(SeeResiduals(reshape(Result.traces_res(n,:),1,1,[]),mcT(1,:).*mcT(2,:)));
        Result.AvgImg=mean(mov,3);
    end

    % Result.spike=zeros(size(Result.traces));
    % tmp=squeeze(Result.traces_res) - movmedian(squeeze(Result.traces_res),250,2); tmp=tmp./get_threshold(tmp,1);
    % Result.spike=(Result.spike+find_spike_bh(tmp,5.5,2))>0;

    save(fullfile(fpath{i},'PC_Result.mat'),'Result','fpath','-v7.3')
      load(fullfile(fpath{i},'PC_Result.mat'),'Result');
backupServer(fpath{i},'BHL18TB_D2','cohen_lab/Lab/Labmembers/Byung Hun Lee/Data','PC_Result.mat')
end
%% Recollect

for i=1:length(fpath)
    Result_tmp=load(fullfile(fpath{i},'PC_Result.mat'));
    PC_Result{i}=Result_tmp.Result;
end

%% Clean up and bleach correction

exclude_frq=[241.7 242]; %monitor
%exclude_frq2=[483.5 484]; %monitor
exclude_frq2=[20 65]; %motion
time_bin=50000; Fs=1000;

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
figure; clf;
for i=9%1:length(fpath)
    clear traces_res_filtered noise noise_intp norm_trace sp_height SpHeight_intp sp_time
    tN=[1:time_bin:size(PC_Result{i}.traces,2)]; tN=[tN size(PC_Result{i}.traces,2)];
    sp_time=zeros(size(PC_Result{i}.traces,1),size(PC_Result{i}.traces,2));
    sp_height=zeros(size(PC_Result{i}.traces,1),size(PC_Result{i}.traces,2));
    PC_Result{i}.subThreshold=[]; PC_Result{i}.theta=[];

    for n=1:size(PC_Result{i}.traces,1)

        mcTrace=squeeze(PC_Result{i}.mcTrace(:,:,n));
        tr=PC_Result{i}.traces(n,:);
        trhi=tr-movmedian(tr,100);
        sphi=find_spike_bh(trhi,6,4);
        sphi_vec=ind2vec(size(tr,2),unique(find(sphi)'+[-2:10]),1);
        blue_vec=ind2vec(size(tr,2),unique(find(PC_Result{i}.Blue)'+[-30:100]),1);
       
        t_fit=setdiff([1:size(tr,2)],[find(sphi_vec) find(blue_vec)]);

        [y_fit_tr t_consts coeffY]  = expfitDM_2(t_fit',-tr(1,t_fit)',[1:size(tr,2)]',[10^5]);

        % lwpass_fit=nan(1,size(tr,2));
        % lwpass_fit(:,t_fit)=sum(tr(1,t_fit),1,'omitnan'); lwpass_fit=movmedian(movprc(lwpass_fit,20000,30,2),30000,2);
        tr=squeeze(SeeResiduals(reshape(tr,1,1,[]),y_fit_tr))';

        % regress out motion frequency
        for t=1:length(tN)-1
            tr(tN(t):tN(t+1))=squeeze(SeeResiduals(reshape(tr(tN(t):tN(t+1)),1,1,[]),mcTrace(:,(tN(t):tN(t+1)))));
            tr(tN(t):tN(t+1))=squeeze(SeeResiduals(reshape(tr(tN(t):tN(t+1)),1,1,[]),mcTrace(:,(tN(t):tN(t+1))).^2));
            tr(tN(t):tN(t+1))=squeeze(SeeResiduals(reshape(tr(tN(t):tN(t+1)),1,1,[]),mcTrace(1,(tN(t):tN(t+1))).*mcTrace(2,(tN(t):tN(t+1)))));
            imcorr_tmp=PC_Result{i}.im_corr{n}(tN(t):tN(t+1));
            imcorr_tmp=movmean(imcorr_tmp,5,'omitnan');
            %tr(tN(t):tN(t+1))=squeeze(SeeResiduals(reshape(tr(tN(t):tN(t+1)),1,1,[]),imcorr_tmp));
        end

        % regress out motion frequency
        traces_res_filtered(n,:) = filtfilt(b, a, tr); %monitor
        %traces_res_filtered(n,:) = filtfilt(b2, a2,traces_res_filtered(n,:)); %motion
        norm_trace(n,:)=traces_res_filtered(n,:);%-movmedian(traces_res_filtered(n,:),500,2);


        for t=1:length(tN)-2
            tr_tmp=norm_trace(n,tN(t):tN(t+1));
            tr_tmp=tr_tmp-movmedian(tr_tmp,300);
            noise(t,n)=get_threshold(tr_tmp,1);
            tr_tmp_norm=tr_tmp./noise(t,n);
            [sp_temp, pks, prom]=find_spike_bh(tr_tmp_norm,4,3);
            sp_time(n,tN(t):tN(t+1))=sp_temp;
        end
        t_fit=find(sp_time(n,:));
        sp_height(n,t_fit)=norm_trace(n,t_fit);
        tx=tN(1:end-1)+time_bin/2;
        %noise_intp(n,:)=movmean(interp1(tx,noise(:,n),[1:size(PC_Result{i}.traces,2)],'linear','extrap'),10000);

        [Ny_fit t_consts coeffY]  = expfitDM_2(tx(~isnan(noise(1:end-1,n)))',noise(~isnan(noise(1:end-1,n)),n),[1:size(PC_Result{i}.traces,2)]',10^7);
        noise_intp(n,:)=Ny_fit;
        %[y_fit t_consts coeffY]  = expfitDM_2(tx(~isnan(sp_height(1:end-1,n)))',sp_height(~isnan(sp_height(1:end-1,n)),n),[1:size(PC_Result{i}.traces,2)]',10^7);
        [y_fit t_consts coeffY]  = expfitDM_2(t_fit',sp_height(n,t_fit)',[1:size(PC_Result{i}.traces,2)]',10^7);
        SpHeight_intp(n,:)=y_fit;

    end
    nexttile([1 1])
    plot(norm_trace')
    hold all
    plot(SpHeight_intp')
    title(num2str(i))

    %norm_trace=norm_trace./noise_intp;
    norm_trace=norm_trace./SpHeight_intp;
    %norm_trace=norm_trace;%./(SpHeight_intp./SpHeight_intp(:,1));
    PC_Result{i}.normTraces=norm_trace./get_threshold(norm_trace,1);
    %PC_Result{i}.spike=find_spike_bh(PC_Result{i}.normTraces-movmedian(PC_Result{i}.normTraces,300,2),4,3);


    % for n=1:size(PC_Result{i}.traces,1)
    %     PC_Result{i}.subThreshold(n,:) = filtfilt(b3, a3, PC_Result{i}.normTraces(n,:));
    %     PC_Result{i}.theta(n,:) = filtfilt(b4, a4, PC_Result{i}.normTraces(n,:));
    % 
    %     %Result{i}.subThreshold(n,:) = movmean(Result{i}.normTraces(n,:),1000);
    % end

end


%%

for g=1:length(fpath)
    g
    [sz_fprntY, sz_fprntX, nNeurons]=size(PC_Result{g}.c_ftprnt);
    sz=size(PC_Result{g}.FOV);
    PC_Result{g}.c_fprnt_merge=zeros(sz(1),sz(2),nNeurons);

    PC_Result{g}.Lap_FR=[]; PC_Result{g}.Lap_sub=[]; PC_Result{g}.Lap_theta=[];
    PC_Result{g}.Lap_F=[]; PC_Result{g}.Lap_V=[];
    for n=1:nNeurons
        %offset=round(PC_Result{g}.centers(n,:)+[-floor(sz_fprntY/2):1:floor(sz_fprntY/2); -floor(sz_fprntX/2):1:floor(sz_fprntX/2)]');
        %PC_Result{g}.c_fprnt_merge(offset(:,2),offset(:,1),n)=PC_Result{g}.c_ftprnt(:,:,n)';
        [PC_Result{g}.Lap_FR(:,:,n) PC_Result{g}.Lap_V]=PlaceTrigger_average(PC_Result{g}.spike(n,:),place_bin,PC_Result{g}.VR,velocity_threshold,115);
        [PC_Result{g}.Lap_F(:,:,n) PC_Result{g}.Lap_V]=PlaceTrigger_average(PC_Result{g}.normTraces(n,:),place_bin,PC_Result{g}.VR,velocity_threshold,115);
        LP_subthreshold=movmean(PC_Result{g}.subThreshold(n,:),1000);
        %LP_subthreshold=Result{g}.subThreshold(n,:);
        [PC_Result{g}.Lap_sub(:,:,n) PC_Result{g}.Lap_V]=PlaceTrigger_average(LP_subthreshold,place_bin,PC_Result{g}.VR,velocity_threshold,115);
        [PC_Result{g}.Lap_theta(:,:,n) PC_Result{g}.Lap_V]=PlaceTrigger_average(abs(PC_Result{g}.theta(n,:)),place_bin,PC_Result{g}.VR,velocity_threshold,115);
    end
    Result=PC_Result{g};
    save(fullfile(fpath{g},'PC_Result.mat'),'Result','fpath','-v7.3')
end

%% Recollect

for i=1:length(fpath)
    Result_tmp=load(fullfile(fpath{i},'PC_Result.mat'));
    PC_Result{i}=Result_tmp.Result;
end


%% Show place fields
save_figto='/Volumes/BHL_WD18TB/YQ601_PlaceCellResults';

f1=figure(1); clf;
for i=1:length(PC_Result)
    [sz_fprntY, sz_fprntX, nNeurons]=size(PC_Result{i}.c_ftprnt);
    StimN=find(mean(PC_Result{i}.spike(:,find(PC_Result{i}.Blue>0)),2)>mean(PC_Result{i}.spike,2)*5)';
    UnstimN=setdiff([1:nNeurons],StimN);
    for n=[StimN UnstimN]
        nexttile([1 1])
        pos_bin=size(PC_Result{i}.Lap_FR,2);
        Lap_tmp=repmat(PC_Result{i}.Lap_FR(:,:,n),1,3);
        Lap_tmp=movmean(Lap_tmp,7,2);
        Lap_tmp=Lap_tmp(:,pos_bin+1:2*pos_bin);
        imagesc(Lap_tmp)
        colormap('turbo')
        if ismember(n,UnstimN)
            title(['Unstimulated ' num2str(i) '-' num2str(n)])
        else
            title([num2str(i) '-' num2str(n)])
        end
    end
end

set(f1, 'Position', [100, 100, 800, 400]);
saveas(f1,fullfile(save_figto ,['PF_map' '.fig']))
%print(f1, fullfile(save_figto ,['PF_map' '.jpg']),'-djpeg', ['-r', num2str(400)]);

f2=figure(2); clf;
cmap=turbo(length(PC_Result));
g=1;
for i=1:length(PC_Result)
    [sz_fprntY, sz_fprntX, nNeurons]=size(PC_Result{i}.c_ftprnt);
    StimN=find(mean(PC_Result{i}.spike(:,find(PC_Result{i}.Blue>0)),2)>mean(PC_Result{i}.spike,2)*5)';
    UnstimN=setdiff([1:nNeurons],StimN);
    for n=[StimN UnstimN]
        tr_rescaled=rescale(PC_Result{i}.normTraces(n,1:end-20));
        spike_time=find(PC_Result{i}.spike(n,1:end-20));

        plot([1:length(tr_rescaled)],tr_rescaled+g,'color',cmap(i,:)); hold all
        plot(spike_time,tr_rescaled(spike_time)+g,'r.')
        g=g+1;
    end
end
save(fullfile(save_figto,['PF_Result_20231126.mat']),'PC_Result','-v7.3')
%% Load Data
load(fullfile(save_figto,['PF_Result_20231126.mat']))

%% Calculate Firing rate change upon stimulation
FRMat_Time=[]; VolMat_Time=[];
FRMat_Pos=[]; VolMat_Pos=[];

StimulationWaveform_indicator=[];
StimulationWaveform_actuator=[];
UnstimulationVoltageMat=[];
UnstimulationSpikeMat=[];
AverageLapN=5;
TimeWindow=[-10000:10000];
PositionWindow=[-12 12]; %VR units
Wvf_window=[-500:1500];
Source_file=[];

g=1; %Each stimulation
for i=1:size(NeuronsToUse,1)
    if ~isnan(NeuronsToUse{i,1})

        NOI=NeuronsToUse{i,1};
        StimLaps=unique(double(PC_Result{i}.Blue>0).*PC_Result{i}.VR(8,:)); StimLaps=StimLaps(2:end);
        BlueLap=zeros(1,max(PC_Result{i}.VR(8,:)));
        World1Lap=unique(PC_Result{i}.VR(8,:).*(PC_Result{i}.VR(2,:)==1)); World1Lap=World1Lap(2:end);
        BlueLap(StimLaps)=1;
        BlueLap=bwlabel(~BlueLap); BlueLap(World1Lap)=NaN;
        LapOfInterest=[];
        %Lap_FR_vec=reshape(PC_Result{i}.Lap_FR(:,:,NeuronInterest)',1,[]);
        
        for NeuronInterest=NOI
        for b=1:max(BlueLap)-1
            LapOfInterest{b,1}=find(BlueLap==b); %before
            if length(LapOfInterest{b,1})>AverageLapN
                LapOfInterest{b,1}=LapOfInterest{b,1}(end-AverageLapN+1:end);
            end
            LapOfInterest{b,2}=find(BlueLap==b+1); %after
            if length(LapOfInterest{b,2})>AverageLapN
                LapOfInterest{b,2}=LapOfInterest{b,2}(1:AverageLapN);
            end

            StimulationWaveform_indicator{g}=[];
            StimulationWaveform_actuator{g}=[];
            Stimulation_Lap{g}=[];
            Stimulation_Pos{g}=[];
            for StimulationLap=LapOfInterest{b,1}(end)+1:LapOfInterest{b,2}(1)-1

                StimOnFrame = find(PC_Result{i}.VR(8,:)==StimulationLap & PC_Result{i}.Blue>0,1);
                StimOnPosition = PC_Result{i}.VR(5,StimOnFrame);
                StimOnFrame_Tau = StimOnFrame + Wvf_window;
                StimulationWaveform_indicator{g}=[StimulationWaveform_indicator{g}; PC_Result{i}.normTraces(NeuronInterest,StimOnFrame_Tau)];
                StimulationWaveform_actuator{g}=[StimulationWaveform_actuator{g}; PC_Result{i}.Blue(StimOnFrame_Tau)];
                Stimulation_Lap{g}=StimulationLap;
                Stimulation_Pos{g}=StimOnPosition;

            end

            for ab=1:2
                FRMat_Time{g,ab}=[];
                VolMat_Time{g,ab}=[];
                FRMat_Pos{g,ab}=[];

                for l=1:length(LapOfInterest{b,ab})
                    loi=LapOfInterest{b,ab}(l);
                    l_frame=find(PC_Result{i}.VR(8,:)==loi);
                    [~, PF_frame]=min(abs(PC_Result{i}.VR(5,l_frame)-StimOnPosition));
                    PF_frame=PF_frame+l_frame(1)-1;

                    if PF_frame+TimeWindow(1)<0  %Truncated in the front
                        repnan=sum((PF_frame+TimeWindow)<1);
                        FRMat_Time{g,ab}=[FRMat_Time{g,ab}; [NaN(1,repnan) PC_Result{i}.spike(NeuronInterest,PF_frame+TimeWindow(repnan+1:end))]];
                        VolMat_Time{g,ab}=[VolMat_Time{g,ab}; [NaN(1,repnan) PC_Result{i}.normTraces(NeuronInterest,PF_frame+TimeWindow(repnan+1:end))]];

                    else if PF_frame+TimeWindow(end)> size(PC_Result{i}.normTraces,2) %Truncated in the end

                            repnan=sum((PF_frame+TimeWindow)> size(PC_Result{i}.normTraces,2));
                            FRMat_Time{g,ab}=[FRMat_Time{g,ab}; [PC_Result{i}.spike(NeuronInterest,PF_frame+TimeWindow(1:end-repnan)) NaN(1,repnan)]];
                            VolMat_Time{g,ab}=[VolMat_Time{g,ab};  [PC_Result{i}.normTraces(NeuronInterest,PF_frame+TimeWindow(1:end-repnan)) NaN(1,repnan)]];

                    else
                        FRMat_Time{g,ab}=[FRMat_Time{g,ab}; PC_Result{i}.spike(NeuronInterest,PF_frame+TimeWindow)];
                        VolMat_Time{g,ab}=[VolMat_Time{g,ab}; PC_Result{i}.normTraces(NeuronInterest,PF_frame+TimeWindow)];
                    end
                    end

                    TakeoutPosition = 115 * (loi-1) * (loi>1) + StimOnPosition + PositionWindow;
                    TakeoutFrame =find(PC_Result{i}.VR(end-1,:) > TakeoutPosition(1) & PC_Result{i}.VR(end-1,:) < TakeoutPosition(2));
                    RunFrame=PC_Result{i}.VR(end,TakeoutFrame) > velocity_threshold;
                    FR_position = sum(PC_Result{i}.spike(NeuronInterest,TakeoutFrame).*RunFrame)/sum(RunFrame)*1000;

                    FRMat_Pos{g,ab}=[FRMat_Pos{g,ab}; FR_position];
                    
                end
            end
            Source_file(g)=i;
            g=g+1;
        end
        end
    end
end

%% Plot Firing rate change for all the cases
t=[-10:0.001:10];
Bin_time=500; %0.5 sec
t_bin=t([1:Bin_time:end-1])+Bin_time/1000/2;
for stim_type=1:3
    figure(stim_type); clf;
        StimTypeSort=find(Stimulation_type==stim_type);
dFRTrace=cellfun(@(x) movsum(x,Bin_time,2,'omitnan')/Bin_time*1000, FRMat_Time(StimTypeSort,:),'UniformOutput',false);
dFRTrace=cellfun(@(x) mean(x(:,Bin_time/2+1:Bin_time:end),1), dFRTrace,'UniformOutput',false);
beforeFR=cell2mat(dFRTrace(:,1)); AfterFR=cell2mat(dFRTrace(:,2));
for s=1:length(StimTypeSort)
nexttile([1 1])
plot(t_bin,beforeFR(s,:)); hold all
plot(t_bin,AfterFR(s,:));
end
end
%% Firing rate change

figure(3); clf;
cmap=distinguishable_colors(3);  cmap_light=cmap*1.2; cmap_light(cmap_light>1)=1; 
Stimulation_type=cell2mat(cellfun(@(x) max(max(bwlabel(x>0.1),[],2)),StimulationWaveform_actuator,'UniformOutput',false));
Stimulation_type(Stimulation_type>2 & Stimulation_type<10)=2; 
Stimulation_type(Stimulation_type>10)=3;
Bin_time=500; %0.5 sec
Mean_time=[-2000 2000];
MeanFR_stimtype=[]; MeanFR_stimtype_Pos=[];

tiledlayout(5,6)
for stim_type=1:3
    nexttile([1 2])
    StimTypeSort=find(Stimulation_type==stim_type);
dFRTrace=cellfun(@(x) movsum(x,Bin_time,2,'omitnan')/Bin_time*1000, FRMat_Time(StimTypeSort,:),'UniformOutput',false);
dFRTrace=cellfun(@(x) mean(x(:,Bin_time/2+1:Bin_time:end),1), dFRTrace,'UniformOutput',false);
beforeFR=cell2mat(dFRTrace(:,1)); AfterFR=cell2mat(dFRTrace(:,2));
deltaFR=AfterFR-beforeFR;
tax=[1:size(beforeFR,2)]*Bin_time/1000;
M=mean(deltaFR,1,'omitnan'); S=std(deltaFR,0,1,'omitnan');
errorbar_shade(tax,M,S,cmap(stim_type,:))
hold all
set(gca,'xtick',tax([1:size(M,2)/4:size(M,2)]),'xticklabel',(tax([1:size(M,2)/4:size(M,2)])-tax(size(M,2)/2+1)))
axis tight
ylabel('\Delta Firing rate (Hz)')
xlabel('Time (s)')
MeanFR_stimtype{stim_type}=cell2mat(cellfun(@(x) mean(sum(x(:,round(end/2)+Mean_time(1):round(end/2)+Mean_time(2)),2,'omitnan')/(Mean_time(2)-Mean_time(1))*1000,1,'omitnan'), FRMat_Time(StimTypeSort,:),'UniformOutput',false));
MeanFR_stimtype_Pos{stim_type}=cell2mat(cellfun(@(x) mean(x,1,'omitnan'),FRMat_Pos(StimTypeSort,:),'UniformOutput',false));
ylim([-3 6])
end

nexttile([2 3])
for stim_type=1:3
plot(MeanFR_stimtype{stim_type}(:,1),MeanFR_stimtype{stim_type}(:,2),'.','color',cmap(stim_type,:),'markersize',20)
hold all
end
plot([0 21],[0 21],'--','color','k','linewidth',2)
title(['Mean Response within ± 2 sec of' char(10) 'Stimulation Onset Position'])
xlabel('Firing rate before stimulation')
ylabel('Firing rate after stimulation')
axis tight

nexttile([2 3])
for stim_type=1:3
plot(MeanFR_stimtype_Pos{stim_type}(:,1),MeanFR_stimtype_Pos{stim_type}(:,2),'.','color',cmap(stim_type,:),'markersize',20)
hold all
end
plot([0 30],[0 30],'--','color','k','linewidth',2)
title(['Mean Response within ± 12 cm of' char(10) 'Stimulation Onset Position'])
xlabel('Firing rate before stimulation')
ylabel('Firing rate after stimulation')
axis tight


nexttile([2 3])
LapBinSize=6;
for stim_type=1:3
StimTypeSort=find(Stimulation_type==stim_type);
ZapLap=cell2mat(Stimulation_Lap(StimTypeSort));
thingsToPlot=MeanFR_stimtype_Pos{stim_type}(:,2)./MeanFR_stimtype_Pos{stim_type}(:,1);
thingsToPlot(thingsToPlot==inf)=NaN;
%plot(ZapLap,MeanFR_stimtype_Pos{stim_type}(:,2)-MeanFR_stimtype_Pos{stim_type}(:,1),'.','color',cmap(stim_type,:),'markersize',20)
plot(ZapLap,thingsToPlot,'.','color',cmap_light(stim_type,:),'markersize',7)
hold all
Average_FR_Ratio=[]; Std_FR_Ratio=[];
for b=1:ceil(max(ZapLap)/LapBinSize)
    bin_list=find(ceil(ZapLap/LapBinSize)==b);
    Average_FR_Ratio(b)=mean(thingsToPlot(bin_list),'omitnan');
    Std_FR_Ratio(b)=std(thingsToPlot(bin_list),'omitnan');
end
ZapLap_bin=[1:ceil(max(ZapLap)/LapBinSize)]*LapBinSize-LapBinSize/2;
errorbar(ZapLap_bin,Average_FR_Ratio,Std_FR_Ratio,'color',cmap(stim_type,:),'linewidth',2,'marker','+','capsize',10)
end
%set(gca,'yscale','log')
plot([0 35],[1 1],'--','color','k','linewidth',2)
xlabel('Stimulation Lap #')
ylabel(['Firing rate change:' char(10) 'Post/pre stimulation'])
axis tight
ylim([0 10])

nexttile([2 3])
PosBinSize=0.25*115;
for stim_type=1:3
StimTypeSort=find(Stimulation_type==stim_type);
ZapPos=cell2mat(Stimulation_Pos(StimTypeSort));
thingsToPlot=MeanFR_stimtype_Pos{stim_type}(:,2)./MeanFR_stimtype_Pos{stim_type}(:,1);
thingsToPlot(thingsToPlot==inf)=NaN;
%plot(ZapLap,MeanFR_stimtype_Pos{stim_type}(:,2)-MeanFR_stimtype_Pos{stim_type}(:,1),'.','color',cmap(stim_type,:),'markersize',20)
plot(ZapPos*2/115,thingsToPlot,'.','color',cmap_light(stim_type,:),'markersize',7)
hold all
Average_FR_Ratio=[]; Std_FR_Ratio=[];
for b=1:ceil(max(ZapPos)/PosBinSize)
    bin_list=find(ceil(ZapPos/PosBinSize)==b);
    Average_FR_Ratio(b)=mean(thingsToPlot(bin_list),'omitnan');
    Std_FR_Ratio(b)=std(thingsToPlot(bin_list),'omitnan');
end
ZapPos_bin=[1:ceil(max(ZapPos)/PosBinSize)]*PosBinSize-PosBinSize/2;
errorbar(ZapPos_bin*2/115,Average_FR_Ratio,Std_FR_Ratio,'color',cmap(stim_type,:),'linewidth',2,'marker','+','capsize',10)
end
%set(gca,'yscale','log')
plot([0 2],[1 1],'--','color','k','linewidth',2)
plot([1.6 1.6],[0 10],'--','color','r','linewidth',2)
xlabel('Stimulation position (m)')
ylabel(['Firing rate change:' char(10) 'Post/pre stimulation'])
axis tight
ylim([0 10])


%% Complex spike detection
Fs=1000; pass_frq=15;
freq_lowhigh=pass_frq/(Fs/2);
[b, a] = butter(4, freq_lowhigh, 'low');

for i=8:length(PC_Result)
    tr_highpass=[];
    for j=1:size(PC_Result{i}.normTraces,1)
tr_highpass(j,:)=PC_Result{i}.normTraces(j,:)-movprc(PC_Result{i}.normTraces(j,:),300,30);
    end
nNeuron=size(PC_Result{i}.normTraces,1);

%tr_pass = cellfun(@(x) filtfilt(b, a, x),mat2cell(tr_highpass,ones(1,nNeuron)),'UniformOutput',false);
%tr_pass=cell2mat(tr_pass);
tr_pass=tr_highpass;
[trans tr_trace]=detect_transient(tr_pass,[3 0.5],PC_Result{i}.spike);
CS_trace=[];
for n=1:nNeuron
CS_ind=find(trans(n).spike_number>2 & trans(n).mean_ISI<18);
CS_trace(n,:)=ismember(tr_trace(n,:),CS_ind);
end
PC_Result{i}.CS_trace=CS_trace;
end

%% Calculate place field for CS and SS

place_bin=150;
velocity_threshold=0.002;

for g=1:length(PC_Result)
    g
    PC_Result{g}.Lap_CS=[]; PC_Result{g}.Lap_SS=[];
    nNeuron=size(PC_Result{g}.normTraces,1);

    for n=1:nNeuron

        CS_trace = PC_Result{g}.CS_trace(n,:).*PC_Result{g}.spike(n,:);
        SS_trace = (PC_Result{g}.CS_trace(n,:) == 0).*PC_Result{g}.spike(n,:);
        [PC_Result{g}.Lap_CS(:,:,n), ~]=PlaceTrigger_average(CS_trace,place_bin,PC_Result{g}.VR,velocity_threshold,115);
        [PC_Result{g}.Lap_SS(:,:,n), ~]=PlaceTrigger_average(SS_trace,place_bin,PC_Result{g}.VR,velocity_threshold,115);
    end
end
    
save(fullfile(save_figto,['PF_Result_20240111.mat']),'PC_Result','-v7.3')

%% Show Complex spike firing map
cmap=[0.5 0.05 0.15;0.1 0.1 0.1]/5;
f1=figure(1); clf;
coarse_bin=3;
for i=1:length(PC_Result)
    [sz_fprntY, sz_fprntX, nNeurons]=size(PC_Result{i}.c_ftprnt);
    StimN=find(mean(PC_Result{i}.spike(:,find(PC_Result{i}.Blue>0)),2)>mean(PC_Result{i}.spike,2)*5)';
    UnstimN=setdiff([1:nNeurons],StimN);
    for n=[StimN UnstimN]
        nexttile([1 1])
        pos_bin=size(PC_Result{i}.Lap_CS,2);
        
        Lap_tmp_CS=repmat(PC_Result{i}.Lap_CS(:,:,n),1,3);
        Lap_tmp_CS=movmean(Lap_tmp_CS,coarse_bin,2);
        Lap_tmp_CS=(Lap_tmp_CS(:,pos_bin+1:2*pos_bin));
        

        Lap_tmp_SS=repmat(PC_Result{i}.Lap_FR(:,:,n),1,3);
        Lap_tmp_SS=movmean(Lap_tmp_SS,coarse_bin,2);
        Lap_tmp_SS=(Lap_tmp_SS(:,pos_bin+1:2*pos_bin));
        
        composite_LR=(squeeze(sum(cat(3,Lap_tmp_CS,Lap_tmp_SS).*reshape(cmap,1,1,[],3),3)));
        %imagesc(reshape(rescale(composite_LR(:)),[],pos_bin,3))
        imagesc(composite_LR)
 
        if ismember(n,UnstimN)
            title(['Unstimulated ' num2str(i) '-' num2str(n)])
        else
            title([num2str(i) '-' num2str(n)])
        end
    end
end

f2=figure(2); clf;
for i=1:length(PC_Result)
    
    [sz_fprntY, sz_fprntX, nNeurons]=size(PC_Result{i}.c_ftprnt);
    StimN=find(mean(PC_Result{i}.spike(:,find(PC_Result{i}.Blue>0)),2)>mean(PC_Result{i}.spike,2)*5)';
    UnstimN=setdiff([1:nNeurons],StimN);
    for n=[StimN UnstimN]
    nexttile([1 1])
    plot(mean(PC_Result{i}.Lap_CS(:,:,n),1,'omitnan'),'color',cmap(1,:)*10); hold all
    plot(mean(PC_Result{i}.Lap_SS(:,:,n),1,'omitnan'),'color',cmap(2,:)*10)
    if ismember(n,UnstimN)
            title(['Unstimulated ' num2str(i) '-' num2str(n)])
        else
            title([num2str(i) '-' num2str(n)])
        end
    end
end

f3=figure(3); clf;
NtoPlot=[3 3;4 3;5 1;6 1;9 1;18 1;15 1];
tiledlayout(4,4)
for n=1:size(NtoPlot,1)

    Lap_tmp_CS=repmat(PC_Result{NtoPlot(n,1)}.Lap_CS(:,:,NtoPlot(n,2)),1,3);
        Lap_tmp_CS=movmean(Lap_tmp_CS,coarse_bin,2);
        Lap_tmp_CS=(Lap_tmp_CS(:,pos_bin+1:2*pos_bin));
        

        Lap_tmp_FR=repmat(PC_Result{NtoPlot(n,1)}.Lap_FR(:,:,NtoPlot(n,2)),1,3);
        Lap_tmp_FR=movmean(Lap_tmp_FR,coarse_bin,2);
        Lap_tmp_FR=(Lap_tmp_FR(:,pos_bin+1:2*pos_bin));
        
        composite_LR=(squeeze(sum(cat(3,Lap_tmp_CS,Lap_tmp_FR).*reshape(cmap,1,1,[],3),3)));

    nexttile(n+(n>4)*4,[1 1])
    imagesc(composite_LR)
    set(gca,'xtick',[0:75:150],'xticklabel',[0:75:150]*2/150)
    xlabel('Position (m)')
    ylabel('VR trial')
    
    nexttile(n+(n>4)*4+4,[1 1])
    % plot(mean(PC_Result{NtoPlot(n,1)}.Lap_CS(:,:,NtoPlot(n,2)),1,'omitnan'),'color',cmap(1,:)*10); hold all
    % plot(mean(PC_Result{NtoPlot(n,1)}.Lap_FR(:,:,NtoPlot(n,2)),1,'omitnan'),'color',cmap(2,:)*10)
    plot(mean(Lap_tmp_CS,1,'omitnan'),'color',cmap(1,:)*10); hold all
    plot(mean(Lap_tmp_FR,1,'omitnan'),'color',cmap(2,:)*10)
    %plot(mean(Lap_tmp_FR-Lap_tmp_CS,1,'omitnan'),'color',[0 0.4 0.8])
    set(gca,'xtick',[0:75:150],'xticklabel',[0:75:150]*2/150)
    xlabel('Position (m)')
    ylabel('Firing rate (Hz)')

   % nexttile(size(NtoPlot,1)*2+n,[1 1])
   % errorbar_shade([1:pos_bin],mean(Lap_tmp_CS./Lap_tmp_FR,1,'omitnan'),std(Lap_tmp_CS./Lap_tmp_FR,0,1,'omitnan'),cmap(1,:)*10)
end
legend({'Complex spike','All spike'})

