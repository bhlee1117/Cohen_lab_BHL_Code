
clear
clc;
cd '/Volumes/BHL18TB_D1/Arranged_Data/Prism_OptopatchResult';
[~, ~, raw] = xlsread(['/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/' ...
    'Prism_OptopatchData_Arrangement.xlsx'], 'Sheet1', 'B5:K104');

save_to='/Volumes/BHL18TB_D1/Arranged_Data/Prism_OptopatchResult';
fpath=raw(:,1);
Mouse=cell2mat(raw(:,2));
NeuronInd=cell2mat(raw(:,5));
CamType=raw(:,3);
StructureData=raw(:,10);
%%
figure;
[~, unqInd] = unique([Mouse NeuronInd] ,'row');
for i=unqInd'
    nexttile([1 1])
    load(fullfile(fpath{i},"output_data.mat"))
    switch char(CamType(i))
        case 'flash'
            sz=double(Device_Data{1, 4}.ROI([2 4]));
        case 'fusion'
            sz=double(Device_Data{1, 3}.ROI([2 4]));
    end
    ref_time=[2000:4000];
    mov_test=double(readBinMov_times([fpath{i} '/frames1.bin'],sz(2),sz(1),[ref_time(1):ref_time(end)]));
    imshow2(mean(mov_test,3),[])
    title(['Mouse #' num2str(Mouse(i)) '-Neuron#' num2str(NeuronInd(i))])
end

%% MC

for i=[87:90]%length(fpath)

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
for i=unqInd([27])'
    nexttile([1 1])
    load(fullfile(fpath{i},"output_data.mat"))
    disp(fpath{i})

    switch char(CamType(i))
        case 'flash'
            sz=double(Device_Data{1, 4}.ROI([2 4]));
        case 'fusion'
            sz=double(Device_Data{1, 3}.ROI([2 4]));
    end

    ref_time=[2000:4000];
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
    save(fullfile(fpath{i},'Result.mat'),'Result','-v7.3')

    SameCellInd=find(Mouse==Mouse(i) & NeuronInd==NeuronInd(i));
    for j=SameCellInd'
        Result.ROIpoly=ROIpoly;

        mov_mc=double(readBinMov_times([fpath{j} '/mc_ShutterReg' num2str(1,'%02d') '.bin'],sz(2),sz(1),[3000:4000]));
        [offsetY offsetX]=calculate_shift(Result.ref_im,mean(mov_mc,3));
        ROIpoly_shift=cellfun(@(x) x+[offsetX offsetY],Result.ROIpoly,'UniformOutput',false);
        Result.ROIpoly=ROIpoly_shift;
        save([fpath{j} '/Result.mat'],'Result','-v7.3')
    end
end

%% Set Footprint
bound=6;
for i=95:100%length(fpath)
    load([fpath{i} '/Result.mat'])
    load(fullfile(fpath{i},"output_data.mat"))
    switch char(CamType(i))
        case 'flash'
            sz=double(Device_Data{1, 4}.ROI([2 4]));
        case 'fusion'
            sz=double(Device_Data{1, 3}.ROI([2 4]));
    end 
    
    ref_time=[2000:4000];
    load(fullfile(fpath{i},['/mcTrace' num2str(1,'%02d') '.mat']));

    mov_mc=double(readBinMov_times([fpath{i} '/mc_ShutterReg' num2str(1,'%02d') '.bin'],sz(2),sz(1),ref_time));

    Result.ref_im=mean(mov_mc,3);
    Result.mc=mcTrace.xymean;
    mov_res= mov_mc-mean(mov_mc,3);
    %mcTrace.xymean=movmean(mcTrace.xymean,3,2);
    % mov_res = SeeResiduals(mov_res,mcTrace.xymean(ref_time,:));
    % mov_res = SeeResiduals(mov_res,mcTrace.xymean(ref_time,:).^2);
    % mov_res = SeeResiduals(mov_res,mcTrace.xymean(ref_time,1).*mcTrace.xymean(ref_time,2));

  n_comp=6;
    mov_filt=imgaussfilt3(mov_res,[2 2 0.1]);
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
    figure; clf;
    show_im=Result.ref_im;
    show_im(show_im<15)=median(show_im(:));
    show_im=show_im-imgaussfilt(show_im,15);
    show_im(show_im<prctile(show_im(:),20))=prctile(show_im(:),20);
    show_footprnt(Result.ftprnt,show_im)
    save([fpath{i} '/Result.mat'],'Result','-v7.3')
    saveas(gca,[char(fpath{i}) '/' 'ftprnt.fig'])
end


%% Signal extraction
bound=5;

for i=[85:88]%length(fpath)]

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
    Result.mc=mcTrace.xymean;
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
    bkg(1,:)=movmedian(get_blueoffTrace(squeeze(mean(mov_mc,[1 2])),[Result.Blue],30),3000,'omitnan');
    mov_res = SeeResiduals(mov_res,Result.mc);
    mov_res = SeeResiduals(mov_res,Result.mc.^2);
    mov_res = SeeResiduals(mov_res,Result.mc(:,1).*Result.mc(:,end));
    mov_res= SeeResiduals(mov_res,bkg,1);


    Result.traces=[-(tovec(mov_res)'*tovec(Result.ftprnt))'];
    Result.im_corr=[sum(mov_mc_vec.*ref_im_vec,1)/(size(mov_mc_vec,1)-1)];  %image correlation
    save([fpath{i} '/Result.mat'],'Result','-v7.3')
end

%%
for i=[90]%length(fpath)]
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

