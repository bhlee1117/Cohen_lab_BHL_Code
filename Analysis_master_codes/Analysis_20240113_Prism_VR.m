 clear
clc;
cd '/Volumes/BHL_WD18TB/PP72_PlaceCellResults';
[~, ~, raw] = xlsread(['/Volumes/BHL_WD18TB/' ...
    'PrismPCdata_Arrangement.xlsx'], 'Sheet1', 'C5:N19');

% [~, ~, NeuronsToUse]=xlsread(['/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/' ...
%     'PlaceCellData_Arrangement.xlsx'], 'Sheet1', 'L8:M46');
% 
% NeuronsToUse=cellfun(@(x) (str2num(num2str(x))),NeuronsToUse,'UniformOutput',false);
ref_ROI=cellfun(@(x) (str2num(num2str(x))),raw(:,11),'UniformOutput',false);
fpath=raw(:,1)';
StructureData=raw(:,10);
save_figto='/Volumes/BHL_WD18TB/PP72_PlaceCellResults';
place_bin=150; time_segment=15000; overlap=200;
%% Motion correction

for f=14:length(fpath)
load(fullfile(fpath{f},"output_data.mat"))
sz=double(Device_Data{1, 3}.ROI([2 4]));
ref_time=[9000:10000];

frm_end=max(Device_Data{1, 2}.Counter_Inputs(1, 1).data);
f_seg=[[1:time_segment:frm_end] frm_end+1]; f_seg(2:end)=f_seg(2:end)-1;

if length(ref_time)>2000
    mov_test=double(readBinMov_times([fpath{f} '/frames1.bin'],sz(2),sz(1),[ref_time(1):ref_time(1)+2000]));
else
    mov_test=double(readBinMov_times([fpath{f} '/frames1.bin'],sz(2),sz(1),[ref_time(1):ref_time(end)]));
end
mov_test=rollingShutter_correction(mov_test,Device_Data{1, 3}.exposuretime,'fusion');
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

    mov=rollingShutter_correction(mov,Device_Data{1, 3}.exposuretime,'fusion');
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

for f=15%:length(fpath)
disp(fpath{f}); Result=[];
DAQ_rate=0.000005;
load([fpath{f} '/output_data.mat'])
sz=double(Device_Data{1, 3}.ROI([2 4]));
ref_time=[2000:3000];
mov_test=double(readBinMov_times([fpath{f} '/mc_ShutterReg' num2str(5,'%02d') '.bin'],sz(2),sz(1),ref_time));
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

load(fullfile(fpath{f},'mcTrace05.mat'));
frm_end=max(Device_Data{1, 2}.Counter_Inputs(1, 1).data);
Result.Blue=Device_Data{1, 2}.buffered_tasks(1, 1).channels(1, 2).data;
Result.Reward=Device_Data{1, 2}.buffered_tasks(1, 3).channels.data;
CamCounter=Device_Data{1, 2}.Counter_Inputs(1, 1).data;
CamTrigger=find(CamCounter(2:end)-CamCounter(1:end-1));

Result.Blue=Result.Blue(CamTrigger); Result.Reward=Result.Reward(CamTrigger);
frm_rate=double((CamTrigger(2)-CamTrigger(1))*DAQ_rate);
mov_mc=double(readBinMov([fpath{f} '/mc_ShutterReg' num2str(5,'%02d') '.bin'],sz(2),sz(1)));
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

n_comp=5;
mov_filt=imgaussfilt3(mov_res,[3 3 0.1]);
movVec=tovec(mov_filt);
Npoly=size(Result.ROIpoly,1);
ftprnt = zeros(size(mov_filt,1)*size(mov_filt,2),Npoly);
clear mask
for p=1:Npoly %each ROIs
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
    coeff = mat2gray(mean(abs(V(:,1:n_comp)).*D(1:n_comp)',2));
    ftprnt(pixelList,p)=coeff;
end

Result.ftprnt=toimg(ftprnt,sz(2),sz(1));
figure; clf;
imshow2(squeeze(sum(toimg(Result.ftprnt,sz(2),sz(1)).*reshape(jet(Npoly),1,1,[],3),3)),[]);
save(fullfile(fpath{f},'PC_Result.mat'),'Result','fpath','-v7.3')

end
%% Signal extraction

for f=14:length(fpath)
load(fullfile(fpath{f},'PC_Result.mat')); 
load([fpath{f} '/output_data.mat'])
sz=double(Device_Data{1, 3}.ROI([2 4])); blueDMDcontour=[];
CamCounter=Device_Data{1, 2}.Counter_Inputs(1, 1).data;
CamTrigger=find(CamCounter(2:end)-CamCounter(1:end-1));
Rfixed = imref2d(repmat(Device_Data{1, 3}.virtualSensorSize,1,2));
inverseTform = invert(Device_Data{1, 6}.tform);
revertedImage = imwarp(double(Device_Data{6}.pattern_stack(:,:,find(sum(Device_Data{6}.pattern_stack,[1 2])>0))), inverseTform,'OutputView',Rfixed);
[blueDMDimg]=imcrop_3d(revertedImage,double(Device_Data{1, 3}.ROI([1 3 2 4]))+[0 0 -1 -1]);
for d=1:size(blueDMDimg,3)
blueDMDcontour{d}=bwboundaries(blueDMDimg(:,:,d));
end

Result.blueDMDimg=blueDMDimg;
Result.blueDMDcontour=blueDMDcontour;
Result.traces=[];
Result.traces_res=[];
Result.mcTrace=[];
Result.im_corr=[];
bound=5;
ref_im_vec=tovec(Result.ref_im(bound:end-bound,bound:end-bound));
frm_end=max(Device_Data{1, 2}.Counter_Inputs(1, 1).data);
f_seg=[[1:time_segment:frm_end] frm_end+1]; f_seg(2:end)=f_seg(2:end)-1;
take_window=repmat([1 time_segment],length(f_seg)-1,1);
take_window(2:end,1)=take_window(2:end,1)+overlap; take_window(1:end-1,2)=take_window(1:end-1,2)+overlap;
take_window(end)=mod(f_seg(end),time_segment);

for j=1:length(f_seg)-1
    j
    mov_mc=double(readBinMov([fpath{f} '/mc_ShutterReg' num2str(j,'%02d') '.bin'],sz(2),sz(1)));
    load([fpath{f} '/mcTrace' num2str(j,'%02d') '.mat']);

        mov_mc=mov_mc(:,:,[take_window(j,1):take_window(j,2)]);
        mc=mcTrace.xymean([take_window(j,1):take_window(j,2)],:);
 
    mov_mc_vec=tovec(mov_mc(bound:end-bound,bound:end-bound,:)); 
    mov_mc_vec=(mov_mc_vec-mean(mov_mc_vec,1))./std(mov_mc_vec,0,1);

    mov_res= mov_mc-mean(mov_mc,3);
    bkg = zeros(2, size(mov_mc,3));
    bkg(1,:) = linspace(-1, 1, size(mov_mc,3));  % linear term
    bkg(2,:) = linspace(-1, 1, size(mov_mc,3)).^2;  % quadratic term
    mov_res = SeeResiduals(mov_res,mc);
    mov_res = SeeResiduals(mov_res,mc.^2);
    mov_res = SeeResiduals(mov_res,mc(:,1).*mc(:,end));
    mov_res= SeeResiduals(mov_res,bkg,1);

    Result.traces=[Result.traces -(tovec(mov_res)'*tovec(Result.ftprnt))'];
    Result.mcTrace=[Result.mcTrace; mc];
    Result.im_corr=[Result.im_corr sum(mov_mc_vec.*ref_im_vec,1)/(size(mov_mc_vec,1)-1)];  %image correlation

end

Result.traces=Result.traces(:,1:length(CamTrigger)-1);
Result.mcTrace=Result.mcTrace(1:length(CamTrigger)-1,:);

save(fullfile(fpath{f},'PC_Result.mat'),'Result','fpath','-v7.3')

end
%% Load Virmen data

for f=13:length(fpath)

load(fullfile(fpath{f},'PC_Result.mat'))
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
t_DAQ=CamTrigger/DAQ_rate;
t_VR = datetime(datetime(VRdata(1,:),'ConvertFrom','datenum'), 'InputFormat', 'yyyy-MM-dd HH:mm:ss.SSS');
t_VR= t_VR-t_VR(1);
t_VR= milliseconds(t_VR)/1000;
t_VR= t_VR*t_DAQ(end)/t_VR(end);
VRdata(1,:)=t_VR;
[Virmen_data_int vel_trace]=virmen_interpolate(VRdata,115,t_DAQ);
Virmen_data_int(end+1,:)=vel_trace;

Result.VR=Virmen_data_int;
Result.Blue=Result.Blue.*Result.VR(12,:);
Result.VR=Result.VR(:,1:end-1);

save(fullfile(fpath{f},'PC_Result.mat'),'Result','fpath','-v7.3')
end
%%
load(fullfile(save_figto,'Result_PC_Prism_20240219.mat'))
for i=13:length(fpath)
    Result_tmp=load(fullfile(fpath{i},'PC_Result.mat'));
    PC_Result{i}=Result_tmp.Result;
end

%% Clean up and norm
exclude_frq=[241.7 242]; %monitor
%exclude_frq2=[483.5 484]; %monitor
exclude_frq2=[25 65.7]; %motion
    time_bin=15000; Fs=1000; ref_trace=[2 4 4 2 3 2 1 1 1 3 1 1 1 1 1]; %2nd trunk is the reliable trace

for f=13:length(fpath)
    nTime=size(PC_Result{f}.traces,2);
    nROI=size(PC_Result{f}.traces,1);
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


clear traces_res_filtered noise noise_intp norm_trace sp_height SpHeight_intp sp_time tr_mc_imcorr tr_mc tr
    tN=[1:time_bin:nTime]; tN=[tN nTime];
    sp_time=zeros(nROI,nTime);
    sp_height=zeros(nROI,nTime);

mcTrace=squeeze(PC_Result{f}.mcTrace)';

for n=1:size(PC_Result{f}.traces,1)
    tr=PC_Result{f}.traces(n,1:nTime);

    % regress out motion frequency
    for t=1:length(tN)-1
      tr(tN(t):tN(t+1))=squeeze(SeeResiduals(reshape(tr(tN(t):tN(t+1)),1,1,[]),mcTrace(:,(tN(t):tN(t+1)))));
                tr(tN(t):tN(t+1))=squeeze(SeeResiduals(reshape(tr(tN(t):tN(t+1)),1,1,[]),mcTrace(:,(tN(t):tN(t+1))).^2));
                tr(tN(t):tN(t+1))=squeeze(SeeResiduals(reshape(tr(tN(t):tN(t+1)),1,1,[]),mcTrace(1,(tN(t):tN(t+1))).*mcTrace(end,(tN(t):tN(t+1)))));

                tr_mc(tN(t):tN(t+1))=tr(tN(t):tN(t+1));

                imcorr_tmp=PC_Result{f}.im_corr(tN(t):tN(t+1));
                imcorr_tmp=movmean(imcorr_tmp,5,'omitnan');
    end

    % regress out motion frequency
    traces_res_filtered(n,:) = filtfilt(b, a, tr);
    tr_mc_bpMonitor=traces_res_filtered(n,:);
    %traces_res_filtered(n,:) = filtfilt(b2, a2, traces_res_filtered(n,:));
    tr_mc_bpMonitorMotion=traces_res_filtered(n,:);

    norm_trace(n,:)=traces_res_filtered(n,:);%-movmedian(traces_res_filtered(n,:),500,2);

    

    if n==ref_trace(f)
    for t=1:length(tN)-1
        tr_tmp=norm_trace(n,tN(t):tN(t+1));
            tr_tmp=tr_tmp-movmedian(tr_tmp,300);
            noise(t,n)=get_threshold(tr_tmp,1);
            tr_tmp_norm=tr_tmp./noise(t,n);
            [sp_temp, pks, prom]=find_spike_bh(tr_tmp_norm,5,4);
            sp_time(n,tN(t):tN(t+1))=sp_temp;
    end

    

    t_fit=find(sp_time(n,:));
    sp_height(n,t_fit)=norm_trace(n,t_fit);
    tx=tN(1:end-1)+time_bin/2;
    %noise_intp(n,:)=movmean(interp1(tx,noise(:,n),[1:size(PC_Result{f}{i}.traces,2)],'linear','extrap'),10000);

    
     [Ny_fit t_consts coeffY]  = expfitDM_2(tx(~isnan(noise(1:end-1,n)))',noise(~isnan(noise(1:end-1,n)),n),[1:nTime]',10^7);
     noise_intp=Ny_fit;
    %[y_fit t_consts coeffY]  = expfitDM_2(tx(~isnan(sp_height(1:end-1,n)))',sp_height(~isnan(sp_height(1:end-1,n)),n),[1:size(Result{i}.traces,2)]',10^7);
    [y_fit t_consts coeffY]  = expfitDM_2(t_fit',sp_height(n,t_fit)',[1:nTime]',10^7);
    SpHeight_intp=y_fit;


    figure; clf;
    ax1=nexttile([1 1]);
    plot(rescale2([PC_Result{f}.traces(n,1:nTime);tr_mc;tr_mc_bpMonitorMotion;mcTrace(1,1:nTime)],2)')%+[1:4])
    ax2=nexttile([1 1]);
    plot([1:nTime],norm_trace(n,:))
    hold all
    plot(find(sp_time(n,:)),norm_trace(n,find(sp_time(n,:))),'r.')
    plot([1:nTime],y_fit([1:nTime]),'k')
    plot([1:nTime],Ny_fit([1:nTime]),'g')

    end
end

   norm_trace=norm_trace./SpHeight_intp';        
    %norm_trace=norm_trace;%./(SpHeight_intp./SpHeight_intp(:,1));        
    PC_Result{f}.normTraces=norm_trace./get_threshold(norm_trace,1);
    PC_Result{f}.spike=find_spike_bh(PC_Result{f}.normTraces-movmedian(PC_Result{f}.normTraces,300,2),5,3);

    
for n=1:size(PC_Result{f}.traces,1)
PC_Result{f}.subThreshold(n,:) = filtfilt(b3, a3, PC_Result{f}.normTraces(n,:));
PC_Result{f}.theta(n,:) = filtfilt(b4, a4, PC_Result{f}.normTraces(n,:));
end

end
% 
 %save(fullfile(save_figto,'Result_PC_Prism_20240225.mat'),'PC_Result','fpath','-v7.3')
% 
%save(fullfile(fpath,'Result_20230928.mat'),'Result','fpath','-v7.3')

%% Structure segmentation
Struct_valid=find(1-cell2mat(cellfun(@(x) sum(isnan(x)), StructureData, 'UniformOutput', false)));

for i=Struct_valid(3)'
StructureStack=mat2gray(double(tiffreadVolume(StructureData{i})));
StructureStack(StructureStack==0)=median(StructureStack(:));
%StructureStack_med=medfilt2_mov(StructureStack,[15 15]);
StructureStack_Gauss=imgaussfilt3(StructureStack,[7 7 0.1]);
%StructureStack_med(StructureStack_med==0)=median(StructureStack_med(:));
%StructureStack=(StructureStack-StructureStack_med)./StructureStack_med;
StructureStack_filt=(StructureStack-StructureStack_Gauss);
StructureStack_filt=mat2gray(StructureStack_filt);
StructureStack_bin=[]; level=[];
level = graythresh(StructureStack_filt);
StructureStack_bin=StructureStack_filt>level*0.94;
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


rot_ang=0;
Structure_ref=(imrotate(StructureStack_final,rot_ang));
ref_img=PC_Result{i}.ref_im; ref_img(ref_img<prctile(ref_img(:),20))=median(ref_img(:)); ref_img=ref_img-prctile(ref_img(:),20);
[RegImg,tformReg]=imReg(ref_img,max(Structure_ref,[],3));
saveastiff(uint16(mat2gray(Structure_ref)*255), [fpath{i} 'Structure.tiff']);

    PC_Result{i}.Structure=max(Structure_ref,[],3);
    PC_Result{i}.Structure_bin=max(imrotate(dendrite_bin,rot_ang),[],3);
    PC_Result{i}.tform=tformReg;
end

q

%% Place averaged
place_bin=150;
velocity_threshold=-0.002;
%ref_ROI={[2],[3 4],[2 3],[2],[3],[2],[1],[1],[1],[3],[1],[1]};
for g=13:15%1:length(fpath)
    g

    PC_Result{g}.Lap_FR=[]; PC_Result{g}.Lap_sub=[]; PC_Result{g}.Lap_theta=[];
    PC_Result{g}.Lap_F=[]; PC_Result{g}.Lap_V=[];
    
    ref_trace=mean(PC_Result{g}.normTraces(ref_ROI{g},:),1);
    PC_Result{g}.ref_spike=find_spike_bh(ref_trace-movmedian(ref_trace,300),5,3);

        [PC_Result{g}.Lap_FR PC_Result{g}.Lap_V]=PlaceTrigger_average(PC_Result{g}.ref_spike,place_bin,PC_Result{g}.VR,velocity_threshold,115);
for n=1:size(PC_Result{g}.normTraces,1)
        [PC_Result{g}.Lap_F(:,:,n) PC_Result{g}.Lap_V]=PlaceTrigger_average(PC_Result{g}.normTraces(n,:),place_bin,PC_Result{g}.VR,velocity_threshold,115);
end
        ref_subthreshold=movmean(mean(PC_Result{g}.subThreshold(ref_ROI{g},:),1),1000);
        
        [PC_Result{g}.Lap_sub PC_Result{g}.Lap_V]=PlaceTrigger_average(ref_subthreshold,place_bin,PC_Result{g}.VR,velocity_threshold,115);
        %[PC_Result{g}.Lap_theta(:,:,n) PC_Result{g}.Lap_V]=PlaceTrigger_average(abs(PC_Result{g}.theta(n,:)),place_bin,PC_Result{g}.VR,velocity_threshold,115);
    
    Result=PC_Result{g};
    save(fullfile(fpath{g},'PC_Result.mat'),'Result','fpath','-v7.3')
end


%% Show place fields
save_figto='/Volumes/BHL_WD18TB/PP72_PlaceCellResults';

f1=figure(1); clf;
for i=1:length(PC_Result)
    
        nexttile([1 1])
        pos_bin=size(PC_Result{i}.Lap_FR,2);
        Lap_tmp=repmat(PC_Result{i}.Lap_FR,1,3);
        Lap_tmp=movmean(Lap_tmp,7,2);
        Lap_tmp=Lap_tmp(:,pos_bin+1:2*pos_bin);
        imagesc(Lap_tmp)
        colormap('turbo')
            title([num2str(i)])
end

set(f1, 'Position', [100, 100, 800, 400]);
saveas(f1,fullfile(save_figto ,['PF_map' '.fig']))
%print(f1, fullfile(save_figto ,['PF_map' '.jpg']),'-djpeg', ['-r', num2str(400)]);

%% realign by footprint's location from soma (1st footprint)
for f=1:length(PC_Result)
centroid_ftprnt = get_coord(PC_Result{f}.ftprnt);
dist_centroid = distance_mat(centroid_ftprnt(1,:),centroid_ftprnt);
[~, PC_Result{f}.dist_order]=sort(dist_centroid,'ascend');
end
%% Spike triggered average
nFront=30; nBack=50;
nTau=[-nFront:nBack];
tiledlayout(6,length(PC_Result))

for i=13:length(PC_Result)
   ax1=nexttile(i,[1 1]); colormap(ax1,'gray');
   imshow2(PC_Result{i}.ref_im,[])
   s_tau=find(PC_Result{i}.ref_spike)'+nTau;
   nFrame=size(PC_Result{i}.normTraces,2); nROI=size(PC_Result{i}.normTraces,1);
   omit_sp=find(min(s_tau,[],2)<1 | max(s_tau,[],2)>nFrame);
   s_tau(omit_sp,:)=[];
   PC_Result{i}.STA=[];
   PC_Result{i}.STA=reshape(PC_Result{i}.normTraces(:,s_tau),nROI,[],nFront+nBack+1);
   ax2=nexttile(length(PC_Result)+i,[5 1]);
   imagesc(rescale2(squeeze(mean(PC_Result{i}.STA(PC_Result{i}.dist_order,:,:),2)),2)); hold all
   colormap(ax2,'turbo')
   l= plot(-rescale2(squeeze(mean(PC_Result{i}.STA(PC_Result{i}.dist_order,:,:),2)),2)'+[1:nROI]+0.5,'LineWidth',2);
   arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(turbo(nROI),2))
end

%%

%% Spike triggered average movie
nFront=5; nBack=15;
nTau=[-nFront:nBack];
for f=13:length(PC_Result)
    f
    load([fpath{f} '/output_data.mat'])
    sz=double(Device_Data{1, 3}.ROI([2 4]));
    frm_end=max(Device_Data{1, 2}.Counter_Inputs(1, 1).data);
    f_seg=[[1:time_segment:frm_end] frm_end+1]; f_seg(2:end)=f_seg(2:end)-1;

    fileInd=ceil(find(PC_Result{f}.ref_spike)/time_segment);
    frameInd=mod(find(PC_Result{f}.ref_spike),time_segment);
    %line up movie segments at simple spike
    STA_movie_align=[];
    for j=unique(fileInd)
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
        mov_res = SeeResiduals(mov_res,mc(:,1).*mc(:,3));
        mov_res= SeeResiduals(mov_res,bkg,1);

        for k=find(fileInd==j)
            try
                STA_movie_align = cat(3,STA_movie_align, tovec(mov_res(:,:,frameInd(k)+nTau)));
            end
        end
    end
save(fullfile(fpath{f},'STA_movie_align.mat'),'STA_movie_align','-v7.3');
PC_Result{f}.STA_movie=toimg(squeeze(mean(STA_movie_align,3)),sz(2),sz(1));
save(fullfile(save_figto,'Result_PC_Prism_20240225.mat'),'PC_Result','fpath','-v7.3')
end

%% Write STA movie
STAmov=[]; MaxSize=500;
for f=1:length(PC_Result)
    movTmp=zeros(size(PC_Result{f}.ref_im,1),MaxSize,size(PC_Result{f}.STA_movie,3));
    bd=(MaxSize-size(PC_Result{f}.ref_im,2))/2;
    movTmp(:,bd+1:MaxSize-bd,:)=PC_Result{f}.STA_movie;
    STAmov=[STAmov; movTmp];
end
writeMov([save_figto,'/STA_movie_total'],-STAmov,[],10,1,[0 30])

%% Generate SNAPT movie
for i=5%length(PC_Result)
mask=max(PC_Result{i}.Structure_bin,[],3)>0.01;
maskSTA=max(-PC_Result{i}.STA_movie,[],3)./PC_Result{i}.ref_im>0.05;
StrImg=max(PC_Result{i}.Structure,[],3);
STAmovie=mat2gray(-PC_Result{i}.STA_movie);
STAmovie=STAmovie-prctile(STAmovie,10,3);
STAmovie=mat2gray(STAmovie(:,:,1:21));
tformReg=PC_Result{i}.tform;
[PC_Result{i}.SNAPT PC_Result{i}.dtimg]=generate_SNAPTmov(mat2gray(STAmovie),mask,StrImg,tformReg);

dtimg_Reg=imwarp(PC_Result{i}.dtimg, tformReg, 'OutputView', imref2d(size(StrImg)));
dtimg_Reg=dtimg_Reg.*mask; dtimg_Reg(dtimg_Reg==0)=NaN;
imshow2(dtimg_Reg-prctile(dtimg_Reg(:),5), [0 3]); title('Timing')
hold all
Maskboundary = cell2mat(bwboundaries(mask));
plot(Maskboundary(:,2),Maskboundary(:,1),'r.','markersize',2); 
colormap('turbo')


[yR xR zR]=size(PC_Result{i}.Structure);
figure(20); clf;
v = VideoWriter([fpath{i} '/SNAPT_movie'],'Uncompressed AVI');

open(v);
subframeT = 0.025; % ms
initialT = -2; % ms
finalT = 2; % ms
times = initialT:subframeT:finalT;

for j = 1:length(times)
    clf;
    %set(gca,'units','pixels','position',[200 0 1000 800])
    imshow(PC_Result{i}.SNAPT(:,:,:,j),[])
    pbaspect([size(double(PC_Result{i}.SNAPT(:,:,:,j)),2) size(double(PC_Result{i}.SNAPT(:,:,:,j)),1) 1]),colormap(gray)
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

%% detect complex spikes

%ref_ROI={[2],[3 4],[3]};
figure(2); clf;
tiledlayout(6,1)
ax1=[];
for g=1:length(fpath)
    g
    ref_trace=mean(PC_Result{g}.normTraces(ref_ROI{g},:),1);
    spike=find_spike_bh(ref_trace-movmedian(ref_trace,300),4,3);
    t=[1:length(ref_trace)];
    spike_time=find(spike); spike_height=ref_trace(spike_time);

tr_sub=get_subthreshold(ref_trace,spike,5,15);
show_t=[350000:550000];
ax1=[ax1 nexttile([1 1])];
plot(t(show_t),ref_trace(show_t)); hold all
plot(t(show_t),tr_sub(show_t));
plot(t(show_t),zeros(length(show_t),1));
legend({'F trace','Subthreshold'})

[trans tr_trace]=detect_transient(tr_sub,[4 0.8],spike);
CS_ind=find(trans.spike_number>2 & trans.mean_ISI<30);
PC_Result{g}.CS_trace=ismember(tr_trace,CS_ind);

ax1=[ax1 nexttile([1 1])];
plot(t(show_t),ref_trace(show_t))
hold all
tr_nan=ref_trace; tr_nan(PC_Result{g}.CS_trace==0)=NaN;
plot(t(show_t),tr_nan(show_t),'r')
legend({'F trace','Complex Spikes'})

clear CS_time
CS_trace_spike=spike.*bwlabel(PC_Result{g}.CS_trace);
for CS=1:max(CS_trace_spike)
PC_Result{g}.CS_time(CS)= find(CS_trace_spike==CS,1,'first');
end

%fileInd=ceil(CS_time/time_segment);
%frameInd=mod(CS_time,time_segment);
end
linkaxes(ax1,'xy')


%%
f=1;
figure; clf;
nBack = -50;
nFront = 100;
tau = nBack:nFront;
fgauss=fspecial('gaussian',[1 3],1);
%renormTrace=rescale2(PC_Result{f}.normTraces(PC_Result{f}.dist_order,:),2);
renormTrace=rescale2(PC_Result{f}.normTraces,2);
for i=41:50%length(PC_Result{f}.CS_time)
nexttile([1 1])
nROI=size(PC_Result{f}.normTraces,1);
trace_seg=renormTrace(:,PC_Result{f}.CS_time(i)+tau);
%trace_seg=rescale2(trace_seg,2);
trace_seg_filter=conv2(trace_seg,fgauss,'same');
imagesc(trace_seg_filter);
hold all
colormap(Aurora)
caxis([0.25 0.9])
line=plot(-trace_seg'+[1:nROI]+0.5);
arrayfun(@(l,c) set(l,'Color',c{:}),line,num2cell(jet(nROI),2))
title(['Complex spike #' num2str(i)])
end

%% Detect dSpike without AP

nFront=100; nBack=50;
nTau=[-nFront:nBack];

%ref_ROI={[2 3], [2 3], [2 3],[1 2 3],[1 2 3],[1 2 3]};
%DOI={[5 6 8 9 10],[4 5 6 7 8 9 10],[4:9],[4:8 10 11 14 15 16],[4 5 6 16 17 14],[4 5 6 12 13 14]};
nROIs=cellfun(@(x) size(x.normTraces,1),PC_Result);
scale=5;
dSP_Mat=[]; sp_dSP_time=[];
motion_frq=[20 40];

for i=4%:length(PC_Result)

    DOI{i}=setdiff([2:nROIs(i)],ref_ROI{i});
    [wvletTr wvletF] = cwt(PC_Result{i}.im_corr,1000); % Compute the CWT
    motionArtTrace=zscore(sum(abs(wvletTr(find(wvletF>(motion_frq(1)) & wvletF<motion_frq(2)),:)),1));

    tr=PC_Result{i}.normTraces(ref_ROI{i},:);
    sp_ref=max(find_spike_bh(tr-movmedian(tr,100,2),5,3),[],1);

    nROI=size(PC_Result{i}.normTraces,1);
    %sp=PC_Result{i}.spike;
    sp=find_spike_bh(PC_Result{i}.normTraces-movmedian(PC_Result{i}.normTraces,100,2),5,3);
    sp(:,motionArtTrace>5)=0; sp_ref(:,motionArtTrace>5)=0;

    [~, shift]=max(reshape(PC_Result{i}.normTraces(1,find(sp_ref)+[-1:0]'),2,[]),[],1);
    shift=shift-2;
    sp_timeROI = find(sp_ref>0)+shift;
    %sp_timeROI = find(sp_ref>0);

    rejectSP=find((sp_timeROI+nTau(1))<1 | (sp_timeROI+nTau(end))>size(tr,2));
    sp_timeROI(rejectSP)=[];

    sp_dSp=sp;
    sp_dSp(:,sp_timeROI'+[0:3])=0;
    

    somSP_mat{i} = reshape(PC_Result{i}.normTraces(:,sp_timeROI' + nTau),nROI,[],length(nTau));
    SPTA= squeeze(mean(somSP_mat{i},2));
    SPTA=SPTA-median(SPTA(:,[1:50]),2); F_Ref=mean(SPTA(:,[nFront+15:nFront+25]),2);
    SPTA=SPTA./F_Ref;
    
    nValidDend=sum(sum(sp_dSp(DOI{i},:),2)>1);
    figure; clf; tiledlayout(5, nValidDend+1)
    nexttile([1 nValidDend+1])
    show_footprnt_contour(PC_Result{i}.ftprnt,PC_Result{i}.ref_im,turbo(nROIs(i)))
    cmapAll=turbo(size(sp,1));

    for d=DOI{i}
    sp_dSP_time{i,d}=find(sum(sp_dSp(d,:),1)>0);
    rejectdSP=find((sp_dSP_time{i,d}+nTau(1))<1 | (sp_dSP_time{i,d}+nTau(end))>size(tr,2));
    sp_dSP_time{i,d}(rejectdSP)=[];

    dSP_Mat{i,d} = reshape(PC_Result{i}.normTraces(:,sp_dSP_time{i,d}' + nTau),nROI,[],length(nTau));
    dSPTA= squeeze(mean(dSP_Mat{i,d},2));
    dSPTA=dSPTA-median(dSPTA(:,[1:50]),2);
    dSPTA=dSPTA./F_Ref;
    cmap=repmat([0.3],size(dSPTA,1),3); cmap(d,:)=cmapAll(d,:);
    %cmap=turbo(size(dSPTA,1));
    if length(sp_dSP_time{i,d})>1 % more than 1 dSpike
    nexttile([4 1])
    plot(nTau,dSPTA'+[1:size(dSPTA,1)]*scale); hold all
    %dl=plot(nTau,rescale2(dSPTA,2)'+[1:size(dSPTA,1)]); hold all
    dl=plot(nTau,dSPTA'+[1:size(dSPTA,1)]*scale); hold all
    plot(nTau([1 end]),repmat(median(dSPTA,2),1,2)'+[1:size(dSPTA,1)]*scale,'color',[0.7 0.7 0.7]);
    
    arrayfun(@(l,c) set(l,'Color',c{:}),dl,num2cell(cmap,2))
    title(['N =' num2str(length(sp_dSP_time{i,d})) ' spikes']); ylim([0 (size(dSPTA,1)+1)*scale]+scale/2)
    set(gca,'YTick',[1:size(dSPTA,1)]*scale,'YTickLabel',[1:size(dSPTA,1)])
    end
    end

    nexttile([4 1])
    %sl=plot(nTau,rescale2(SPTA,2)'+[1:size(SPTA,1)]); hold all
    sl=plot(nTau,SPTA'+[1:size(SPTA,1)]*scale); hold all
    plot(nTau([1 end]),repmat(median(rescale2(SPTA,2),2),1,2)'+[1:size(SPTA,1)]*scale,'color',[0.7 0.7 0.7]);
    cmap_som=repmat([0.3],size(SPTA,1),3); cmap_som([1 DOI{i}],:)=cmapAll([1 DOI{i}],:);
    arrayfun(@(l,c) set(l,'Color',c{:}),sl,num2cell(cmapAll,2))
    title(['N =' num2str(length(sp_timeROI)) ' spikes']); ylim([0 (size(dSPTA,1)+1)*scale]+scale/2)
    set(gca,'YTick',[1:size(dSPTA,1)]*scale,'YTickLabel',[1:size(dSPTA,1)])
end

%% Calculate correlation image

for i=1%:3
    load([fpath{i} '/output_data.mat'])
    sz=double(Device_Data{1, 3}.ROI([2 4]));
    frm_end=max(Device_Data{1, 2}.Counter_Inputs(1, 1).data);
    f_seg=[[1:time_segment:frm_end] frm_end+1]; f_seg(2:end)=f_seg(2:end)-1;

    fileInd=ceil(sp_dSP_time/time_segment);
    frameInd=mod(sp_dSP_time,time_segment);
    %line up movie segments at simple spike
    dSTA_movie_align=[]; used_s=[];
    for j=unique(fileInd)
        mov_mc=double(readBinMov([fpath{i} '/mc_ShutterReg' num2str(j,'%02d') '.bin'],sz(2),sz(1)));
        load([fpath{i} '/mcTrace' num2str(j,'%02d') '.mat']);

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
        mov_res = SeeResiduals(mov_res,mc(:,1).*mc(:,3));
        mov_res= SeeResiduals(mov_res,bkg,1);

        for k=find(fileInd==j)
            try
                dSTA_movie_align = cat(3,dSTA_movie_align, tovec(mov_res(:,:,frameInd(k)+nTau)));
                used_s=[used_s; (fileInd(j)-1)*time_segment+frameInd(k)];
            end
        end
    end

    dSP_CorrImg=[];
    for s=1:size(dSTA_movie_align,3)
        dST_mov=dSTA_movie_align(:,:,s);
        dST_mov=Gaussfilt2_mov(toimg(dST_mov-movmedian(dST_mov,10,2),size(PC_Result{i}.ref_im,1),size(PC_Result{i}.ref_im,2)),3);
        dST_mov=tovec(dST_mov);
        dST_tr=squeeze(dSP_Mat(:,s,:));
        dST_tr=dST_tr-movmedian(dST_tr,10,2);

        % for roi = 1:size(dST_tr, 1)
        %     for p = 1:size(dST_mov, 1)
        %         dSP_CorrImg(roi, p,s) = corr(dST_tr(roi,:)', dST_mov(p,:)');
        %     end
        % end
        denominator=sqrt(sum((dST_tr-mean(dST_tr,2)).^2,2)*sum((dST_mov-mean(dST_mov,2)).^2,2)');
        dSP_CorrImg(:,:,s)=(dST_tr-mean(dST_tr,2))*(dST_mov-mean(dST_mov,2))'./denominator;
    end
    PC_Result{i}.dSP_CorrImg=dSP_CorrImg;
    PC_Result{i}.dSTA_movie_align=dSTA_movie_align;
end
save(fullfile(save_figto,'Result_PC_Prism_20240204.mat'),'PC_Result','fpath','-v7.3')

%% STA in Place field

i=10; PF_bin=[60 80]; Lap_segment={[3:14],[16:27]};
nFront=30; nBack=50;
nTau=[-nFront:nBack];
nROI= size(PC_Result{i}.normTraces,1);
cmapAll=turbo(nROI);
figure(10); clf;
nexttile([1 1]); show_footprnt_contour(PC_Result{i}.ftprnt(:,:,PC_Result{i}.dist_order),PC_Result{i}.ref_im,turbo(nROI));
nexttile([1 1]); imagesc(PC_Result{i}.Lap_FR);
nexttile([1 1]); imagesc(PC_Result{i}.Lap_sub);
for n=PC_Result{i}.dist_order
    nexttile([1 1]); imagesc(PC_Result{i}.Lap_F(:,:,n));
end
colormap('turbo')

PFpos=PF_bin*115/150;
pos_time=find_PosTime(PC_Result{i}.VR,PFpos);
tr=PC_Result{i}.normTraces(ref_ROI{i},:);
sp_ref=max(find_spike_bh(tr-movmedian(tr,100,2),5,3),[],1);
[~, shift]=max(reshape(PC_Result{i}.normTraces(1,find(sp_ref)+[-1:0]'),2,[]),[],1);
shift=shift-2;
sp_timeROI = find(sp_ref>0)+shift;
STA_PFmat=[];
STAtrace=squeeze(mean(PC_Result{i}.STA-prctile(PC_Result{i}.STA,40,3),2));
F_ref=mean(STAtrace(:,11+[10:15]),2);
for l=1:size(pos_time,1)
    
sp_PF=sp_timeROI(find(ismember(sp_timeROI,[pos_time(l,1):pos_time(l,2)])));
STA_PFmat{l}=reshape(PC_Result{i}.normTraces(:,sp_PF' + nTau),nROI,[],length(nTau));
STA_PFmat{l}=(STA_PFmat{l}-prctile(STA_PFmat{l},40,3))./F_ref;

end
SpikePFmat=cell2mat(STA_PFmat);
%SpikePFmat=SpikePFmat-prctile(SpikePFmat,5,3);
%SpikePFmat=SpikePFmat-median(SpikePFmat,3);
%SpikePFmat=SpikePFmat./mean(mean(SpikePFmat(:,:,nFront+[15:20]),3),2);
STA_PF=cellfun(@(x) squeeze(mean(x,2)),STA_PFmat, 'UniformOutput', false);
STA_PF=cell2mat(reshape(STA_PF,1,1,[]));

figure(11); clf;
nexttile([1 1])
pfl=plot(squeeze(mean(SpikePFmat(PC_Result{i}.dist_order,:,:),2))');
arrayfun(@(l,c) set(l,'Color',c{:}),pfl,num2cell(cmapAll,2))
axis tight
for n=PC_Result{i}.dist_order   
    nexttile([1 1])
    show_lap=find(squeeze(sum(isnan(STA_PF),[1 2]))==0);
    imagesc(squeeze(STA_PF(n,nFront+[-25:25],show_lap))')
    set(gca,'ytick',[1:length(show_lap)],'YTickLabel',show_lap)
end
colormap('turbo')

figure(12); clf;
STA_PF_LapSeg=[];
    for lseg=1:length(Lap_segment)
    STA_PF_LapSeg(:,:,lseg)=squeeze(mean(cell2mat(STA_PFmat(Lap_segment{lseg})),2));
    end
g=1;
for n=PC_Result{i}.dist_order   
    nexttile([1 1])
    lsg=plot(rescale2(squeeze(STA_PF_LapSeg(n,:,:)),1)); hold all
    plot([nFront+[-29:51]],rescale(STAtrace(n,:)),'k')
    arrayfun(@(l,c) set(l,'Color',c{:}),lsg,num2cell(distinguishable_colors(length(Lap_segment)),2))
    title(g); g=g+1;
end


%% dSTA movie

Sz_st=size(PC_Result{i}.DendriteBin);
figure(111); clf;
imshow2(imwarp(PC_Result{i}.ref_im,PC_Result{i}.tform,'OutputView', imref2d(Sz_st([1 2]))),[])
[~, roi]=imcrop; close(figure(111));
dSTAtrace=squeeze(mean(dSP_Mat,2));
dSTAmovie=toimg(mean(dSTA_movie_align,3),size(PC_Result{i}.ref_im,1),size(PC_Result{i}.ref_im,2));

dSTAmovie_Reg=imwarp(dSTAmovie,PC_Result{i}.tform,'OutputView', imref2d(Sz_st([1 2]))).*mat2gray(max(PC_Result{i}.DendriteBin,[],3)>0.01);
dSTAmovie_Reg=imgaussfilt3(dSTAmovie_Reg,[3 3 0.1]);
writeMov_wTrace([save_figto '/BHLm78_dSTAmovie'],-dSTAmovie(5:end-5,5:end-5,:),[],10,1,[],[],dSTAtrace')
writeMov_wTrace([save_figto '/BHLm78_dSTAmovieReg'],-imcrop_3d(dSTAmovie_Reg,roi),[],10,1,[],[],dSTAtrace')
%% Show correlation image of dSpike
Sz_st=size(PC_Result{i}.DendriteBin);
distalCorrImg=-toimg(squeeze(PC_Result{i}.dSP_CorrImg(6,:,:)),size(PC_Result{i}.ref_im,1),size(PC_Result{i}.ref_im,2));
distalCorrImg_Reg=imwarp(distalCorrImg,PC_Result{i}.tform,'OutputView', imref2d(Sz_st([1 2]))).*(max(PC_Result{i}.DendriteBin,[],3)>0.01);
distalCorrImg_Reg(distalCorrImg_Reg==0)=NaN;
moviefixsc(distalCorrImg_Reg)
colormap(turbo)

figure(2); clf;
for s=1:size(distalCorrImg_Reg,3)
    nexttile([1 1])
    plot(squeeze(dSP_Mat(6,s,:)))
    nexttile([1 1])
    imagesc(imcrop(distalCorrImg_Reg(:,:,s),roi),[-0.7 0.7]); axis equal tight off
    colormap(turbo)
end

figure(3); clf; cmap=turbo(7); M=[];
for n=[2]
distalCorrImg=-imgaussfilt3(toimg(squeeze(dSP_CorrImg(n,:,:)),size(PC_Result{i}.ref_im,1),size(PC_Result{i}.ref_im,2)),[3 3 0.1]);
distalCorrImg_Reg=imwarp(distalCorrImg,PC_Result{i}.tform,'OutputView', imref2d(Sz_st([1 2]))).*(max(PC_Result{i}.DendriteBin,[],3)>0.01);
distalCorrImg_Reg(distalCorrImg_Reg==0)=NaN;
Mean_corr=squeeze(mean(tovec(distalCorrImg_Reg),1,'omitnan'));
STD_corr=squeeze(std(tovec(distalCorrImg_Reg),0,1,'omitnan'));
%errorbar_shade(sp_dSP_time/1000,Mean_corr,STD_corr,cmap(n,:)); hold all
errorbar_shade([1:length(sp_dSP_time)],Mean_corr,STD_corr,cmap(n,:)); hold all
%plot(Median_corr,'color',cmap(n,:),'marker','.'); hold all
M(n,:)=[mean(Mean_corr(1:16)) mean(Mean_corr(end-15:end))];
[~, p(n)]=ttest2(tovec(tovec(distalCorrImg_Reg(:,:,10:20))),tovec(tovec(distalCorrImg_Reg(:,:,end-10:end))));
end
axis tight
xlabel('Spike order')
ylabel('Mean correlation coefficient')

% PCA/ICA on the correlation img
figure(4); 
Sz_st=size(PC_Result{i}.DendriteBin);
distalCorrImg=-toimg(squeeze(dSP_CorrImg(6,:,:)),size(PC_Result{i}.ref_im,1),size(PC_Result{i}.ref_im,2));
distalCorrImg_Reg=imwarp(distalCorrImg,PC_Result{i}.tform,'OutputView', imref2d(Sz_st([1 2]))).*(max(PC_Result{i}.DendriteBin,[],3)>0.01);
%distalCorrImg_Reg(distalCorrImg_Reg==0)=NaN;

datDS = imresize(imcrop_3d(distalCorrImg_Reg,roi), 0.4, 'bilinear', 'Antialiasing',true);

MovVec = tovec(datDS);
covMat = MovVec*MovVec';
[V, D] = eig(covMat);
D = diag(D);
D = D(end:-1:1);
V = V(:,end:-1:1);
vSign = sign(max(V) - max(-V));  % make the largest value always positive
V = V.*vSign;
eigTraces = V'*MovVec;

figure(8); clf
nKeep = 20;
for i=1:nKeep
    plot(rescale(eigTraces(i,:)-movmedian(eigTraces(i,:),30))+i-0.5)
    hold all
end
set(gca, 'YDir','reverse')
sz=size(distalCorrImg_Reg);
eigImgs = zeros(sz(1), sz(2), nKeep);
for j = 1:nKeep
    eigImgs(:,:,j) = mean(distalCorrImg_Reg.*reshape(eigTraces(j,:), [1, 1,sz(3)]),3);
end
figure; clf;
for j = 1:nKeep
    nexttile([1 1]);
    imshow2(eigImgs(3:end-3,3:end-3,j), []);
    title(num2str(j))
end

keep_ind=[1:5];

[ics, mixmat, sepmat] = sorted_ica(eigTraces(keep_ind,:)',length(keep_ind));
figure(9); clf
stackplot(ics)

eigImgsVec = tovec(eigImgs(:,:,keep_ind));
footPrintsVec = eigImgsVec*sepmat';
footPrints = toimg(footPrintsVec, [sz(1), sz(2)]);
figure(10); clf
for j = 1:size(ics,2)
    nexttile([1 1]);
    imshow2(footPrints(3:end-3,3:end-3,j), []);
    title(num2str(j))
end


figure(5);
    tr=PC_Result{i}.normTraces(ref_ROI{i},:);
    sp_ref=sum(find_spike_bh(tr-movmedian(tr,300,2),4,2),1);
    nROI=size(PC_Result{i}.normTraces,1);
    sp=PC_Result{i}.spike;
    SS=sp_ref;
    CS=sp_ref.*PC_Result{i}.CS_trace;

    sp_timeROI = find(sp_ref>0);
    sp_dSp=sp;
    sp_dSp(:,sp_timeROI'+[0:3])=0;

    sp_dSP_time=find(sum(sp_dSp(6,:),1)>0);
    peri_spiketime=[]; peri_CStime=[];
    for s=1:length(sp_dSP_time)
lagTime=find(SS)-sp_dSP_time(s);
post_sp=find(lagTime>0,1);
pre_sp=post_sp-1;
peri_spiketime(s,:)=[lagTime(pre_sp) lagTime(post_sp)];

lagTimeCS=find(CS)-sp_dSP_time(s);
post_sp=find(lagTimeCS>0,1);
if ~isempty(post_sp)
pre_sp=post_sp-1;
peri_CStime(s,:)=[lagTimeCS(pre_sp) lagTimeCS(post_sp)];
else
peri_CStime(s,:)=[NaN NaN];
end
    end
clf;
nexttile([1 1])
    [periSp_hist, bin]=histcounts(abs(peri_spiketime(:,1)),100); hold all
    plot(bin(1:end-1)+(bin(2)-bin(1))/2,periSp_hist,'color',[0.6 0.6 0.6])
    [periCS_hist, ~]=histcounts(abs(peri_CStime(:,1)),bin); hold all
    plot(bin(1:end-1)+(bin(2)-bin(1))/2,periCS_hist,'color',[1 0.2 0.2])
    set(gca,'xscale','log','Xdir','reverse')
    title('Nearest somatic spike before dSpike')
    xlabel('Time (ms)')
    xlim([1 10000])
nexttile([1 1])
    [periSp_hist, bin]=histcounts(abs(peri_spiketime(:,2)),100); hold all
    plot(bin(1:end-1)+(bin(2)-bin(1))/2,periSp_hist,'color',[0.6 0.6 0.6])
    [periCS_hist, ~]=histcounts(abs(peri_CStime(:,2)),bin); hold all
    plot(bin(1:end-1)+(bin(2)-bin(1))/2,periCS_hist,'color',[1 0.2 0.2])
    set(gca,'xscale','log')
    title('Nearest somatic spike after dSpike')
    xlabel('Time (ms)')
    legend('Total spike','Complex spike')
    xlim([1 10000])

%%


%%



t_ref=[6600:6650];
avgImg=mean(mov_mc(bound:end-bound,bound:end-bound,:),3);
mask=(avgImg-medfilt2(avgImg,[15 15]))>80;
mov_seg=mov_res(bound:end-bound,bound:end-bound,t_ref);
mov_seg=mov_seg-mean(mov_seg,3);
movSegfilt = spatialfilt(mov_seg, 7, 3);
movSegfilt=movSegfilt-movmedian(movSegfilt,10,3);
[movSegfiltPCA, eigVecs, eigVals] = pcafilt(movSegfilt, 20);
eigImgs = toimg(tovec(mov_seg)*eigVecs(:,1:20), size(mov_seg,1), size(mov_seg,2));
figure(13); clf
for j = 1:20;
    nexttile([1 1])
    imshow2(eigImgs(:,:,j), []);
end;
clf;
for j=1:14
trace=rescale(PC_Result{i}.normTrace(PC_Result{i}.dist_order(j),t_ref));    
trace=trace-movmedian(trace,10);
corrImg=-toimg(corr(tovec(movSegfilt)',trace'),size(movSegfilt,1),size(movSegfilt,2));
corrImg=imgaussfilt(corrImg,3).*mask;
nexttile([1 1]);
%imagesc(im_merge(cat(3,avgImg,mat2gray(corrImg)/10),[0.5 0.5 0.5; 0.5 0 0]))
imshow2(corrImg,[0 0.6])
title(num2str(j))
end
colormap(turbo)


%%
f=1; tau=[0:3];
sp_list=find(PC_Result{f}.ref_spike);
shortISI_list=find((sp_list(2:end)-sp_list(1:end-1))<3);
rmv_list=[];
for s=shortISI_list
 [~, arg]=max(PC_Result{f}.normTraces(2,sp_list(s:s+1)));
 if arg==1
     rmv_list=[rmv_list s+1];
 else
     rmv_list=[rmv_list s];
 end
end
sp_list(rmv_list)=[];


% histogram(sp_list(2:end)-sp_list(1:end-1),[0:1:500])
% set(gca,'xscale','log')

Spike_amp_mat=squeeze(max(reshape(PC_Result{f}.normTraces(:,sp_list+tau'),[],length(tau),length(sp_list)),[],2));
figure(4); clf;
Bin_trace=zeros(size(PC_Result{f}.normTraces,1),size(PC_Result{f}.normTraces,2));
%for n=1:size(Spike_amp_mat,1)
%     nexttile([1 1])
%     [trace_hist{n}.distr trace_hist{n}.bin]=histcounts(PC_Result{f}.normTraces(n,:),1000);
%     [trace_hist{n}.spikeH] = histcounts(Spike_amp_mat(n,:),trace_hist{n}.bin);
%     plot(trace_hist{n}.bin(2:end),rescale(trace_hist{n}.distr))
% hold all
% plot(trace_hist{n}.bin(2:end),rescale(trace_hist{n}.spikeH),'color',[0.5 0.2 0.8])
% threshold(n)=prctile(Spike_amp_mat(n,:),10);
% plot([threshold(n) threshold(n)],[0 1],'r','LineWidth',2)
% Prob_thres(n,1)=sum(PC_Result{f}.normTraces(n,:)>threshold(n))/size(PC_Result{f}.normTraces,2);
% title('Probability upper than threshold: ',num2str(Prob_thres(n)))
% Bin_trace(n,:)=PC_Result{f}.normTraces(n,:)>threshold(n);
%end
figure(2); clf;
ax1=nexttile([1 1]);
imagesc(PC_Result{f}.normTraces./median(Spike_amp_mat,2))
colormap(turbo)
ax2=nexttile([1 1]);
Tracesum=movsum(sum(PC_Result{f}.normTraces./median(Spike_amp_mat,2),1),3,2);
plot(Tracesum); hold all
linkaxes([ax1 ax2],'x');
bAP_frame=sp_list+tau'; bAP_frame=bAP_frame(:);
thr=prctile(Tracesum(bAP_frame),15);
dSp_frame=find(Tracesum>thr); dSp_frame=setdiff(dSp_frame,bAP_frame);
plot(dSp_frame,Tracesum(dSp_frame),'ko')
plot(sp_list,Tracesum(sp_list),'r.','markersize',10)
f1 = gcf;
f1.WindowScrollWheelFcn = @scrollWheel;

figure(2); clf;

[Total_dist trsum_bin]=histcounts(Tracesum,1000);
[bAPsum_dist]=histcounts(Tracesum(bAP_frame),trsum_bin);
[omitbAPsum_dist]=histcounts(Tracesum(setdiff([1:length(Tracesum)],bAP_frame)),trsum_bin);

plot(trsum_bin(2:end),rescale(Total_dist)); hold all
plot(trsum_bin(2:end),rescale(bAPsum_dist))
plot(trsum_bin(2:end),rescale(omitbAPsum_dist))
line([thr thr],[0 1],'color','r','linewidth',2)


%% See Raw movie of candidate frame
f=1; cand_Frame=285156; 
nFront=20; nBack=50;
nTau=[-nFront:nBack];

fileInd=ceil(cand_Frame/time_segment);
frameInd=mod(cand_Frame,time_segment);

load([fpath{f} '/output_data.mat'])
sz=double(Device_Data{1, 3}.ROI([2 4]));
frm_end=max(Device_Data{1, 2}.Counter_Inputs(1, 1).data);
f_seg=[[1:time_segment:frm_end] frm_end+1]; f_seg(2:end)=f_seg(2:end)-1;
for j=unique(fileInd)
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
    mov_res = SeeResiduals(mov_res,mc(:,1).*mc(:,3));
    mov_res= SeeResiduals(mov_res,bkg,1);

    for k=find(fileInd==j)
    try 
        mov_tmp = mov_res(:,:,frameInd(k)+nTau);
    end
    end
end

mov_tmpFilt=-imgaussfilt3(mov_tmp,[3 3 0.1]).*mat2gray(medfilt2(mean(mov_mc,3),[1 1]));
figure; moviefixsc(mov_tmpFilt,[-5 30])

figure;
kymo_tmp=-polyLineKymo2(mov_tmp,25,25);
clf;
ax1=nexttile([1 1]);
imagesc(rescale2(kymo_tmp,1)')
ax2=nexttile([1 1]);
ls=plot(rescale2(kymo_tmp,1)+[1:size(kymo_tmp,2)]+0.5); axis tight;
arrayfun(@(l,c) set(l,'Color',c{:}),ls,num2cell(jet(size(kymo_tmp,2)),2))
linkaxes([ax1 ax2],'x')
%% Complex spike average vs Simple spike average
nBack = -30;
nFront = 200;
tau = nBack:nFront;
nTau = length(tau);
nCS=length(CS_time);

CS_tau = CS_time' + tau;
CS_spikeMat = reshape(Result.normTraces(:,CS_tau), 6, nCS , nTau);

SS_trace_spike=spike.*(~CS_trace);
SS_time = find(SS_trace_spike);
SS_time(end) = [];
SS_time([find(diff(SS_time) < 45), find(diff(SS_time) < 45)+1]) = [];

nSS = length(SS_time);

SS_tau = SS_time' + tau;
SS_spikeMat = reshape(Result.normTraces(:,SS_tau), 6, nSS, nTau);
    
        figure(6); clf; 
        ax1=[]; ax2=[];
        tiledlayout(10,6)
        for i=1:nCS
            ax1=[ax1 nexttile([1 1])];
            imagesc(rescale2(squeeze(CS_spikeMat(:,i,:)),2))
            colormap('turbo')
            title(i)
        end
        linkaxes(ax1,'xy')
        
        figure(5); clf; 
        tiledlayout(10,6)
        for i=1:60
            ax2=[ax2 nexttile([1 1])];
            imagesc(rescale2(squeeze(SS_spikeMat(:,i,:)),2))
            colormap('turbo')
            title(i)
        end
        linkaxes(ax2,'xy')



%% PCA/ICA
fi=1;
load(fullfile(fpath{fi},'PC_Result.mat'))

f=5;   bound=5;
load([fpath{fi} '/output_data.mat'])
sz=double(Device_Data{1, 3}.ROI([2 4]));
mov_mc=double(readBinMov([fpath{fi} '/mc_ShutterReg' num2str(f,'%02d') '.bin'],sz(2),sz(1)));
load([fpath{fi} '/mcTrace' num2str(f,'%02d') '.mat']);

mov_mc=mov_mc(:,:,1:time_segment); 
mc=mcTrace.xymean(1:time_segment,:);

mov_res= mov_mc-mean(mov_mc,3);
bkg = zeros(2, size(mov_mc,3));
bkg(1,:) = linspace(-1, 1, size(mov_mc,3));  % linear term
bkg(2,:) = linspace(-1, 1, size(mov_mc,3)).^2;  % quadratic term
mov_res = SeeResiduals(mov_res,mc);
mov_res = SeeResiduals(mov_res,mc.^2);
mov_res = SeeResiduals(mov_res,mc(:,1).*mc(:,3));
mov_res= SeeResiduals(mov_res,bkg,1);

nFrames2=size(mov_res,3);
fGauss=fspecial('gaussian',10,3);
mov_res = imfilter(mov_res,fGauss,'same');
datDS = imresize(mov_res(bound:end-bound,bound:end-bound,:), 0.4, 'bilinear', 'Antialiasing',true);

MovVec = tovec(datDS);
covMat = MovVec*MovVec';
[V, D] = eig(covMat);
D = diag(D);
D = D(end:-1:1);
V = V(:,end:-1:1);
vSign = sign(max(V) - max(-V));  % make the largest value always positive
V = V.*vSign;
eigTraces = V'*MovVec;

figure(8); clf
nKeep = 20;
for i=1:nKeep
    plot(rescale(eigTraces(i,:)-movmedian(eigTraces(i,:),30))+i-0.5)
    hold all
end
set(gca, 'YDir','reverse')

eigImgs = zeros(sz(2), sz(1), nKeep);
for j = 1:nKeep
    eigImgs(:,:,j) = mean(mov_res.*reshape(eigTraces(j,:), [1, 1, nFrames2]),3);
end
figure; clf;
for j = 1:nKeep
    nexttile([1 1]);
    imshow2(eigImgs(3:end-3,3:end-3,j), []);
    title(num2str(j))
end

keep_ind=[1:7];

[ics, mixmat, sepmat] = sorted_ica(eigTraces(keep_ind,:)',length(keep_ind));
figure(9); clf
stackplot(ics)

eigImgsVec = tovec(eigImgs(:,:,keep_ind));
footPrintsVec = eigImgsVec*sepmat';
footPrints = toimg(footPrintsVec, [sz(2), sz(1)]);
figure(10); clf
for j = 1:size(ics,2)
    nexttile([1 1]);
    imshow2(footPrints(3:end-3,3:end-3,j), []);
    title(num2str(j))
end

% MovVecInv=pinv(MovVec);
% Vics=(ics'*MovVecInv)';
keep_ind_ics=[1 3 4 5 8 9];
% coeffsICA=MovVec'*Vics(:,keep_ind_ics);
% ReconMovICA=toimg((coeffsICA*Vics(:,keep_ind_ics)')',size(datDS,1),size(datDS,2));

coeffs=MovVec'*(V(:,keep_ind));
ReconMov=toimg(coeffs*(V(:,keep_ind)*D(keep_ind))',size(datDS,1),size(datDS,2));


ReconMovICA=toimg(footPrintsVec(:,keep_ind_ics)*ics(:,keep_ind_ics)',sz(2),sz(1));
%ReconMov=toimg(tovec(mov_res).*mean(footPrintsVec(:,keep_ind_ics),2),sz(2),sz(1));
%ReconMov_filt=-imgaussfilt3(ReconMov,[2 2 0.5]);
figure; writeMov_wTrace([fpath{1},'/dFMov_mov34PCA'],-ReconMov(bound:end-bound,bound:end-bound,:),[4000:6250],10,1,[-5 35],[],mean(Result.traces([4 5 6 7],f_seg(f):f_seg(f+1)-1),1))
figure; writeMov_wTrace([fpath{1},'/dFMov_mov34ICA'],-(ReconMovICA(bound:end-bound,bound:end-bound,:)),[4000:6250],10,1,[],[],mean(Result.traces([4 5 6 7],f_seg(f):f_seg(f+1)-1),1))

%% Polyline
f=34;   bound=5;
load([fpath '/output_data.mat'])
sz=double(Device_Data{1, 3}.ROI([2 4]));
mov_mc=double(readBinMov([fpath '/mc_ShutterReg' num2str(f,'%02d') '.bin'],sz(2),sz(1)));
load([fpath '/mcTrace' num2str(f,'%02d') '.mat']);

mov_mc=mov_mc(:,:,1:time_segment); 
mc=mcTrace.xymean(1:time_segment,:);

mov_res= mov_mc-mean(mov_mc,3);
bkg = zeros(2, size(mov_mc,3));
bkg(1,:) = linspace(-1, 1, size(mov_mc,3));  % linear term
bkg(2,:) = linspace(-1, 1, size(mov_mc,3)).^2;  % quadratic term
mov_res = SeeResiduals(mov_res,mc);
mov_res = SeeResiduals(mov_res,mc.^2);
mov_res = SeeResiduals(mov_res,mc(:,1).*mc(:,3));
mov_res= SeeResiduals(mov_res,bkg,1);

nFrames2=size(mov_res,3);
fGauss=fspecial('gaussian',10,3);
mov_res = imfilter(mov_res,fGauss,'same');

[kymo, kymoROI]=polyLineKymo2(-mov_res,40,20,Result.ref_im);
figure;
clf;
imagesc(kymo')
colormap(turbo)
caxis([-1 50])


%%
tr=Result.normTraces(2,:);
bwCStrace=bwlabel(CS_trace);
for CS=1:max(bwCStrace)
CS_time_tmp= find(bwCStrace==CS);
tr_sub_tmp=tr_sub(CS_time_tmp);
slope=tr_sub_tmp(2:end)-tr_sub_tmp(1:end-1);
[~, max_slope]=max(slope);
CS_time_slope(CS)=max_slope+CS_time_tmp(1);
end

nBack = -200;
nFront = 200;
tau = nBack:nFront;
nTau = length(tau);
nCS=length(CS_time);

CS_tau = CS_time_slope' + tau;
CS_spikeMat = reshape(Result.normTraces(:,CS_tau), 6, nCS , nTau);

figure(4); clf; 
        ax1=[]; ax2=[];
        tiledlayout(10,6)
        for i=1:nCS
            ax1=[ax1 nexttile([1 1])];
            imagesc(rescale2(squeeze(CS_spikeMat(:,i,:)),2))
            colormap('turbo')
            title(i)
        end
        linkaxes(ax1,'xy')
%%
clear sp_relamp_soma sp_relamp_dd sp_amp_soma sp_amp_dd
shift_window=[0:2]; CS_baseline_wnd=[-100:200];
tr_soma=Result.normTraces(2,:);
tr_dd=Result.normTraces(6,:);

for c=1:length(CS_time)
    t_tmp=CS_time(c)+CS_baseline_wnd;
    tr_soma_seg=tr_soma(t_tmp);

    tr_dd_seg=tr_dd(t_tmp);
    sp_soma=find_spike_bh(tr_soma_seg,5,3.5);
    sp_dd=zeros(1,length(t_tmp));
    tr_soma_pass=tr_soma_seg;
    tr_dd_pass=tr_dd_seg;
    se = strel('square', 5); % 0 degree means horizontal
    for s=find(sp_soma)
        [~, shift]=max(tr_dd_seg(s+shift_window));
        sp_dd(s+shift-1)=1;
    
    tmp=zeros(1,length(t_tmp)); tmp_d=zeros(1,length(t_tmp));
    tmp(s)=1; tmp_d(s+shift-1)=1;
    sp_soma_di = imdilate(tmp, se);
    sp_dd_di = imdilate(tmp_d, se);

    tr_soma_pass(find(sp_soma_di))=NaN;
    valid_point=find(~isnan(tr_soma_pass));
    tr_soma_pass=interp1(valid_point,tr_soma_pass(valid_point),[1:length(t_tmp)],'linear');

    tr_dd_pass(find(sp_dd_di))=NaN;
    valid_point=find(~isnan(tr_dd_pass));
    tr_dd_pass=interp1(valid_point,tr_dd_pass(valid_point),[1:length(t_tmp)],'linear');

    end
    tr_dd_pass=movmean(tr_dd_pass,10);
    tr_soma_pass=movmean(tr_soma_pass,10);

%     figure; 
%     cmap=[1:sum(sp_soma)];
%     plot(tr_soma_seg); hold all
%     plot([1 length(t_tmp)],[prctile(tr_soma_seg,30) prctile(tr_soma_seg,30)])
%     plot(tr_soma_pass);
%     scatter(find(sp_soma),tr_soma_seg(find(sp_soma)),50,cmap,'filled');
%     colormap(winter)
% 
%     plot(tr_dd_seg+15);
%     plot(tr_dd_pass+15);
%     plot([1 length(t_tmp)],[prctile(tr_dd_seg,30) prctile(tr_dd_seg,30)]+15)   
%     scatter(find(sp_dd),tr_dd_seg(find(sp_dd))+15,50,cmap,'filled');

sp_amp_soma{c}=tr_soma_seg(find(sp_soma))-prctile(tr_soma_seg,30);
sp_amp_dd{c}=tr_dd_seg(find(sp_dd))-prctile(tr_dd_seg,30);

sp_relamp_soma{c}=tr_soma_seg(find(sp_soma))-tr_soma_pass(find(sp_soma));
sp_relamp_dd{c}=tr_dd_seg(find(sp_dd))-tr_dd_pass(find(sp_dd));
end

%% N th spike, amplitude
clear Sp_order_somdd Sp_ratio_somdd
figure; cmap=distinguishable_colors(7); cmap_light=cmap+0.5; cmap_light(cmap_light>1)=1;
for n=1:7
for c=1:length(sp_amp_soma)
if length(sp_amp_soma{c})<n
Sp_order_somdd{n}(:,c)= NaN(2,1);    
Sp_ratio_somdd{n}(:,c) =NaN;
else
Sp_order_somdd{n}(:,c) = [sp_amp_soma{c}(n); sp_amp_dd{c}(n)];
Sp_ratio_somdd{n}(:,c) = [sp_amp_dd{c}(n)/sp_amp_soma{c}(n)];
end
end
    %plot(Sp_order_somdd{n}(1,:),Sp_order_somdd{n}(2,:),'.','color',cmap_light(n,:),'markersize',10)
hold all
end
plot([1:7],cell2mat(Sp_ratio_somdd')','color',[0.7 0.7 0.7])
hold all

M=cellfun(@(x) mean(x,2,'omitnan'),Sp_ratio_somdd,'UniformOutput',false);
S=cellfun(@(x) std(x,0,2,'omitnan'),Sp_ratio_somdd,'UniformOutput',false);
N=cellfun(@(x) sum(~isnan(x(1,:))),Sp_ratio_somdd,'UniformOutput',false);
errorbar([1:7],cell2mat(M),cell2mat(S),'color',[0 0 0])

%% N th spike, amplitude scatter plot
figure; cmap=jet(7); cmap_light=cmap+0.5; cmap_light(cmap_light>1)=1;
M=cellfun(@(x) mean(x,2,'omitnan'),Sp_order_somdd,'UniformOutput',false);
S=cellfun(@(x) std(x,0,2,'omitnan'),Sp_order_somdd,'UniformOutput',false);
N=cellfun(@(x) sum(~isnan(x(1,:))),Sp_order_somdd,'UniformOutput',false);

for n=1:7
plot(Sp_order_somdd{n}(1,:),Sp_order_somdd{n}(2,:),'.','color',cmap_light(n,:),'markersize',10); 
hold all
errorbar(M{n}(1),M{n}(2),S{n}(2),S{n}(2),S{n}(1),S{n}(1),'color',cmap(n,:),'marker','*','markersize',20)
end


%% N th spike, relative amplitude
clear Sp_relorder_somdd Sp_relratio_somdd
figure; cmap=distinguishable_colors(7); cmap_light=cmap+0.5; cmap_light(cmap_light>1)=1;
for n=1:7
for c=1:length(sp_amp_soma)
if length(sp_amp_soma{c})<n
Sp_relorder_somdd{n}(:,c)= NaN(2,1);    
Sp_relratio_somdd{n}(:,c) =NaN;
else
Sp_relorder_somdd{n}(:,c) = [sp_relamp_soma{c}(n); sp_relamp_dd{c}(n)];
Sp_relratio_somdd{n}(:,c) = [sp_relamp_dd{c}(n)/sp_relamp_soma{c}(n)];
end
end
    %plot(Sp_order_somdd{n}(1,:),Sp_order_somdd{n}(2,:),'.','color',cmap_light(n,:),'markersize',10)
hold all
end
plot([1:7],cell2mat(Sp_relratio_somdd')','color',[0.7 0.7 0.7])
hold all

M=cellfun(@(x) mean(x,2,'omitnan'),Sp_relratio_somdd,'UniformOutput',false);
S=cellfun(@(x) std(x,0,2,'omitnan'),Sp_relratio_somdd,'UniformOutput',false);
N=cellfun(@(x) sum(~isnan(x(1,:))),Sp_relratio_somdd,'UniformOutput',false);
errorbar([1:7],cell2mat(M),cell2mat(S),'color',[0 0 0])
xlabel('N^t^h spike of CS')
ylabel('Spike local amplitude ratio (Distal dend/Soma)')

%%

figure; cmap=jet(7); cmap_light=cmap+0.5; cmap_light(cmap_light>1)=1;
M=cellfun(@(x) mean(x,2,'omitnan'),Sp_relorder_somdd,'UniformOutput',false);
S=cellfun(@(x) std(x,0,2,'omitnan'),Sp_relorder_somdd,'UniformOutput',false);
N=cellfun(@(x) sum(~isnan(x(1,:))),Sp_relorder_somdd,'UniformOutput',false);

for n=1:7
plot(Sp_relorder_somdd{n}(1,:),Sp_relorder_somdd{n}(2,:),'.','color',cmap_light(n,:),'markersize',10); 
hold all
errorbar(M{n}(1),M{n}(2),S{n}(2),S{n}(2),S{n}(1),S{n}(1),'color',cmap(n,:),'marker','*','markersize',20)
end

c = colorbar;
c.Ticks = [0:6]/6;
c.TickLabels=num2str([1:7]');

%% line up movie segments at complex spike
wnd=[-15:250];
CS_movie_align=[];
for j=unique(fileInd)
    mov_mc=double(readBinMov([fpath '/mc_ShutterReg' num2str(j,'%02d') '.bin'],sz(2),sz(1)));
    load([fpath '/mcTrace' num2str(j,'%02d') '.mat']);

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
    mov_res = SeeResiduals(mov_res,mc(:,1).*mc(:,3));
    mov_res= SeeResiduals(mov_res,bkg,1);

    for k=find(fileInd==j)
    try 
        CS_movie_align = cat(3,CS_movie_align, tovec(mov_res(:,:,frameInd(k)+wnd)));
    end
    end
end


%%
CS_movie_align=vm(CS_movie_align);
CS_movie_align.transpose.savebin([fpath '/CS_movie_align.bin'])
%%
bound=5;
CS_movie_align=double(CS_movie_align);
CS_movie_align=reshape(CS_movie_align,sz(2),sz(1),length(wnd),[]);        
subCS_movie_align=CS_movie_align()
CSVec = tovec(CS_movie_align);
covMat = CSVec*CSVec';
[V, D] = eig(covMat);
D = diag(D);
D = D(end:-1:1);
V = V(:,end:-1:1);
vSign = sign(max(V) - max(-V));  % make the largest value always positive
V = V.*vSign;
eigTraces = V'*MovVec;

figure(8); clf
nKeep = 20;
for i=1:nKeep
    plot(rescale(eigTraces(i,:)-movmedian(eigTraces(i,:),30))+i-0.5)
    hold all
end
set(gca, 'YDir','reverse')

eigImgs = zeros(sz(2), sz(1), nKeep);
for j = 1:nKeep
    eigImgs(:,:,j) = mean(mov_res.*reshape(eigTraces(j,:), [1, 1, nFrames2]),3);
end
figure; clf;
for j = 1:nKeep
    nexttile([1 1]);
    imshow2(eigImgs(3:end-3,3:end-3,j), []);
    title(num2str(j))
end

%%

[Lap_FR Lap_V]=PlaceTrigger_average(sum(Result.spike(noi,:)),50,Result.VR,-0.1,115);
[LickFR]=PlaceTrigger_average(Result.VR(9,:),50,Result.VR,-0.1,115);

%%
load(fullfile(fpath,'Result_20231004.mat'))

f=34;   bound=6;
load([fpath '/output_data.mat'])
sz=double(Device_Data{1, 3}.ROI([2 4]));
mov_mc=double(readBinMov([fpath '/mc_ShutterReg' num2str(f,'%02d') '.bin'],sz(2),sz(1)));
load([fpath '/mcTrace' num2str(f,'%02d') '.mat']);

mov_mc=mov_mc(:,:,1:time_segment); 
mc=mcTrace.xymean(1:time_segment,:);

mov_res= mov_mc-mean(mov_mc,3);
bkg = zeros(2, size(mov_mc,3));
bkg(1,:) = linspace(-1, 1, size(mov_mc,3));  % linear term
bkg(2,:) = linspace(-1, 1, size(mov_mc,3)).^2;  % quadratic term
mov_res = SeeResiduals(mov_res,mc);
mov_res = SeeResiduals(mov_res,mc.^2);
mov_res = SeeResiduals(mov_res,mc(:,1).*mc(:,3));
mov_res= SeeResiduals(mov_res,bkg,1);

nFrames2=size(mov_res,3);
mov_res_filter=zeros(size(mov_res));
 fGauss=fspecial('gaussian',20,7);
 mov_res_filter = imfilter(mov_res(bound:end-bound,bound:end-bound,:),fGauss,'same');
 F=imgaussfilt(Result.ref_im,5);
avgImg=Result.ref_im-imgaussfilt(Result.ref_im,20);
%H=ones(8,8);
%     for z=1:nFrames2
%     %mov_res_filter (:,:,z) = medfilt2(mov_res(:,:,z),[15 15]);                                  
%     mov_res_filter (:,:,z) = filter2(H,mov_res(:,:,z));                                  
%     end

tr=Result.traces([1 2 3 5 4 6],1:end-1);
tr= (tr - min(tr,[],2)) ./ (max(tr,[],2) - min(tr,[],2));

mov_res_filter= (mov_res_filter)./F(bound:end-bound,bound:end-bound,:).*mat2gray(avgImg(bound:end-bound,bound:end-bound,:));
figure; writeMov_wKymo([fpath,'/mov_Res'],-mov_res_filter,[4100:4350],10,1,[-0.01 0.05],[],tr(:,f_seg(f):f_seg(f+1)-1))
figure; writeMov_wKymo([fpath,'/mov_Res2'],-mov_res_filter,[4840:4940],10,1,[-0.01 0.05],[],tr(:,f_seg(f):f_seg(f+1)-1))

datDS = imresize(mov_res(bound:end-bound,bound:end-bound,:), 0.4, 'bilinear', 'Antialiasing',true);

%%
tr=Result.normTraces(2,:); t=[1:length(tr)];
spike=find_spike_bh(tr-movmedian(tr,300,2),4,3);
spike_time=find(spike); spike_height=tr(spike_time);

pass_frq=[15];
Fs=1000;

freq_lowhigh=pass_frq/(Fs/2);
[b, a] = butter(4, freq_lowhigh, 'low');
tr_pass = filtfilt(b, a, tr);
figure(2); clf; show_t=[750000:850000];
plot(t(show_t),tr(show_t)); hold all
plot(t(show_t),tr_pass(show_t));
plot(t(show_t),zeros(length(show_t),1));

[trans tr_trace]=detect_transient(tr_pass,[3 0.5],spike);
CS_ind=find(trans.spike_number>2 & trans.meanISI<40);
CS_trace=ismember(tr_trace,CS_ind);
figure(3); clf;
plot(t,tr)
hold all
tr_nan=tr; tr_nan(CS_trace==0)=NaN;
plot(t,tr_nan,'r')

CS_trace_spike=spike.*bwlabel(CS_trace);
for CS=1:max(CS_trace_spike)
    CS_time(CS)= find(CS_trace_spike==CS,1,'first');
    tmp = find(CS_trace_spike==CS,3,'first');
    CS2_time(CS) = tmp(2);
    CS3_time(CS) = tmp(3);
end

nBack = -20;
nFront = 150;
tau = nBack:nFront;
nTau = length(tau);

CS_tau = CS_time' + tau;
CS_spikeMat = reshape(Result.traces(:,CS_tau), 6, 122, nTau);
CS_SubthMat = reshape(tr_pass(CS_tau), 1, 122, nTau);
for j = 1:6;
    subplot(2,3,j)
    imshow2(squeeze(CS_spikeMat(j,:,:)), [])
end;

Ca_height=max(squeeze(CS_SubthMat),[],2);
[~, Ca_heightIdx] = sort(Ca_height); 

fileInd=ceil(CS_time(Ca_heightIdx)/time_segment);
frameInd=mod(CS_time(Ca_heightIdx),time_segment);

trsc=Result.traces([1 2 3 5 4 6],1:end-1);
trsc= (trsc - min(trsc,[],2)) ./ (max(trsc,[],2) - min(trsc,[],2));

    figure;
%    tiledlayout(10,1)
for i=[108 113]

f=fileInd(i);   bound=6;
load([fpath '/output_data.mat'])
sz=double(Device_Data{1, 3}.ROI([2 4]));
mov_mc=double(readBinMov([fpath '/mc_ShutterReg' num2str(f,'%02d') '.bin'],sz(2),sz(1)));
load([fpath '/mcTrace' num2str(f,'%02d') '.mat']);

mov_mc=mov_mc(:,:,1:time_segment); 
mc=mcTrace.xymean(1:time_segment,:);

mov_res= mov_mc-mean(mov_mc,3);
bkg = zeros(2, size(mov_mc,3));
bkg(1,:) = linspace(-1, 1, size(mov_mc,3));  % linear term
bkg(2,:) = linspace(-1, 1, size(mov_mc,3)).^2;  % quadratic term
mov_res = SeeResiduals(mov_res,mc);
mov_res = SeeResiduals(mov_res,mc.^2);
mov_res = SeeResiduals(mov_res,mc(:,1).*mc(:,3));
mov_res= SeeResiduals(mov_res,bkg,1);

nFrames2=size(mov_res,3);
 fGauss=fspecial('gaussian',20,7);
 mov_res_filter = imfilter(mov_res(bound:end-bound,bound:end-bound,:),fGauss,'same');
avgImg=Result.ref_im-imgaussfilt(Result.ref_im,6);
normImg=imgaussfilt(std(mov_res(bound:end-bound,bound:end-bound,:),0,3),5);
%refcorrImg=mat2gray(avgImg(bound:end-bound,bound:end-bound,:));
refcorrImg=mat2gray(std(mov_res(bound:end-bound,bound:end-bound,:),0,3));

%[~, refcorrImg] = GenerateDMD_dendrite(Result.ref_im(bound:end-bound,bound:end-bound),5,'Edge');
%mov_res_filter= (mov_res_filter)./F(bound:end-bound,bound:end-bound,:).*mat2gray(avgImg(bound:end-bound,bound:end-bound,:));
mov_res_filter= (mov_res_filter)./normImg.*refcorrImg;
%nexttile([1 1])
%imagesc(trsc(:,f_seg(f)+[frameInd(i)+tau]))
figure; writeMov_wKymo([fpath,'/mov_Res', num2str(i)],-mov_res_filter,[frameInd(i)+tau],10,1,[-0.1 0.7],[],trsc(:,f_seg(f):f_seg(f+1)-1))
end

%%
tr=Result.normTraces(2,:); t=[1:length(tr)];
spike=find_spike_bh(tr-movmedian(tr,300,2),4,3);
spike_time=find(spike); spike_height=tr(spike_time);

pass_frq=[15];
Fs=1000;

freq_lowhigh=pass_frq/(Fs/2);
[b, a] = butter(4, freq_lowhigh, 'low');
tr_pass = filtfilt(b, a, tr);
figure(2); clf; show_t=[750000:850000];
plot(t(show_t),tr(show_t)); hold all
plot(t(show_t),tr_pass(show_t));
plot(t(show_t),zeros(length(show_t),1));

[trans tr_trace]=detect_transient(tr_pass,[3 0.5],spike);
CS_ind=find(trans.spike_number>2 & trans.meanISI<40);
CS_trace=ismember(tr_trace,CS_ind);
figure(3); clf;
plot(t,tr)
hold all
tr_nan=tr; tr_nan(CS_trace==0)=NaN;
plot(t,tr_nan,'r')

CS_trace_spike=spike.*bwlabel(CS_trace);
for CS=1:max(CS_trace_spike)
    CS_time(CS)= find(CS_trace_spike==CS,1,'first');
    tmp = find(CS_trace_spike==CS,3,'first');
    CS2_time(CS) = tmp(2);
    CS3_time(CS) = tmp(3);
end

nBack = -20;
nFront = 20;
tau = nBack:nFront;
nTau = length(tau);

CS_tau = CS_time' + tau;
CS_spikeMat = reshape(Result.traces(:,CS_tau), 6, 122, nTau);

SS_trace_spike=spike.*(~CS_trace);
SS_time = find(SS_trace_spike);
SS_time(end) = [];
SS_time([find(diff(SS_time) < 45), find(diff(SS_time) < 45)+1]) = [];

plot(diff(SS_time))
nSS = length(SS_time);

SS_tau = SS_time' + tau;
SS_spikeMat = reshape(Result.traces(:,SS_tau), 6, nSS, nTau);

fileInd=ceil(SS_time/time_segment);
frameInd=mod(SS_time,time_segment);
%line up movie segments at simple spike
wnd=[-20:20];
SS_movie_align=[];
for j=unique(fileInd)
    mov_mc=double(readBinMov([fpath '/mc_ShutterReg' num2str(j,'%02d') '.bin'],sz(2),sz(1)));
    load([fpath '/mcTrace' num2str(j,'%02d') '.mat']);

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
    mov_res = SeeResiduals(mov_res,mc(:,1).*mc(:,3));
    mov_res= SeeResiduals(mov_res,bkg,1);

    for k=find(fileInd==j)
    try 
        SS_movie_align = cat(3,SS_movie_align, tovec(mov_res(:,:,frameInd(k)+wnd)));
    end
    end
end
%%
SS_movie_align=vm(SS_movie_align);
SS_movie_align.transpose.savebin([fpath '/SS_movie_align.bin'])
%%
SS_movie_align=double(SS_movie_align);
STA_movie_SS=-toimg(mean(SS_movie_align,3),sz(2),sz(1));
STA_trace_SS=squeeze(mean(SS_spikeMat(2,:,:),2));
interpolationFactor = 10; % Doubling the resolution in the third dimension

originalSize = size(STA_movie_SS);
interpolatedSize = [originalSize(1:2), originalSize(3) * interpolationFactor];

[Xq, Yq, Zq] = meshgrid(1:originalSize(2), 1:originalSize(1), linspace(1, originalSize(3), interpolatedSize(3)));

interpSpikeTA_mov = interp3(STA_movie_SS, Xq, Yq, Zq, 'linear');
%interpSpikeTA_mov=imgaussfilt3(interpSpikeTA_mov,[0.5 0.5 1]);
interSTA=interp1([1:length(STA_trace_SS)],STA_trace_SS,[1:originalSize(3) * interpolationFactor]/interpolationFactor);
 writeMov_wTrace([fpath ,'/SpikeTA_SS'],interpSpikeTA_mov,[],50,10,[-2 35],[],interSTA)


