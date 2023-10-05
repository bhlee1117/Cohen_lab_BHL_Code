clear
clc;
cd '/Volumes/BHL_WD18TB/20230917_BHLm077_78_vrPrism'
[fpath] = uigetfile_n_dir();
time_segment=15000;
foi=1;
%% Motion correction

load(fullfile(fpath{foi},"output_data.mat"))
sz=double(Device_Data{1, 3}.ROI([2 4]));
ref_time=[8000:9100];

frm_end=max(Device_Data{1, 2}.Counter_Inputs(1, 1).data);
f_seg=[[1:time_segment:frm_end] frm_end+1]; f_seg(2:end)=f_seg(2:end)-1;

if length(ref_time)>1000
    mov_test=double(readBinMov_times([fpath{foi} '/frames1.bin'],sz(2),sz(1),[ref_time(1):ref_time(1)+1000]));
else
    mov_test=double(readBinMov_times([fpath{foi} '/frames1.bin'],sz(2),sz(1),[ref_time(1):ref_time(end)]));
end

options_rigid = NoRMCorreSetParms('d1',size(mov_test,1),'d2',size(mov_test,2),'bin_width',200,'max_shift',30,'us_fac',50,'init_batch',200);
tic; [mov_test,shifts1,template1,options_rigid] = normcorre(mov_test,options_rigid); toc
mov_test=vm(mov_test);
mov_test = single(mov_test)./single(max(mov_test.data(:)));
mov_test = movmean(mov_test,10,3);
mov_ref = squeeze(median(mov_test,3));

for j=1:length(f_seg)-1
    try
        mov=double(readBinMov_times([fpath{foi} '/frames1.bin'],sz(2),sz(1),[f_seg(j):f_seg(j+1)+10]));
    catch % when the image ends
        mov=double(readBinMov_times([fpath{foi} '/frames1.bin'],sz(2),sz(1),[f_seg(j):f_seg(j+1)]));
    end

    mov=rollingShutter_correction(mov,Device_Data{1, 3}.exposuretime,'fusion');
    mov=vm(mov(:,:,2:end));
    if j==1
        mov=mov(:,:,[1 1:size(mov,3)]);
    end

    [mov_mc,xyField]=optical_flow_motion_correction_LBH(mov,mov_ref,'normcorre');

    ave_im=mean(mov_mc,3);
    mov_mc=vm(mov_mc);
    mov_mc.transpose.savebin([fpath{foi} '/mc_ShutterReg' num2str(j,'%02d') '.bin'])

    %        mcTrace = squeeze(mean(xyField,[1 2])); %optic flow
    mcTrace=xyField; % Normcorre
    save([fpath{foi} '/mcTrace' num2str(j,'%02d') '.mat'],'mcTrace','ave_im')

    %  clear mov_mc mov
end

%% Set ROIs
disp(fpath{foi})
DAQ_rate=0.000005;
load([fpath{foi} '/output_data.mat'])
sz=double(Device_Data{1, 3}.ROI([2 4]));
ref_time=[2000:3000];
mov_test=double(readBinMov_times([fpath{foi} '/mc_ShutterReg' num2str(50,'%02d') '.bin'],sz(2),sz(1),ref_time));
avgImg=mean(mov_test,3);
figure(3); clf;
imshow2(avgImg,[])
g=1; coord=[];
while g
    [x y]=ginput(1);
    if isempty(x)
        g=0;
    end
    coord=[coord; [x y]];
    hold all
    plot(x,y,'ro')
end
close(figure(3));
Result.centers=coord;

load(fullfile(fpath{foi},'mcTrace50.mat'));
frm_end=max(Device_Data{1, 2}.Counter_Inputs(1, 1).data);
Result.Blue=Device_Data{1, 2}.buffered_tasks(1, 1).channels(1, 2).data;
Result.Reward=Device_Data{1, 2}.buffered_tasks(1, 3).channels.data;
CamCounter=Device_Data{1, 2}.Counter_Inputs(1, 1).data;
CamTrigger=find(CamCounter(2:end)-CamCounter(1:end-1));

Result.Blue=Result.Blue(CamTrigger); Result.Reward=Result.Reward(CamTrigger);
frm_rate=double((CamTrigger(2)-CamTrigger(1))*DAQ_rate);
sz=double(Device_Data{1, 3}.ROI([2 4]));
mov_mc=double(readBinMov([fpath{foi} '/mc_ShutterReg' num2str(50,'%02d') '.bin'],sz(2),sz(1)));

Result.ref_im=mean(mov_mc,3);

mov_res= mov_mc-mean(mov_mc,3);
mov_res = SeeResiduals(mov_res,mcTrace.xymean);
mov_res = SeeResiduals(mov_res,mcTrace.xymean.^2);
mov_res = SeeResiduals(mov_res,mcTrace.xymean(:,1).*mcTrace.xymean(:,2));

% calculate footprint
if frm_end<10000; end_frame=size(mov_mc,3); else end_frame=time_segment; end
Result.c_ftprnt=mask_footprint(Result.centers,movmean(mov_res(:,:,1000:end),10,3),[],10);
for n=1:size(Result.c_ftprnt,3)
    Result.c_ftprnt(:,:,n)=imgaussfilt(Result.c_ftprnt(:,:,n),0.6);
end
N=size(Result.c_ftprnt,3);
Result.coord=get_coord(Result.c_ftprnt);
figure;  show_footprnt(Result.c_ftprnt,Result.ref_im)

%% Signal extraction
Result.traces=[];
Result.traces_res=[];
Result.mcTrace=[];
Result.im_corr=[];
bound=10;
ref_im_vec=tovec(Result.ref_im(bound:end-bound,bound:end-bound));

for j=1:length(f_seg)-1
    mov_mc=double(readBinMov([fpath{foi} '/mc_ShutterReg' num2str(j,'%02d') '.bin'],sz(2),sz(1)));
    load([fpath{foi} '/mcTrace' num2str(j,'%02d') '.mat']);

    try mov_mc=mov_mc(:,:,1:end_frame); catch mov_mc=mov_mc; end
    try mc=mcTrace.xymean(1:end_frame,:); catch mc=mcTrace.xymean; end
    mov_mc_vec=tovec(mov_mc(bound:end-bound,bound:end-bound,:));
    mov_mc_vec=(mov_mc_vec-mean(mov_mc_vec,1))./std(mov_mc_vec,0,1);

    mov_res= mov_mc-mean(mov_mc,3);
    bkg = zeros(2, size(mov_mc,3));
    bkg(1,:) = linspace(-1, 1, size(mov_mc,3));  % linear term
    bkg(2,:) = linspace(-1, 1, size(mov_mc,3)).^2;  % quadratic term
    mov_res = SeeResiduals(mov_res,mc);
    mov_res = SeeResiduals(mov_res,mc.^2);
    mov_res = SeeResiduals(mov_res,mc(:,1).*mc(:,2));
    mov_res= SeeResiduals(mov_res,bkg,1);

    Result.traces=[Result.traces -(tovec(mov_res)'*tovec(Result.c_ftprnt))'];
    Result.mcTrace=[Result.mcTrace; mcTrace.xymean];
    Result.im_corr=[Result.im_corr sum(mov_mc_vec.*ref_im_vec,1)/(size(mov_mc_vec,1)-1)];

end

Result.traces=Result.traces(:,1:length(CamTrigger));
Result.mcTrace=Result.mcTrace(1:length(CamTrigger),:);

save(fullfile(fpath{foi},'Result_20230928.mat'),'Result','fpath','-v7.3')

%% Load Virmen data

fid =   fopen(fullfile(fpath{1},'2309171127_BHLm078_optopatch_VR_virmenLog.data'));
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



%% Clean up and bleach correction
exclude_frq=[20 70];
time_bin=10000; Fs=1000;

freq_lowhigh=exclude_frq/(Fs/2);
[b, a] = butter(4, freq_lowhigh, 'stop');

clear traces_res_filtered noise noise_intp norm_trace
tN=[1:time_bin:size(Result.traces,2)]; tN=[tN size(Result.traces,2)];

mcTrace=squeeze(Result.mcTrace)';

for n=1:size(Result.traces,1)
    tr=Result.traces(n,:);

    % regress out motion frequency
    for t=1:length(tN)-1
        tr(tN(t):tN(t+1))=squeeze(SeeResiduals(reshape(tr(tN(t):tN(t+1)),1,1,[]),mcTrace(:,(tN(t):tN(t+1)))));
        tr(tN(t):tN(t+1))=squeeze(SeeResiduals(reshape(tr(tN(t):tN(t+1)),1,1,[]),mcTrace(:,(tN(t):tN(t+1))).^2));
        tr(tN(t):tN(t+1))=squeeze(SeeResiduals(reshape(tr(tN(t):tN(t+1)),1,1,[]),mcTrace(1,(tN(t):tN(t+1))).*mcTrace(2,(tN(t):tN(t+1)))));
    end

    % regress out motion frequency
    traces_res_filtered(n,:) = filtfilt(b, a, tr);
    norm_trace(n,:)=traces_res_filtered(n,:);%-movmedian(traces_res_filtered(n,:),500,2);

    for t=1:length(tN)-1
        tr_tmp=norm_trace(n,tN(t):tN(t+1));
        tr_tmp=tr_tmp-movmedian(tr_tmp,100);
        noise(t,n)=get_threshold(tr_tmp,1);
        [sp_temp, pks, prom]=find_spike_bh(tr_tmp,4,3);
        sp_height(t,n)=median(pks);%/noise(t,n);
    end
    tx=tN(1:end-1)+time_bin/2;
    noise_intp(n,:)=interp1(tx,noise(:,n),[1:size(Result.traces,2)],'linear','extrap');
    [y_fit t_consts coeffY]  = expfitDM_2(tx(~isnan(sp_height(:,n)))',sp_height(~isnan(sp_height(:,n)),n),[1:size(Result.traces,2)]',10^7);
    SpHeight_intp(n,:)=y_fit;
end

norm_trace=norm_trace./SpHeight_intp;
Result.normTraces=norm_trace;
Result.spike=find_spike_bh(Result.normTraces-movmedian(Result.normTraces,300,2),5,3);

%save(fullfile(fpath{foi},'Result_20230928.mat'),'Result','fpath','-v7.3')
%% Generate representative dF movie
%load(fullfile(fpath{1},'Result_20230928.mat'))
f=34;   bound=5;
load([fpath{foi} '/output_data.mat'])
sz=double(Device_Data{1, 3}.ROI([2 4]));
mov_mc=double(readBinMov([fpath{foi} '/mc_ShutterReg' num2str(f,'%02d') '.bin'],sz(2),sz(1)));
load([fpath{foi} '/mcTrace' num2str(f,'%02d') '.mat']);

try mov_mc=mov_mc(:,:,1:end_frame); catch mov_mc=mov_mc; end
try mc=mcTrace.xymean(1:end_frame,:); catch mc=mcTrace.xymean; end

mov_res= mov_mc-mean(mov_mc,3);
bkg = zeros(2, size(mov_mc,3));
bkg(1,:) = linspace(-1, 1, size(mov_mc,3));  % linear term
bkg(2,:) = linspace(-1, 1, size(mov_mc,3)).^2;  % quadratic term
mov_res = SeeResiduals(mov_res,mc);
mov_res = SeeResiduals(mov_res,mc.^2);
mov_res = SeeResiduals(mov_res,mc(:,1).*mc(:,2));
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
    plot(rescale(eigTraces(i,:)')+i-0.5)
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
ReconMov=toimg((coeffs*(V(:,keep_ind)*D(keep_ind))',size(datDS,1),size(datDS,2));


ReconMovICA=toimg(footPrintsVec(:,keep_ind_ics)*ics(:,keep_ind_ics)',sz(2),sz(1));
%ReconMov=toimg(tovec(mov_res).*mean(footPrintsVec(:,keep_ind_ics),2),sz(2),sz(1));
%ReconMov_filt=-imgaussfilt3(ReconMov,[2 2 0.5]);
figure; writeMov_wTrace([fpath{1},'/dFMov_mov34PCA'],-ReconMov(bound:end-bound,bound:end-bound,:),[4000:6250],10,1,[-5 35],[],mean(Result.traces([4 5 6 7],f_seg(f):f_seg(f+1)-1),1))
figure; writeMov_wTrace([fpath{1},'/dFMov_mov34ICA'],-(ReconMovICA(bound:end-bound,bound:end-bound,:)),[4000:6250],10,1,[],[],mean(Result.traces([4 5 6 7],f_seg(f):f_seg(f+1)-1),1))

%%
f=34;   bound=5;
load([fpath{foi} '/output_data.mat'])
sz=double(Device_Data{1, 3}.ROI([2 4]));
mov_mc=double(readBinMov([fpath{foi} '/mc_ShutterReg' num2str(f,'%02d') '.bin'],sz(2),sz(1)));
load([fpath{foi} '/mcTrace' num2str(f,'%02d') '.mat']);

try mov_mc=mov_mc(:,:,1:end_frame); catch mov_mc=mov_mc; end
try mc=mcTrace.xymean(1:end_frame,:); catch mc=mcTrace.xymean; end

mov_res= mov_mc-mean(mov_mc,3);
bkg = zeros(2, size(mov_mc,3));
bkg(1,:) = linspace(-1, 1, size(mov_mc,3));  % linear term
bkg(2,:) = linspace(-1, 1, size(mov_mc,3)).^2;  % quadratic term
mov_res = SeeResiduals(mov_res,mc);
mov_res = SeeResiduals(mov_res,mc.^2);
mov_res = SeeResiduals(mov_res,mc(:,1).*mc(:,2));
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

%% trace normalization using spike height
noi=[3 4 5 7];
tr=mean(Result.normTraces(noi,:)); t=[1:length(tr)];
spike=find_spike_bh(tr-movmedian(tr,300,2),4,3);
spike_time=find(spike); spike_height=tr(spike_time);

figure; clf;

[SpH_ini bin_ini]=histcounts(spike_height(spike_time<600000 & spike_time>40000),100,'normalization','pdf');
[SpH_end bin_end]=histcounts(spike_height(spike_time>800000),100,'normalization','pdf');
plot(mean([bin_ini(1:end-1); bin_ini(2:end)]),SpH_ini); hold all
plot(mean([bin_end(1:end-1); bin_end(2:end)]),SpH_end);

pass_frq=[1 15];
Fs=1000;

freq_lowhigh=pass_frq/(Fs/2);
[b, a] = butter(4, freq_lowhigh, 'bandpass');
tr_pass = filtfilt(b, a, tr);
figure(2); clf;
plot(t,tr,t,tr_pass,t,zeros(size(t,1)));

[trans tr_trace]=detect_transient(tr_pass,[2 1]);
figure(3); clf;
plot(t,tr)
hold all
tr_nan=tr; tr_nan(tr_trace==0)=NaN;
plot(t,tr_nan,'r')

[CS_list CS]=find_CS(tr,spike,15,5);
figure(3); clf;
plot(t,tr,t(find(CS>0)),tr(find(CS>0)),'color','r')
    
%%

[Lap_FR Lap_V]=PlaceTrigger_average(sum(Result.spike(noi,:)),50,Result.VR,-0.1,115);
[LickFR]=PlaceTrigger_average(Result.VR(9,:),50,Result.VR,-0.1,115);

