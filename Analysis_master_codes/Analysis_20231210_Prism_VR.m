clear
clc;
cd '/Volumes/BHL_WD18TB/20230917_BHLm077_78_vrPrism/112801BHLm078_optopatch_VR'
fpath = '/Volumes/BHL_WD18TB/20230917_BHLm077_78_vrPrism/112801BHLm078_optopatch_VR';
time_segment=15000;

%% Motion correction

load(fullfile(fpath,"output_data.mat"))
sz=double(Device_Data{1, 3}.ROI([2 4]));
ref_time=[9000:10000]; overlap=200;

frm_end=max(Device_Data{1, 2}.Counter_Inputs(1, 1).data);
f_seg=[[1:time_segment:frm_end] frm_end+1]; f_seg(2:end)=f_seg(2:end)-1;

if length(ref_time)>2000
    mov_test=double(readBinMov_times([fpath '/frames1.bin'],sz(2),sz(1),[ref_time(1):ref_time(1)+2000]));
else
    mov_test=double(readBinMov_times([fpath '/frames1.bin'],sz(2),sz(1),[ref_time(1):ref_time(end)]));
end
mov_test=rollingShutter_correction(mov_test,Device_Data{1, 3}.exposuretime,'fusion');
mov_test=mov_test(:,:,2:end);
[mov_test_mc,xyField]=optical_flow_motion_correction_LBH(mov_test,mean(mov_test,3),'normcorre');
mov_test=vm(mov_test);
mov_test = single(mov_test)./single(max(mov_test.data(:)));
mov_test = movmean(mov_test,10,3);
mov_ref = squeeze(median(mov_test,3));

for j=1%:length(f_seg)-1
    try
        mov=double(readBinMov_times([fpath '/frames1.bin'],sz(2),sz(1),[f_seg(j):f_seg(j+1)+overlap]));
    catch % when the image ends
        mov=double(readBinMov_times([fpath '/frames1.bin'],sz(2),sz(1),[f_seg(j):f_seg(j+1)]));
    end

    mov=rollingShutter_correction(mov,Device_Data{1, 3}.exposuretime,'fusion');
    mov=vm(mov(:,:,2:end));
    if j==1
        mov=mov(:,:,[1 1:size(mov,3)]);
    end

    [mov_mc,xyField]=optical_flow_motion_correction_LBH(mov,mov_ref,'normcorre');

    ave_im=mean(mov_mc,3);
    mov_mc=vm(mov_mc);
    mov_mc.transpose.savebin([fpath '/mc_ShutterReg' num2str(j,'%02d') '.bin'])

    %        mcTrace = squeeze(mean(xyField,[1 2])); %optic flow
    mcTrace=xyField; % Normcorre
    save([fpath '/mcTrace' num2str(j,'%02d') '.mat'],'mcTrace','ave_im')

    %  clear mov_mc mov
end


%%
disp(fpath)
DAQ_rate=0.000005;
load([fpath '/output_data.mat'])
sz=double(Device_Data{1, 3}.ROI([2 4]));
ref_time=[2000:3000];
mov_test=double(readBinMov_times([fpath '/mc_ShutterReg' num2str(49,'%02d') '.bin'],sz(2),sz(1),ref_time));
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

load(fullfile(fpath,'mcTrace49.mat'));
frm_end=max(Device_Data{1, 2}.Counter_Inputs(1, 1).data);
Result.Blue=Device_Data{1, 2}.buffered_tasks(1, 1).channels(1, 2).data;
Result.Reward=Device_Data{1, 2}.buffered_tasks(1, 3).channels.data;
CamCounter=Device_Data{1, 2}.Counter_Inputs(1, 1).data;
CamTrigger=find(CamCounter(2:end)-CamCounter(1:end-1));

Result.Blue=Result.Blue(CamTrigger); Result.Reward=Result.Reward(CamTrigger);
frm_rate=double((CamTrigger(2)-CamTrigger(1))*DAQ_rate);
mov_mc=double(readBinMov([fpath '/mc_ShutterReg' num2str(49,'%02d') '.bin'],sz(2),sz(1)));
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
figure(4); clf;
imshow2(squeeze(sum(toimg(Result.ftprnt,sz(2),sz(1)).*reshape(jet(Npoly),1,1,[],3),3)),[]);

%% Signal extraction
Result.traces=[];
Result.traces_res=[];
Result.mcTrace=[];
Result.im_corr=[];
bound=5;
ref_im_vec=tovec(Result.ref_im(bound:end-bound,bound:end-bound));
f_seg=[[1:time_segment:frm_end] frm_end+1]; f_seg(2:end)=f_seg(2:end)-1;
take_window=repmat([1 time_segment],length(f_seg)-1,1);
take_window(2:end,1)=take_window(2:end,1)+overlap; take_window(1:end-1,2)=take_window(1:end-1,2)+overlap;
take_window(end)=mod(f_seg(end),time_segment);


for j=1:length(f_seg)-1
    j
    mov_mc=double(readBinMov([fpath '/mc_ShutterReg' num2str(j,'%02d') '.bin'],sz(2),sz(1)));
    load([fpath '/mcTrace' num2str(j,'%02d') '.mat']);

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
    mov_res = SeeResiduals(mov_res,mc(:,1).*mc(:,3));
    mov_res= SeeResiduals(mov_res,bkg,1);

    Result.traces=[Result.traces -(tovec(mov_res)'*tovec(Result.ftprnt))'];
    Result.mcTrace=[Result.mcTrace; mc];
    Result.im_corr=[Result.im_corr sum(mov_mc_vec.*ref_im_vec,1)/(size(mov_mc_vec,1)-1)];  %image correlation

end

Result.traces=Result.traces(:,1:length(CamTrigger)-1);
Result.mcTrace=Result.mcTrace(1:length(CamTrigger)-1,:);

save(fullfile(fpath,'Result_20231207.mat'),'Result','fpath','-v7.3')

%% Load Virmen data

fid =   fopen(fullfile(fpath,'2309171127_BHLm078_optopatch_VR_virmenLog.data'));
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

save(fullfile(fpath,'Result_20231207.mat'),'Result','fpath','-v7.3')
%%
load(fullfile(fpath,'Result_20231207.mat'))
%% Clean up and norm
exclude_frq=[241.7 242]; %monitor
%exclude_frq2=[483.5 484]; %monitor
exclude_frq2=[25 65.7]; %motion
    time_bin=15000; Fs=1000; ref_trace=2; %2nd trunk is the reliable trace

    nTime=size(Result.traces,2);
    nROI=size(Result.traces,1);
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


    figure(2); clf;
    [p f]=fft_simple(Result.traces(2,:),1000);
    ax1=nexttile([1 1]);
    plot(f(1,:),p')
    set(gca,'YScale','log')
    
    [p f]=fft_simple(Result.mcTrace,1000);
    ax2=nexttile([1 1]);
    plot(f(1,:),p')
    set(gca,'YScale','log')
    linkaxes([ax1 ax2],'x')


clear traces_res_filtered noise noise_intp norm_trace sp_height SpHeight_intp sp_time tr_mc_imcorr tr_mc tr
    tN=[1:time_bin:nTime]; tN=[tN nTime];
    sp_time=zeros(nROI,nTime);
    sp_height=zeros(nROI,nTime);

mcTrace=squeeze(Result.mcTrace)';

for n=1:size(Result.traces,1)
    tr=Result.traces(n,1:nTime);

    % regress out motion frequency
    for t=1:length(tN)-1
      tr(tN(t):tN(t+1))=squeeze(SeeResiduals(reshape(tr(tN(t):tN(t+1)),1,1,[]),mcTrace(:,(tN(t):tN(t+1)))));
                tr(tN(t):tN(t+1))=squeeze(SeeResiduals(reshape(tr(tN(t):tN(t+1)),1,1,[]),mcTrace(:,(tN(t):tN(t+1))).^2));
                tr(tN(t):tN(t+1))=squeeze(SeeResiduals(reshape(tr(tN(t):tN(t+1)),1,1,[]),mcTrace(1,(tN(t):tN(t+1))).*mcTrace(3,(tN(t):tN(t+1)))));

                tr_mc(tN(t):tN(t+1))=tr(tN(t):tN(t+1));

                imcorr_tmp=Result.im_corr(tN(t):tN(t+1));
                imcorr_tmp=movmean(imcorr_tmp,5,'omitnan');
    end

    % regress out motion frequency
    traces_res_filtered(n,:) = filtfilt(b, a, tr);
    tr_mc_bpMonitor=traces_res_filtered(n,:);
    traces_res_filtered(n,:) = filtfilt(b2, a2, traces_res_filtered(n,:));
    tr_mc_bpMonitorMotion=traces_res_filtered(n,:);

    norm_trace(n,:)=traces_res_filtered(n,:);%-movmedian(traces_res_filtered(n,:),500,2);

    

    if n==ref_trace
    for t=1:length(tN)-1
        tr_tmp=norm_trace(n,tN(t):tN(t+1));
            tr_tmp=tr_tmp-movmedian(tr_tmp,300);
            noise(t,n)=get_threshold(tr_tmp,1);
            tr_tmp_norm=tr_tmp./noise(t,n);
            [sp_temp, pks, prom]=find_spike_bh(tr_tmp_norm,4,4);
            sp_time(n,tN(t):tN(t+1))=sp_temp;
    end

    

    t_fit=find(sp_time(n,:));
    sp_height(n,t_fit)=norm_trace(n,t_fit);
    tx=tN(1:end-1)+time_bin/2;
    %noise_intp(n,:)=movmean(interp1(tx,noise(:,n),[1:size(Result{i}.traces,2)],'linear','extrap'),10000);

    
     [Ny_fit t_consts coeffY]  = expfitDM_2(tx(~isnan(noise(1:end-1,n)))',noise(~isnan(noise(1:end-1,n)),n),[1:nTime]',10^7);
     noise_intp=Ny_fit;
    %[y_fit t_consts coeffY]  = expfitDM_2(tx(~isnan(sp_height(1:end-1,n)))',sp_height(~isnan(sp_height(1:end-1,n)),n),[1:size(Result{i}.traces,2)]',10^7);
    [y_fit t_consts coeffY]  = expfitDM_2(t_fit',sp_height(n,t_fit)',[1:nTime]',10^7);
    SpHeight_intp=y_fit;


    figure(3); clf;
    ax1=nexttile([1 1]);
    plot(rescale2([Result.traces(n,1:nTime);tr_mc;tr_mc_bpMonitorMotion;mcTrace(1,1:nTime)],2)')%+[1:4])
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
    Result.normTraces=norm_trace./get_threshold(norm_trace,1);
    Result.spike=find_spike_bh(Result.normTraces-movmedian(Result.normTraces,300,2),5,3);

    
for n=1:size(Result.traces,1)
Result.subThreshold(n,:) = filtfilt(b3, a3, Result.normTraces(n,:));
Result.theta(n,:) = filtfilt(b4, a4, Result.normTraces(n,:));
end
% 
 save(fullfile(fpath,'Result_20231028.mat'),'Result','fpath','-v7.3')
% 
% %save(fullfile(fpath,'Result_20230928.mat'),'Result','fpath','-v7.3')

%% Show segment

t_interest=[401935:406680];
tr_soma=rescale(Result.normTraces([2],1:end-10)); tr_dd=rescale(Result.normTraces([6],1:end-10));
figure; clf;    
plot(tr_soma(t_interest),tr_dd(t_interest))

%% PCA/ICA
load(fullfile(fpath,'Result_20231004.mat'))

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

%% detect complex spikes
noi=[2];
tr=Result.normTraces(noi,:); t=[1:length(tr)];
spike=find_spike_bh(tr-movmedian(tr,300,2),4,3);
spike_time=find(spike); spike_height=tr(spike_time);

tr_sub=get_subthreshold(tr,Result.spike(2,:),5,15);
figure(2); clf; show_t=[310000:320000];
nexttile([1 1])
plot(t(show_t),tr(show_t)); hold all
plot(t(show_t),tr_sub(show_t));
plot(t(show_t),zeros(length(show_t),1));
legend({'F trace','Subthreshold'})

[trans tr_trace]=detect_transient(tr_sub,[6 0.5],spike);
CS_ind=find(trans.spike_number>2 & trans.mean_ISI<40);
CS_trace=ismember(tr_trace,CS_ind);
nexttile([1 1])
plot(t,tr)
hold all
tr_nan=tr; tr_nan(CS_trace==0)=NaN;
plot(t,tr_nan,'r')
legend({'F trace','Complex Spikes'})

clear CS_time
CS_trace_spike=spike.*bwlabel(CS_trace);
for CS=1:max(CS_trace_spike)
CS_time(CS)= find(CS_trace_spike==CS,1,'first');
end

fileInd=ceil(CS_time/time_segment);
frameInd=mod(CS_time,time_segment);

%% Show Complex spikes
figure;
tiledlayout(1,2)
ax2=nexttile([1 1]);
line=plot(rescale2(reshape(Result.traces(2,CS_time'+[-100:200])',length(CS_time),[]),2)'+[1:length(CS_time)]);
arrayfun(@(l,c) set(l,'Color',c{:}),line,num2cell(jet(57),2))
title('Original Trace')
ax1=nexttile([1 1]);
line=plot(rescale2(reshape(Result.normTraces(2,CS_time'+[-100:200])',length(CS_time),[]),2)'+[1:length(CS_time)]);
arrayfun(@(l,c) set(l,'Color',c{:}),line,num2cell(jet(57),2))
title('Filtered Trace')
linkaxes([ax1 ax2],'xy')


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


