clear
clc;
cd '/Volumes/BHL_WD18TB/20230917_BHLm077_78_vrPrism/111838BHLm078_optopatch_pulse'
fpath = '/Volumes/BHL_WD18TB/20230917_BHLm077_78_vrPrism/111838BHLm078_optopatch_pulse';

%% Motion correction

load(fullfile(fpath,"output_data.mat"))
sz=double(Device_Data{1, 3}.ROI([2 4]));
ref_time=[9000:10000]; overlap=200;

frm_end=max(Device_Data{1, 2}.Counter_Inputs(1, 1).data);

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
CamCounter=Device_Data{1, 2}.Counter_Inputs(1, 1).data;
CamTrigger=find(CamCounter(2:end)-CamCounter(1:end-1));

mov=double(readBinMov([fpath '/frames1.bin'],sz(2),sz(1)));
mov=mov(:,:,1:length(CamTrigger));
mov=rollingShutter_correction(mov,Device_Data{1, 3}.exposuretime,'fusion');
mov=vm(mov(:,:,2:end));
[mov_mc,xyField]=optical_flow_motion_correction_LBH(mov,mov_ref,'normcorre');

ave_im=mean(mov_mc,3);
mov_mc=vm(mov_mc);
mcTrace=xyField;
mov_mc.transpose.savebin([fpath '/mc_ShutterReg' num2str(1,'%02d') '.bin'])
save([fpath '/mcTrace' '.mat'],'mcTrace','ave_im')


%%
disp(fpath)
DAQ_rate=0.000005;
load([fpath '/output_data.mat'])
sz=double(Device_Data{1, 3}.ROI([2 4]));
ref_time=[2000:3000];
mov_test=double(readBinMov_times([fpath '/mc_ShutterReg' num2str(1,'%02d') '.bin'],sz(2),sz(1),ref_time));
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

load(fullfile(fpath,'mcTrace.mat'));
frm_end=max(Device_Data{1, 2}.Counter_Inputs(1, 1).data);
Result.Blue=Device_Data{1, 2}.buffered_tasks(1, 1).channels(1, 2).data;
Result.Reward=Device_Data{1, 2}.buffered_tasks(1, 3).channels.data;
CamCounter=Device_Data{1, 2}.Counter_Inputs(1, 1).data;
CamTrigger=find(CamCounter(2:end)-CamCounter(1:end-1));

Result.Blue=Result.Blue(CamTrigger); Result.Reward=Result.Reward(CamTrigger);
Result.Blue=Result.Blue(2:end);
frm_rate=double((CamTrigger(2)-CamTrigger(1))*DAQ_rate);
mov_mc=double(readBinMov([fpath '/mc_ShutterReg' num2str(1,'%02d') '.bin'],sz(2),sz(1)));
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

    mov_mc=double(readBinMov([fpath '/mc_ShutterReg' num2str(1,'%02d') '.bin'],sz(2),sz(1)));
    load([fpath '/mcTrace' '.mat']);

        mov_mc=mov_mc;
        mc=mcTrace.xymean;
 
    mov_mc_vec=tovec(mov_mc(bound:end-bound,bound:end-bound,:)); 
    mov_mc_vec=(mov_mc_vec-mean(mov_mc_vec,1))./std(mov_mc_vec,0,1);

    

    mov_res= mov_mc-mean(mov_mc,3);
    bkg = zeros(2, size(mov_mc,3));
    bkg(1,:) = linspace(-1, 1, size(mov_mc,3));  % linear term
    bkg(2,:) = linspace(-1, 1, size(mov_mc,3)).^2;  % quadratic term
    bkg(3,:)=movmedian(get_blueoffTrace(squeeze(mean(mov_mc,[1 2])),Result.Blue,30),1000);

    mov_res = SeeResiduals(mov_res,mc);
    mov_res = SeeResiduals(mov_res,mc.^2);
    mov_res = SeeResiduals(mov_res,mc(:,1).*mc(:,3));
    mov_res= SeeResiduals(mov_res,bkg,1);

    Result.traces=[Result.traces -(tovec(mov_res)'*tovec(Result.ftprnt))'];
    Result.mcTrace=[Result.mcTrace; mc];
    Result.im_corr=[Result.im_corr sum(mov_mc_vec.*ref_im_vec,1)/(size(mov_mc_vec,1)-1)];  %image correlation

Result.traces=Result.traces(:,1:length(CamTrigger)-1);
Result.mcTrace=Result.mcTrace(1:length(CamTrigger)-1,:);

save(fullfile(fpath,'Result_optopatch_20231028.mat'),'Result','fpath','-v7.3')

%% Clean up and norm
exclude_frq=[241.7 242]; %monitor
%exclude_frq2=[483.5 484]; %monitor
exclude_frq2=[25 65.7]; %motion
    time_bin=1000; Fs=1000; ref_trace=2; %2nd trunk is the reliable trace

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
 save(fullfile(fpath,'Result_optopatch_20231028.mat'),'Result','fpath','-v7.3')
%%
load(fullfile(fpath,'Result_optopatch_20231028.mat'))
%%
noi=[2]; nTau=[-30:60];
tr=Result.normTraces(noi,:); t=[1:length(tr)];
spike=find_spike_bh(tr-movmedian(tr,300,2),3,2);
Blue=Result.Blue;
blueOff = Blue == 0;
blueOff2 = imerode(blueOff, [ones(1,20), zeros(1, 20)]);
Blue_di=~blueOff2;

bwBlue_di=bwlabel(Blue_di);
sp_pulse=[];
for b=1:max(bwBlue_di)
t_tmp=find(bwBlue_di==b);
sp_pulse=[sp_pulse find(spike(t_tmp),1,'first')+t_tmp(1)-1];
end     
Pulse_Tau=sp_pulse(1:end-1)'+nTau;
Pulse_spikeMat=reshape(Result.normTraces(:,Pulse_Tau),4,length(sp_pulse)-1,[]);
Pulse_STA=squeeze(mean(Pulse_spikeMat,2));
Pulse_STA=Pulse_STA-prctile(Pulse_STA, 5, 2);

%% STA movie Short Pulse

load(fullfile(fpath,"output_data.mat"))
frm_end=max(Device_Data{1, 2}.Counter_Inputs(1, 1).data);
f_seg=[[1:time_segment:frm_end] frm_end+1]; f_seg(2:end)=f_seg(2:end)-1;
take_window=repmat([1 time_segment],length(f_seg)-1,1);
take_window(2:end,1)=take_window(2:end,1)+overlap; take_window(1:end-1,2)=take_window(1:end-1,2)+overlap;
take_window(end)=mod(f_seg(end),time_segment);
sz=double(Device_Data{1, 3}.ROI([2 4]));
    mov_mc=double(readBinMov([fpath '/mc_ShutterReg' num2str(1,'%02d') '.bin'],sz(2),sz(1)));
    load([fpath '/mcTrace' '.mat']);

        mov_mc=mov_mc;
        mc=mcTrace.xymean;
 
    mov_mc_vec=tovec(mov_mc(bound:end-bound,bound:end-bound,:)); 
    mov_mc_vec=(mov_mc_vec-mean(mov_mc_vec,1))./std(mov_mc_vec,0,1);

    

    mov_res= mov_mc-mean(mov_mc,3);
    bkg = zeros(2, size(mov_mc,3));
    bkg(1,:) = linspace(-1, 1, size(mov_mc,3));  % linear term
    bkg(2,:) = linspace(-1, 1, size(mov_mc,3)).^2;  % quadratic term
    bkg(3,:)=movmedian(get_blueoffTrace(squeeze(mean(mov_mc,[1 2])),Result.Blue,30),1000);

    mov_res = SeeResiduals(mov_res,mc);
    mov_res = SeeResiduals(mov_res,mc.^2);
    mov_res = SeeResiduals(mov_res,mc(:,1).*mc(:,3));
    mov_res= SeeResiduals(mov_res,bkg,1);


%% Spike height Analysis

Result.spike=find_spike_bh(Result.normTraces-movmedian(Result.normTraces,50,2),2.5,0.5);
Result.spike(1,:) = remove_closeSpike(Result.spike(1,:),Result.normTraces(1,:),4);
Blue=Result.Blue>0; delay_max=3; delay_front=2;

blueOff = Blue == 0;
blueOff2 = imerode(blueOff, [ones(1,20), zeros(1, 20)]);
Blue_di=~blueOff2;

bwBlue=bwlabel(Blue_di);
nback=20; nfront=20; nSpike=4;
g=0;

nROI=size(Result.normTraces,1);
spike_by_blue=Blue_di.*Result.spike(1,:);
spike_list=NaN(nROI,nSpike,max(bwBlue));
spike_trace=zeros(size(Result.spike,2),nROI);
DelayMat=NaN(nROI,nSpike,max(bwBlue));
AmpMat=NaN(nROI,nSpike,max(bwBlue));
LocAmpMat=NaN(nROI,nSpike,max(bwBlue));

for b=1:max(bwBlue)
t_tmp=find(bwBlue==b);
soma_spike=find(spike_by_blue(t_tmp));
trace_segment=Result.normTraces(:,t_tmp);

soma_spike(find(soma_spike<1 | soma_spike>length(t_tmp)-delay_max))=[];

dend_delay=[];
for s=1:length(soma_spike)
[~, dend_delay(s,:)]=max(trace_segment(:,soma_spike(s)+[-delay_front:delay_max]),[],2);
end

if length(soma_spike)<nSpike
spike_list(:,1:length(soma_spike),b)=dend_delay' + soma_spike - delay_front - 1 + t_tmp(1)-1;    
%DelayMat(:,1:length(soma_spike),b)=(dend_delay-min(dend_delay,[],2))';
if ~isempty(dend_delay)
DelayMat(:,1:length(soma_spike),b)=(dend_delay-dend_delay(:,1))';
end
else
spike_list(:,1:nSpike,b)=dend_delay(1:nSpike,:)' + soma_spike(1:nSpike) - delay_front -1 + t_tmp(1)-1;
DelayMat(:,1:nSpike,b)=(dend_delay(1:nSpike,:)-dend_delay(1:nSpike,1))';
end

end

sp_list_vec=reshape(spike_list,nROI,[]);
sp_list_vec(:,isnan(sum(sp_list_vec,1)))=[];
spike_trace((sp_list_vec+size(Result.spike,2)*[0:nROI-1]'))=1;
spike_trace=spike_trace';
%show_traces_spikes(Result.normTraces,spike_trace,double(Blue_di))



tr_seg_subthres=zeros(nROI,size(Result.normTraces,2));
for n=1:nROI
tr_seg_subthres(n,:)=get_subthreshold(Result.normTraces(n,:),spike_trace(n,:),round(5+n*1.3),4);
end

tr_local=Result.normTraces-tr_seg_subthres;
%show_traces_spikes_subth(Result.normTraces,spike_trace,tr_seg_subthres,double(Blue_di))

trace_transpose=Result.normTraces';
LocTrace_transpose=tr_local';

figure(2); clf;
rscale_trace=rescale2(trace_transpose,1)+[1:nROI];
plot(rscale_trace)
hold all
cmap2=jet(nSpike);
for sp=1:nSpike
sp_tmp=squeeze(spike_list(:,sp,:));
sp_tmp(:,isnan(sum(sp_tmp,1)))=[];
for n=1:nROI
    plot(sp_tmp(n,:),rscale_trace(sp_tmp(n,:),n),'.','color',cmap2(sp,:),'markersize',10)
end
end
yyaxis right
plot(bwBlue)
   
for b=1:max(bwBlue)
sp_list_tmp=spike_list(:,:,b)+size(Result.normTraces,2)*[0:nROI-1]';
ValidIndice=~isnan(sp_list_tmp);
AmpMat(:,1:sum(ValidIndice(:))/nROI,b)=reshape(trace_transpose(sp_list_tmp(ValidIndice)),nROI,[]);
LocAmpMat(:,1:sum(ValidIndice(:))/nROI,b)=reshape(LocTrace_transpose(sp_list_tmp(ValidIndice)),nROI,[]);
end

figure(4); clf; 
tiledlayout(3,3)
cmap=distinguishable_colors(nROI);
cmap_light=cmap*1.8; cmap_light(cmap_light>1)=1;

nexttile([1 1])
% for i=1:nROI
% plot(squeeze(DelayMat(i,:,:)),'color',cmap_light(i,:),'linewidth',0.5); hold all;
% end

l=errorbar(repmat([1:nSpike]',1,nROI),mean(DelayMat,3,'omitnan')',std(DelayMat,0,3,'omitnan')','linewidth',2,'marker','o');
%l=plot(mean(DelayMat,3,'omitnan')','linewidth',2);
arrayfun(@(l,c) set(l,'Color',c{:}),l',num2cell(cmap,2))
xlabel('N^t^h spike')   
ylabel('Delay (ms)')


nexttile([1 1])
% for i=1:nROI
% plot(squeeze(AmpMat(i,:,:)),'color',cmap_light(i,:),'linewidth',0.5); hold all;
% end
%
l=errorbar(repmat([1:nSpike]',1,nROI),mean(AmpMat./mean(AmpMat,2,'omitnan'),3,'omitnan')',std(AmpMat./mean(AmpMat,2,'omitnan'),0,3,'omitnan')','linewidth',2,'marker','o');
%l=plot(mean(AmpMat,3,'omitnan')','linewidth',2);
arrayfun(@(l,c) set(l,'Color',c{:}),l',num2cell(cmap,2))
xlabel('N^t^h spike')   
ylabel('Normalized amplitude')

nexttile([1 1])
% for i=1:nROI
% plot(squeeze(LocAmpMat(i,:,:)),'color',cmap_light(i,:),'linewidth',0.5); hold all;
% end

l=errorbar(repmat([1:nSpike]',1,nROI),mean(LocAmpMat./mean(LocAmpMat,2,'omitnan'),3,'omitnan')',std(LocAmpMat./mean(LocAmpMat,2,'omitnan'),0,3,'omitnan')','linewidth',2,'marker','o');
%l=plot(mean(LocAmpMat,3,'omitnan')','linewidth',2);
arrayfun(@(l,c) set(l,'Color',c{:}),l',num2cell(cmap,2))
xlabel('N^t^h spike')   
ylabel('Normalized local amplitude')

nexttile([2 3])
imagesc(rescale2(Result.normTraces,2))
colormap('turbo')

