clear
clc;
cd '/Volumes/BHL_WD18TB/20230831_PP72_Prism/BHLm078_Optopatch/190937Prism_PP72_mouse2_488_soma_step_ramp_2';
fpath = '/Volumes/BHL_WD18TB/20230831_PP72_Prism/BHLm078_Optopatch/190937Prism_PP72_mouse2_488_soma_step_ramp_2';
%addpath(genpath('/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Code/2020RigControl'))
%% Motion correction

load(fullfile(fpath,"output_data.mat"))
sz=double(Device_Data{1, 4}.ROI([2 4]));
ref_time=[9000:10000]; overlap=200;
time_segment=6000;

frm_end=max(Device_Data{1, 2}.Counter_Inputs(1, 1).data);
f_seg=[[1:time_segment:frm_end] frm_end+1]; f_seg(2:end)=f_seg(2:end)-1;

if length(ref_time)>2000
    mov_test=double(readBinMov_times([fpath '/frames1.bin'],sz(2),sz(1),[ref_time(1):ref_time(1)+2000]));
else
    mov_test=double(readBinMov_times([fpath '/frames1.bin'],sz(2),sz(1),[ref_time(1):ref_time(end)]));
end
mov_test=rollingShutter_correction(mov_test,Device_Data{1, 4}.exposuretime,'flash');
mov_test=mov_test(:,:,2:end);
[mov_test_mc,xyField]=optical_flow_motion_correction_LBH(mov_test,mean(mov_test,3),'normcorre');
mov_test=vm(mov_test);
mov_test = single(mov_test)./single(max(mov_test.data(:)));
mov_test = movmean(mov_test,10,3);
mov_ref = squeeze(median(mov_test,3));

for j=1:length(f_seg)-1
    try
        mov=double(readBinMov_times([fpath '/frames1.bin'],sz(2),sz(1),[f_seg(j):f_seg(j+1)+overlap]));
    catch % when the image ends
        mov=double(readBinMov_times([fpath '/frames1.bin'],sz(2),sz(1),[f_seg(j):f_seg(j+1)]));
    end

    mov=rollingShutter_correction(mov,Device_Data{1, 4}.exposuretime,'flash');
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
sz=double(Device_Data{1, 4}.ROI([2 4]));
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

load(fullfile(fpath,['/mcTrace' num2str(1,'%02d') '.mat']));
frm_end=max(Device_Data{1, 2}.Counter_Inputs(1, 1).data);
Result.Blue=Device_Data{1, 2}.buffered_tasks(1, 1).channels(1, 2).data;
CamCounter=Device_Data{1, 2}.Counter_Inputs(1, 1).data;
CamTrigger=find(CamCounter(2:end)-CamCounter(1:end-1));

Result.Blue=Result.Blue(CamTrigger);
Result.Blue=Result.Blue(2:end);
frm_rate=double((CamTrigger(2)-CamTrigger(1))*DAQ_rate);
mov_mc=double(readBinMov([fpath '/mc_ShutterReg' num2str(1,'%02d') '.bin'],sz(2),sz(1)));
Result.ref_im=mean(mov_mc,3);
% [clickyROI, original_trace]=clicky(mov_mc);

mov_res= mov_mc-mean(mov_mc,3);
mcTrace.xymean=movmean(mcTrace.xymean,3,2);
mov_res = SeeResiduals(mov_res,mcTrace.xymean);
mov_res = SeeResiduals(mov_res,mcTrace.xymean.^2);
mov_res = SeeResiduals(mov_res,mcTrace.xymean(:,1).*mcTrace.xymean(:,3));

% [SeeRes_trace]=apply_clicky(clickyROI,mov_res);
% figure(3); clf;
% plot(rescale(original_trace)); hold all;
% plot(rescale(SeeRes_trace)+1);

n_comp=5;
mov_filt=imgaussfilt3(mov_res(:,:,3000:5000),[3 3 0.1]);
movVec=tovec(mov_filt);
Npoly=size(Result.ROIpoly,1);
ftprnt = zeros(size(mov_filt,1)*size(mov_filt,2),Npoly);

for p=1:Npoly %each ROIs
    p
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
f1=figure(4); clf;
imshow2(squeeze(sum(toimg(Result.ftprnt,sz(2),sz(1)).*reshape(jet(Npoly),1,1,[],3),3)),[]);
saveas(f1,fullfile(fpath ,['footprint.fig']))
print(f1, fullfile(fpath ,['footprint.jpg']),'-djpeg', ['-r', num2str(600)]);
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
%     bkg(1,:) = linspace(-1, 1, size(mov_mc,3));  % linear term
%     bkg(2,:) = linspace(-1, 1, size(mov_mc,3)).^2;  % quadratic term
    mov_res = SeeResiduals(mov_res,mc);
    mov_res = SeeResiduals(mov_res,mc.^2);
    mov_res = SeeResiduals(mov_res,mc(:,1).*mc(:,3));
%     mov_res= SeeResiduals(mov_res,bkg,1);

    Result.traces=[Result.traces -(tovec(mov_res)'*tovec(Result.ftprnt))'];
    Result.mcTrace=[Result.mcTrace; mc];
    Result.im_corr=[Result.im_corr sum(mov_mc_vec.*ref_im_vec,1)/(size(mov_mc_vec,1)-1)];  %image correlation

end

Result.traces=Result.traces(:,1:length(CamTrigger)-1);
Result.mcTrace=Result.mcTrace(1:length(CamTrigger)-1,:);

save(fullfile(fpath,'Result_Optopatch_20231029.mat'),'Result','fpath','-v7.3')

%% Clean up and norm
exclude_frq=[182 184]; %monitor
%exclude_frq2=[483.5 484]; %monitor
    time_bin=1000; Fs=1000; ref_trace=2; %2nd trunk is the reliable trace

    nTime=12200;
    nROI=size(Result.traces,1);
freq_lowhigh=exclude_frq/(Fs/2);
[b, a] = butter(4, freq_lowhigh, 'stop');

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
    bkg=movmedian(get_blueoffTrace(tr,Result.Blue(1:nTime),30),1000)';
    tr_res= tr-bkg;

    % regress out motion frequency
    traces_res_filtered(n,:) = filtfilt(b, a, tr_res);
    norm_trace(n,:)=traces_res_filtered(n,:);%-movmedian(traces_res_filtered(n,:),500,2);

end
      
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

%% Spike height Analysis

Result.spike=find_spike_bh(Result.normTraces-movmedian(Result.normTraces,100,2),3.5,0);

Blue=Result.Blue; boi=[23]; delay_max=4; boi2=boi-18;
bwBlue=bwlabel(Blue);
nback=100; nfront=100;
figure(2); clf;
g=0;
for b=boi
        t_tmp=find(bwBlue==b);
        t_tmp=[t_tmp(1)-nback:t_tmp(end)+nfront];
nROI=size(Result.normTraces,1);
tr_seg=Result.normTraces(:,t_tmp);
som_spike=find(Result.spike(1,t_tmp));
Dend_tr=Result.normTraces(2:end,t_tmp);


tr_seg_subthres=zeros(length(t_tmp),nROI)';
for n=1:nROI
%tr_seg_subthres(n,:)=get_subthreshold(tr_seg(n,:),sp_tr(:,n)',round(7+n*1.3),4);
[upperEnv, tr_seg_subthres(n,:)] = envelope(tr_seg(n,:),round(6+n*1.3),'peak');
end


dend_delay=[];
Hp_tr=Dend_tr - tr_seg_subthres(2:end,:);
for s=1:length(som_spike)
%[~, dend_delay(s,:)]=max(Dend_tr(:,som_spike(s)+[0:delay_max]),[],2);
[~, dend_delay(s,:)]=max(Hp_tr(:,som_spike(s)+[-1:delay_max]),[],2);
end

dend_spike=dend_delay' + som_spike - 2;
sp_list=[som_spike; dend_spike];
sp_tr=zeros(length(t_tmp),nROI);
sp_tr(sp_list+length(t_tmp)*[0:nROI-1]')=1;



spike_information{g+1}.speed=dend_delay-2;
tr_seg=tr_seg';
tr_seg_local=tr_seg-tr_seg_subthres';
spike_information{g+1}.amplitude=tr_seg(sp_list+size(tr_seg,1)*[0:nROI-1]');
spike_information{g+1}.local_amp=tr_seg_local(sp_list+size(tr_seg,1)*[0:nROI-1]');


show_traces_spikes_subth(tr_seg',sp_tr',tr_seg_subthres,tr_seg_subthres(1,:))
% 
% 
% plot(tr_seg+[1:4]*10)
% hold all
% plot(tr_seg_subthres'+[1:4]*10,'k')
% 
% 
% nexttile([1 2])
% plot(Blue(t_tmp))
% ylim([0 1.6])
% axis tight off
% 
% nexttile(12*g+5,[2 2])
% imagesc(rescale2(Result.normTraces(:,t_tmp),2))
% colormap('turbo')
% axis tight off
% 
% nexttile(12*g+[3],[3 2])
% l=plot(rescale2(Result.normTraces(:,t_tmp),2)'-[1:4]);
% arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(distinguishable_colors(4),2))
% axis tight
g=g+1;
end

figure(3); clf;

    delays=[];
amplitudes=[]; amplitudes_ratio=[];
loc_amps=[]; loc_amps_ratio=[];
g=1; N_spike=9;


cmap_delay=turbo(4)+0.3; cmap_delay(cmap_delay>1)=1;
cmap_delay=cmap_delay(2:end,:);
cmap2=turbo(4); cmap2=cmap2(2:end,:);
cmap=turbo(4)+0.3; cmap(cmap>1)=1;
nexttile([1 1])
g=1;
for b=1:length(spike_information)
    delays(:,:,g)=spike_information{b}.speed(1:N_spike,:);
    l=plot([1:N_spike],spike_information{b}.speed(1:N_spike,:));  hold all
    arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(cmap_delay,2))
    g=g+1;
end
l=plot([1:N_spike],mean(delays,3),'linewidth',2);
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(cmap2,2))
xlabel('N^t^h spike')
ylabel('Delay (ms)')


nexttile([1 1])
g=1;
for b=1:length(spike_information)
    amplitudes(:,:,g)=spike_information{b}.amplitude(:,1:N_spike);
    l=plot([1:N_spike],spike_information{b}.amplitude(:,1:N_spike));  hold all
    arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(cmap,2))
    g=g+1;
end
l=plot([1:N_spike],mean(amplitudes,3)','linewidth',2);
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(turbo(4),2))
xlabel('N^t^h spike')
ylabel('Amplitude')

nexttile([1 1])
g=1;
for b=1:length(spike_information)
    loc_amps(:,:,g)=spike_information{b}.local_amp(:,1:N_spike);
    l=plot([1:N_spike],spike_information{b}.local_amp(:,1:N_spike));  hold all
    arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(cmap,2))
    g=g+1;
end
l=plot([1:N_spike],mean(loc_amps,3)','linewidth',2);
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(turbo(4),2))
xlabel('N^t^h spike')
ylabel('Local amplitude')

nexttile([1 1])
g=1;
for b=1:length(spike_information)
    amplitudes_ratio(:,:,g)=spike_information{b}.amplitude(:,1:N_spike)./mean(spike_information{b}.amplitude(:,N_spike-5:N_spike),2);
    l=plot([1:N_spike],amplitudes_ratio(:,:,g));  hold all
    arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(cmap,2))
    g=g+1;
end
l=plot([1:N_spike],mean(amplitudes_ratio,3)','linewidth',2);
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(turbo(4),2))
xlabel('N^t^h spike')
ylabel('Norm. Amplitude')

nexttile([1 1])
g=1;
for b=1:length(spike_information)
    loc_amps_ratio(:,:,g)=spike_information{b}.local_amp(:,1:N_spike)./mean(spike_information{b}.local_amp(:,N_spike-5:N_spike),2);
    l=plot([1:N_spike],loc_amps_ratio(:,:,g));  hold all
    arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(cmap,2))
    g=g+1;
end
l=plot([1:N_spike],mean(loc_amps_ratio,3)','linewidth',2);
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(turbo(4),2))
xlabel('N^t^h spike')       
ylabel('Norm. local amplitude')     


%%
noi=[2]; nTau=[-30:60];
tr=Result.normTraces(noi,:); t=[1:length(tr)];
spike=find_spike_bh(tr-movmedian(tr,300,2),3,2);
Blue=Result.Blue;
se = strel('square', [5]); % 0 degree means horizontal
Blue_di = imdilate(Blue, se);

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

