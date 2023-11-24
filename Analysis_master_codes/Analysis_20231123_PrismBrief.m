clear
clc;
cd '/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/20231118_BHLm85_87_88_89_107_109_Prism'
fpath = uigetfile_n_dir;
%%
i=1;
disp(fpath{i})
DAQ_rate=0.000005;
load([fpath{i} '/output_data.mat'])
sz=double(Device_Data{1, 3}.ROI([2 4]));
mov=double(readBinMov([fpath{i} '/frames1.bin'],sz(2),sz(1)));
mov=mov(:,:,1:end-1); avgImg=mean(mov,3);

CamCounter=Device_Data{1, 2}.Counter_Inputs(1, 1).data;
CamTrigger=find(CamCounter(2:end)-CamCounter(1:end-1));
Blue=Device_Data{1, 2}.buffered_tasks(1, 1).channels(1, 2).data;
Blue=Blue(CamTrigger);

bkg(1,:)=movmedian(get_blueoffTrace(squeeze(mean(mov,[1 2])),Blue,30),3000);
mov=mov-avgImg;
mov_res= SeeResiduals(mov,bkg,1);

[~, intens]=clicky(mov_res,avgImg);

%%
noi=[1]; nTau=[-30:60];
tr=-intens'; t=[1:length(tr)];
spike=find_spike_bh(tr-movmedian(tr,300,2),5,4);
blueOff = Blue == 0;
blueOff2 = imerode(blueOff, [ones(1,20), zeros(1, 20)]);
Blue_di=~blueOff2;

f1=figure(2); clf;
tiledlayout(1,3)
nexttile([1 2])
plot(tr); hold all
plot(find(spike),tr(find(spike)),'r.')
bwBlue_di=bwlabel(Blue_di);
sp_pulse=[];
for b=1:max(bwBlue_di)
t_tmp=find(bwBlue_di==b);
if sum(spike(t_tmp))>0
sp_pulse=[sp_pulse find(spike(t_tmp),1,'first')+t_tmp(1)-1];
end
end
BlueBi=Blue>0;
%sp_pulse=find((BlueBi(2:end)-BlueBi(1:end-1))==1);
Pulse_Tau=sp_pulse(1:end-1)'+nTau;
Pulse_spikeMat=reshape(tr(Pulse_Tau),1,length(sp_pulse)-1,[]);
Pulse_STA=squeeze(mean(Pulse_spikeMat,2));
Pulse_STA=Pulse_STA-prctile(Pulse_STA, 5, 1);
nexttile([1 1])
plot(rescale(Pulse_STA))
set(f1, 'Position', [100, 100, 800, 400]);
saveas(f1,fullfile(fpath{1} ,['BriefSTA' '.fig']))
print(f1, fullfile(fpath{1} ,['BriefSTA' '.jpg']),'-djpeg', ['-r', num2str(400)]);

close all
TriggerMov=double(reshape(mov_res(:,:,Pulse_Tau),sz(2),sz(1),[],length(nTau)));
STAMov=-squeeze(mean(TriggerMov,3));
ref_im_filt=medfilt2(avgImg-prctile(avgImg,5),[5 5]);
ref_im_filt(ref_im_filt<0)=0;
STAMov_filt=STAMov.*mat2gray(ref_im_filt); STAMov_filt=STAMov_filt-mean(STAMov_filt(:,:,1:15),3);
for i=1:size(STAMov_filt,3)
STAMov_filt(:,:,i)=medfilt2(STAMov_filt(:,:,i),[5 5]);
end
writeMov_wTrace([fpath{1} ,'/STAmov'],mat2gray(STAMov_filt),[],10,1,[0.1 0.7],[],rescale(Pulse_STA))

%% Clean up and norm

   figure(2); clf;
    [p f]=fft_simple(Result{1}.traces,1000);
    ax1=nexttile([1 1]);
    plot(f(1,:),p')
    set(gca,'YScale','log')

    [p f]=fft_simple(Result{2}.traces,1000);
    ax2=nexttile([1 1]);
    plot(f(1,:),p')
    set(gca,'YScale','log')
    linkaxes([ax1 ax2],'x')


exclude_frq=[182 184]; %monitor
exclude_frq2=[186 188]; %monitor
time_bin=1000; Fs=1000; %2nd trunk is the reliable trace

    nROI=size(Result{1}.traces,1);

freq_lowhigh=exclude_frq/(Fs/2);
[b, a] = butter(4, freq_lowhigh, 'stop');

 freq_lowhigh2=exclude_frq2/(Fs/2);
 [b2, a2] = butter(4, freq_lowhigh2, 'stop');

 for n=1:nROI
Result{1}.normTrace(n,:)=filtfilt(b,a,Result{1}.traces(n,:));
Result{2}.normTrace(n,:)=filtfilt(b2,a2,Result{2}.traces(n,:));
 end

 for i=1:2
    Result{i}.normTrace=Result{i}.normTrace./get_threshold(Result{i}.normTrace,1);
    Result{i}.spike=find_spike_bh(Result{i}.normTrace-movmedian(Result{i}.normTrace,300,2),5,3);
 end
% 
 save('Result_SomaStim_20231114.mat','Result','fpath','-v7.3')
%%
 load('Result_SomaStim_20231114.mat')
%%

noi=[1]; nTau=[-30:60];
tr=Result{2}.normTrace(noi,:); t=[1:length(tr)];
spike=find_spike_bh(tr-movmedian(tr,300,2),4,1.5);
Blue=Result{2}.Blue;
blueOff = Blue == 0;
blueOff2 = imerode(blueOff, [ones(1,20), zeros(1, 20)]);
Blue_di=~blueOff2;

figure(2); clf;
plot(tr); hold all
plot(find(spike),tr(find(spike)),'r.')
bwBlue_di=bwlabel(Blue_di);
sp_pulse=[];
for b=1:max(bwBlue_di)
t_tmp=find(bwBlue_di==b);
if sum(spike(t_tmp))>2
sp_pulse=[sp_pulse find(spike(t_tmp),1,'first')+t_tmp(1)-1];
end
end
Pulse_Tau=sp_pulse(1:end-1)'+nTau;
Pulse_spikeMat=reshape(Result{2}.normTrace(:,Pulse_Tau),4,length(sp_pulse)-1,[]);
Pulse_STA=squeeze(mean(Pulse_spikeMat,2));
Pulse_STA=Pulse_STA-prctile(Pulse_STA, 5, 2);

%% STA movie Short Pulse
TriggerMov=double(reshape(mov_res(:,:,Pulse_Tau),sz(2),sz(1),[],length(nTau)));
STAMov=-squeeze(mean(TriggerMov,3));
ref_im_filt=medfilt2(Result{2}.ref_im-prctile(Result{2}.ref_im,5),[5 5]);
ref_im_filt(ref_im_filt<0)=0;
STAMov_filt=STAMov.*mat2gray(ref_im_filt); STAMov_filt=STAMov_filt-mean(STAMov_filt(:,:,1:15),3);
for i=1:size(STAMov_filt,3)
STAMov_filt(:,:,i)=medfilt2(STAMov_filt(:,:,i),[5 5]);
end
writeMov_wTrace([fpath{2} ,'/STAmov'],STAMov_filt,[],10,1,[-0.1 2.5],[],rescale2(Pulse_STA,2)')
%% Spike height Analysis

Result{1}.spike=find_spike_bh(Result{1}.normTrace-movmedian(Result{1}.normTrace,50,2),2.5,0.5);
Result{1}.spike(1,:) = remove_closeSpike(Result{1}.spike(1,:),Result{1}.normTrace(1,:),4);
Blue=Result{1}.Blue>0; delay_max=3; delay_front=2;

blueOff = Blue == 0;
blueOff2 = imerode(blueOff, [ones(1,20), zeros(1, 20)]);
Blue_di=~blueOff2;

bwBlue=bwlabel(Blue_di);
nback=20; nfront=20; nSpike=5;
g=0;

nROI=size(Result{1}.normTrace,1);
spike_by_blue=Blue_di.*Result{1}.spike(1,:);
spike_list=NaN(nROI,nSpike,max(bwBlue));
spike_trace=zeros(size(Result{1}.spike,2),nROI);
DelayMat=NaN(nROI,nSpike,max(bwBlue));
AmpMat=NaN(nROI,nSpike,max(bwBlue));
LocAmpMat=NaN(nROI,nSpike,max(bwBlue));

for b=1:max(bwBlue)
t_tmp=find(bwBlue==b);
soma_spike=find(spike_by_blue(t_tmp));
trace_segment=Result{1}.normTrace(:,t_tmp);

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
spike_trace((sp_list_vec+size(Result{1}.spike,2)*[0:nROI-1]'))=1;
spike_trace=spike_trace';
%show_traces_spikes(Result.normTrace,spike_trace,double(Blue_di))



tr_seg_subthres=zeros(nROI,size(Result{1}.normTrace,2));
for n=1:nROI
tr_seg_subthres(n,:)=get_subthreshold(Result{1}.normTrace(n,:),spike_trace(n,:),round(5+n*1.3),4);
end

tr_local=Result{1}.normTrace-tr_seg_subthres;
%show_traces_spikes_subth(Result.normTrace,spike_trace,tr_seg_subthres,double(Blue_di))

trace_transpose=Result{1}.normTrace';
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
sp_list_tmp=spike_list(:,:,b)+size(Result{1}.normTrace,2)*[0:nROI-1]';
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
imagesc(rescale2(Result{2}.normTrace,2))
colormap('turbo')


%% Spike height Analysis

Result{1}.spike=find_spike_bh(Result{1}.normTrace-movmedian(Result{1}.normTrace,100,2),2.5,0.5);
Result{1}.spike(1,:) = remove_closeSpike(Result{1}.spike(1,:),Result{1}.normTrace(1,:),4);
Blue=Result{1}.Blue; boi=[22]; delay_max=3; boi2=boi-18;
bwBlue=bwlabel(Blue); N_spike=10;
nback=10; nfront=50;
f1=figure(2); clf;
g=0;
for b=boi
        t_tmp=find(bwBlue==b);
        t_tmp=[t_tmp(1)-nback:t_tmp(end)+nfront];
nROI=size(Result{1}.normTrace,1);
tr_seg=Result{1}.normTrace(:,t_tmp);
som_spike=find(Result{1}.spike(1,t_tmp));
Dend_tr=Result{1}.normTrace(2:end,t_tmp);


dend_delay=[];
Hp_tr=Dend_tr - tr_seg_subthres(2:end,:);
for s=1:length(som_spike)
[~, dend_delay(s,:)]=max(Dend_tr(:,som_spike(s)+[-1:delay_max]),[],2);
%[~, dend_delay(s,:)]=max(Hp_tr(:,som_spike(s)+[-1:delay_max]),[],2);
end

dend_spike=dend_delay' + som_spike - 2;
sp_list=[som_spike; dend_spike];
sp_tr=zeros(length(t_tmp),nROI);
sp_tr(sp_list+length(t_tmp)*[0:nROI-1]')=1;

tr_seg_subthres=zeros(length(t_tmp),nROI)';
for n=1:nROI
tr_seg_subthres(n,:)=get_subthreshold(tr_seg(n,:),sp_tr(:,n)',round(7+n*1.5),5);
%[upperEnv, tr_seg_subthres(n,:)] = envelope(tr_seg(n,:),round(6+n*1.3),'peak');
end


spike_information{g+1}.speed=dend_delay-2;
tr_seg=tr_seg';
tr_seg_local=tr_seg-tr_seg_subthres';
spike_information{g+1}.amplitude=tr_seg(sp_list+size(tr_seg,1)*[0:nROI-1]');
spike_information{g+1}.local_amp=tr_seg_local(sp_list+size(tr_seg,1)*[0:nROI-1]');


show_traces_spikes_subth(tr_seg',sp_tr',tr_seg_subthres,double(Blue(t_tmp)))
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
% imagesc(rescale2(Result{1}.normTrace(:,t_tmp),2))
% colormap('turbo')
% axis tight off
% 
% nexttile(12*g+[3],[3 2])
% l=plot(rescale2(Result{1}.normTrace(:,t_tmp),2)'-[1:4]);
% arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(distinguishable_colors(4),2))
% axis tight
g=g+1;
end

set(f1, 'Position', [100, 100, 600, 400]);
saveas(f1,fullfile(fpath{1} ,[num2str(boi2) 'th_trace.fig']))
print(f1, fullfile(fpath{1} ,[num2str(boi2) 'th_trace.jpg']),'-djpeg', ['-r', num2str(300)]);

f2=figure(3); clf;
tiledlayout(1,5)
    delays=[];
amplitudes=[]; amplitudes_ratio=[];
loc_amps=[]; loc_amps_ratio=[];
g=1; 


cmap_delay=distinguishable_colors(4)+0.3; cmap_delay(cmap_delay>1)=1;
cmap_delay=cmap_delay(2:end,:);
cmap2=distinguishable_colors(4); cmap2=cmap2(2:end,:);
cmap=distinguishable_colors(4)+0.3; cmap(cmap>1)=1;
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
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(distinguishable_colors(4),2))
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
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(distinguishable_colors(4),2))
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
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(distinguishable_colors(4),2))
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
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(distinguishable_colors(4),2))
xlabel('N^t^h spike')       
ylabel('Norm. local amplitude')  

set(f2, 'Position', [100, 100, 1300, 300]);
saveas(f2,fullfile(fpath{1},[num2str(boi2) 'th_Pulse.fig']))
print(f2, fullfile(fpath{1},[num2str(boi2) 'th_Pulse.jpg']),'-djpeg', ['-r', num2str(300)]);