clear
clc;
cd '/Volumes/BHL_WD18TB/20230831_PP72_Prism/BHLm078_Optopatch'
fpath{1} = '/Volumes/BHL_WD18TB/20230831_PP72_Prism/BHLm078_Optopatch/191502Prism_PP72_mouse2_488_dend_step_ramp_2';
fpath{2} = '/Volumes/BHL_WD18TB/20230831_PP72_Prism/BHLm078_Optopatch/191940Prism_PP72_mouse2_488_dend_step_30ms';

%% Motion correction

for i=1:2

load(fullfile(fpath{i},"output_data.mat"))
sz=double(Device_Data{1, 4}.ROI([2 4]));
ref_time=[9000:10000]; overlap=200;
time_segment=7000;

frm_end=max(Device_Data{1, 2}.Counter_Inputs(1, 1).data);
f_seg=[[1:time_segment:frm_end] frm_end+1]; f_seg(2:end)=f_seg(2:end)-1;

if length(ref_time)>2000
    mov_test=double(readBinMov_times([fpath{i} '/frames1.bin'],sz(2),sz(1),[ref_time(1):ref_time(1)+2000]));
else
    mov_test=double(readBinMov_times([fpath{i} '/frames1.bin'],sz(2),sz(1),[ref_time(1):ref_time(end)]));
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
        mov=double(readBinMov_times([fpath{i} '/frames1.bin'],sz(2),sz(1),[f_seg(j):f_seg(j+1)+overlap]));
    catch % when the image ends
        mov=double(readBinMov_times([fpath{i} '/frames1.bin'],sz(2),sz(1),[f_seg(j):f_seg(j+1)]));
    end

    mov=rollingShutter_correction(mov,Device_Data{1, 4}.exposuretime,'flash');
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

end

%%
for i=1:2
disp(fpath{i})
DAQ_rate=0.000005;
load([fpath{i} '/output_data.mat'])
sz=double(Device_Data{1, 4}.ROI([2 4]));
ref_time=[2000:3000];
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
Result{i}.ROIpoly=ROIpoly;

load(fullfile(fpath{i},['/mcTrace' num2str(1,'%02d') '.mat']));
frm_end=max(Device_Data{1, 2}.Counter_Inputs(1, 1).data);

mov_mc=double(readBinMov([fpath{i} '/mc_ShutterReg' num2str(1,'%02d') '.bin'],sz(2),sz(1)));
Result{i}.ref_im=mean(mov_mc,3);
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
Npoly=size(Result{i}.ROIpoly,1);
ftprnt = zeros(size(mov_filt,1)*size(mov_filt,2),Npoly);

for p=1:Npoly %each ROIs
    mask(:,:,p) = poly2mask(Result{i}.ROIpoly{p}(:,1), Result{i}.ROIpoly{p}(:,2), sz(2), sz(1));
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

Result{i}.ftprnt=toimg(ftprnt,sz(2),sz(1));
figure(4); clf;
imshow2(squeeze(sum(toimg(Result{i}.ftprnt,sz(2),sz(1)).*reshape(jet(Npoly),1,1,[],3),3)),[]);

end
%% Signal extraction

f_seg=[[1:time_segment:frm_end] frm_end+1]; f_seg(2:end)=f_seg(2:end)-1;
take_window=repmat([1 time_segment],length(f_seg)-1,1);
take_window(2:end,1)=take_window(2:end,1)+overlap; take_window(1:end-1,2)=take_window(1:end-1,2)+overlap;
take_window(end)=mod(f_seg(end),time_segment);


for i=1:2
Result{i}.traces=[];
Result{i}.traces_res=[];
Result{i}.mcTrace=[];
Result{i}.im_corr=[];
bound=5;
ref_im_vec=tovec(Result{i}.ref_im(bound:end-bound,bound:end-bound));

load(fullfile(fpath{i},"output_data.mat"))
sz=double(Device_Data{1, 4}.ROI([2 4]));
frm_end=max(Device_Data{1, 2}.Counter_Inputs(1, 1).data);
f_seg=[[1:time_segment:frm_end] frm_end+1]; f_seg(2:end)=f_seg(2:end)-1;
take_window=repmat([1 time_segment],length(f_seg)-1,1);
take_window(2:end,1)=take_window(2:end,1)+overlap; take_window(1:end-1,2)=take_window(1:end-1,2)+overlap;
take_window(end)=mod(f_seg(end),time_segment);
mc=[]; mov_mc=[];

CamCounter=Device_Data{1, 2}.Counter_Inputs(1, 1).data;
CamTrigger=find(CamCounter(2:end)-CamCounter(1:end-1));

Result{i}.Blue=Device_Data{1, 2}.buffered_tasks(1, 1).channels(1, 2).data;
Result{i}.Blue=Result{i}.Blue(CamTrigger);

for j=1:length(f_seg)-1
    j
    mov=double(readBinMov([fpath{i} '/mc_ShutterReg' num2str(j,'%02d') '.bin'],sz(2),sz(1)));
    load([fpath{i} '/mcTrace' num2str(j,'%02d') '.mat']);

        mov=mov(:,:,[take_window(j,1):take_window(j,2)]);
        mc=[mc; mcTrace.xymean([take_window(j,1):take_window(j,2)],:)];

        mov_mc(:,:,end+1:end+size(mov,3))=mov;
end
mov_mc=mov_mc(:,:,2:end);

    mov_mc_vec=tovec(mov_mc(bound:end-bound,bound:end-bound,:)); 
    mov_mc_vec=(mov_mc_vec-mean(mov_mc_vec,1))./std(mov_mc_vec,0,1);

    mov_res= mov_mc-mean(mov_mc,3);
    bkg = zeros(1, size(mov_mc,3));
%     bkg(1,:) = linspace(-1, 1, size(mov_mc,3));  % linear term
%     bkg(2,:) = linspace(-1, 1, size(mov_mc,3)).^2;  % quadratic term
    bkg(1,:)=movmedian(get_blueoffTrace(squeeze(mean(mov_mc,[1 2])),Result{i}.Blue,30),3000);
    mov_res = SeeResiduals(mov_res,mc);
    mov_res = SeeResiduals(mov_res,mc.^2);
    mov_res = SeeResiduals(mov_res,mc(:,1).*mc(:,3));
    mov_res= SeeResiduals(mov_res,bkg,1);

    Result{i}.traces=[Result{i}.traces -(tovec(mov_res)'*tovec(Result{i}.ftprnt))'];
    Result{i}.mcTrace=[Result{i}.mcTrace; mc];
    Result{i}.im_corr=[Result{i}.im_corr sum(mov_mc_vec.*ref_im_vec,1)/(size(mov_mc_vec,1)-1)];  %image correlation


Result{i}.traces=Result{i}.traces(:,1:length(CamTrigger)-1);
Result{i}.mcTrace=Result{i}.mcTrace(1:length(CamTrigger)-1,:);

end
save('Result_DendStim_20240107.mat','Result','fpath','-v7.3')
%% Clean up and norm

   % figure(2); clf;
   %  [p f]=fft_simple(Result{1}.traces,1000);
   %  ax1=nexttile([1 1]);
   %  plot(f(1,:),p')
   %  set(gca,'YScale','log')
   % 
   %  [p f]=fft_simple(Result{2}.traces,1000);
   %  ax2=nexttile([1 1]);
   %  plot(f(1,:),p')
   %  set(gca,'YScale','log')
   %  linkaxes([ax1 ax2],'x')


exclude_frq=[180 182]; %monitor
exclude_frq2=[173 175]; %monitor
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
 save('Result_DendStim_20240107.mat','Result','fpath','-v7.3')
%%
 load('Result_DendStim_20240107.mat')
%%

noi=[1]; nTau=[-30:60];
tr=Result{2}.normTrace(noi,:); t=[1:length(tr)];
spike=find_spike_bh(tr-movmedian(tr,300,2),4,2);
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
plot(rescale2(Pulse_STA,2)')
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

Result{2}.spike=find_spike_bh(Result{2}.normTrace-movmedian(Result{2}.normTrace,50,2),2.5,0.5);
Result{2}.spike(1,:) = remove_closeSpike(Result{2}.spike(1,:),Result{2}.normTrace(1,:),4);
Blue=Result{2}.Blue>0; delay_max=3; delay_front=2;

blueOff = Blue(1:end-1) == 0;
blueOff2 = imerode(blueOff, [ones(1,20), zeros(1, 20)]);
Blue_di=~blueOff2;

bwBlue=bwlabel(Blue_di);
nback=20; nfront=20; nSpike=4;
g=0;

nROI=size(Result{2}.normTrace,1);
spike_by_blue=Blue_di.*Result{2}.spike(1,:);
spike_list=NaN(nROI,nSpike,max(bwBlue));
spike_trace=zeros(size(Result{2}.spike,2),nROI);
DelayMat=NaN(nROI,nSpike,max(bwBlue));
AmpMat=NaN(nROI,nSpike,max(bwBlue));
LocAmpMat=NaN(nROI,nSpike,max(bwBlue));

for b=1:max(bwBlue)
t_tmp=find(bwBlue==b);
soma_spike=find(spike_by_blue(t_tmp));
trace_segment=Result{2}.normTrace(:,t_tmp);

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
spike_trace((sp_list_vec+size(Result{2}.spike,2)*[0:nROI-1]'))=1;
spike_trace=spike_trace';
%show_traces_spikes(Result{2}.normTrace,spike_trace,double(Blue_di))



tr_seg_subthres=zeros(nROI,size(Result{2}.normTrace,2));
for n=1:nROI
tr_seg_subthres(n,:)=get_subthreshold(Result{2}.normTrace(n,:),spike_trace(n,:),round(5+n*1.3),4);
end

tr_local=Result{2}.normTrace-tr_seg_subthres;
%show_traces_spikes_subth(Result{2}.normTrace,spike_trace,tr_seg_subthres,double(Blue_di))

trace_transpose=Result{2}.normTrace';
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
sp_list_tmp=spike_list(:,:,b)+size(Result{2}.normTrace,2)*[0:nROI-1]';
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


%% Putative dSpike

i=1; segment_dsp=[9100:9200]; segment_CS=[6600:6700];
load(fullfile(fpath{i},"output_data.mat"))
sz=double(Device_Data{1, 4}.ROI([2 4]));
frm_end=max(Device_Data{1, 2}.Counter_Inputs(1, 1).data);
f_seg=[[1:time_segment:frm_end] frm_end+1]; f_seg(2:end)=f_seg(2:end)-1;
take_window=repmat([1 time_segment],length(f_seg)-1,1);
take_window(2:end,1)=take_window(2:end,1)+overlap; take_window(1:end-1,2)=take_window(1:end-1,2)+overlap;
take_window(end)=mod(f_seg(end),time_segment);
mc=[]; mov_mc=[];

CamCounter=Device_Data{1, 2}.Counter_Inputs(1, 1).data;
CamTrigger=find(CamCounter(2:end)-CamCounter(1:end-1));

Result{i}.Blue=Device_Data{1, 2}.buffered_tasks(1, 1).channels(1, 2).data;
Result{i}.Blue=Result{i}.Blue(CamTrigger);

for j=1:length(f_seg)-1
    j
    mov=double(readBinMov([fpath{i} '/mc_ShutterReg' num2str(j,'%02d') '.bin'],sz(2),sz(1)));
    load([fpath{i} '/mcTrace' num2str(j,'%02d') '.mat']);

        mov=mov(:,:,[take_window(j,1):take_window(j,2)]);
        mc=[mc; mcTrace.xymean([take_window(j,1):take_window(j,2)],:)];

        mov_mc(:,:,end+1:end+size(mov,3))=mov;
end
mov_mc=mov_mc(:,:,2:end);

    mov_mc_vec=tovec(mov_mc(bound:end-bound,bound:end-bound,:)); 
    mov_mc_vec=(mov_mc_vec-mean(mov_mc_vec,1))./std(mov_mc_vec,0,1);

    mov_res= mov_mc-mean(mov_mc,3);
    bkg = zeros(1, size(mov_mc,3));
%     bkg(1,:) = linspace(-1, 1, size(mov_mc,3));  % linear term
%     bkg(2,:) = linspace(-1, 1, size(mov_mc,3)).^2;  % quadratic term
    bkg(1,:)=movmedian(get_blueoffTrace(squeeze(mean(mov_mc,[1 2])),Result{i}.Blue,30),3000);
    mov_res = SeeResiduals(mov_res,mc);
    mov_res = SeeResiduals(mov_res,mc.^2);
    mov_res = SeeResiduals(mov_res,mc(:,1).*mc(:,3));
    mov_res= SeeResiduals(mov_res,bkg,1);

% dSpike Case

    clear mov_res_filt_dsp
for z=1:length(segment_dsp)
mov_res_filt_dsp(:,:,z)=imgaussfilt(mov_res(:,:,segment_dsp(z)),2);
end

nFrames2=size(mov_res_filt_dsp,3);
datDS = imresize(mov_res_filt_dsp(10:end-10,10:end-10,:), 0.3, 'bilinear', 'Antialiasing',true);
datDS = datDS - medfilt1(datDS, 20, [], 3);

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
for i=1:10
plot(rescale(eigTraces(i,:)')+i-0.5)
hold all
end

nKeep = 10;
eigImgs = zeros(sz(2), sz(1), nKeep);
for j = 1:nKeep;
    eigImgs(:,:,j) = mean(mov_res_filt_dsp.*reshape(eigTraces(j,:), [1, 1, nFrames2]),3);
end;
figure; clf;
for j = 1:nKeep;
    nexttile([1 1]);
    imshow2(eigImgs(3:end-3,3:end-3,j), []);
    title(num2str(j))
end;

keep_ind=[1:3];
[ics, mixmat, sepmat] = sorted_ica(eigTraces(keep_ind,:)',length(keep_ind));
figure(9); clf
stackplot(ics)

eigImgsVec = tovec(eigImgs(:,:,keep_ind));
footPrintsVec = eigImgsVec*sepmat';
footPrints = toimg(footPrintsVec, [sz(2), sz(1)]);
figure(10); clf
for j = 1:length(keep_ind);
    nexttile([1 1]);
    imshow2(footPrints(3:end-3,3:end-3,j), []);
    title(num2str(j))
end;

dsp_img=mat2gray(mean(mov_res_filt_dsp(10:end-10,10:end-10,:).*reshape(ics(:,1),1,1,[]),3));

% CS case
clear mov_res_filt_CS
for z=1:length(segment_CS)
mov_res_filt_CS(:,:,z)=imgaussfilt(mov_res(:,:,segment_CS(z)),2);
end

nFrames2=size(mov_res_filt_CS,3);
datDS = imresize(mov_res_filt_CS(10:end-10,10:end-10,:), 0.3, 'bilinear', 'Antialiasing',true);
datDS = datDS - medfilt1(datDS, 20, [], 3);

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
for i=1:10
plot(rescale(eigTraces(i,:)')+i-0.5)
hold all
end

nKeep = 10;
eigImgs = zeros(sz(2), sz(1), nKeep);
for j = 1:nKeep;
    eigImgs(:,:,j) = mean(mov_res_filt_CS.*reshape(eigTraces(j,:), [1, 1, nFrames2]),3);
end;
figure; clf;
for j = 1:nKeep;
    nexttile([1 1]);
    imshow2(eigImgs(3:end-3,3:end-3,j), []);
    title(num2str(j))
end;

keep_ind=[1:4];
[ics, mixmat, sepmat] = sorted_ica(eigTraces(keep_ind,:)',length(keep_ind));
figure(9); clf
stackplot(ics)

eigImgsVec = tovec(eigImgs(:,:,keep_ind));
footPrintsVec = eigImgsVec*sepmat';
footPrints = toimg(footPrintsVec, [sz(2), sz(1)]);
figure(10); clf
for j = 1:length(keep_ind);
    nexttile([1 1]);
    imshow2(footPrints(3:end-3,3:end-3,j), []);
    title(num2str(j))
end;

keep_ind_ics=[1 2];

CS_img=1-mat2gray(squeeze(sum(mean(mov_res_filt_CS(10:end-10,10:end-10,:).*reshape(ics(:,keep_ind_ics),1,1,[],length(keep_ind_ics)),3),4)));

avgImg=mat2gray(medfilt2(Result{1}.ref_im,[3 3]));
avgImg=avgImg(10:end-10,10:end-10);
imshow2(avgImg,[])

figure; clf;

nexttile([1 1])
dsp_img_filt=medfilt2(dsp_img.*avgImg,[5 5]);
imshow2(dsp_img_filt,[]); colormap("turbo")
axis equal
nexttile([1 1])
CS_img_filt=medfilt2(CS_img.*avgImg,[5 5]);
imshow2(CS_img_filt,[]); colormap("turbo")
axis equal
nexttile([1 1])
imshow2(imfuse(dsp_img_filt,CS_img_filt),[])

nexttile([1 1])
[kymo_trace_dSp, kymoROI]=polyLineKymo2(-imgaussfilt3(mov_res_filt_dsp(10:end-10,10:end-10,:).*avgImg,[1.5 1.5 0.1]),50,35,avgImg);
[kymo_trace_CS]=polyLineKymo2(-imgaussfilt3(mov_res_filt_CS(10:end-10,10:end-10,:).*avgImg,[1.5 1.5 0.1]),50,35,avgImg);

imagesc([(kymo_trace_dSp)' (kymo_trace_CS)'])
colormap('turbo'); axis off tight

%%

