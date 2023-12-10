clear
clc;
cd '/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/20231123_BHLm112_Prism'
fpath = uigetfile_n_dir;
%% Motion correction

for i=13:length(fpath)

    load(fullfile(fpath{i},"output_data.mat"))
    sz=double(Device_Data{1, 3}.ROI([2 4]));
    ref_time=[9000:10000]; overlap=200;
    time_segment=16000;

    frm_end=max(Device_Data{1, 2}.Counter_Inputs(1, 1).data);
    f_seg=[[1:time_segment:frm_end] frm_end+1]; f_seg(2:end)=f_seg(2:end)-1;

    if length(ref_time)>2000
        mov_test=double(readBinMov_times([fpath{i} '/frames1.bin'],sz(2),sz(1),[ref_time(1):ref_time(1)+2000]));
    else
        mov_test=double(readBinMov_times([fpath{i} '/frames1.bin'],sz(2),sz(1),[ref_time(1):ref_time(end)]));
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
            mov=double(readBinMov_times([fpath{i} '/frames1.bin'],sz(2),sz(1),[f_seg(j):f_seg(j+1)+overlap]));
        catch % when the image ends
            mov=double(readBinMov_times([fpath{i} '/frames1.bin'],sz(2),sz(1),[f_seg(j):f_seg(j+1)]));
        end

        mov=rollingShutter_correction(mov,Device_Data{1, 3}.exposuretime,'fusion');
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

%% PCA analysis
i=13; bound=5;

load(fullfile(fpath{i},"output_data.mat"))
mov_mc=double(readBinMov([fpath{i} '/mc_ShutterReg' num2str(1,'%02d') '.bin'],sz(2),sz(1)));
load([fpath{i} '/mcTrace' num2str(1,'%02d') '.mat']);
mc=[mcTrace.xymean];

mov_mc_vec=tovec(mov_mc(bound:end-bound,bound:end-bound,:));
mov_mc_vec=(mov_mc_vec-mean(mov_mc_vec,1))./std(mov_mc_vec,0,1);

CamCounter=Device_Data{1, 2}.Counter_Inputs(1, 1).data;
CamTrigger=find(CamCounter(2:end)-CamCounter(1:end-1));
Blue=Device_Data{1, 2}.buffered_tasks(1, 1).channels(1, 2).data;
Blue=Blue(CamTrigger);

mov_res= mov_mc-mean(mov_mc,3);
bkg = zeros(1, size(mov_mc,3));
bkg(1,:)=movmedian(get_blueoffTrace(squeeze(mean(mov_mc,[1 2])),Blue,30),3000);
mov_res = SeeResiduals(mov_res,mc);
mov_res = SeeResiduals(mov_res,mc.^2);
mov_res = SeeResiduals(mov_res,mc(:,1).*mc(:,2));
mov_res= SeeResiduals(mov_res,bkg,1);

nFrames2=size(mov_res,3);
datDS = imresize(mov_res(bound:end-bound,bound:end-bound,:), 0.3, 'bilinear', 'Antialiasing',true);
% datDS = imfilter(dat2, fspecial('average', [3, 3]), 'replicate');
% datDS = datDS(3:end-3,22:289,:);
% datDS = datDS(1:3:end,1:3:end,:);
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
for i=1:20
    plot(rescale(eigTraces(i,:)')+i-0.5)
    hold all
end

nKeep = 20;
eigImgs = zeros(sz(2), sz(1), nKeep);
for j = 1:nKeep;
    eigImgs(:,:,j) = mean(mov_res.*reshape(eigTraces(j,:), [1, 1, nFrames2]),3);
end;
figure; clf;
for j = 1:nKeep;
    nexttile([1 1]);
    imshow2(eigImgs(3:end-3,3:end-3,j), []);
    title(num2str(j))
end;


keep_ind=[1:20];
[ics, mixmat, sepmat] = sorted_ica(eigTraces(keep_ind,:)',length(keep_ind));
figure(9); clf
stackplot(ics)

eigImgsVec = tovec(eigImgs(:,:,keep_ind));
footPrintsVec = eigImgsVec*sepmat';
footPrints = toimg(footPrintsVec, [sz(2), sz(1)]);
figure(10); clf
for j = 1:size(footPrints,3);
    nexttile([1 1]);
    imshow2(footPrints(3:end-3,3:end-3,j), []);
    title(num2str(j))
end;

keep_ind_ics=[1 3];
ReconMov=toimg(tovec(mov_res).*mean(footPrintsVec(:,keep_ind_ics),2),sz(2),sz(1));

%% ROI extract
for i=[13]
    disp(fpath{i})
    DAQ_rate=0.000005;
    load([fpath{i} '/output_data.mat'])
    sz=double(Device_Data{1, 3}.ROI([2 4]));
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

    frm_rate=double((CamTrigger(2)-CamTrigger(1))*DAQ_rate);
    mov_mc=double(readBinMov([fpath{i} '/mc_ShutterReg' num2str(1,'%02d') '.bin'],sz(2),sz(1)));
    Result{i}.ref_im=mean(mov_mc,3);
    %[clickyROI, original_trace]=clicky(mov_mc);

    mov_res= mov_mc-mean(mov_mc,3);
    try
        mov_res = SeeResiduals(mov_res,mcTrace.xymean(:,1).*mcTrace.xymean(:,3));
    end

    %[SeeRes_trace]=apply_clicky(clickyROI,mov_res);
    % figure(3); clf;
    % plot(rescale(original_trace)); hold all;
    % plot(rescale(SeeRes_trace)+1);

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
%%
for i=2:12
    Result{i}.ftprnt=Result{1}.ftprnt;
    Result{i}.ROIpoly=Result{1}.ROIpoly;
    Result{i}.ref_im=Result{1}.ref_im;
end

for i=14:15
    Result{i}.ftprnt=Result{13}.ftprnt;
    Result{i}.ROIpoly=Result{13}.ROIpoly;
    Result{i}.ref_im=Result{13}.ref_im;
end

save('Result_BHLm112_20231126.mat','Result','fpath','-v7.3')
%%
for i=1:15
    i
    Result{i}.traces=[];
    Result{i}.traces_res=[];
    Result{i}.mcTrace=[];
    Result{i}.im_corr=[];
    bound=5;
    ref_im_vec=tovec(Result{i}.ref_im(bound:end-bound,bound:end-bound));

    load(fullfile(fpath{i},"output_data.mat"))
    sz=double(Device_Data{1, 3}.ROI([2 4]));
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
    bkg(1,:)=movmedian(get_blueoffTrace(squeeze(mean(mov_mc,[1 2])),Result{i}.Blue,30),3000,'omitnan');
    mov_res = SeeResiduals(mov_res,mc);
    mov_res = SeeResiduals(mov_res,mc.^2);
    try
        mov_res = SeeResiduals(mov_res,mc(:,1).*mc(:,3));
    end
    mov_res= SeeResiduals(mov_res,bkg,1);

    Result{i}.traces=[Result{i}.traces -(tovec(mov_res)'*tovec(Result{i}.ftprnt))'];
    Result{i}.mcTrace=[Result{i}.mcTrace; mc];
    Result{i}.im_corr=[Result{i}.im_corr sum(mov_mc_vec.*ref_im_vec,1)/(size(mov_mc_vec,1)-1)];  %image correlation


    Result{i}.traces=Result{i}.traces(:,1:length(CamTrigger));
    Result{i}.mcTrace=Result{i}.mcTrace(1:length(CamTrigger),:);

    nROI=size(Result{i}.ftprnt,3);
    Result{i}.traces_res=squeeze(SeeResiduals(reshape(Result{i}.traces,1,nROI,[]),Result{i}.mcTrace));
    Result{i}.traces_res=squeeze(SeeResiduals(reshape(Result{i}.traces_res,1,nROI,[]),Result{i}.mcTrace.^2'));
    Result{i}.traces_res=squeeze(SeeResiduals(reshape(Result{i}.traces_res,1,nROI,[]),bkg'));

end

save('Result_BHLm112_20231126.mat','Result','fpath','-v7.3')

%%

for i=1:length(Result)
    Result{i}.normTrace=Result{i}.traces_res./get_threshold(Result{i}.traces_res,1);
    Result{i}.spike=find_spike_bh(Result{i}.normTrace-movmedian(Result{i}.normTrace,300,2),5,3);
end
%
save('Result_BHLm112_20231126.mat','Result','fpath','-v7.3')

%% show traces
i=6; cmap=distinguishable_colors(3);


tr=rescale2(Result{i}.normTrace,2);
show_traces_spikes(tr,Result{i}.spike,Result{i}.Blue)
figure(2); clf;
imagesc(tr); hold all
colormap('turbo')
l=plot(-tr'+[1:3]+0.5,'linewidth',1);
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(cmap,2))
plot(-rescale(Result{i}.Blue)/2+4.1,'color',[0 0.6 1],'LineWidth',1)
axis tight off


%% Show traces
i=12; cmap=distinguishable_colors(3);
DMDPatt=mod(cumsum(Result{i}.DMDtrigger),2)+1;
tr=rescale2(Result{i}.normTrace,2);
show_traces_spikes(tr,Result{i}.spike,[Result{i}.Blue; DMDPatt])
figure(2); clf;
imagesc(tr); hold all
colormap('turbo')
l=plot(-tr'+[1:3]+0.5,'linewidth',1);
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(cmap,2))
plot(-rescale(Result{i}.Blue)/2+4.1,'color',[0 0.6 1],'LineWidth',1)
plot(-DMDPatt/2+4.6,'r')
axis tight off


%% make a STA
i=8;
noi=[2]; nTau=[-30:60]; nROI=size(Result{i}.normTrace,1);
tr=Result{i}.normTrace(noi,:); t=[1:length(tr)];
spike=find_spike_bh(tr-movmedian(tr,300,2),4,1.5);
tr=rescale(tr);
Blue=rescale(Result{i}.Blue);
blueOff = Blue == 0;
blueOff2 = imerode(blueOff, [ones(1,20), zeros(1, 20)]);
Blue_di=~blueOff2;

f1=figure(2); clf;
tiledlayout(1,3)
nexttile([1 2])
plot(tr); hold all
plot(find(spike),tr(find(spike)),'r.')
plot(Blue_di)
bwBlue_di=bwlabel(Blue_di);
sp_pulse=[];
for b=1:max(bwBlue_di)
    t_tmp=find(bwBlue_di==b);
    if sum(spike(t_tmp))>0
        sp_pulse=[sp_pulse find(spike(t_tmp),1,'first')+t_tmp(1)-1];
    end
end
Pulse_Tau=sp_pulse(1:end-1)'+nTau;
Pulse_spikeMat=reshape(Result{i}.normTrace(:,Pulse_Tau),nROI,length(sp_pulse)-1,[]);
Pulse_STA=squeeze(mean(Pulse_spikeMat,2));
Pulse_STA=Pulse_STA-prctile(Pulse_STA, 5, 2);
plot(sp_pulse,tr(sp_pulse),'marker','o','color',[1 0 1],'LineStyle','none')
nexttile([1 1])
plot(rescale2(Pulse_STA,2)')
legend({'Soma','Trunk','Distal'})
saveas(f1,fullfile(fpath{i} ,['STA_trace.fig']))

load(fullfile(fpath{i},"output_data.mat"))
sz=double(Device_Data{1, 3}.ROI([2 4])); mov_mc=[]; mc=[];
mov=double(readBinMov([fpath{i} '/mc_ShutterReg' num2str(1,'%02d') '.bin'],sz(2),sz(1)));
load([fpath{i} '/mcTrace' num2str(1,'%02d') '.mat']);
mc=[mc; mcTrace.xymean];
mov_mc(:,:,end+1:end+size(mov,3))=mov;
mov_mc=mov_mc(:,:,2:end);
mov_res= mov_mc-mean(mov_mc,3);
bkg = zeros(1, size(mov_mc,3));
bkg(1,:)=movmedian(get_blueoffTrace(squeeze(mean(mov_mc,[1 2])),Result{i}.Blue,30),3000,'omitnan');
mov_res = SeeResiduals(mov_res,mc);
mov_res = SeeResiduals(mov_res,mc.^2);
try
    mov_res = SeeResiduals(mov_res,mc(:,1).*mc(:,3));
end
mov_res= SeeResiduals(mov_res,bkg,1);

TriggerMov=double(reshape(mov_res(:,:,Pulse_Tau),sz(2),sz(1),[],length(nTau)));
STAMov=-squeeze(mean(TriggerMov,3));
ref_im_filt=medfilt2(Result{i}.ref_im-prctile(Result{i}.ref_im,5),[5 5]);
ref_im_filt(ref_im_filt<0)=0;
STAMov_filt=STAMov.*mat2gray(ref_im_filt); STAMov_filt=STAMov_filt-mean(STAMov_filt(:,:,1:15),3);
for j=1:size(STAMov_filt,3)
    STAMov_filt(:,:,j)=medfilt2(STAMov_filt(:,:,j),[5 5]);
end
writeMov_wTrace([fpath{i} ,'/STAmov'],STAMov_filt,[],10,1,[1 15],[],rescale2(Pulse_STA,2)')
