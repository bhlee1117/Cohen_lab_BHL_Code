clear; clc;
f1='/Volumes/cohen_lab/Lab/Labmembers/Pojeong Park/Data/230717_carbachol_without_picrotoxin/Cell1_carbachol/174347PP070_P15_spont_carbachol_50uM';
load(fullfile(f1,"output_data.mat"))
sz=double(Device_Data{1, 4}.ROI([2 4]));
mov_mc=double(readBinMov([f1 '/frames1.bin'],sz(2),sz(1)));
mov_mc=mov_mc(:,:,50:13580);
Camtrigger=find([0 Device_Data{1, 2}.Counter_Inputs.data(2:end)-Device_Data{1, 2}.Counter_Inputs.data(1:end-1)]);
%Blue=Blue(Camtrigger);
%Blue=Blue(50:13580);
[roi, tr_soma] = clicky(mov_mc);
%%
%[bkg blue_off]=get_blueoffTrace(squeeze(mean(mov_mc,[1 2])),[Blue],30);
bkg=movmedian(squeeze(mean(mov_mc,[1 2])),50,'omitnan');
% [y_fit t_consts coeffY]  = expfitDM_2(find(blue_off)',squeeze(mean(mov_mc(:,:,blue_off),[1 2])), ...
%     [1:size(mov_mc,3)]',2000);
t_fit=setdiff([1:size(mov_mc,3)],[9600:11500]);
[y_fit t_consts coeffY]  = expfitDM_2(t_fit',bkg(t_fit), ...
    [1:size(mov_mc,3)]',[10000 1000 100]);
figure(2); clf;
plot(bkg); hold all
plot(y_fit);
mov_res=SeeResiduals(mov_mc,y_fit);
[tr_raw]=apply_clicky(roi,-mov_res);
%%
tr=tr_raw'./get_threshold(tr_raw',1);
tr_soma_hi=(tr-movmedian(tr,50));
sp=find_spike_bh(tr_soma_hi./get_threshold(tr_soma_hi,1),5,2);

tr_sub=mean(tr,1)-movprc(mean(tr,1),1500,20,2);
tr_sub=get_subthreshold(tr_sub,sp,5,10);

[trans tr_trace]=detect_transient2(tr_sub,[7 1.5],sp,15);
transcand=cell2mat(cellfun(@(x) length(x)>2,trans.ISI,'UniformOutput',false));
meanISI_frnt=cellfun(@(x) mean(x(1:2)),trans.ISI(transcand));
meanISI_first3=zeros(1,length(trans.length));
meanISI_first3(transcand)=meanISI_frnt;

%CS_ind=find(trans.spike_number>2 & trans.mean_ISI<15);
CS_ind=find(trans.spike_number>2 & meanISI_first3<30);
CS_trace=ismember(tr_trace,CS_ind);
CS_spike=sp.*bwlabel(CS_trace);
[~, CS_spike_time]=unique(CS_spike);

SS_s=find(sp.*(1-CS_trace));
CS_s=find(sp.*CS_trace,1,'first');
nTau={[-20:20],[-100:1000]};

STA_SSmovie=reshape(mov_res(:,:,SS_s'+nTau{1}),sz(2),sz(1),[],length(nTau{1}));
STA_CSmovie=reshape(mov_res(:,:,CS_s'+nTau{2}),sz(2),sz(1),[],length(nTau{2}));

[SS_kymo,ROI]=polyLineKymo3(-squeeze(mean(STA_SSmovie,3)),30,30);
[CS_kymo]=apply_clicky(ROI,-squeeze(mean(STA_CSmovie,3)));
SS_kymo=SS_kymo-mean(SS_kymo(1:5,:),1); 
CS_kymo=CS_kymo-mean(CS_kymo(1:5,:),1); 
F_ref=mean(SS_kymo(-nTau{1}(1)+[7:13],:))';

%%
figure(4); clf; cax=[-3 15];
tiledlayout(3,2)
nexttile([1 1])
imagesc(SS_kymo'./F_ref,cax)
nexttile([1 1])
imagesc(CS_kymo'./F_ref,cax)

nexttile([2 1])
Z=imresize(SS_kymo./F_ref',3);
X=repmat(imresize(nTau{1},3),size(Z,2)/3,1)';
Y=repmat([1:size(Z,2)],size(Z,1),1);
C=grs2rgb(Z,jet,cax(1),cax(2));
l=surf(X,Y,Z,C); colormap(turbo);
l.FaceAlpha=0.8; l.EdgeColor='none';
xlabel('Peri-spike time (ms)')
ylabel('ROIs')
zlabel('\DeltaF/F_r_e_f')
nexttile([2 1])
Z=imresize(CS_kymo./F_ref',3);
X=repmat(imresize(nTau{2},3),size(Z,2)/3,1)';
Y=repmat([1:size(Z,2)],size(Z,1),1);
C=grs2rgb(Z,jet,cax(1),cax(2));
l=surf(X,Y,Z,C); colormap(turbo);
l.FaceAlpha=0.8; l.EdgeColor='none';
xlabel('Peri-spike time (ms)')
ylabel('ROIs')
zlabel('\DeltaF/F_r_e_f')

%%

clear; clc;
f1='/Volumes/cohen_lab/Lab/Labmembers/Pojeong Park/Data/230206_PP071_P15_electrode_PTX(50uM)_and_Mg2+(1mM)/Cell2/194824PP071_P15_Cell2_soma_30ms+EC_30V_interval_-30ms';
f2='/Volumes/cohen_lab/Lab/Labmembers/Pojeong Park/Data/230206_PP071_P15_electrode_PTX(50uM)_and_Mg2+(1mM)/Cell2/193627PP071_P15_Cell2_soma_30ms+EC_30V_interval_0ms';
f3='/Volumes/cohen_lab/Lab/Labmembers/Pojeong Park/Data/230206_PP071_P15_electrode_PTX(50uM)_and_Mg2+(1mM)/Cell2/200224PP071_P15_Cell2_soma_30ms+EC_30V_interval_-20 ms';

load(fullfile(f1,"output_data.mat"))
sz=double(Device_Data{1, 4}.ROI([2 4]));
mov_mc=double(readBinMov([f1 '/frames1.bin'],sz(2),sz(1)));
mov_mc=mov_mc(:,:,50:3870);
DAQ_rate=1/Device_Data{1, 2}.Counter_Inputs.rate;
Camtrigger=find([0 Device_Data{1, 2}.Counter_Inputs.data(2:end)-Device_Data{1, 2}.Counter_Inputs.data(1:end-1)]);
Blue=Device_Data{1, 2}.buffered_tasks(1, 1).channels(1, 2).data;
Blue=Blue(Camtrigger);
Blue=Blue(50:3870);
EFS=Device_Data{1, 2}.buffered_tasks(1, 1).channels(1, 3).data;
EFS=EFS(Camtrigger);
EFS=EFS(50:3870);
Camtrigger=Camtrigger(50:3870);
[tr_soma, roi] = polyLineKymo3(mov_mc,30,40);

load(fullfile(f2,"output_data.mat"))
mov_mc2=double(readBinMov([f2 '/frames1.bin'],sz(2),sz(1)));
mov_mc2=mov_mc2(:,:,50:3870);
DAQ_rate=1/Device_Data{1, 2}.Counter_Inputs.rate;
Camtrigger2=find([0 Device_Data{1, 2}.Counter_Inputs.data(2:end)-Device_Data{1, 2}.Counter_Inputs.data(1:end-1)]);
Blue2=Device_Data{1, 2}.buffered_tasks(1, 1).channels(1, 2).data;
Blue2=Blue2(Camtrigger2);
Blue2=Blue2(50:3870);
EFS2=Device_Data{1, 2}.buffered_tasks(1, 1).channels(1, 3).data;
EFS2=EFS2(Camtrigger2);
EFS2=EFS2(50:3870);
Camtrigger2=Camtrigger2(50:3870);

load(fullfile(f3,"output_data.mat"))
mov_mc3=double(readBinMov([f3 '/frames1.bin'],sz(2),sz(1)));
mov_mc3=mov_mc3(:,:,50:3870);
DAQ_rate=1/Device_Data{1, 2}.Counter_Inputs.rate;
Camtrigger3=find([0 Device_Data{1, 2}.Counter_Inputs.data(2:end)-Device_Data{1, 2}.Counter_Inputs.data(1:end-1)]);
Blue3=Device_Data{1, 2}.buffered_tasks(1, 1).channels(1, 2).data;
Blue3=Blue3(Camtrigger3);
Blue3=Blue3(50:3870);
EFS3=Device_Data{1, 2}.buffered_tasks(1, 1).channels(1, 3).data;
EFS3=EFS3(Camtrigger3);
EFS3=EFS3(50:3870);
Camtrigger3=Camtrigger3(50:3870);
%%
%[bkg blue_off]=get_blueoffTrace(squeeze(mean(mov_mc,[1 2])),[Blue],30);
bkg=movmedian(squeeze(mean(mov_mc,[1 2])),50,'omitnan');
bkg2=movmedian(squeeze(mean(mov_mc2,[1 2])),50,'omitnan');
bkg3=movmedian(squeeze(mean(mov_mc3,[1 2])),50,'omitnan');
% [y_fit t_consts coeffY]  = expfitDM_2(find(blue_off)',squeeze(mean(mov_mc(:,:,blue_off),[1 2])), ...
%     [1:size(mov_mc,3)]',2000);
t_fit=setdiff([1:size(mov_mc,3)],[1500:2400]);
[y_fit t_consts coeffY]  = expfitDM_2(t_fit',bkg(t_fit), ...
    [1:size(mov_mc,3)]',[1000 100]);
[y_fit2 t_consts coeffY]  = expfitDM_2(t_fit',bkg2(t_fit), ...
    [1:size(mov_mc,3)]',[1000 100]);
[y_fit3 t_consts coeffY]  = expfitDM_2(t_fit',bkg3(t_fit), ...
    [1:size(mov_mc,3)]',[1000 100]);
figure(2); clf;
nexttile([1 1])
plot(bkg); hold all
plot(y_fit);
nexttile([1 1])
plot(bkg2); hold all
plot(y_fit2);
nexttile([1 1])
plot(bkg3); hold all
plot(y_fit3);

mov_res=SeeResiduals(mov_mc,y_fit);
mov_res2=SeeResiduals(mov_mc2,y_fit2);
mov_res3=SeeResiduals(mov_mc3,y_fit3);
[tr_raw]=apply_clicky(roi,-mov_res);
[tr_raw2]=apply_clicky(roi,-mov_res2);
[tr_raw3]=apply_clicky(roi,-mov_res3);
%%
tr=tr_raw'./get_threshold(tr_raw',1);
tr_soma_hi=(tr(1,:)-movmedian(tr(1,:),50));
sp=find_spike_bh(tr_soma_hi./get_threshold(tr_soma_hi,1),5,2);
F_ref=mean(tr(:,923:930),2);
tr=tr./F_ref;

tr2=tr_raw2'./get_threshold(tr_raw2',1);
tr2=tr2./F_ref;

tr3=tr_raw3'./get_threshold(tr_raw3',1);
tr3=tr3./F_ref;
%%
figure(4); clf; cax=[-0.5 4];
tiledlayout(6,1); ax1=[];

ax1=[ax1 nexttile([1 1])];
imagesc(Camtrigger*DAQ_rate,[1:length(roi)],tr,cax)
colormap("turbo")
ax1=[ax1 nexttile([1 1])];
plot(Camtrigger*DAQ_rate,Blue,'color',[0 0.5 1]); hold all
plot(Camtrigger*DAQ_rate,EFS,'color',[0 0 0])
xlabel('Time (s)')
legend({'Mod488','EFS'})
title(f1,'Interpreter','none')

ax1=[ax1 nexttile([1 1])];
imagesc(Camtrigger2*DAQ_rate,[1:length(roi)],tr2,cax)
colormap("turbo")
ax1=[ax1 nexttile([1 1])];
plot(Camtrigger2*DAQ_rate,Blue2,'color',[0 0.5 1]); hold all
plot(Camtrigger2*DAQ_rate,EFS2,'color',[0 0 0])
legend({'Mod488','EFS'})
title(f2,'Interpreter','none')

ax1=[ax1 nexttile([1 1])];
imagesc(Camtrigger3*DAQ_rate,[1:length(roi)],tr3,cax)
colormap("turbo")
ax1=[ax1 nexttile([1 1])];
plot(Camtrigger3*DAQ_rate,Blue3,'color',[0 0.5 1]); hold all
plot(Camtrigger3*DAQ_rate,EFS3,'color',[0 0 0])
linkaxes([ax1],'x')
legend({'Mod488','EFS'})
title(f3 ,'Interpreter','none')

%%
tr_soma_hi=(tr(1,:)-movmedian(tr(1,:),50));
sp=find_spike_bh(tr_soma_hi./get_threshold(tr_soma_hi,1),5,2);
sp_list=find(sp);
%SS_tmplt=tr(:,sp_list(1:3)+[-3:5]);
SS_tmplt=squeeze(mean(reshape(tr(:,sp_list(1:3)'+[-3:5]),size(tr,1),[],9),2));

SPconvMat=conv2(SS_tmplt,sp);
SPconvMat=SPconvMat(:,4:end-5);

tr_soma_hi2=(tr2(1,:)-movmedian(tr2(1,:),50));
sp2=find_spike_bh(tr_soma_hi2./get_threshold(tr_soma_hi2,1),5,2);
sp_list=find(sp2);
SS_tmplt2=squeeze(mean(reshape(tr2(:,sp_list(1:3)'+[-3:5]),size(tr2,1),[],9),2));
SPconvMat2=conv2(SS_tmplt2,sp2);
SPconvMat2=SPconvMat2(:,4:end-5);

tr_soma_hi3=(tr3(1,:)-movmedian(tr3(1,:),50));
sp3=find_spike_bh(tr_soma_hi3./get_threshold(tr_soma_hi3,1),5,2);
sp_list=find(sp3);
SS_tmplt3=squeeze(mean(reshape(tr3(:,sp_list(1:3)'+[-3:5]),size(tr3,1),[],9),2));
SPconvMat3=conv2(SS_tmplt3,sp3);
SPconvMat3=SPconvMat3(:,4:end-5);

figure(5); clf; cax=[-0.5 4];
tiledlayout(6,1); ax1=[];
ax1=[ax1 nexttile([1 1])];
imagesc(Camtrigger*DAQ_rate,[1:length(roi)],tr-SPconvMat,cax)
colormap("turbo")
ax1=[ax1 nexttile([1 1])];
plot(Camtrigger*DAQ_rate,Blue,'color',[0 0.5 1]); hold all
plot(Camtrigger*DAQ_rate,EFS,'color',[0 0 0])
plot(Camtrigger*DAQ_rate,sp,'color',[1 0 0])
xlabel('Time (s)')
legend({'Mod488','EFS','Somatic spike'})
title(f1,'Interpreter','none')
linkaxes(ax1,'x')

ax1=[ax1 nexttile([1 1])];
imagesc(Camtrigger*DAQ_rate,[1:length(roi)],tr2-SPconvMat2,cax)
colormap("turbo")
ax1=[ax1 nexttile([1 1])];
plot(Camtrigger*DAQ_rate,Blue2,'color',[0 0.5 1]); hold all
plot(Camtrigger*DAQ_rate,EFS2,'color',[0 0 0])
plot(Camtrigger*DAQ_rate,sp2,'color',[1 0 0])
xlabel('Time (s)')
legend({'Mod488','EFS','Somatic spike'})
title(f1,'Interpreter','none')
linkaxes(ax1,'x')

ax1=[ax1 nexttile([1 1])];
imagesc(Camtrigger*DAQ_rate,[1:length(roi)],tr3-SPconvMat3,cax)
colormap("turbo")
ax1=[ax1 nexttile([1 1])];
plot(Camtrigger*DAQ_rate,Blue3,'color',[0 0.5 1]); hold all
plot(Camtrigger*DAQ_rate,EFS3,'color',[0 0 0])
plot(Camtrigger*DAQ_rate,sp3,'color',[1 0 0])
xlabel('Time (s)')
legend({'Mod488','EFS','Somatic spike'})
title(f1,'Interpreter','none')
linkaxes(ax1,'x')
%%
clear; clc;
f1='/Volumes/cohen_lab/Lab/Labmembers/Pojeong Park/Data/230126_PP070_P17_electrode_PTX(50uM)_and_Mg2+(1mM)/Cell2/225407PP070_P17_Cell2_soma_30ms+EC_20V';

load(fullfile(f1,"output_data.mat"))
sz=double(Device_Data{1, 4}.ROI([2 4]));
mov_mc=double(readBinMov([f1 '/frames1.bin'],sz(2),sz(1)));
mov_mc=mov_mc(:,:,50:3000);
DAQ_rate=1/Device_Data{1, 2}.Counter_Inputs.rate;
Camtrigger=find([0 Device_Data{1, 2}.Counter_Inputs.data(2:end)-Device_Data{1, 2}.Counter_Inputs.data(1:end-1)]);
Blue=Device_Data{1, 2}.buffered_tasks(1, 1).channels(1, 2).data;
Blue=Blue(Camtrigger);
Blue=Blue(50:3000);
EFS=Device_Data{1, 2}.buffered_tasks(1, 2).channels(1, 2).data;
EFS=movsum(EFS,Camtrigger(2)-Camtrigger(1));
EFS=EFS(Camtrigger);
EFS=EFS(50:3000);
Camtrigger=Camtrigger(50:3000);
[tr_soma, roi] = polyLineKymo3(mov_mc,30,40);

%%
%[bkg blue_off]=get_blueoffTrace(squeeze(mean(mov_mc,[1 2])),[Blue],30);
bkg=movmedian(squeeze(mean(mov_mc,[1 2])),50,'omitnan');
% [y_fit t_consts coeffY]  = expfitDM_2(find(blue_off)',squeeze(mean(mov_mc(:,:,blue_off),[1 2])), ...
%     [1:size(mov_mc,3)]',2000);
t_fit=setdiff([1:size(mov_mc,3)],[1400:1600]);
[y_fit t_consts coeffY]  = expfitDM_2(t_fit',bkg(t_fit), ...
    [1:size(mov_mc,3)]',[1000 100]);
figure(2); clf;
nexttile([1 1])
plot(bkg); hold all
plot(y_fit);

mov_res=SeeResiduals(mov_mc,y_fit);
[tr_raw]=apply_clicky(roi,-mov_res);
%%
tr=tr_raw'./get_threshold(tr_raw',1);
tr_soma_hi=(tr-movmedian(tr,50));
sp=find_spike_bh(tr_soma_hi./get_threshold(tr_soma_hi,1),5,2);
F_ref=mean(tr(:,1186:1190),2);
tr=tr./F_ref;
%%
figure(4); clf; cax=[-0.5 4];
tiledlayout(2,1)
ax1=nexttile([1 1]);
imagesc(Camtrigger*DAQ_rate,[1:length(roi)],tr,cax)
colormap("turbo")
ax2=nexttile([1 1]);
plot(Camtrigger*DAQ_rate,Blue,'color',[0 0.5 1]); hold all
plot(Camtrigger*DAQ_rate,EFS,'color',[0 0 0])
linkaxes([ax1,ax2],'x')
xlabel('Time (s)')
legend({'Mod488','EFS'})
title(f1,'Interpreter','none')

%%
clear; clc;
f1='/Volumes/cohen_lab/Lab/Labmembers/Pojeong Park/Data/230202_PP070_P17_electrode_PTX(50uM)_and_Mg2+(1mM)/Cell2/202348PP070_P16_Cell2_soma_30ms+EC_30V';

load(fullfile(f1,"output_data.mat"))
sz=double(Device_Data{1, 4}.ROI([2 4]));
mov_mc=double(readBinMov([f1 '/frames1.bin'],sz(2),sz(1)));
mov_mc=mov_mc(:,:,50:3000);
DAQ_rate=1/Device_Data{1, 2}.Counter_Inputs.rate;
Camtrigger=find([0 Device_Data{1, 2}.Counter_Inputs.data(2:end)-Device_Data{1, 2}.Counter_Inputs.data(1:end-1)]);
Blue=Device_Data{1, 2}.buffered_tasks(1, 1).channels(1, 2).data;
Blue=Blue(Camtrigger);
Blue=Blue(50:3000);
EFS=Device_Data{1, 2}.buffered_tasks(1, 1).channels(1, 3).data;
EFS=EFS(Camtrigger);
EFS=EFS(50:3000);
Camtrigger=Camtrigger(50:3000);
[tr_soma, roi] = polyLineKymo3(mov_mc,30,40);

%%
%[bkg blue_off]=get_blueoffTrace(squeeze(mean(mov_mc,[1 2])),[Blue],30);
bkg=movmedian(squeeze(mean(mov_mc,[1 2])),50,'omitnan');
% [y_fit t_consts coeffY]  = expfitDM_2(find(blue_off)',squeeze(mean(mov_mc(:,:,blue_off),[1 2])), ...
%     [1:size(mov_mc,3)]',2000);
t_fit=setdiff([1:size(mov_mc,3)],[400:650 1200:1550]);
[y_fit t_consts coeffY]  = expfitDM_2(t_fit',bkg(t_fit), ...
    [1:size(mov_mc,3)]',[1000 100]);
figure(2); clf;
nexttile([1 1])
plot(bkg); hold all
plot(y_fit);

mov_res=SeeResiduals(mov_mc,y_fit);
[tr_raw]=apply_clicky(roi,-mov_res);
%%
tr=tr_raw'./get_threshold(tr_raw',1);
tr_soma_hi=(tr-movmedian(tr,50));
sp=find_spike_bh(tr_soma_hi./get_threshold(tr_soma_hi,1),5,2);
F_ref=mean(tr(:,1263:1270),2);
tr=tr./F_ref;
%%
figure(4); clf; cax=[-0.5 4];
tiledlayout(2,1)
ax1=nexttile([1 1]);
imagesc(Camtrigger*DAQ_rate,[1:length(roi)],tr,cax)
colormap("turbo")
ax2=nexttile([1 1]);
plot(Camtrigger*DAQ_rate,Blue,'color',[0 0.5 1]); hold all
plot(Camtrigger*DAQ_rate,EFS,'color',[0 0 0])
linkaxes([ax1,ax2],'x')
xlabel('Time (s)')
legend({'Mod488','EFS'})
title(f1,'Interpreter','none')
