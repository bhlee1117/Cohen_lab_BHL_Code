
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
j=84; 
load(fullfile(fpath{j},"output_data.mat"))
load([fpath{j} '/mcTrace' num2str(1,'%02d') '.mat']);
CamCounter=Device_Data{1, 2}.Counter_Inputs(1, 1).data;
CamTrigger=find(CamCounter(2:end)-CamCounter(1:end-1));
Blue=Device_Data{1, 2}.buffered_tasks(1, 1).channels(1, 2).data;
Blue=Blue(CamTrigger);
mc=mcTrace.xymean;
sz=double(Device_Data{1, 3}.ROI([2 4]));  
mov_mc=double(readBinMov([fpath{j} '/mc_ShutterReg' num2str(1,'%02d') '.bin'],sz(2),sz(1)));
bkg(1,:)=movmedian(get_blueoffTrace(squeeze(mean(mov_mc,[1 2])),[Blue],100),3000,'omitnan');
bkg(2,:) = linspace(-1, 1, size(mov_mc,3));  % linear term
bkg(3,:) = linspace(-1, 1, size(mov_mc,3)).^2;  % quadratic term
mov_res= mov_mc-median(mov_mc,3);
mov_res = SeeResiduals(mov_res,mc);
mov_res = SeeResiduals(mov_res,mc.^2);
mov_res = SeeResiduals(mov_res,mc(:,1).*mc(:,end));
mov_res= SeeResiduals(mov_res,bkg,1);

Rfixed = imref2d(repmat(Device_Data{1, 3}.virtualSensorSize,1,2));
inverseTform = invert(Device_Data{1, 6}.tform);
revertedImage = imwarp(double(Device_Data{1, 6}.Target), inverseTform,'OutputView',Rfixed);
[blueDMDimg blueDMDcontour]=imcrop(revertedImage,double(Device_Data{1, 3}.ROI([1 3 2 4]))+[0 0 -1 -1]);
figure; imagesc(im_merge(cat(3,mean(mov_mc,3),blueDMDimg),[1 1 1;0 0.5 1]))
axis equal tight off
%%
j=86; 
load(fullfile(fpath{j},"output_data.mat"))
load([fpath{j} '/mcTrace' num2str(1,'%02d') '.mat']);
CamCounter=Device_Data{1, 2}.Counter_Inputs(1, 1).data;
CamTrigger=find(CamCounter(2:end)-CamCounter(1:end-1));
Blue=Device_Data{1, 2}.buffered_tasks(1, 1).channels(1, 2).data;
Blue=Blue(CamTrigger);
mc=mcTrace.xymean;
sz=double(Device_Data{1, 3}.ROI([2 4]));  
mov_mc=double(readBinMov([fpath{j} '/mc_ShutterReg' num2str(1,'%02d') '.bin'],sz(2),sz(1)));
bkg(1,:)=movmedian(get_blueoffTrace(squeeze(mean(mov_mc,[1 2])),[Blue],100),3000,'omitnan');
bkg(2,:) = linspace(-1, 1, size(mov_mc,3));  % linear term
bkg(3,:) = linspace(-1, 1, size(mov_mc,3)).^2;  % quadratic term
mov_res_DD= mov_mc-mean(mov_mc,3);
mov_res_DD = SeeResiduals(mov_res_DD,mc);
mov_res_DD = SeeResiduals(mov_res_DD,mc.^2);
mov_res_DD = SeeResiduals(mov_res_DD,mc(:,1).*mc(:,end));
mov_res_DD= SeeResiduals(mov_res_DD,bkg,1);

Rfixed = imref2d(repmat(Device_Data{1, 3}.virtualSensorSize,1,2));
inverseTform = invert(Device_Data{1, 6}.tform);
revertedImage = imwarp(double(Device_Data{1, 6}.Target), inverseTform,'OutputView',Rfixed);
[blueDMDimg blueDMDcontour]=imcrop(revertedImage,double(Device_Data{1, 3}.ROI([1 3 2 4]))+[0 0 -1 -1]);
figure; imagesc(im_merge(cat(3,mean(mov_mc,3),blueDMDimg),[1 1 1;0 0.5 1]))
axis equal tight off
%%
[roi, tr_som]=clicky(mov_res);
[tr_DD]=apply_clicky(roi,mov_res_DD);
%%
tr_som=tr_som-median(tr_som(end-100:end,:));
tr_DD=tr_DD-median(tr_DD(end-100:end,:));
figure(4); clf;
tiledlayout(length(roi),1)
ax1=[];
for r=1:length(roi)
    ax1=[ax1 nexttile([1 1])];
    plot(-tr_som(:,r));
    hold all
    plot(-tr_DD(:,r)+80)
    title(['ROI#', num2str(r)])
end
linkaxes(ax1,'x')

%%
[kymoSom, kymoROI]=polyLineKymo2(mov_res,15,15);
[kymoDD]=apply_polyLineKymo2(mov_res_DD,15,15,mean(mov_res_DD,3),kymoROI);
[F0]=apply_polyLineKymo2(mean(mov_mc,3),15,15,mean(mov_res_DD,3),kymoROI);
%%

figure(4); clf; ax2=[];
tiledlayout(3,1)
ax2=[ax2 nexttile([1 1])];
imagesc(-kymoSom'./F0')
colormap(turbo)
ax2=[ax2 nexttile([1 1])];
imagesc(-kymoDD'./F0')
colormap(turbo)
ax2=[ax2 nexttile([1 1])];
plot(Blue)
linkaxes(ax2,'x')

%%
figure(6); clf; ax2=[];
tiledlayout(3,1)
rois=[1:3:size(kymoSom,2)];
ax2=[ax2 nexttile([1 1])];
l=plot(rescale2(-kymoSom(:,rois)./F0(rois),1)+[1:length(rois)]*0.6);
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(turbo(length(rois)),2))
ax2=[ax2 nexttile([1 1])];
l=plot(rescale2(-kymoDD(:,rois)./F0(rois),1)+[1:length(rois)]*0.6);
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(turbo(length(rois)),2))
ax2=[ax2 nexttile([1 1])];
plot(Blue)
linkaxes(ax2,'xy')

%%
  som_tr=-tr_som(:,1)'; som_tr=som_tr./get_threshold(som_tr,1);
  som_tr_dd=-tr_som(:,end)';
  dd_tr=-tr_DD(:,1)'; dd_tr=dd_tr./get_threshold(dd_tr,1);
  dd_tr_dd=-tr_DD(:,end)';
  sp_soma=find_spike_bh(som_tr+movmedian(som_tr,50,2),5,2);
  sp_dd=find_spike_bh(dd_tr+movmedian(dd_tr,50,2),5,2);
 
  tr_sub_soma=som_tr-movprc(som_tr,200,20,2);
  tr_sub_soma=get_subthreshold(tr_sub_soma,sp_soma,5,10);

    tr_sub_DD=dd_tr-movprc(dd_tr,200,20,2);
  tr_sub_DD=get_subthreshold(tr_sub_DD,sp_dd,5,10);

    [trans tr_trace]=detect_transient2(tr_sub_soma,[4 1.5],sp_soma,15);
    CS_ind=find(trans.spike_number>1 & trans.mean_ISI<20);
    CS_trace=ismember(tr_trace,CS_ind);
    CS_spike=sp_soma.*bwlabel(CS_trace);

        [transDD tr_traceDD]=detect_transient2(tr_sub_DD,[4 1.5],sp_dd,15);
    CS_indDD=find(transDD.spike_number>1 & transDD.mean_ISI<20);
    CS_traceDD=ismember(tr_traceDD,CS_indDD);
    CS_spikeDD=sp_dd.*bwlabel(CS_traceDD);

    %%
    Area_som_dd=[]; Area_DD_dd=[];
    Amp_som_dd=[];  Amp_DD_dd=[];
    for b=1:max(bwlabel(CS_trace))
        bind=find(bwlabel(CS_trace)==b);
        Area_som_dd(b)=sum(som_tr_dd(bind));
        Area_som(b)=sum(som_tr(bind));
        Amp_som_dd(b)=max(som_tr_dd(bind));
    end

    for b=1:max(bwlabel(CS_traceDD))
        bind=find(bwlabel(CS_traceDD)==b);
        Area_DD_dd(b)=sum(dd_tr_dd(bind));
        Area_DD(b)=sum(dd_tr(bind));
        Amp_DD_dd(b)=max(som_tr_dd(bind));
    end

figure(7); clf;
nexttile([1 1])
plot_errorbar2([1 2],{trans.mean_ISI(CS_ind)',transDD.mean_ISI(CS_indDD)'},'ttest2',20,'Mean ISI (ms)',{'S','D'},1)
nexttile([1 1])
%plot_errorbar2([1 2],{trans.int(CS_ind)',transDD.int(CS_indDD)'},'ttest2',[],'Area (soma)',{'S','D'},1)
plot_errorbar2([1 2],{Area_som',Area_DD'},'ttest2',[1000],'Area (soma)',{'S','D'},1)
nexttile([1 1])
plot_errorbar2([1 2],{Area_som_dd',Area_DD_dd'},'ttest2',[],'Area (distal dendrite)',{'S','D'},1)
nexttile([1 1])
plot_errorbar2([1 2],{trans.spike_number(CS_ind)',transDD.spike_number(CS_indDD)'},'ttest2',[],'# of spike',{'S','D'},1)
nexttile([1 1])
plot_errorbar2([1 2],{Amp_som_dd',Amp_DD_dd'},'ttest2',[],'Amplitude (distal dendrite)',{'S','D'},1)
nexttile([1 1])
plot_errorbar2([1 2],{trans.amp(CS_ind)',transDD.amp(CS_indDD)'},'ttest2',[20],'Amplitude (soma)',{'S','D'},1)
nexttile([1 1])
plot_errorbar2([1 2],{trans.length(CS_ind)',transDD.length(CS_indDD)'},'ttest2',[],'Length (frame)',{'S','D'},1)

%%
nTau=[-200:150];
 for b=1:max(bwlabel(CS_trace))
        C_start(b)=find(bwlabel(CS_trace)==b,1);        
    end
for b=1:max(bwlabel(CS_traceDD))
        C_startDD(b)=find(bwlabel(CS_traceDD)==b,1);        
end
somMat=som_tr(C_start'+nTau);
somMatDD=som_tr_dd(C_start'+nTau);

DDMat=dd_tr(C_startDD'+nTau);
DDMatDD=dd_tr_dd(C_startDD'+nTau);
figure(8); clf;
nexttile([1 1])
plot(nTau,mean(somMat,1)); hold all
plot(nTau,mean(DDMat,1)); hold all

nexttile([1 1])
plot(nTau,mean(somMatDD,1)); hold all
plot(nTau,mean(DDMatDD,1)); hold all
