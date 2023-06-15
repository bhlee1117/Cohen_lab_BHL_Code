% Simple analysis on AAV expression sample and plot, in house YQ201
% 2023/02/21, Byung Hun Lee
clear

addpath('/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Code')
[fpath] = uigetfile_n_dir(); %select folders to analyze

%% Correct Motion

for i=1%:length(fpath)

mov=vm(fpath{i}); %load movie
load(fullfile(fpath{i},'settings.mat'));
mov_test=mov(:,:,150:250);

try mov_test = single(mov_test)./single(max(mov_test.data(:))); 
catch disp('change to vm')
mov_test=vm(mov_test); mov_test = single(mov_test)./single(max(mov_test.data(:))); end
%/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Code/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Code
mov_test = movmean(mov_test,10,3);
mov_ref = squeeze(median(mov_test,3));
[mov_mc,xyField]=optical_flow_motion_correction_LBH((mov),mov_ref,'optic_flow');

clear mov
mov_mc=vm(mov_mc);
mov_mc.transpose.savebin([fpath{i} '/mc.bin']) %save movie

mcTrace = squeeze(mean(xyField,[1 2]));
save([fpath{i} '/mcTrace.mat'],'mcTrace') %save motion trace

clear mov_mc
delete([fpath{i} '/Sq_camera.bin'])
end
%% Find Cells
for i=1%:length(fpath)
    load(fullfile(fpath{i},'settings.mat'));
    %load(fullfile(fpath{i},'mcTrace.mat'));
    Sz = importdata([fpath{i} '/experimental_parameters.txt']);
    sz1=Sz.data(1); sz2=Sz.data(2);
    mov_mc=double(readBinMov_times([fpath{i} '/mc.bin'],sz2,sz1,[600:900])); %get Movie part

    im_G=imgaussfilt(mean(mov_mc,3),2);
    [centers radii]=Cell_segment_circle_2x(im_G,0.95);
    Result{i}.centers=cell_detection_manual(mean(mov_mc,3),centers,[]);
end

%% Extract Voltage signals

for i=1:length(fpath)
    load(fullfile(fpath{i},'settings.mat'));
    load(fullfile(fpath{i},'mcTrace.mat'));
    Sz = importdata([fpath{i} '/experimental_parameters.txt']);
    sz1=Sz.data(1); sz2=Sz.data(2);
    mov_mc=double(readBinMov([fpath{i} '/mc.bin'],sz2,sz1));
    a=find(DAQ_waves.amplitude(1,[1:500])); frm_rate=1/((a(2)-a(1)-2)*1e-5);
    Result{i}.Blue=DAQ_waves.amplitude(4,round([1:size(mov_mc,3)]*(1/frm_rate+0.00002)/1e-5));

    %mov_res=SeeResiduals(mov_mc,squeeze(mean(movmean(mov_mc,200,3),[1 2])));
    mov_res= mov_mc-mean(mov_mc,3); %regress out motion
    mov_res = SeeResiduals(mov_res,mcTrace);
    mov_res = SeeResiduals(mov_res,mcTrace.^2);
    mov_res = SeeResiduals(mov_res,mcTrace(:,1).*mcTrace(:,2));

    Result{i}.c_ftprnt=mask_footprint(Result{i}.centers,movmean(mov_res(:,:,1000:end),10,3),[],6); %find footprints
    Result{i}.ref_im=mean(mov_mc,3);
    show_footprnt(Result{i}.c_ftprnt,mov_mc)
    Result{i}.coord=get_coord(Result{i}.c_ftprnt);

    Result{i}.traces=-(tovec(mov_mc)'*tovec(Result{i}.c_ftprnt))'; %weighted product
    Result{i}.traces_bin=-(tovec(mov_res)'*tovec(Result{i}.c_ftprnt>0))'; %binary

    Result{i}.traces_hi=squeeze(Result{i}.traces) - movmedian(squeeze(Result{i}.traces),150,2);
    Result{i}.traces_bin_hi=squeeze(Result{i}.traces_bin) - movmedian(squeeze(Result{i}.traces_bin),200,2);
    %traces_hi=squeeze(traces) - mean(squeeze(traces),2);
    Result{i}.noise=get_threshold(Result{i}.traces,1);
    %traces_hi=traces_hi./range(traces_hi,2);
    Result{i}.traces_hi=Result{i}.traces_hi./Result{i}.noise;
    Result{i}.traces_bin_hi=Result{i}.traces_bin_hi./get_threshold(Result{i}.traces_bin_hi,1);
    ns=size(Result{i}.traces,1); Result{i}.mcTrace=mcTrace;

    Result{i}.traces_res=squeeze(SeeResiduals(reshape(Result{i}.traces,ns,1,[]),Result{i}.mcTrace')); %regress out motion
    Result{i}.traces_res=squeeze(SeeResiduals(reshape(Result{i}.traces_res,ns,1,[]),Result{i}.mcTrace.^2'));
    Result{i}.traces_res=squeeze(SeeResiduals(reshape(Result{i}.traces_res,ns,1,[]),(Result{i}.mcTrace(:,1).*Result{i}.mcTrace(:,2))'));

    Result{i}.spike=[];
    tmp=squeeze(Result{i}.traces_res) - movmedian(squeeze(Result{i}.traces_res),4,2); tmp=tmp./get_threshold(tmp,1,2000);
    Result{i}.spike=find_spike_bh(tmp,4,3); %single spike
    tmp=squeeze(Result{i}.traces_res) - movmedian(squeeze(Result{i}.traces_res),100,2); tmp=tmp./get_threshold(tmp,1,2000);
    Result{i}.spike=(Result{i}.spike+find_spike_bh(tmp,3,1.5))>0; %complex spikes
end

%% Plot traces

figure; i=1;
scale=10; t=[1:size(Result{i}.traces_hi,2)]/frm_rate;
tr=Result{i}.traces-median(Result{i}.traces,2); fprnt=Result{i}.c_ftprnt;
%tr=tr./range(tr,2);
tr=tr./get_threshold(Result{i}.traces,1);
%order=[1:size(tr,1)];
order=find(sum(Result{i}.spike,2)>0)';
%order=[47 setdiff([1:size(tr,1)],47)];

tiledlayout(10,4)

ax1 = nexttile([2 2]);
colr = max(colormap(jet(size(fprnt,3))),0);
imshow2(squeeze(sum(fprnt.*reshape(colr,1,1,[],3),3)),[]); hold all;
text(Result{i}.coord(:,1)',Result{i}.coord(:,2)',num2str([1:size(fprnt,3)]'),'color','w')

ax4 = nexttile([2 2]);
imshow2(Result{i}.ref_im,[]); colormap('gray')
linkaxes([ax1 ax4],'xy')

ax2 = nexttile([6 4]);

lines=plot(t,tr(order,:)'+[1:size(order,2)]*scale);
arrayfun(@(l,c) set(l,'Color',c{:}),lines,num2cell(jet(size(order,2)),2))
axis tight
hold all
S=tr; S(~Result{i}.spike)=NaN;
plot(t,S(order,:)'+[1:size(order,2)]*scale,'r.')
set(gca,'ytick',[1:size(order,2)]*scale,'yticklabel',order)
%title(fpath_mac2window(fpath{i}),'Interpreter', 'none')

ax3 = nexttile([2 4]);
plot(t,Result{i}.Blue)
%show_multiROI_waveform(DAQ_waves,3);

linkaxes([ax2 ax3],'x')
saveas(gcf,[fpath{i} '/voltage_trace_plot'])
save([fpath{1} '/result.mat'],'Result','fpath')

%%
i=1;
mov_res= mov_mc-movmean(mov_mc,150,3);
tr=Result{i}.traces-movmedian(Result{i}.traces,150,2);
sz=size(mov_res);
corr_im_proj=zeros(sz(1),sz(2),size(Result{i}.traces,1));
coord=round(Result{i}.coord);
mov_res_vec=tovec(mov_res);
for n=1:size(Result{i}.traces,1)
    n
    %     try
    %     corr_im=normxcorr2(tr(n,:),tovec(mov_res(coord(n,2)-crop_xy(2):coord(n,2)+crop_xy(2),...
    %                                              coord(n,1)-crop_xy(1):coord(n,1)+crop_xy(1),:)));
    %     corr_im_proj(coord(n,2)-crop_xy(2):coord(n,2)+crop_xy(2),coord(n,1)-crop_xy(1):coord(n,1)+crop_xy(1),n)=...
    %           max(reshape(corr_im,2*crop_xy(2)+1,2*crop_xy(1)+1,[]),[],3);
    %     end
    xx=zeros(size(mov_res_vec,1),41);
    for j=1:size(mov_res_vec,1)
        xx(j,:)=xcorr(tr(n,4000:end),-mov_res_vec(j,4000:end),20,'normalized');
    end
    corr_im{n}=xx;
    [corr_im_proj(:,:,n) arg(:,:,n)]=max(reshape(corr_im{n},sz(1),sz(2),[]),[],3);
end

%%

test=corr_im_proj;
for j=1:32
    t=test(:,:,j);
    imG=imgaussfilt(corr_im_proj(:,:,j),3);
    t(imG<0.07)=0;
    t=(t-min(t(:)))./(max(t(:))-min(t(:)));
    test(:,:,j)=t;
end
show_footprnt(reshape(tovec(corr_im_proj)./median(tovec(corr_im_proj),2),304,512),mov_mc)


%% Write movie
writeMov('DMD2_stim',-movmean(mov_res(:,:,4000:end),30,3),[],200,1,[-5 70],[],Blue(4000:end))
%writeMov('DMD3_stim',-movmean(mov_res(:,:,8701:17400),30,3),[],200,1.25,[-5 70],[],Blue(8701:17400))
%writeMov('DMD23_stim',-movmean(mov_res(:,:,17401:end),30,3),[],200,1.25,[-5 70],[],Blue(17401:end))


