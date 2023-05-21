%FF experiment, 20230206
addpath('Y:\Lab\Labmembers\Byung Hun Lee\Code')
addpath('Y:\Lab\Labmembers\Byung Hun Lee\Code\Codes(~2021)\uigetfile_n_dir')
[fpath] = uigetfile_n_dir();

%% motion correction
for i=1:length(fpath)
mov=vm(fpath{i});
load(fullfile(fpath{i},'settings.mat'));
mov_test=mov(:,:,150:250);
try mov_test = single(mov_test)./single(max(mov_test.data(:)));
catch disp('change to vm')

mov_test=vm(mov_test); mov_test = single(mov_test)./single(max(mov_test.data(:))); end
mov_test = movmean(mov_test,10,3);
mov_ref = squeeze(median(mov_test,3));
[mov_mc,flow_xy_out]= optical_flow_motion_correction(mov,mov_ref);
%[mov_mc,xyField]=optical_flow_motion_correction_LBH((mov),mov_ref,'optic_flow');
clear mov
mov_mc=vm(mov_mc);
mov_mc.transpose.savebin([fpath{i} '/mc.bin'])
mcTrace = squeeze(mean(xyField,[1 2]));
save([fpath{i} '/mcTrace.mat'],'mcTrace')
clear mov_mc
delete([fpath{i} '/Sq_camera.bin'])
end

%% Signal extraction

for i=1:length(fpath)
load(fullfile(fpath{i},'settings.mat')); load(fullfile(fpath{i},'mcTrace.mat'));
Sz = importdata([fpath{i} '/experimental_parameters.txt']); sz1=Sz.data(1); sz2=Sz.data(2);
mov_mc=double(readBinMov([fpath{i} '/mc.bin'],sz2,sz1));
a=find(DAQ_waves.amplitude(1,[1:500])); frm_rate=1/((a(2)-a(1)-2)*1e-5);

if ~isempty(dmd_mask_sequence_rois)
[~, Result{i}.Blue]=show_multiROI_waveform(DAQ_waves,size(dmd_mask_sequence_rois,2),frm_rate,false);
else
Result{i}.Blue=DAQ_waves.amplitude(4,round([1:size(mov_mc,3)]*(1/frm_rate+0.00002)/1e-5));
end

mov_res= mov_mc-mean(mov_mc,3);
mov_res = SeeResiduals(mov_res,mcTrace); mov_res = SeeResiduals(mov_res,mcTrace.^2);
mov_res = SeeResiduals(mov_res,mcTrace(:,1).*mcTrace(:,2));

im_G=imgaussfilt(mean(mov_mc,3),2);
[centers radii]=Cell_segment_circle_2x(im_G,0.95);
Result{i}.centers=cell_detection_manual(mean(mov_mc,3),centers);
Result{i}.c_ftprnt=mask_footprint(Result{i}.centers,movmean(mov_res(:,:,1000:end),10,3),[],6);
Result{i}.ref_im=mean(mov_mc,3);
Result{i}.coord=get_coord(Result{i}.c_ftprnt);

if ~isempty(dmd_mask_sequence_rois)
dmd_mask_sequence_rois{1}=dmd_mask_sequence_rois{1};
dmd_mask_sequence_rois{2}=dmd_mask_sequence_rois{1};

for j=1:size(dmd_mask_sequence_rois,2)
roi=dmd_pat2roi(dmd_mask_sequence_rois{j},size(Result{i}.ref_im));
Result{i}.clist_p{j}=find(roi(sub2ind(size(Result{i}.ref_im),round(Result{i}.coord(:,2)),round(Result{i}.coord(:,1)))));
end
end
Result{i}.frm_rate=frm_rate;
Result{i}.DMDPatt=dmd_mask_sequence_rois;
Result{i}.traces=-(tovec(mov_mc)'*tovec(Result{i}.c_ftprnt))';
Result{i}.traces_bin=-(tovec(mov_mc)'*tovec(Result{i}.c_ftprnt>0))';
Result{i}.traces_hi=squeeze(Result{i}.traces) - movmedian(squeeze(Result{i}.traces),150,2);

Result{i}.noise=get_threshold(Result{i}.traces,1);
Result{i}.traces_hi=Result{i}.traces_hi./Result{i}.noise;

tmp=squeeze(Result{i}.traces) - movmedian(squeeze(Result{i}.traces),20,2); tmp=tmp./get_threshold(tmp,1);
Result{i}.spike=find_spike_bh(tmp,4,1.5);
tmp=squeeze(Result{i}.traces) - movmedian(squeeze(Result{i}.traces),100,2); tmp=tmp./get_threshold(tmp,1);
Result{i}.spike=(Result{i}.spike+find_spike_bh(tmp,4,1.5))>0;
end
save([fpath{i} '/Result.mat'],'Result','fpath')
%% Plot traces
figure;
show_footprnt(Result{1}.c_ftprnt,mean(mov_mc,3))
show_traces_spikes(Result{1}.traces,Result{1}.spike,Result{1}.Blue)
save(gca,fullfile(fpath{1},'Voltage_trace_plot.fig'))
%% PCA/ICA
n_comp=15;
[CAs ic_sub]=assembly_PCAICA(Result{1},30*1e-3,2,n_comp,1);

order=[];
for i=1:n_comp; order=[order; find(CAs(:,i))];
end
S_trace=movsum(Result{1}.spike(order,:),30/1.25,2);
imshow2(corrcoef(S_trace'),[])
hold all
for i=1:n_comp
    line([sum(CAs(:,1:i),[1 2])+0.5 sum(CAs(:,1:i),[1 2])+0.5],[0.5 112.5],'color','r')
    line([0.5 112.5],[sum(CAs(:,1:i),[1 2])+0.5 sum(CAs(:,1:i),[1 2])+0.5],'color','r')
end
%% Calculate Cross-correlation Patterns
[Patt]=Assembly_cross_corr(Result{1},100,0.25,15);
%STA=make_STA(mov_res,Patt{14}.Corr_frm(:,2)'-50,[-100:150]);
%% Number of cells stimulaion

c_ftprnt=Result{1}.c_ftprnt;
Patt_ind=2;
L = hDMDApp.dmd.device.height;
W = hDMDApp.dmd.device.width;
cam_offset = [hDMDApp.HorizontalOffsetEditField.Value hDMDApp.VerticalOffsetEditField.Value];
ref_im = hDMDApp.ref_im;
sq_size = 6;

coord=get_coord(c_ftprnt>0);
dist_threshold=[50 50 50 10 0];

for i=1%:5 % percentage
    %list=cluster.cluster_N{s};
    N_neuron_patt=length(Patt{Patt_ind}.maxtime);
    list=Patt{Patt_ind}.order(N_neuron_patt-10:N_neuron_patt);
    %list=Patt{Patt_ind}.order(1:4*i);
centers={coord,coord(list,:)};
stim_list(i,list)=1;
file_name_prefix = ['BHLm008_cluster' num2str(Patt_ind) '_Bp1'];
%hDragonflyApp.appendedfoldernameEditField.Value = [file_name_prefix '_r' num2str(4*i)];
hDragonflyApp.appendedfoldernameEditField.Value = [file_name_prefix '_last10neuron'];
[mask_sequence]=generate_dmd_sequence(centers,L,W,cam_offset,hDMDApp,5);
try
    hDragonflyApp.RunSynchronizedAQ([])
end
pause(5)
end

%%
L = hDMDApp.dmd.device.height;
W = hDMDApp.dmd.device.width;
cam_offset = [hDMDApp.HorizontalOffsetEditField.Value hDMDApp.VerticalOffsetEditField.Value];
ref_im = hDMDApp.ref_im;
sq_size = 7;
compete_ICA=[1 4;4 1;1 15;1 3;3 1;3 5;11 12];
coord=get_coord(c_ftprnt>0);
dist_threshold=[50 50 50 10 0];

for i=1:size(compete_ICA,1) % percentage
    %list=cluster.cluster_N{s};
    N_neuron_patt=length(Patt{Patt_ind}.maxtime);
    list1=find(CAs(:,compete_ICA(i,1)));
    list2=find(CAs(:,compete_ICA(i,2)));
    %list=Patt{Patt_ind}.order(1:4*i);
centers={coord,coord(list1,:),coord(list2,:)};
stim_list(i,list)=1;
file_name_prefix = ['BHLm017_' num2str(Patt_ind) '_B4'];
hDragonflyApp.appendedfoldernameEditField.Value = [file_name_prefix '_PCAICA_compete' num2str(compete_ICA(i,1)) 'v' num2str(compete_ICA(i,2))];
%hDragonflyApp.appendedfoldernameEditField.Value = [file_name_prefix '_last10neuron'];
[mask_sequence]=generate_dmd_sequence(centers,L,W,cam_offset,hDMDApp,sq_size);
try
    hDragonflyApp.RunSynchronizedAQ([])
end
pause(5)
end
save(['protocol_PCAICA_stim_20230206.mat'],'CAs','compete_ICA')
%save([file_name_prefix '_cluster_sequence.mat'],'CAs','ic_sub','perm_list','Result','coord')
