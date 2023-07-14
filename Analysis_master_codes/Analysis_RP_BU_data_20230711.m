% Analysis on AAV expression sample and plot, in house YQ201
% 2023/07/11, Byung Hun Lee

%clear
[fpath] = uigetfile_n_dir;
%%
for i=1:length(fpath)
load([fpath{i} '/output_data.mat'])
sz=double(Device_Data{1, 4}.ROI([2 4]));
mov=double(readBinMov([fpath{i} '/frames1.bin'],sz(2),sz(1)));
%mov=double(readBinMov_times([fpath{i} '/frames1.bin'],sz(2),sz(1),[1:30000]));
CamTrigger=find(Device_Data{1, 2}.Counter_Inputs.data(2:end)-Device_Data{1, 2}.Counter_Inputs.data(1:end-1));
mov_test=mov(:,:,150:250);
try mov_test = single(mov_test)./single(max(mov_test.data(:)));
catch disp('change to vm')

mov_test=vm(mov_test); mov_test = single(mov_test)./single(max(mov_test.data(:))); end
mov_test = movmean(mov_test,10,3);
mov_ref = squeeze(median(mov_test,3));
[mov_mc,xyField]=optical_flow_motion_correction_LBH(vm(mov(:,:,1:length(CamTrigger))),mov_ref,'normcorre');
mov_mc=vm(mov_mc);
mov_mc.transpose.savebin([fpath{i} '/mc.bin'])
%mcTrace = squeeze(mean(xyField,[1 2]));

shifts_r = squeeze(cat(3,xyField(:).shifts));
shifts_nr = cat(ndims(xyField(1).shifts)+1,xyField(:).shifts);
shifts_nr = reshape(shifts_nr,[],ndims(mov_mc)-1,size(mov_mc,3));
%shifts_x = squeeze(shifts_nr(:,2,:))';
%shifts_y = squeeze(shifts_nr(:,1,:))';
mcTrace=squeeze(shifts_nr)';

save([fpath{i} '/mcTrace.mat'],'mcTrace')
end
