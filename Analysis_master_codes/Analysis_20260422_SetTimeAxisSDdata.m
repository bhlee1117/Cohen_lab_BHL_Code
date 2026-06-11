%% Load file path
clear
clc;
[~, ~, raw] = xlsread(['/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/Prism_V2+Glu_Data_Arrangement.xlsx'], 'Sheet1', 'C5:AA31');

fpath=raw(:,1)';
V2moviemaxTime=15000;
GlumoviemaxTime=5000;
foi=[3 5 8 9 10 12 13];
%%
for f=foi

Devicedata_filename=fullfile(fpath{f},'output_data.mat');
load(Devicedata_filename);

o_laser = 200001;
DAQ_rate=Device_Data{1, 2}.Counter_Inputs(1, 1).rate;

% voltage camera clock
cam1_vsyn = Device_Data{1, 2}.Counter_Inputs(1, 1).data;
start_idx = min(find(cam1_vsyn ==max (cam1_vsyn)));
end_idx = length(Device_Data{1, 2}.buffered_tasks(1, 2).channels(1, 1).data);
CamTrigger1_DAQax=find((cam1_vsyn(2:end)-cam1_vsyn(1:end-1))>0);
segment_size=unique(diff(CamTrigger1_DAQax));
last_val = cam1_vsyn(start_idx - 1);
n_to_add = end_idx - start_idx + 1;
n_segments = ceil(n_to_add / segment_size);
added_part = repelem((last_val + 1 : last_val + n_segments)', segment_size);
added_part = added_part(1:n_to_add);
cam1_vsyn(start_idx:end_idx) = added_part;

% Prepare
VmovTimesegments=[(cam1_vsyn(o_laser)+2):V2moviemaxTime:cam1_vsyn(end)];
VmovTimesegments(end+1)=cam1_vsyn(end);
nFrame2analyze=VmovTimesegments(end)-VmovTimesegments(1);
CamTrigger1_DAQaxVec=ind2vec(length(cam1_vsyn),CamTrigger1_DAQax(1):segment_size:cam1_vsyn(end)*segment_size,1);
CamTrigger1_DAQaxVec=CamTrigger1_DAQaxVec(o_laser:end);
FirstFrameDAQax=find(CamTrigger1_DAQaxVec>0,1);
t_vol=[FirstFrameDAQax:segment_size:(FirstFrameDAQax+segment_size*nFrame2analyze)]/DAQ_rate;

clear Result
load(fullfile(fpath{f}, 'Volt_Result.mat'), 'Result');
Result.t_ax=t_vol;
save(fullfile(fpath{f}, 'Volt_Result.mat'), 'Result', '-v7.3');


% Glutamate image motion correction

cam2_vsyn = Device_Data{1, 2}.Counter_Inputs(1, 2).data;
cam2_vsyn_trig = find (diff(cam2_vsyn) == 1)+1;
segment_size2 = cam2_vsyn_trig (10) - cam2_vsyn_trig(9);
start_idx = min(find(cam2_vsyn ==max (cam2_vsyn)));
end_idx = length(Device_Data{1, 2}.buffered_tasks(1, 2).channels(1, 2).data);
last_val = cam2_vsyn(start_idx - 1);
n_to_add = end_idx - start_idx + 1;
n_segments = ceil(n_to_add / segment_size2);
added_part = repelem((last_val + 1 : last_val + n_segments)', segment_size2);
added_part = added_part(1:n_to_add);
cam2_vsyn(start_idx:end_idx) = added_part;

mod488 = 200001;
mod607 = 200001;
cam2_trig = find (diff(cam2_vsyn) ==1)+1; %DAQ frame
%cam2_trig = cam1_vsyn (cam2_trig); %Cam1 frame
t_start= cam2_trig (cam2_vsyn(mod488)+2)-mod488+1;
if cam2_vsyn(end)<GlumoviemaxTime
    GluemovTimesegments=[(cam2_vsyn(mod488)+2) cam2_vsyn(end)];
else
    GluemovTimesegments=[(cam2_vsyn(mod488)+2):GlumoviemaxTime:cam2_vsyn(end)];
    GluemovTimesegments(end+1)=cam2_vsyn(end);
end
nFrame2analyze2=GluemovTimesegments(end)-GluemovTimesegments(1);
t_glu=[t_start:segment_size2:(t_start+segment_size2*nFrame2analyze2)]/DAQ_rate;

clear Result GluResult
GluResult=importdata(fullfile(fpath{f}, 'Glu_Result.mat'));
GluResult.t_ax=t_glu;
save(fullfile(fpath{f}, 'Glu_Result.mat'), 'GluResult', '-v7.3');
disp(['Saved,' fullfile(fpath{f}, 'Glu_Result.mat')]);
end