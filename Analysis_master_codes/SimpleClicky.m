clear;
clc;
fpath='/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/20221207_BHLm007_Treadmill/173155_BHLm007_Sqr_Train_5Hz_Bp6';
cd(fpath)
file = dir('*.txt');              % Find txt file in folder
txt = fileread(file(1).name);     % Read the first txt file found
nums = sscanf(txt, '%*[^0-9]%d%*[^0-9]%d');

horizontal = nums(1);
vertical   = nums(2);
pixelsize=6.5; % µm

fprintf('Horizontal pixel: %d\n', horizontal);
fprintf('Vertical pixel: %d\n', vertical);
%%
load("settings.mat");
a=find(DAQ_waves.amplitude(1,[1:500])); frm_rate=1/((a(2)-a(1)-2)*1e-5); 
BlueAOTF=DAQ_waves.amplitude(6,round([1:size(mov,3)]*(1/frm_rate+0.00002)/1e-5));
mov=readBinMov('Sq_camera.bin',vertical,horizontal);
[roi,tr]=clicky_faster(mov);

timeaxis=[1:size(mov,3)]/frm_rate;
