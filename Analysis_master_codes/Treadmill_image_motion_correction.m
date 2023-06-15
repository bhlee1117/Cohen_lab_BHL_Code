clear
clc;

[fpath] = uigetfile_n_dir();
%%
for i=1:length(fpath)

load(fullfile(fpath{i},"settings.mat"))
Sz = importdata([fpath{i} '/experimental_parameters.txt']);
sz1=Sz.data(1); sz2=Sz.data(2);

% [a,b]=system(sprintf('GetFileInfo "%s"',fullfile(fpath{i},'settings.mat'))); s=strfind(b,'modified:')+10; crdat=b(s:s+18);
%     image_start=datestr(datenum(crdat));
%     dt=datetime(image_start(end-8:end),'InputFormat','HH:mm:ss'); dt.Format = 'HH:mm:ss';
%     fnm=dir(fullfile(fpath{i}, '*.csv'));
% 
%     a=find(DAQ_waves.amplitude(1,[1:500])); frm_rate=1/((a(2)-a(1)-2)*1e-5);
%     frm_end=sum(DAQ_waves.amplitude(1,:)); f_seg=[[1:10000:frm_end] frm_end+1];
%     [Arduino_data reward_pos lap_dist]=match_treadmill_DAQ(fullfile(fpath{i},fnm(1).name),1e-5,DAQ_data,(1/frm_rate+0.00002),2); % time, treadmill, Reward, Run
%      Arduino_data(:,1)=Arduino_data(:,1)*(1/frm_rate*1000)/(1/frm_rate*1000+0.02);
%    aa=find(bwlabel(Arduino_data(:,4))==1);
aa=[55000:56100];
frm_end=sum(DAQ_waves.amplitude(1,:)); f_seg=[[1:10000:frm_end] frm_end+1];
if length(aa)>1000
    mov_test=double(vm([fpath{i}],[aa(1):aa(1)+1000]));
else
    mov_test=double(vm([fpath{i}],[aa(1):aa(end)]));
end
options_rigid = NoRMCorreSetParms('d1',size(mov_test,1),'d2',size(mov_test,2),'bin_width',200,'max_shift',30,'us_fac',50,'init_batch',200);
tic; [mov_test,shifts1,template1,options_rigid] = normcorre(mov_test,options_rigid); toc
mov_test=vm(mov_test);
mov_test = single(mov_test)./single(max(mov_test.data(:)));
mov_test = movmean(mov_test,10,3);
mov_ref = squeeze(median(mov_test,3));

for j=3:length(f_seg)-1
    g=1;
    try 
    mov=vm([fpath{i}],[f_seg(j):f_seg(j+1)+10]); 
    catch % when the image ends
    g=2;
    mov=vm([fpath{i}],[f_seg(j):f_seg(j+1)-1]);    
    end
    [mov_mc,xyField]=optical_flow_motion_correction_LBH(mov,mov_ref,'optic_flow');

     options_rigid = NoRMCorreSetParms('d1',size(mov_mc,1),'d2',size(mov_mc,2),'bin_width',200,'max_shift',30,'us_fac',50,'init_batch',200);
     tic; [mov_mc,shifts1,template1,options_rigid] = normcorre(mov_mc,options_rigid); toc

    ave_im=mean(mov_mc,3);
    mov_mc=vm(mov_mc);
    mov_mc.transpose.savebin([fpath{i} '/mc' num2str(j,'%02d') '.bin'])
    
    mcTrace = squeeze(mean(xyField,[1 2]));
    save([fpath{i} '/mcTrace' num2str(j,'%02d') '.mat'],'mcTrace','ave_im')

  %  clear mov_mc mov
end
end


