clear; clc;
fpath = '/Volumes/cohen_lab/Lab/Labmembers/Xiang Wu/Data/20250131/180409XW_C008_VR1';

cd(fpath);

%%
t_seg=[5000:20000];
load(fullfile(fpath,"output_data.mat"))
load(fullfile(fpath,'cal_mcTrace.mat'));
load(fullfile(fpath,'cal_new_roi.mat'));


mov_mc=double(readBinMov_times(fullfile(fpath,'calcium_mc.bin'),nRow_cal_new,nCol_cal_new,t_seg));

[~, tr]=clicky(mov_mc); % calcium trace
mov_res=SeeResiduals(mov_mc,movmedian(tr,200));

artifactImage=toimg(tovec(mov_mc-mean(mov_mc,3))*movmedian(tr,200),size(mov_mc,1),size(mov_mc,2));
imshow2(artifactImage,[])