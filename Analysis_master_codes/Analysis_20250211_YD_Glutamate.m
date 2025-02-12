clear; clc;
fpath='/Volumes/cohen_lab/Lab/Labmembers/Yangdong Wang/Data/NorthwestScreening/20250204/131724YD039_2';
cd(fpath)
%%
load(fullfile(fpath,'output_data.mat'))
nrows = Device_Data{1, 2}.ROI(4);
ncols = Device_Data{1, 2}.ROI(2);
mov = readBinMov(fullfile(fpath,'frames1.bin'), nrows, ncols);

CamCounter=Device_Data{1, 4}.Counter_Inputs(1, 1).data;
CamTrigger=find(CamCounter(2:end)-CamCounter(1:end-1));

frm_rate=1/((CamTrigger(2)-CamTrigger(1))/Device_Data{1, 2}.hsync_rate);


%% photobleaching correction & cut front and back

meanF=squeeze(mean(mov,[1 2]));
meanF=meanF(200:end-100);
[y_fit]=expfitDM_2([1:length(meanF)]',meanF,[1:length(meanF)]',1000);
mov_res=SeeResiduals(double(mov(:,:,200:end-100)),y_fit);
% movhi=movmean(mov-movmedian(mov,200,3),5,3);
% movhi=movhi(:,:,200:end);
%% 2d bandpass filter (spatial filter)
filtered_image=bandpass_2d(mov_res,10,20);
filtered_image=imgaussfilt(filtered_image,2);
%% convolution with response function (temporal filter)
load('iGluSnFR3_response_fun.mat','response_fnt')
conv_mov=conv2(tovec(filtered_image),response_fnt');
conv_mov=toimg(conv_mov,size(mov_res,1),size(mov_res,2));
Lag1_cov=mean(conv_mov(:,:,2:end).*conv_mov(:,:,1:end-1),3);
%%
figure; 
maxConvMov=grs2rgb(max(conv_mov,[],3),colormap('hot'),11000,30000);
ColorImg = grs2rgb(double(mean(mov,3)),colormap(gray))+maxConvMov;

nexttile([1 1])
imshow2(ColorImg,[])


