clear
clc;
cd '/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/FromBackup/PP72_PlaceCellResults';
[~, ~, raw] = xlsread(['/Volumes/cohen_lab/Lab/Labmembers' ...
    '/Byung Hun Lee/Data/PrismPCdata_Arrangement.xlsx'], 'Sheet1', 'C5:P27');

% [~, ~, NeuronsToUse]=xlsread(['/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/' ...
%     'PlaceCellData_Arrangement.xlsx'], 'Sheet1', 'L8:M46');
%
% NeuronsToUse=cellfun(@(x) (str2num(num2str(x))),NeuronsToUse,'UniformOutput',false);
ref_ROI=cellfun(@(x) (str2num(num2str(x))),raw(:,10),'UniformOutput',false);
basal_ROI=cellfun(@(x) (str2num(num2str(x))),raw(:,11),'UniformOutput',false);
apical_ROI=cellfun(@(x) (str2num(num2str(x))),raw(:,12),'UniformOutput',false);
fpath=raw(:,1)';
StructureData=raw(:,10);
EndFrame=cell2mat(raw(:,13));
ifmotionReject=cell2mat(raw(:,14));
save_figto='/Volumes/BHL_WD18TB/PP72_PlaceCellResults';
place_bin=150; time_segment=15000; overlap=200;
alignedMovFN = {'STA_Mat_SS','STA_Mat_CS','STA_Mat_dSP'};
bound=6;
%%
f=20; 

load(fullfile(fpath{f},'PC_Result.mat'),'Result')
load(fullfile(fpath{f},"output_data.mat"))
sz=double(Device_Data{1, 3}.ROI([2 4]));
j=20;
mov_mc=double(readBinMov([fpath{f} '/mc_ShutterReg' num2str(j,'%02d') '.bin'],sz(2),sz(1)));
load([fpath{f} '/mcTrace' num2str(j,'%02d') '.mat']);

mc=mcTrace.xymean;
bkg = zeros(2, size(mov_mc,3));
bkg(1,:) = linspace(-1, 1, size(mov_mc,3));  % linear term
bkg(2,:) = linspace(-1, 1, size(mov_mc,3)).^2;  % quadratic term

mov_res=mov_mc-mean(mov_mc,3);
mov_res = SeeResiduals(mov_res,mc);
mov_res = SeeResiduals(mov_res,mc.^2);
mov_res = SeeResiduals(mov_res,mc(:,1).*mc(:,end));
mov_res= SeeResiduals(mov_res,bkg,1);

%%
bound=15;
mov_int=imresize(mov_res(bound:end-bound,bound:end-bound,1000:5000),1/2);
mov_int_filt=imgaussfilt(mov_int-imgaussfilt(mov_int,1),5);
mov_int_filt_vec=tovec(mov_int_filt);
mov_int_filt_vec=movmedian(mov_int_filt_vec,30,2);
mov_int_filt2=toimg(mov_int_filt_vec,size(mov_int_filt,1),size(mov_int_filt,2));


conv_result=[];
for k=1:size(mov_int_filt2,3)
conv_result(:,:,k) = conv2(mov_int_filt2(:,:,k), z, 'same');
end
moviefixsc(conv_result)

%%
filtered_frame=[];
for k=1:size(mov_int_filt2,3)
filtered_frame(:,:,k) = bandpass_2d(double(mov_res(:,:,1000+k)), 5,15);
end
moviefixsc(filtered_frame)

%%
bound=10;
mov_int=imresize(mov_res(bound:end-bound,bound:end-bound,:),1/2);
mov_int_filt=imgaussfilt3(mov_int,[8 8 5])-imgaussfilt3(mov_int,[3 3 5]);
figure(3); clf;
moviefixsc([mat2gray(mov_int); mat2gray(mov_int_filt)])

z2=mov_int_filt(3:43,45:85,1864);
N = 20;         % Grid size
x = linspace(-N, N, 2*N+1);
y = linspace(-N, N, 2*N+1);
[X, Y] = meshgrid(x, y);
r = sqrt(X.^2 + Y.^2);

z = ( sin(r/6) ./ ( r/6)).^2;
z(r == 0) = 1;  % Avoid division by zero at the center
z = -z / max(z(:));
imagesc(imfuse(z,z2))
match_img=[]; match_data=[];
for k=1:size(mov_int_filt,3)
[ match_data{k}, match_img(:,:,k)] = matchPattern(mov_int_filt(:,:,k), z2, 0.5, 2);
end
match_data=cellfun(@(x) x+[21 21 0],match_data,'UniformOutput',false);
filteredPtl=trackParticle(match_data,5,50);
traveldist=cell2mat(cellfun(@(x) sum(sqrt(sum((x(2:end,1:2)-x(1:end-1,1:2)).^2,2)),1),filteredPtl,'UniformOutput',false));
valid_ptl=traveldist>30;

clf;
imshow2(Result.ref_im,[])
hold all
for p=find(valid_ptl)
plot(filteredPtl{p}(:,2)*2,filteredPtl{p}(:,1)*2)
end
%%
ref=imresize(mov_res(:,:,1864),1/2);
ref1=imgaussfilt(ref,8);
ref2=imgaussfilt(ref,3);
figure(3); clf;
nexttile([1 1])
imshow2(ref1,[])
nexttile([1 1])
imshow2(ref2,[])
nexttile([1 1])
imshow2(ref1-ref2,[])

