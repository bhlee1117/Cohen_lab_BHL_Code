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
fpath=raw(:,1)';
StructureData=raw(:,10);
EndFrame=cell2mat(raw(:,13));
ifmotionReject=cell2mat(raw(:,14));
save_figto='/Volumes/BHL_WD18TB/PP72_PlaceCellResults';
place_bin=150; time_segment=15000; overlap=200;
alignedMovFN = {'STA_Mat_SS','STA_Mat_CS','STA_Mat_dSP'};
%%
bound=7;
f=23;
load(fullfile(fpath{f},'PC_Result.mat'))
load(fullfile(fpath{f},"output_data.mat"))
sz=double(Device_Data{1, 3}.ROI([2 4]));
frm_end=EndFrame(f);
f_seg=[[1:time_segment:frm_end] frm_end+1]; f_seg(2:end)=f_seg(2:end)-1;
take_window=repmat([1 time_segment],length(f_seg)-1,1);
take_window(2:end,1)=take_window(2:end,1)+overlap; take_window(1:end-1,2)=take_window(1:end-1,2)+overlap;
take_window(end)=mod(f_seg(end),time_segment);
take_window(take_window==0)=time_segment;
Blue_on_Seg=unique(ceil(find(Result.Blue)/time_segment));

j=3;
mov_mc=double(readBinMov([fpath{f} '/mc_ShutterReg' num2str(j,'%02d') '.bin'],sz(2),sz(1)));
load([fpath{f} '/mcTrace' num2str(j,'%02d') '.mat']);

mov_mc=mov_mc(:,:,[take_window(j,1):take_window(j,2)]);
mc=mcTrace.xymean([take_window(j,1):take_window(j,2)],:);

if ismember(j,Blue_on_Seg)
    isfirst_seg=double(j==1);
    [~, blueomitTr]=get_blueoffTrace(squeeze(mean(mov_mc,[1 2])), ...
        Result.Blue(f_seg(j)+take_window(j,1)-isfirst_seg:f_seg(j)+take_window(j,2)-isfirst_seg),30);
    bkg = zeros(1, size(mov_mc,3));
    bkg(1,:)=movmedian(blueomitTr,4000,'omitnan');
else
    bkg = zeros(2, size(mov_mc,3));
    bkg(1,:) = linspace(-1, 1, size(mov_mc,3));  % linear term
    bkg(2,:) = linspace(-1, 1, size(mov_mc,3)).^2;  % quadratic term
end

mov_res=mov_mc;
mov_res= mov_mc-median(mov_mc,3);
mov_res = SeeResiduals(mov_res,mc);
mov_res = SeeResiduals(mov_res,mc.^2);
mov_res = SeeResiduals(mov_res,mc(:,1).*mc(:,end));
mov_res= SeeResiduals(mov_res,bkg,1);
%%
frameInterest=[5000:10000];
mov_res_mask=mov_res(:,:,frameInterest).*double(max(Result.bvMask,[],3)==0);
subMov=tovec(imresize(imgaussfilt3(mov_res_mask(bound:end-bound,bound:end-bound,:),[1 1 0.1]),1/4));
subMov=subMov-mean(subMov,2);
covMat=subMov'*subMov;

[V, D] = eig(covMat);
D = diag(D);
D = D(end:-1:1);
V = V(:,end:-1:1);
vSign = sign(max(V) - max(-V));  % make the largest value always positive
V = V.*vSign;

nPCs=28;
eigImg=toimg(tovec(mov_res(:,:,frameInterest))*V(:,1:nPCs),size(mov_res,1),size(mov_res,2));
figure(4); clf;
for n=1:nPCs
    nexttile([1 1])
    imshow2(eigImg(:,:,n),[])
    title(['PC #', num2str(n), ' Fraction : ' num2str(D(n)/sum(D),2)])
end

[V_ics, mixmat, sepmat]=sorted_ica(V(:,1:nPCs),10);
icsImg=toimg(tovec(mov_res(:,:,frameInterest))*V_ics,size(mov_res,1),size(mov_res,2));
figure(5); clf;
for n=1:size(V_ics,2)
    nexttile([1 1])
    %show_footprnt_contour(Result.bvMask,icsImg(:,:,n))
    imshow2(icsImg(:,:,n),[])
    title(['ICS #', num2str(n)])
end
colormap('gray')