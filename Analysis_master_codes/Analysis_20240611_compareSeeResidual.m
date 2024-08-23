clear
clc;
cd '/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/FromBackup/PP72_PlaceCellResults';
[~, ~, raw] = xlsread(['/Volumes/cohen_lab/Lab/Labmembers' ...
    '/Byung Hun Lee/Data/PrismPCdata_Arrangement.xlsx'], 'Sheet1', 'C5:P23');

% [~, ~, NeuronsToUse]=xlsread(['/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/' ...
%     'PlaceCellData_Arrangement.xlsx'], 'Sheet1', 'L8:M46');
%
% NeuronsToUse=cellfun(@(x) (str2num(num2str(x))),NeuronsToUse,'UniformOutput',false);
ref_ROI=cellfun(@(x) (str2num(num2str(x))),raw(:,10),'UniformOutput',false);
fpath=raw(:,1)';
StructureData=raw(:,10);
EndFrame=cell2mat(raw(:,13));
ifmotionReject=cell2mat(raw(:,14));
place_bin=150; time_segment=15000; overlap=200;
f=18; j=2;

frm_end=EndFrame(f);
f_seg=[[1:time_segment:frm_end] frm_end+1]; f_seg(2:end)=f_seg(2:end)-1;
take_window=repmat([1 time_segment],length(f_seg)-1,1);
take_window(2:end,1)=take_window(2:end,1)+overlap; take_window(1:end-1,2)=take_window(1:end-1,2)+overlap;
take_window(end)=mod(f_seg(end),time_segment);


load(fullfile(fpath{f},'PC_Result.mat'),'Result')
load(fullfile(fpath{f},"output_data.mat"))
sz=double(Device_Data{1, 3}.ROI([2 4]));
mov_mc=double(readBinMov([fpath{f} '/mc_ShutterReg' num2str(j,'%02d') '.bin'],sz(2),sz(1)));
mov_mc=mov_mc(:,:,take_window(j,1):take_window(j,2));
[roi tr]=clicky(mov_mc);
load([fpath{f} '/mcTrace' num2str(j,'%02d') '.mat']);
SilentPeriod=ones(1,size(Result.normTraces,2));
SilentPeriod(find(max(Result.spike,[],1))'+[-10:100])=NaN;
SilentPeriod_seg=SilentPeriod(take_window(j,1)+time_segment*(j-1)-overlap:take_window(j,2)+time_segment*(j-1)-overlap);

mcTrace.xymean=mcTrace.xymean(take_window(j,1):take_window(j,2),:);
mc=mcTrace.xymean;
%mc= movmean(mcTrace.xymean-movmedian(mcTrace.xymean,500,1),3,1);
meanF=squeeze(mean(mov_mc,[1 2]));
t_fit=find(~isnan(SilentPeriod_seg));
[y_fit t_consts coeffY]  = expfitDM_2(t_fit',meanF(t_fit),[1:size(mov_mc,3)]',1000);

bkg = zeros(1, size(mov_mc,3));
bkg(1,:) = y_fit;  % linear term

bkg2 = zeros(2, size(mov_mc,3));
bkg2(1,:) = linspace(-1, 1, size(mov_mc,3));  % linear term
bkg2(2,:) = linspace(-1, 1, size(mov_mc,3)).^2;  % quadratic term

mov_res= mov_mc-median(mov_mc,3);
mov_res2 = SeeResiduals(mov_res,mc);
mov_res2 = SeeResiduals(mov_res2,mc.^2);
mov_res2 = SeeResiduals(mov_res2,mc(:,1).*mc(:,end));
mov_res3= SeeResiduals(mov_res2,bkg,1);
mov_res4= SeeResiduals(mov_res,mcTrace.xymean);
mov_res4 = SeeResiduals(mov_res4,mcTrace.xymean.^2);
mov_res4 = SeeResiduals(mov_res4,mcTrace.xymean(:,1).*mcTrace.xymean(:,end));
mov_res4= SeeResiduals(mov_res4,bkg2,1);
mov_res5= SeeResiduals(mov_res2,bkg2,1);

[tr]=apply_clicky(roi,mov_res,'no');
[tr2]=apply_clicky(roi,mov_res2,'no');
[tr3]=apply_clicky(roi,mov_res3,'no');
[tr4]=apply_clicky(roi,mov_res4,'no');
[tr5]=apply_clicky(roi,mov_res5,'no');

%%
figure(4); clf;
plot(tr(find(SilentPeriod_seg)),mc(find(SilentPeriod_seg),1),'.'); hold all
plot(tr(find(SilentPeriod_seg)),mc(find(SilentPeriod_seg),3),'.'); hold all
xlabel('Fluorescence');
ylabel('Motion');
%%
figure(5); clf; ax1=[]; cmap=distinguishable_colors(6);
tiledlayout(length(roi)+1,1)
for n=1:length(roi)
ax1=[ax1 nexttile([1 1])];
plot(-tr(:,n)+median(tr(1:1000,n)),'Color',cmap(1,:)); hold all
plot(-tr2(:,n)+median(tr2(1:1000,n)),'Color',cmap(2,:))
plot(-tr3(:,n)+median(tr3(1:1000,n)),'Color',cmap(3,:))
plot(-tr4(:,n)+median(tr4(1:1000,n)),'Color',cmap(4,:))
plot(-tr5(:,n)+median(tr5(1:1000,n)),'Color',cmap(6,:))
plot([1 size(tr,1)],[0 0],'r')
legend({'Raw','HP MC','HP MC+regress out exp.','MC+regress out linear/quadratic','HP MC+regress out linear/quadratic'})
end
ax1=[ax1 nexttile([1 1])];
plot(mcTrace.xymean)
linkaxes(ax1,'x')