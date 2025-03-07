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
save_figto='/Volumes/BHL_WD18TB/PP72_PlaceCellResults';
place_bin=150; time_segment=15000; overlap=200;
alignedMovFN = {'STA_Mat_SS','STA_Mat_CS','STA_Mat_dSP'};

%%

f=[18];
disp(fpath{f});
DAQ_rate=0.000005;
load([fpath{f} '/output_data.mat'])
load(fullfile(fpath{f},'PC_Result.mat'),'Result')
sz=double(Device_Data{1, 3}.ROI([2 4]));

load(fullfile(fpath{f},'mcTrace15.mat'));
frm_end=max(Device_Data{1, 2}.Counter_Inputs(1, 1).data);
CamCounter=Device_Data{1, 2}.Counter_Inputs(1, 1).data;
CamTrigger=find(CamCounter(2:end)-CamCounter(1:end-1));
frm_rate=double((CamTrigger(2)-CamTrigger(1))*DAQ_rate);
mov_mc=double(readBinMov([fpath{f} '/mc_ShutterReg' num2str(15,'%02d') '.bin'],sz(2),sz(1)));
Result.ref_im=mean(mov_mc,3);

mov_res= mov_mc-mean(mov_mc,3);
mcTrace.xymean=movmean(mcTrace.xymean,3,2);
mov_res = SeeResiduals(mov_res,mcTrace.xymean);
mov_res = SeeResiduals(mov_res,mcTrace.xymean.^2);
mov_res = SeeResiduals(mov_res,mcTrace.xymean(:,1).*mcTrace.xymean(:,3));

mov_filt=imgaussfilt3(mov_res,[3 3 0.1]);
movVec=tovec(mov_filt);
%%
n_comp=5;
Npoly=size(Result.ROIpoly,1);
ftprnt = zeros(size(mov_filt,1)*size(mov_filt,2),Npoly);
clear mask

for p=2 %each ROIs
    mask(:,:,p) = poly2mask(Result.ROIpoly{p}(:,1), Result.ROIpoly{p}(:,2), sz(2), sz(1));
    pixelList=find(tovec(squeeze(mask(:,:,p))));
    subMov = movVec(pixelList,:);
    covMat = subMov*subMov';  % PCA within each region
    [V, D] = eig(covMat);
    D = diag(D);
    D = D(end:-1:1);
    V = V(:,end:-1:1);
    vSign = sign(max(V) - max(-V));  % make the largest value always positive
    V = V.*vSign;
    eigTrace=subMov'*V(:,1:n_comp);
    figure(3); clf;
    nexttile([1 1])
    stackplot(eigTrace);
    nexttile([1 1])
    plot(D(1:5)./sum(D));

    figure(4); clf;
    for n=1:n_comp
        eigImg=zeros(size(mov_filt,1)*size(mov_filt,2),1);
        nexttile([1 1])
        eigImg(pixelList,1)=V(:,n);
        imshow2(toimg(eigImg,size(mov_filt,1),size(mov_filt,2)),[])
    end

    % coeff = mat2gray(mean(abs(V(:,1:n_comp)).*D(1:n_comp)',2));
    % ftprnt(pixelList,p)=coeff;
end


%%
for p=1:size(Result.ftprnt,3)
    Mov_roi{p}=[];
end
tic;
f=[18];%length(fpath)
load(fullfile(fpath{f},'PC_Result.mat'),'Result');
load([fpath{f} '/output_data.mat'])
sz=double(Device_Data{1, 3}.ROI([2 4])); blueDMDcontour=[];
CamCounter=Device_Data{1, 2}.Counter_Inputs(1, 1).data;
CamTrigger=find(CamCounter(2:end)-CamCounter(1:end-1));

frm_end=EndFrame(f);
f_seg=[[1:time_segment:frm_end] frm_end+1]; f_seg(2:end)=f_seg(2:end)-1;
take_window=repmat([1 time_segment],length(f_seg)-1,1);
take_window(2:end,1)=take_window(2:end,1)+overlap; take_window(1:end-1,2)=take_window(1:end-1,2)+overlap;
take_window(end)=mod(f_seg(end),time_segment);
Blue_on_Seg=unique(ceil(find(Result.Blue)/time_segment));

for j=1:3%length(f_seg)-1
    j
    mov_mc=double(readBinMov([fpath{f} '/mc_ShutterReg' num2str(j,'%02d') '.bin'],sz(2),sz(1)));
    load([fpath{f} '/mcTrace' num2str(j,'%02d') '.mat']);

    mov_mc=mov_mc(:,:,[take_window(j,1):take_window(j,2)]);
    %mc= movmean(mcTrace.xymean-movmedian(mcTrace.xymean,500,1),3,1);
    mc=mcTrace.xymean([take_window(j,1):take_window(j,2)],:);

    mov_mc_filt=imgaussfilt3(mov_mc,[2 2 0.1]);

    if ismember(j,Blue_on_Seg)
        [~, blueomitTr]=get_blueoffTrace(squeeze(mean(mov_mc,[1 2])),Result.Blue((j-1)*time_segment+1+take_window(j,1):j*time_segment+take_window(j,1)),30);
        bkg = zeros(1, size(mov_mc,3));
        bkg(1,:)=movmedian(blueomitTr,4000,'omitnan');
    else
        bkg = zeros(2, size(mov_mc,3));
        bkg(1,:) = linspace(-1, 1, size(mov_mc,3));  % linear term
        bkg(2,:) = linspace(-1, 1, size(mov_mc,3)).^2;  % quadratic term
    end

    mov_res= mov_mc-mean(mov_mc,3);
    mc=movmean(mc,3,2);
    mov_res = SeeResiduals(mov_res,mc);
    mov_res = SeeResiduals(mov_res,mc.^2);
    mov_res = SeeResiduals(mov_res,mc(:,1).*mc(:,3));

    mov_filt=imgaussfilt3(mov_res,[2 2 0.1]);
    movVec=tovec(mov_filt);

    for p=1:size(Result.ftprnt,3)
        mask(:,:,p) = poly2mask(Result.ROIpoly{p}(:,1), Result.ROIpoly{p}(:,2), sz(2), sz(1));
        pixelList=find(tovec(squeeze(mask(:,:,p))));
        Mov_roi{p} = [Mov_roi{p} movVec(pixelList,:)];
    end
end
toc;


%%
p=12; n_comp=10;
mask(:,:,p) = poly2mask(Result.ROIpoly{p}(:,1), Result.ROIpoly{p}(:,2), sz(2), sz(1));
pixelList=find(tovec(squeeze(mask(:,:,p))));
covMat =(Mov_roi{p}-mean(Mov_roi{p},2))*(Mov_roi{p}-mean(Mov_roi{p},2))';
[V, D] = eig(covMat);
D = diag(D);
D = D(end:-1:1);
V = V(:,end:-1:1);
vSign = sign(max(V) - max(-V));  % make the largest value always positive
V = V.*vSign;
eigTrace=Mov_roi{p}'*V(:,1:n_comp);
figure(3); clf;
nexttile([1 1])
%plot(rescale2(eigTrace,1)+[1:n_comp]);
stackplot(eigTrace)
nexttile([1 1])
plot(D(1:n_comp)./sum(D));
xlabel('# components')

figure(4); clf; ax2=[];
for n=1:n_comp
    eigImg=NaN(size(mov_filt,1)*size(mov_filt,2),1);
    ax2=[ax2 nexttile([1 1])];
    eigImg(pixelList,1)=V(:,n);
    eigImg=toimg(eigImg,size(mov_filt,1),size(mov_filt,2));
    imshow2(im_merge(cat(3,Result.ref_im,eigImg),[1 1 1;1 0 0]),[])
    title([num2str(n) ', Fraction: ' num2str(D(n)/sum(D),2)])
end
linkaxes(ax2,'xy')
% reconMov=zeros(sz(1)*sz(2),size(Mov_roi{p},2));
% reconMov(pixelList,:)=(eigTrace(:,[1 2])*V(:,[1 2])')';
% reconMov=reshape(reconMov,sz(2),sz(1),[]);
coeff=(Mov_roi{p}-mean(Mov_roi{p},2))*mean(eigTrace(:,[1 2])*V(:,[1 2])',2);
coeffImg=zeros(size(mov_filt,1)*size(mov_filt,2),1);
coeffImg(pixelList,1)=coeff;
nexttile([1 1])
coeffImg=toimg(coeffImg,size(mov_filt,1),size(mov_filt,2));
    imshow2(im_merge(cat(3,Result.ref_im,coeffImg),[1 1 1;1 0 0]),[])

