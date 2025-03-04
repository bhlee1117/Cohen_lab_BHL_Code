clear
cd '/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/20230805_iGluSnFR3_BU/181643488_iGluSnFR3_100Hz'

load(['output_data.mat']);
load('mcTrace.mat');
sz=double(Device_Data{1, 4}.ROI([2 4]));
mov_mc=double(readBinMov(['mc.bin'],sz(2),sz(1)));

DAQdT=1e-5;

Blue=Device_Data{1, 2}.buffered_tasks(1, 1).channels(1, 2).data;
camTrig=find(Device_Data{1, 2}.Counter_Inputs.data(2:end)-Device_Data{1, 2}.Counter_Inputs.data(1:end-1));
t=camTrig*DAQdT;
nFrames=length(t);
%%
avgImg=mean(mov_mc,3);
imshow2(avgImg, [])

bkg = zeros(2, nFrames);  
bkg(1,:) = linspace(-1, 1, nFrames);  % linear term
bkg(2,:) = linspace(-1, 1, nFrames).^2;  % quadratic term

[mov_res, ~] = SeeResiduals(mov_mc, bkg, 1);  % remove the background and photobleaching
mov_res=SeeResiduals(mov_res,mcTrace);
mov_res=SeeResiduals(mov_res,mcTrace.^2);
clicky(mov_res, avgImg)

croproi=[20, 50, 220, 240];
datDS= mov_res(croproi(2):croproi(2)+croproi(4),croproi(1):croproi(1)+croproi(3),:);
avgImgCrop= avgImg(croproi(2):croproi(2)+croproi(4),croproi(1):croproi(1)+croproi(3));
%datDS = imresize(mov_res(50:end-30,20:280,:), 1, 'bilinear', 'Antialiasing',true);
% datDS = imfilter(dat2, fspecial('average', [3, 3]), 'replicate');
% datDS = datDS(3:end-3,22:289,:);
% datDS = datDS(1:3:end,1:3:end,:);
%datDS = datDS - medfilt1(datDS, 20, [], 3);
clicky(datDS,avgImgCrop)

MovVec = tovec(datDS);
covMat = MovVec*MovVec';
[V, D] = eig(covMat);
D = diag(D);
D = D(end:-1:1);
V = V(:,end:-1:1);
vSign = sign(max(V) - max(-V));  % make the largest value always positive
V = V.*vSign;
eigTraces = V'*MovVec;
figure(8); clf
stackplot(eigTraces(1:30,:)')

nKeep = 30;
eigImgs = zeros(sz(2), sz(1), nKeep);
for j = 1:nKeep;
    eigImgs(:,:,j) = mean(mov_res.*reshape(eigTraces(j,:), [1, 1, nFrames]),3);
end;
figure(9); clf;
for j = 1:nKeep;
    nexttile([1 1]);
    imshow2(eigImgs(3:end-3,3:end-3,j),[]);
    title(num2str(j))
end;


[ics, mixmat, sepmat] = sorted_ica(eigTraces(1:nKeep,blueOn2)',nKeep);
figure(9); clf
stackplot(ics)

eigImgsVec = tovec(eigImgs);
footPrintsVec = eigImgsVec*sepmat';
footPrints = toimg(footPrintsVec, [sz(2), sz(1)]);
figure(10); clf
for j = 1:nKeep;
    nexttile([1 1]);
    imshow2(footPrints(3:end-3,3:end-3,j), []);
    title(num2str(j))
end;

   % ics = pcs*sepmat'

%% Correlation

[roiC1, Ctrace]=clicky(mov_res,avgImg); %Cell of interest
Ctrace_hi=Ctrace-movmedian(Ctrace,400);
mov_res_vec=tovec(mov_res-movmean(mov_res,400,3));
maxlag=1;
for i=1:maxlag
corrImg(:,:,i) = toimg(mov_res_vec(:,i:end-(maxlag-i))*Ctrace_hi(1:end-maxlag+1), sz(2), sz(1));
end
corrImg=toimg(mov_res_vec*Ctrace_hi, sz(2), sz(1));

figure;
corrImg(corrImg<0)=0;
nexttile([1 1]); imshow2(max(corrImg,[],3),[])
nexttile([1 1]); imshow2(imfuse(max(corrImg,[],3),avgImg),[])

%% PCA

MovFilt = imgaussfilt(mov_res, 3);
MovVec=tovec(imresize(MovFilt,1/4));
covMat = MovVec*MovVec';
[V, D] = eig(covMat);
D = diag(D);
D = D(end:-1:1);
V = V(:,end:-1:1);
vSign = sign(max(V) - max(-V));  % make the largest value always positive
V = V.*vSign;

eigImgs = toimg(V,sz(2)/4,sz(1)/4);

for i=1:50; nexttile([1 1]); imshow2(eigImgs(:,:,i),[]); end
stackplot(MovVec'*V(:,1:10))


footprint = mat2gray(toimg(mean(abs(V(:,4)).*D(4)',2), sz(2)/4, sz(1)/4));
submov_res_vec=tovec(imresize(mov_mc,1/4));
ReconMov=toimg(submov_res_vec.*tovec(footprint),sz(2)/4,sz(1)/4);
Recon_corr_C1 = toimg(imresize(ReconMov,4).*corrImg, sz(2), sz(1));

%% STA

C1trace=apply_clicky(roiC1, imresize(ReconMov,4),'no');
C1trace_hi=-(C1trace-movmedian(C1trace,20))';
Spike=find_spike_bh(C1trace_hi/get_threshold(C1trace_hi,1),4.5,4);
plot([1:nFrames],-C1trace,find(Spike),-C1trace(find(Spike)),'r.')

STAmov=zeros(sz(1)*sz(2),36);
wind=[-5:30]; g=1;
for s=find(Spike)
    try
STAmov=STAmov+mov_res_vec(:,s+wind);
g=g+1;
    end
end
STAmov=toimg(STAmov,sz(2),sz(1))/g;
