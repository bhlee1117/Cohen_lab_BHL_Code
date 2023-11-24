% Analysis on AAV expression sample and plot, iGluSnFR3 in vivo imaging

clear
cd('/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/20230819_iGluSnFR_Prism_Cannula/203829iGluSnFR3_Cannula_awake_od0')
fpath='/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/20230819_iGluSnFR_Prism_Cannula/203829iGluSnFR3_Cannula_awake_od0';

%%
i=1; edgebound=10;
ROI=[130 100 250 250];

load([fpath '/output_data.mat'])
sz=double(Device_Data{1, 4}.ROI([2 4]));
mov_mc=double(readBinMov([fpath '/mc.bin'],sz(2),sz(1)));
mov_mc=mov_mc(ROI(1):ROI(1)+ROI(3),ROI(2):ROI(2)+ROI(4),:);
mov_mc=mov_mc(edgebound:end-edgebound,edgebound:end-edgebound,1:1000);
avgImg=mean(mov_mc,3);
nFrames=size(mov_mc,3);
bkg = zeros(1, nFrames);
%bkg(1,:) = linspace(-1, 1, nFrames);  % linear term
%bkg(1,:) = linspace(-1, 1, nFrames).^2;  % quadratic term
bkg(1,:) = movmean(squeeze(mean(mov_mc,[1 2])),10);  % quadratic term
mov_mc= SeeResiduals(mov_mc,bkg);
MovFilt=imgaussfilt3(mov_mc,[2 2 1]);
SpineImg=mean(MovFilt(:,:,2:end).*MovFilt(:,:,1:end-1),3);

G = fspecial('gaussian', [6*2 6*2], 2);
SpineImg = imfilter(SpineImg, G, 'replicate');
localMax = imregionalmax(SpineImg);
[ys, xs] = find(localMax);
intensity=SpineImg(find(localMax));

sz_crop=size(SpineImg);
spine_valid=find(intensity>100);
ys=ys(spine_valid); xs=xs(spine_valid);
spine_list=[ys xs intensity(spine_valid)];
plot(spine_list(:,2),spine_list(:,1),'ro')

radius=3;
[x, y] = meshgrid(1:sz_crop(2), 1:sz_crop(2));
for i=1:size(spine_list,1)
circleMask(:,:,i) = (x - spine_list(i,2)).^2 + (y - spine_list(i,1)).^2 <= 2^2;
end

SpineTrace=tovec(MovFilt)'*tovec(circleMask);

spoi=[1:5 8 10:12 14:16 18];
show_footprnt(circleMask(:,:,spoi),avgImg)
SpineTrace=SpineTrace-movmedian(SpineTrace,100,1);
SpineTrace=rescale2(SpineTrace,1);
figure;
plot(SpineTrace(:,spoi)+[1:length(spoi)])


SpineImg=SpineImg-imgaussfilt(SpineImg,5);
fGauss=fspecial('gaussian',11,2)-fspecial('gaussian',11,5);
[ match_data, match_img ] = matchPattern( SpineImg, fGauss, 0.1, 1);



MovFilt=imgaussfilt3(mov_mc,[2 2 1]);
resize_factor=0.5;

MovFilt_resz=imresize(MovFilt,resize_factor);
d=size(MovFilt_resz);
MovVec=tovec(MovFilt_resz);
covMat = MovVec*MovVec';
[V, D] = eig(covMat);   
D = diag(D);
D = D(end:-1:1);
V = V(:,end:-1:1);
vSign = sign(max(V) - max(-V));  % make the largest value always positive
V = V.*vSign;

eigImgs = toimg(V,d(1),d(2));
figure;
for ind=1:10; nexttile([1 1]); imshow2(eigImgs(:,:,ind),[]); title(num2str(ind)); end
figure;
nexttile([1 1]);
for ind=1:10
plot(rescale(MovVec'*V(:,ind))+ind)
hold all
end
%%
pca_ind=[1 2 3];
sz2=size(mov_mc);
footprint = mat2gray(toimg(mean(abs(V(:,pca_ind)).*D(pca_ind)',2), d(1), d(2)));
reconMov=mat2gray(MovFilt.*imresize(footprint,sz2([2 1])));
dF_recon=reconMov-mean(reconMov,3);
moviefixsc(dF_recon,[0 0.1])

