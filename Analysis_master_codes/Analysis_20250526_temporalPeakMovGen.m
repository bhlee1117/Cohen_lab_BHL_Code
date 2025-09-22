nTauPeak=[150 150];
load([fpath{f} '/output_data.mat'])
load(fullfile(fpath{f},'PC_Result.mat'),'Result');
ftp=Ftprnts{f}(:,:,Dist_order{f}(noi_dist{f}));
sz=double(Device_Data{1, 3}.ROI([2 4]));
Mov_PeakTA=double(readBinMov([fpath{f} '/PeakTriggered_movie.bin'],sz(2)*sz(1),301));

rawPattern = double(Device_Data{1, 4}.Target);
Rfixed = imref2d([Device_Data{1, 3}.virtualSensorSize Device_Data{1, 3}.virtualSensorSize]);
inverseTform = invert(Device_Data{1, 4}.tform);
revertedImage = imwarp(rawPattern, inverseTform, 'OutputView', Rfixed);
ROI = double(Device_Data{1, 3}.ROI([1 3 2 4]));
OrangeDMD = imcrop(revertedImage, ROI + [0 0 -1 -1]);

if isfield(Result,'Structure')
[~, expandROI]=get_ROI(im_merge(cat(3,OrangeDMD,Result.Structure),[0 1 1;1 0 0]),[]);
OrangeDMD(max(expandROI,[],3)>0)=1;
else
[~, expandROI]=get_ROI(im_merge(cat(3,OrangeDMD,Result.ref_im),[0 1 1;1 0 0]),[]);
OrangeDMD(max(expandROI,[],3)>0)=1;
end

cellMask=point2img(Result.SWC(:,[1 2]),5,size(Result.ref_im));

[nROI, nTime]=size(Result.traces_bvMask);
for n=[29]%:size(Mov_PeakTA,3)
if 0%range(tovec(Mov_PeakTA(:,:,n)))<1
else
PeakMov=toimg(Mov_PeakTA(:,:,n),sz(2),sz(1));
PeakMov=imgaussfilt(PeakMov,2);
PeakMov=movmean(PeakMov,10,3);
cax_sub=[prctile(PeakMov(:),5) prctile(PeakMov(:),97)];
colorSTA=grs2rgb(tovec(PeakMov),colormap(turbo),cax_sub(1),cax_sub(2));
colorSTA=reshape(colorSTA,sz(2),sz(1),sum(nTauPeak)+1,[]);
colorSTA=permute(colorSTA,[1 2 4 3]);
if isfield(Result,'Structure')
PeakMov_sub_Struc=colorSTA.*mat2gray(Result.Structure)*3.*cellMask;
else
ref_im_hi=Result.ref_im-medfilt2(Result.ref_im,[40 40]);
level = graythresh(ref_im_hi);
ref_im_hi(ref_im_hi<level)=0;
%PeakMov_sub_Struc=colorSTA.*mat2gray(Result.ref_im-100)*3;
PeakMov_sub_Struc=colorSTA.*mat2gray(ref_im_hi)*5.*cellMask;
end
ftp_boundary=bwboundaries(ftp(:,:,n));
figure(161); clf;
writeMov4d(fullfile(fpath{f},['PeakTriggered_movie_' num2str(n,2)]), ...
PeakMov_sub_Struc(:,:,:,:),[-nTauPeak(1):nTauPeak(2)],10,1,[-0.2 1.2],ftp_boundary{1}(:,[2 1]))
end
end