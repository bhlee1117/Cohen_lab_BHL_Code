clear
clc;
cd '/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/20231123_BHLm112_Prism'
fpath = uigetfile_n_dir;
save_to='/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/20231123_BHLm112_Prism';
%%
for i=1:length(fpath)
load(fullfile(fpath{i},'output.mat'))
nRow = output{2}(1);
nCol = output{2}(2);
[mov1,nframe] = readBinMov(fullfile(fpath{i},'frames.bin'),nCol,nRow);

mov2 = zeros(nRow,nCol,size(mov1,3));

% %Clean signal
for ii = 1:nframe
    mov2(:,:,ii) = medfilt2(mov1(:,:,ii),[3,3]);
end
%rotate and flip to get right direction
% Perform HiLo

    H=hilospeckle(mov2,mov2);%, mov2(:,:,1:2:end));
    HiloMov= H.HiloFinal;
    HiloMov = double(HiloMov);
    
    imwrite(uint16(HiloMov(:,:,1)),fullfile(fpath{i},'HiLoMov.tiff'),'tiff')
    for ii=2:nframe
    imwrite(uint16(HiloMov(:,:,ii)),fullfile(fpath{i},'HiLoMov.tiff'),'tiff','WriteMode','append')
    end
end


%%
save(fullfile(save_to,['SegmentHilo.mat']),'Segment_Otsu')
%% Load fused
load(fullfile(save_to,['SegmentHilo.mat']))
for i=1:100
HiLoImg(:,:,i)=imread("HiLo_Fused.tif",i);
end

%%
bwSeg=bwlabeln(Segment_Otsu);
segments = regionprops3(bwSeg,'Volume','EquivDiameter');
segments = table2array(segments);
%bwlist=find(arrayfun(@(x) x.Volume>4000, segments) & arrayfun(@(x) x.EquivDiameter>10, segments));
bwlist = find(segments(:,1)>7000 & segments(:,2)>30);

se = strel('sphere',2);
Segment_Otsu_di= imdilate(double(ismember(bwSeg,bwlist)),se);
Segment_Otsu_di= imgaussfilt3(Segment_Otsu_di,2);

FiltHiLoImg = double(HiLoImg).* Segment_Otsu_di;
figure(2); clf;
imshow2(imfuse(mat2gray(max(FiltHiLoImg,[],3)),mat2gray(max(HiLoImg,[],3))),[])

    %%
    % imwrite(uint16(FiltHiLoImg(:,:,1)),fullfile(save_to,'FiltHiLoMov.tiff'),'tiff')
    % for ii=2:size(FiltHiLoImg,3)
    % imwrite(uint16(FiltHiLoImg(:,:,ii)),fullfile(save_to,'FiltHiLoMov.tiff'),'tiff','WriteMode','append')
    % end
    writeMov('FiltHiloImg',FiltHiLoImg,[],20,1,[0 15000])