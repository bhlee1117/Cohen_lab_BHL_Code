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
load(fullfile(fpath{20},'SS_STAmovie.mat'))
nTau={[-30:20],[-50:100],[-30:20]}; %SS, CS, dSP
for f=[20]
    load(fullfile(fpath{f},'PC_Result.mat'))
   
    STAmovie=mat2gray(-SS_STA);
    STAmovie=STAmovie-prctile(STAmovie,10,3);
    STAmovie=mat2gray(STAmovie(:,:,-nTau{1}(1)-5:-nTau{1}+10));

    SkelDend = Skeletonize_dendrite(max(STAmovie,[],3),10,0.1,10);
    ref_im=imgaussfilt(max(STAmovie,[],3)-imgaussfilt(max(STAmovie,[],3),20),0.5);
    mask=SkelDend;%max(Result.Structure_bin,[],3)>0.01;
    thres=prctile(ref_im(~SkelDend),20);
    ref_im(ref_im<thres)=0;
    StrImg=ref_im.*SkelDend;
   imshow2(StrImg,[])

    [Result.SNAPT Result.dtimg]=generate_SNAPTmov(mat2gray(STAmovie),mask,StrImg);

figure(20); clf;
v = VideoWriter([fpath{f} '/SNAPT_movie'],'MPEG-4');

open(v);
subframeT = 0.025; % ms
initialT = -2; % ms
finalT = 2; % ms
times = initialT:subframeT:finalT;

for j = 1:length(times)
    clf;
    set(gca,'units','pixels','position',[50 50 800 400])
    imshow(Result.SNAPT(:,:,:,j),[])
    pbaspect([size(double(Result.SNAPT(:,:,:,j)),2) size(double(Result.SNAPT(:,:,:,j)),1) 1]),colormap(gray)
    axis off
    text(2,15,[num2str(times(j)+0.9,'%2.3f') ' ms'], 'FontSize', 20, 'color', [0.99 0.99 0.99])% the value 1. is to adjust timing by eyes
    pause(0.1)
    set(gcf,'color','w')    % Sets background to white
    frame = getframe(gcf);
    writeVideo(v,frame);
    pause(0.1);
end;
close(v);

end
