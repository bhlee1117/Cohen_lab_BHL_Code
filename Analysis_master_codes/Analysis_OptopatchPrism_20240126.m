
clear
clc;
cd '/Volumes/BHL_WD18TB/PP72_PlaceCellResults';
[~, ~, raw] = xlsread(['/Volumes/BHL_WD18TB/' ...
    'Prism_OptopatchData_Arrangement.xlsx'], 'Sheet1', 'B5:K86');

save_to='/Volumes/BHL_WD18TB/Arranged_Data/Prism_OptopatchResult';
fpath=raw(:,1);
Mouse=cell2mat(raw(:,2));
NeuronInd=cell2mat(raw(:,5));
CamType=raw(:,3);
StructureData=raw(:,10);
%%
figure;
[~, unqInd] = unique([Mouse NeuronInd] ,'row');
for i=unqInd'
nexttile([1 1])
load(fullfile(fpath{i},"output_data.mat"))
switch char(CamType(i))
    case 'flash'
sz=double(Device_Data{1, 4}.ROI([2 4]));
    case 'fusion'
sz=double(Device_Data{1, 3}.ROI([2 4]));        
end
ref_time=[2000:4000];
mov_test=double(readBinMov_times([fpath{i} '/frames1.bin'],sz(2),sz(1),[ref_time(1):ref_time(end)]));
imshow2(mean(mov_test,3),[])
title(['Mouse #' num2str(Mouse(i)) '-Neuron#' num2str(NeuronInd(i))])
end

%% MC

for i=[74:82]%length(fpath)

load(fullfile(fpath{i},"output_data.mat"))

ref_time=[6000:7000]; overlap=200;
time_segment=25000;

frm_end=max(Device_Data{1, 2}.Counter_Inputs(1, 1).data);
f_seg=[[1:time_segment:frm_end] frm_end+1]; f_seg(2:end)=f_seg(2:end)-1;


switch char(CamType(i))
    case 'flash'
sz=double(Device_Data{1, 4}.ROI([2 4]));
    case 'fusion'
sz=double(Device_Data{1, 3}.ROI([2 4]));        
end

mov_test=double(readBinMov_times([fpath{i} '/frames1.bin'],sz(2),sz(1),[ref_time(1):ref_time(end)]));
mov_test=mov_test(:,:,2:end);
[mov_test_mc,xyField]=optical_flow_motion_correction_LBH(mov_test,mean(mov_test,3),'normcorre');
mov_test=vm(mov_test);
mov_test = single(mov_test)./single(max(mov_test.data(:)));
mov_test = movmean(mov_test,10,3);
mov_ref = squeeze(median(mov_test,3));

for j=1:length(f_seg)-1

    mov=double(readBinMov_times([fpath{i} '/frames1.bin'],sz(2),sz(1),[f_seg(j):f_seg(j+1)]));


    switch char(CamType(i))
    case 'flash'
mov=rollingShutter_correction(mov,Device_Data{1, 4}.exposuretime,'flash');
    case 'fusion'
mov=rollingShutter_correction(mov,Device_Data{1, 3}.exposuretime,'fusion');
    end
   
    mov=vm(mov(:,:,2:end));
    if j==1
        mov=mov(:,:,[1 1:size(mov,3)]);
    end

    [mov_mc,xyField]=optical_flow_motion_correction_LBH(mov,mov_ref,'normcorre');

    ave_im=mean(mov_mc,3);
    mov_mc=vm(mov_mc);
    mov_mc.transpose.savebin([fpath{i} '/mc_ShutterReg' num2str(j,'%02d') '.bin'])

    %        mcTrace = squeeze(mean(xyField,[1 2])); %optic flow
    mcTrace=xyField; % Normcorre
    save([fpath{i} '/mcTrace' num2str(j,'%02d') '.mat'],'mcTrace','ave_im')

    %  clear mov_mc mov
end
sd_mov=std(double(mov),0,3); sd_mov_mc=std(double(mov_mc),0,3);
figure; clf;
nexttile([1 1]); imshow2(sd_mov,[]); title('before mc')
nexttile([1 1]); imshow2(sd_mov_mc,[]); title('after mc')
nexttile([1 1]); imshow2(imfuse(sd_mov,sd_mov_mc),[]); 
title(fpath{i},'Interpreter',  'none')
saveas(gca,[char(fpath{i}) '/' 'MC_result.fig'])
end

%% ROI setting

[a, unqInd] = unique([Mouse NeuronInd] ,'row');
for i=unqInd([2 5 6 7])'
nexttile([1 1])
load(fullfile(fpath{i},"output_data.mat"))
disp(fpath{i})

switch char(CamType(i))
    case 'flash'
sz=double(Device_Data{1, 4}.ROI([2 4]));
    case 'fusion'
sz=double(Device_Data{1, 3}.ROI([2 4]));        
end

ref_time=[2000:4000];
mov_test=double(readBinMov_times([fpath{i} '/mc_ShutterReg' num2str(1,'%02d') '.bin'],sz(2),sz(1),ref_time));
avgImg=mean(mov_test,3);
figure(3); clf;
imshow2(avgImg,[])
g=1; ROIpoly=[];
while g
    h = drawpolygon('Color','r');
    if size(h.Position,1)==1 %no more ROI
        g=0;
    else
        ROIpoly=[ROIpoly; {h.Position}];
    hold all
    plot(h.Position(:,1),h.Position(:,2))
    end
end
close(figure(3));
Result.ROIpoly=ROIpoly;
Result.ref_im=avgImg;
save(fullfile(fpath{i},'Result_20240210.mat'),'Result','-v7.3')

SameCellInd=find(Mouse==Mouse(i) & NeuronInd==NeuronInd(i));
for j=SameCellInd'
Result.ROIpoly=ROIpoly;

mov_mc=double(readBinMov_times([fpath{j} '/mc_ShutterReg' num2str(1,'%02d') '.bin'],sz(2),sz(1),[3000:4000]));
[offsetY offsetX]=calculate_shift(Result.ref_im,mean(mov_mc,3));
ROIpoly_shift=cellfun(@(x) x+[offsetX offsetY],Result.ROIpoly,'UniformOutput',false);
Result.ROIpoly=ROIpoly_shift;
save([fpath{j} '/Result_20240210.mat'],'Result','-v7.3')
end
end

%% Set Footprint
bound=8;
for i=74:length(fpath)
    load([fpath{i} '/Result_20240210.mat'])
    load(fullfile(fpath{i},"output_data.mat"))
switch char(CamType(i))
    case 'flash'
sz=double(Device_Data{1, 4}.ROI([2 4]));
    case 'fusion'
sz=double(Device_Data{1, 3}.ROI([2 4]));        
end
ref_time=[2000:6000];
load(fullfile(fpath{i},['/mcTrace' num2str(1,'%02d') '.mat']));

mov_mc=double(readBinMov_times([fpath{i} '/mc_ShutterReg' num2str(1,'%02d') '.bin'],sz(2),sz(1),ref_time));

Result.ref_im=mean(mov_mc,3);
Result.mc=mcTrace.xymean;
mov_res= mov_mc-mean(mov_mc,3);
%mcTrace.xymean=movmean(mcTrace.xymean,3,2);
% mov_res = SeeResiduals(mov_res,mcTrace.xymean(ref_time,:));
% mov_res = SeeResiduals(mov_res,mcTrace.xymean(ref_time,:).^2);
% mov_res = SeeResiduals(mov_res,mcTrace.xymean(ref_time,1).*mcTrace.xymean(ref_time,2));

n_comp=5;
mov_filt=imgaussfilt3(mov_res,[3 3 0.1]);
movVec=tovec(mov_filt);
Npoly=size(Result.ROIpoly,1);
ftprnt = zeros(size(mov_filt,1)*size(mov_filt,2),Npoly);
clear mask
for p=1:Npoly %each ROIs
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
    coeff = mat2gray(mean(abs(V(:,1:n_comp)).*D(1:n_comp)',2));
    ftprnt(pixelList,p)=coeff;
end

Result.ftprnt=toimg(ftprnt,sz(2),sz(1));
Result.ftprnt(1:bound,:,:)=0;
Result.ftprnt(:,1:bound,:)=0;
Result.ftprnt(end-bound:end,:,:)=0;
Result.ftprnt(:,end-bound:end,:)=0;
figure; clf;
show_im=Result.ref_im;
show_im(show_im<15)=median(show_im(:));
show_im=show_im-imgaussfilt(show_im,15);
show_im(show_im<prctile(show_im(:),20))=prctile(show_im(:),20);
show_footprnt(Result.ftprnt,show_im)
save([fpath{i} '/Result_20240210.mat'],'Result','-v7.3')
saveas(gca,[char(fpath{i}) '/' 'ftprnt.fig'])
end


%% Signal extraction
bound=5;

for i=[74:length(fpath)]

    load([fpath{i} '/Result_20240210.mat'])
    load(fullfile(fpath{i},"output_data.mat"))
    load([fpath{i} '/mcTrace' num2str(1,'%02d') '.mat']);
switch char(CamType(i))
    case 'flash'
sz=double(Device_Data{1, 4}.ROI([2 4]));
    case 'fusion'
sz=double(Device_Data{1, 3}.ROI([2 4]));        
end

Result.traces=[];
Result.mc=mcTrace.xymean;
Result.im_corr=[];

ref_im_vec=tovec(Result.ref_im(bound:end-bound,bound:end-bound));

CamCounter=Device_Data{1, 2}.Counter_Inputs(1, 1).data;
CamTrigger=find(CamCounter(2:end)-CamCounter(1:end-1));

Result.Blue=Device_Data{1, 2}.buffered_tasks(1, 1).channels(1, 2).data;
Result.Blue=Result.Blue(CamTrigger);
Rfixed = imref2d(Device_Data{1, 3}.ROI([4 2]));
inverseTform = invert(Device_Data{1, 6}.tform);
revertedImage = imwarp(double(Device_Data{1, 6}.Target), inverseTform,'OutputView',Rfixed);
[blueDMDimg blueDMDcontour]=imcrop(revertedImage,double(Device_Data{1, 3}.ROI([1 3 2 4]))+[0 0 -1 -1]);
Result.BlueDMD=blueDMDcontour;
Result.BlueDMDimg=blueDMDimg;


mov_mc=double(readBinMov_times([fpath{i} '/mc_ShutterReg' num2str(1,'%02d') '.bin'],sz(2),sz(1),[1:length(CamTrigger)]));

mov_mc_vec=tovec(mov_mc(bound:end-bound,bound:end-bound,:));
mov_mc_vec=(mov_mc_vec-mean(mov_mc_vec,1))./std(mov_mc_vec,0,1);

mov_res= mov_mc-mean(mov_mc,3);
bkg = zeros(1, size(mov_mc,3));

%     bkg(1,:) = linspace(-1, 1, size(mov_mc,3));  % linear term
%     bkg(2,:) = linspace(-1, 1, size(mov_mc,3)).^2;  % quadratic term
bkg(1,:)=movmedian(get_blueoffTrace(squeeze(mean(mov_mc,[1 2])),[Result.Blue],30),3000,'omitnan');
mov_res = SeeResiduals(mov_res,Result.mc);
mov_res = SeeResiduals(mov_res,Result.mc.^2);
mov_res = SeeResiduals(mov_res,Result.mc(:,1).*Result.mc(:,end));
mov_res= SeeResiduals(mov_res,bkg,1);


    Result.traces=[-(tovec(mov_res)'*tovec(Result.ftprnt))'];
    Result.im_corr=[sum(mov_mc_vec.*ref_im_vec,1)/(size(mov_mc_vec,1)-1)];  %image correlation
save([fpath{i} '/Result_20240210.mat'],'Result','-v7.3')
end

%% Collect all data
%OP_Result=[];

for i=1:length(fpath)
 files = dir([fpath{i} '/Result*.mat']);
 if isempty(files)
    error('No Result files');
 end

 load(fullfile(fpath{i},files(end).name),'Result')
OP_Result{i}=Result;
OP_Result{i}.mouse = Mouse(i);
OP_Result{i}.Neuron = NeuronInd(i);
end

save([save_to 'OP_Result_20240212'],"OP_Result",'fpath','-v7.3')

%% STAs
nTau=[-10:20];

for i=1:length(OP_Result)
OP_Result{i}.normTrace=OP_Result{i}.traces./get_threshold(OP_Result{i}.traces,1);
OP_Result{i}.spike=find_spike_bh(OP_Result{i}.normTrace-movmedian(OP_Result{i}.normTrace,300,2),5,3);

Blue=OP_Result{i}.Blue;
blueOff = Blue == 0;
blueOff2 = imerode(blueOff, [ones(1,20), zeros(1, 20)]);
Blue_di=~blueOff2;
bwBlue_di=bwlabel(Blue_di);

ref_ROI=find(sum(OP_Result{i}.spike(1:5,:).*bwBlue_di,2)==max(sum(OP_Result{i}.spike(1:5,:).*bwBlue_di,2)),1);
nROI=size(OP_Result{i}.normTrace,1);
tr=OP_Result{i}.normTrace(ref_ROI,:); t=[1:length(tr)];
spike=OP_Result{i}.spike(ref_ROI,:);   


sp_pulse=[];
for b=1:max(bwBlue_di)
t_tmp=find(bwBlue_di==b);
if max(bwBlue_di)<15
sp_pulse=[sp_pulse find(spike(t_tmp))+t_tmp(1)-1];
else
sp_pulse=[sp_pulse find(spike(t_tmp),1,'first')+t_tmp(1)-1];
end
end

sTau=sp_pulse(1:end-1)'+nTau;
spikeMat=reshape(OP_Result{i}.normTrace(:,sTau),nROI,length(sp_pulse)-1,[]);
STAtrace=squeeze(mean(spikeMat,2));
STAtrace=STAtrace-prctile(STAtrace, 5, 2);
OP_Result{i}.spikeMat = spikeMat;
OP_Result{i}.STAtrace = STAtrace;

if max(bwBlue_di)>15

load(fullfile(fpath{i},"output_data.mat"))
load([fpath{i} '/mcTrace' num2str(1,'%02d') '.mat']);
switch char(CamType(i))
    case 'flash'
sz=double(Device_Data{1, 4}.ROI([2 4]));
    case 'fusion'
sz=double(Device_Data{1, 3}.ROI([2 4]));        
end

CamCounter=Device_Data{1, 2}.Counter_Inputs(1, 1).data;
CamTrigger=find(CamCounter(2:end)-CamCounter(1:end-1));

mov_mc=double(readBinMov_times([fpath{i} '/mc_ShutterReg' num2str(1,'%02d') '.bin'],sz(2),sz(1),[1:length(CamTrigger)]));

mov_res= mov_mc-mean(mov_mc,3);
bkg = zeros(1, size(mov_mc,3));

bkg(1,:)=movmedian(get_blueoffTrace(squeeze(mean(mov_mc,[1 2])),Blue,30),3000,'omitnan');
mov_res = SeeResiduals(mov_res,OP_Result{i}.mc);
mov_res = SeeResiduals(mov_res,OP_Result{i}.mc.^2);
mov_res = SeeResiduals(mov_res,OP_Result{i}.mc(:,1).*OP_Result{i}.mc(:,end));
mov_res= SeeResiduals(mov_res,bkg,1);

OP_Result{i}.STAmovie=squeeze(mean(reshape(mov_res(:,:,sTau),sz(2),sz(1),[],length(nTau)),3));
end
end

save([save_to '/OP_Result_20240212'],"OP_Result",'fpath','-v7.3')


%% Structure Segmentation
Struct_valid=find(1-cell2mat(cellfun(@(x) sum(isnan(x)), StructureData, 'UniformOutput', false)));
[a, unqInd] = unique([Mouse(Struct_valid) NeuronInd(Struct_valid)] ,'row');

for i=Struct_valid(unqInd(9))'
StructureStack=mat2gray(double(tiffreadVolume(StructureData{i})));
StructureStack(StructureStack==0)=median(StructureStack(:));
%StructureStack_med=medfilt2_mov(StructureStack,[15 15]);
StructureStack_Gauss=imgaussfilt3(StructureStack,[7 7 0.1]);
%StructureStack_med(StructureStack_med==0)=median(StructureStack_med(:));
%StructureStack=(StructureStack-StructureStack_med)./StructureStack_med;
StructureStack_filt=(StructureStack-StructureStack_Gauss);
StructureStack_filt=mat2gray(StructureStack_filt);
StructureStack_bin=[]; level=[];
level = graythresh(StructureStack_filt);
StructureStack_bin=StructureStack_filt>level*1.12;
moviefixsc(StructureStack_bin)

se = strel('sphere', 1);
StructureStack_bin = imdilate(StructureStack_bin, se);
bwSeg=bwlabeln(StructureStack_bin);
segments = regionprops3(bwSeg,'Volume','EquivDiameter');
segments = table2array(segments);
%bwlist=find(arrayfun(@(x) x.Volume>4000, segments) & arrayfun(@(x) x.EquivDiameter>10, segments));
bwlist = find(segments(:,1)>7000 & segments(:,2)>10);

se = strel('sphere',1);
dendrite_bin=double(ismember(bwSeg,bwlist));
dendrite_bin= imdilate(dendrite_bin,se);
dendrite_bin= imgaussfilt3(dendrite_bin,2);

figure(3); clf;
imshow2(max(dendrite_bin,[],3),[])
g=1; ROIrmv=[];
while g
    h = drawpolygon('Color','r');
    if size(h.Position,1)==1 %no more ROI
        g=0;
    else
        ROIrmv=[ROIrmv; {h.Position}];
    hold all
    plot(h.Position(:,1),h.Position(:,2))
    end
end
ROIrmvmask=roi2mask(ROIrmv,size(dendrite_bin,1),size(dendrite_bin,2));
close(figure(3));
dendrite_bin(repmat(ROIrmvmask,1,1,size(dendrite_bin,3)))=0;
figure(3); clf;
imshow2(max(dendrite_bin,[],3),[])

StructureStack_final = double(StructureStack).* dendrite_bin;
figure(2); clf;
imshow2(imfuse(mat2gray(max(StructureStack_final,[],3)),mat2gray(max(StructureStack,[],3))),[])


rot_ang=0;
Structure_ref=(imrotate(StructureStack_final,rot_ang));
ref_img=OP_Result{i}.ref_im; ref_img(ref_img<prctile(ref_img(:),20))=median(ref_img(:)); ref_img=ref_img-prctile(ref_img(:),20);
[RegImg,tformReg]=imReg(ref_img,max(Structure_ref,[],3));
saveastiff(uint16(mat2gray(Structure_ref)*255), [StructureData{i} '_Structure.tiff']);

SameCellInd=find(Mouse==Mouse(i) & NeuronInd==NeuronInd(i));
for j=SameCellInd'
    OP_Result{j}.Structure=max(Structure_ref,[],3);
    OP_Result{j}.Structure_bin=max(imrotate(dendrite_bin,rot_ang),[],3);
    OP_Result{j}.tform=tformReg;
end
end


%%
save([save_to '/OP_Result_20240212'],"OP_Result",'fpath','-v7.3')
%%


%%
for i=[31:33]
mask=max(OP_Result{i}.Structure_bin,[],3)>0.01;
maskSTA=max(-OP_Result{i}.STAmovie,[],3)./OP_Result{i}.ref_im>0.05;
StrImg=max(OP_Result{i}.Structure,[],3);
STAmovie=mat2gray(-OP_Result{i}.STAmovie);
STAmovie=STAmovie-prctile(STAmovie,10,3);
STAmovie=mat2gray(STAmovie(:,:,6:25));

[OP_Result{i}.SNAPT OP_Result{i}.dtimg]=generate_SNAPTmov(mat2gray(STAmovie),mask,StrImg,tformReg);
[yR xR zR]=size(OP_Result{i}.Structure);
bluePatt = bwboundaries(imwarp(OP_Result{i}.BlueDMDimg,OP_Result{i}.tform,'OutputView', imref2d([yR xR])));
figure(20); clf;
v = VideoWriter([fpath{i} '/SNAPT_movie'],'MPEG-4');

open(v);
subframeT = 0.025; % ms
initialT = -2; % ms
finalT = 2; % ms
times = initialT:subframeT:finalT;

for j = 1:length(times)
    imshow(OP_Result{i}.SNAPT(:,:,:,j),[])
    pbaspect([size(double(OP_Result{i}.SNAPT(:,:,:,j)),2) size(double(OP_Result{i}.SNAPT(:,:,:,j)),1) 1]),colormap(gray)
    hold all
    plot(bluePatt{1}(:,2),bluePatt{1}(:,1),'color',[0 0.6 1],'linewidth',2)
    axis off
    text(2,20,[num2str(times(j)+1.7) 'ms'], 'FontSize', 20, 'color', [0.99 0.99 0.99])% the value 1. is to adjust timing by eyes       
     pause(0.1)
    set(gcf,'color','w')    % Sets background to white
    frame = getframe(gcf);
    writeVideo(v,frame);
end;
close(v);
end
save([save_to 'OP_Result_20240127'],"OP_Result",'fpath','-v7.3')

%% realign by footprint's location from soma (1st footprint)
for f=1:length(OP_Result)
centroid_ftprnt = get_coord(OP_Result{f}.ftprnt);
dist_centroid = distance_mat(centroid_ftprnt(1,:),centroid_ftprnt);
[~, OP_Result{f}.dist_order]=sort(dist_centroid,'ascend');
end
%% Show dSpikes

for i=1:length(OP_Result)

    tr=OP_Result{i}.normTrace;
    sp=OP_Result{i}.spike;
    ref_ROI=find(sum(sp(1:5,:),2)==max(sum(sp(1:5,:),2)),1);
    sp_timeROI=find(sp(ref_ROI,:));
    rejectSP=find(min(sp_timeROI'+[-2:0])<1 | max(sp_timeROI'+[-2:0])>size(tr,2));
    sp_timeROI(rejectSP)=[];
    [~, shiftSoma]=max(reshape(tr(1,sp_timeROI'+[-2:0]),length(sp_timeROI),[]),[],2);
    sp_timeSoma=unique([sp_timeROI-(3-shiftSoma') find(sp(1,:))]);
    sp_dSp=sp;
    sp_dSp(:,sp_timeSoma'+[0:3])=0;


t_ref=[6600:6650];
avgImg=mean(mov_mc(bound:end-bound,bound:end-bound,:),3);
mask=(avgImg-medfilt2(avgImg,[15 15]))>80;
mov_seg=mov_res(bound:end-bound,bound:end-bound,t_ref);
mov_seg=mov_seg-mean(mov_seg,3);
movSegfilt = spatialfilt(mov_seg, 7, 3);
movSegfilt=movSegfilt-movmedian(movSegfilt,10,3);
[movSegfiltPCA, eigVecs, eigVals] = pcafilt(movSegfilt, 20);
eigImgs = toimg(tovec(mov_seg)*eigVecs(:,1:20), size(mov_seg,1), size(mov_seg,2));
figure(13); clf
for j = 1:20;
    nexttile([1 1])
    imshow2(eigImgs(:,:,j), []);
end;
clf;
for j=1:14
trace=rescale(OP_Result{i}.normTrace(OP_Result{i}.dist_order(j),t_ref));    
trace=trace-movmedian(trace,10);
corrImg=-toimg(corr(tovec(movSegfilt)',trace'),size(movSegfilt,1),size(movSegfilt,2));
corrImg=imgaussfilt(corrImg,3).*mask;
nexttile([1 1]);
%imagesc(im_merge(cat(3,avgImg,mat2gray(corrImg)/10),[0.5 0.5 0.5; 0.5 0 0]))
imshow2(corrImg,[0 0.6])
title(num2str(j))
end
colormap(turbo)

end
    




%%
i=38; bound=5;

    load(fullfile(fpath{i},"output_data.mat"))
    load([fpath{i} '/mcTrace' num2str(1,'%02d') '.mat']);
switch char(CamType(i))
    case 'flash'
sz=double(Device_Data{1, 4}.ROI([2 4]));
    case 'fusion'
sz=double(Device_Data{1, 3}.ROI([2 4]));        
end

CamCounter=Device_Data{1, 2}.Counter_Inputs(1, 1).data;
CamTrigger=find(CamCounter(2:end)-CamCounter(1:end-1));
mov_mc=double(readBinMov_times([fpath{i} '/mc_ShutterReg' num2str(1,'%02d') '.bin'],sz(2),sz(1),[1:length(CamTrigger)]));

mov_res= mov_mc-mean(mov_mc,3);
bkg = zeros(1, size(mov_mc,3));
bkg(1,:)=movmedian(get_blueoffTrace(squeeze(mean(mov_mc,[1 2])),[OP_Result{i}.Blue],30),3000,'omitnan');
mov_res = SeeResiduals(mov_res,OP_Result{i}.mc);
mov_res = SeeResiduals(mov_res,OP_Result{i}.mc.^2);
mov_res = SeeResiduals(mov_res,OP_Result{i}.mc(:,1).*OP_Result{i}.mc(:,end));
mov_res= SeeResiduals(mov_res,bkg,1);

t_ref=[6600:6650];
avgImg=mean(mov_mc(bound:end-bound,bound:end-bound,:),3);
mask=(avgImg-medfilt2(avgImg,[15 15]))>80;
mov_seg=mov_res(bound:end-bound,bound:end-bound,t_ref);
mov_seg=mov_seg-mean(mov_seg,3);
movSegfilt = spatialfilt(mov_seg, 7, 3);
movSegfilt=movSegfilt-movmedian(movSegfilt,10,3);
[movSegfiltPCA, eigVecs, eigVals] = pcafilt(movSegfilt, 20);
eigImgs = toimg(tovec(mov_seg)*eigVecs(:,1:20), size(mov_seg,1), size(mov_seg,2));
figure(13); clf
for j = 1:20;
    nexttile([1 1])
    imshow2(eigImgs(:,:,j), []);
end;
clf;
for j=1:14
trace=rescale(OP_Result{i}.normTrace(OP_Result{i}.dist_order(j),t_ref));    
trace=trace-movmedian(trace,10);
corrImg=-toimg(corr(tovec(movSegfilt)',trace'),size(movSegfilt,1),size(movSegfilt,2));
corrImg=imgaussfilt(corrImg,3).*mask;
nexttile([1 1]);
%imagesc(im_merge(cat(3,avgImg,mat2gray(corrImg)/10),[0.5 0.5 0.5; 0.5 0 0]))
imshow2(corrImg,[0 0.6])
title(num2str(j))
end
colormap(turbo)
