
clear
clc;
cd '/Volumes/BHL18TB_D2/Arranged_Data/Prism_OptopatchResult';
[~, ~, raw] = xlsread(['/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/' ...
    'Prism_OptopatchData_Arrangement.xlsx'], 'Sheet1', 'B5:K171');

save_to='/Volumes/BHL18TB_D2/Arranged_Data/Prism_OptopatchResult';
fpath=raw(:,1);
Mouse=cell2mat(raw(:,2));
NeuronInd=cell2mat(raw(:,5));
CamType=raw(:,3);
StructureData=raw(:,10);

%% MC

for i=[161:166]%length(fpath)

    load(fullfile(fpath{i},"output_data.mat"))

    ref_time=[6000:7000]; overlap=200;
    time_segment=25000;

    frm_end=max(Device_Data{1, 2}.Counter_Inputs(1, 1).data);
    f_seg=[[1:time_segment:frm_end] frm_end+1]; f_seg(2:end)=f_seg(2:end)-1;
    CamDAQ_rate=Device_Data{1, 2}.Counter_Inputs.rate;
    CamTrig=Device_Data{1, 2}.Counter_Inputs.data;
    CamTrig2=find(CamTrig(2:end)-CamTrig(1:end-1)>0);
    Frm_rate=(CamTrig2(2)-CamTrig2(1))/CamDAQ_rate;
    Blue=Device_Data{1, 2}.buffered_tasks(1, 1).channels(1, 2).data;
    Result.Blue=Blue(CamTrig2);

    for cam=1:2
    sz=double(Device_Data{1, cam+2}.ROI([2 4]));
    mov_test=double(readBinMov_times([fpath{i} '/frames' num2str(cam) '.bin'],sz(2),sz(1),[ref_time(1):ref_time(end)]));
    mov_test=mov_test(:,:,2:end);
    [mov_test_mc,xyField]=optical_flow_motion_correction_LBH(mov_test,mean(mov_test,3),'normcorre');
    mov_test=vm(mov_test);
    mov_test = single(mov_test)./single(max(mov_test.data(:)));
    mov_test = movmean(mov_test,10,3);
    mov_ref = squeeze(median(mov_test,3));

    for j=1:length(f_seg)-1

        mov=double(readBinMov_times([fpath{i} '/frames' num2str(cam) '.bin'],sz(2),sz(1),[f_seg(j):f_seg(j+1)]));
        mov=rollingShutter_correction(mov,1/Frm_rate,'fusion');
        mov=vm(mov(:,:,2:end));

        if j==1
            mov=mov(:,:,[1 1:size(mov,3)]);
        end

        [mov_mc,xyField]=optical_flow_motion_correction_LBH(mov,mov_ref,'normcorre');

        ave_im=mean(mov_mc,3);
        mov_mc=vm(mov_mc);
        mov_mc.transpose.savebin([fpath{i} '/mc_ShutterReg_cam' num2str(cam) '_' num2str(j,'%02d') '.bin'])

        %        mcTrace = squeeze(mean(xyField,[1 2])); %optic flow
        mcTrace=xyField; % Normcorre
        save([fpath{i} '/mcTrace_cam' num2str(cam) '_' num2str(j,'%02d') '.mat'],'mcTrace','ave_im')
    end
    end
    sd_mov=std(double(mov),0,3); sd_mov_mc=std(double(mov_mc),0,3);
    figure; clf;
    nexttile([1 1]); imshow2(sd_mov,[]); title('before mc')
    nexttile([1 1]); imshow2(sd_mov_mc,[]); title('after mc')
    nexttile([1 1]); imshow2(imfuse(sd_mov,sd_mov_mc),[]);
    title(fpath{i},'Interpreter',  'none')
    saveas(gca,[char(fpath{i}) '/' 'MC_result.fig'])
end

%% load movies
i=[166]; clear mov_mc
 load(fullfile(fpath{i},"output_data.mat"))

    CamDAQ_rate=Device_Data{1, 2}.Counter_Inputs.rate;
    CamTrig=Device_Data{1, 2}.Counter_Inputs.data;
    CamTrig2=find(CamTrig(2:end)-CamTrig(1:end-1)>0);
    Frm_rate=(CamTrig2(2)-CamTrig2(1))/CamDAQ_rate;
    Blue=Device_Data{1, 2}.buffered_tasks(1, 1).channels(1, 2).data;
    Blue=Blue(CamTrig2);

    for cam=1:2
    sz=double(Device_Data{1, cam+2}.ROI([2 4]));
    mov_mc{cam}=double(readBinMov([fpath{i} '/mc_ShutterReg_cam' num2str(cam) '_' num2str(1,'%02d') '.bin'],sz(2),sz(1)));
    end

    Result.ref_im=cellfun(@(x) mean(x,3),mov_mc,'UniformOutput',false);
    [~, bluedmd]=get_blueDMDPatt(Device_Data,'normal');

    figure;
imshow2(mean(mov_mc{1},3),[])
hold all
plot(bluedmd{1}(:,2),bluedmd{1}(:,1),'color',[0 0.6 1])
%% Clicky
clear STAmovie mov_res
[~, trOrg]=clicky(mov_mc{1});
tr=-(trOrg-movmedian(trOrg,50))';
trOrg=-trOrg';
tr=tr./get_threshold(tr,1);
sp=find_spike_bh(tr(1,:),5,3);
bwblue=bwlabel(Blue);
sp_time=[];
for b=1:max(bwblue)
    spt=find(bwblue==b & sp,1,'first');
    if ~isempty(spt)
        sp_time=[sp_time spt];
    end
end
STAtrace=mean(tr(sp_time(2:end)'+[-30:20]),1,'omitnan');
for cam=1:2
    mov_res{cam}=mov_mc{cam}-mean(mov_mc{cam},3);
    S=reshape(mov_res{cam}(:,:,sp_time'+[-30:20]),size(mov_mc{cam},1),size(mov_mc{cam},2),[],51);
    STAmovie{cam}=squeeze(mean(S,3));
    STAmovie{cam}=STAmovie{cam}-mean(STAmovie{cam}(:,:,1:10),3);
end
showmov=-[STAmovie{1} flipud(fliplr(STAmovie{2}))];
figure;
moviefixsc(showmov);
figure;
STAtr=squeeze(mean(reshape(trOrg(:,sp_time'+[-30:20]),2,[],51),2,'omitnan'))';
plot(STAtr-mean(STAtr(1:10,:)))

%% Structure registration.
[~, unqInd] = unique([Mouse NeuronInd] ,'row');
Struct_valid=find(1-cell2mat(cellfun(@(x) sum(isnan(x)), StructureData, 'UniformOutput', false)));

for i=unqInd([45])'
    StructureStack=mat2gray(double(tiffreadVolume(StructureData{i})));
    StructureStack(StructureStack==0)=median(StructureStack(:));
    %StructureStack=StructureStack(:,:,1:118);
    %StructureStack_med=medfilt2_mov(StructureStack,[15 15]);
    illumination_field=imgaussfilt(max(StructureStack,[],3),50);
    StructureStack=StructureStack./illumination_field;
    StructureStack_Gauss=imgaussfilt3(StructureStack,[6 6 0.1]);
    %StructureStack_med(StructureStack_med==0)=median(StructureStack_med(:));
    %StructureStack=(StructureStack-StructureStack_med)./StructureStack_med;
    StructureStack_filt=(StructureStack-StructureStack_Gauss);
    StructureStack_filt=mat2gray(StructureStack_filt);
    StructureStack_bin=[]; level=[];
    level = graythresh(StructureStack_filt);
    StructureStack_bin=StructureStack_filt>level*1.2;
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

    rot_ang=90;
    Structure_ref=(imrotate(StructureStack_final,rot_ang));

    for cam=2%:2
    ref_img=Result.ref_im{cam}; ref_img(ref_img<prctile(ref_img(:),20))=median(ref_img(:)); ref_img=ref_img-prctile(ref_img(:),20);
    [RegImg,tformReg]=imReg(max(Structure_ref,[],3),ref_img);

    Result.Structure{cam}=RegImg;
    Result.Structure_bin{cam}=imwarp(max(imrotate(dendrite_bin,rot_ang),[],3), tformReg, 'OutputView', imref2d(size(ref_img)));
    Result.tform{cam}=tformReg;
    end
    save(fullfile(fpath{i},'Result.mat'),'Result','fpath','-v7.3')
end

%%
[ROI]=clicky(mov_mc{1}); ROImask=[];
for r=1:length(ROI) 
    ROImask(:,:,r)=roi2mask({ROI{r}},size(mov_mc{1},1),size(mov_mc{1},2));
end

tform_invcam1=invert(Result.tform{1});
tform_cam1to2=affine2d(tform_invcam1.T*Result.tform{2}.T);

warpRefIm2 = imrotate(imwarp(Result.ref_im{1}, tform_cam1to2, 'OutputView', imref2d(size(Result.ref_im{2}))),180);
ROImaskCam2= imwarp(ROImask, tform_cam1to2, 'OutputView', imref2d(size(Result.ref_im{2})));

Cam1_tr=tovec(mov_mc{1})'*tovec(ROImask);
Cam2_tr=tovec(mov_mc{2})'*tovec(ROImaskCam2);

Cam1_tr_hi=Cam1_tr-movmedian(Cam1_tr,100,1);
Cam2_tr_hi=Cam2_tr-movmedian(Cam2_tr,100,1);
%%

% StructureStack_final=mat2gray(double(tiffreadVolume(fullfile(fpath{i},['Structure_final.tiff']))));
% Structure_ref=(imrotate(StructureStack_final,rot_ang));

warpRefIm=[];
for cam=1:2
warpRefIm(:,:,cam) = imwarp(Result.ref_im{cam}, invert(Result.tform{cam}), 'OutputView', imref2d(size(max(Structure_ref,[],3))));
end

figure(20); clf;
tiledlayout(4,3)
ROIcrop=[363 733;309 1376];
nexttile([1 1]); 
imshow2(warpRefIm(ROIcrop(1,1):ROIcrop(1,2),ROIcrop(2,1):ROIcrop(2,2),1),[])
title('Camera 1')
nexttile([1 1]); 
imshow2(warpRefIm(ROIcrop(1,1):ROIcrop(1,2),ROIcrop(2,1):ROIcrop(2,2),2),[])
title('Camera 2')
nexttile([1 1]); 
imshow2(max(Structure_ref(ROIcrop(1,1):ROIcrop(1,2),ROIcrop(2,1):ROIcrop(2,2),:),[],3),[0 4])
title('Structure reconstruction',[])

show_time=[100:1120];
Cam1_tr_seg=-((Cam1_tr_hi(show_time,:)-median(Cam1_tr_hi(show_time,:),1))./median(Cam1_tr(show_time,:),1));
Cam2_tr_seg=-((Cam2_tr_hi(show_time,:)-median(Cam2_tr_hi(show_time,:),1))./median(Cam2_tr(show_time,:),1));

nexttile([1 1]);
show_footprnt_contour(double(ROImask),Result.ref_im{1},[1 0 0; 0 0.6 1])
ylabel('Camera 1','Rotation',0)
ax2=nexttile([1 2]);
l=plot([1:length(show_time)],Cam1_tr_seg);
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell([1 0 0; 0 0.6 1],2)); axis tight
xlabel('Time (ms)')
ylabel('\DeltaF/F')

nexttile([1 1]);
show_footprnt_contour(imrotate(ROImaskCam2,180),imrotate(Result.ref_im{2},180),[1 0 0; 0 0.6 1])
ylabel('Camera 2','Rotation',0)
ax1=nexttile([1 2]);
l=plot([1:length(show_time)],Cam2_tr_seg);
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell([1 0 0; 0 0.6 1],2));
xlabel('Time (ms)')
ylabel('\DeltaF/F')

ax3=nexttile(11,[1 2]);
plot(Result.Blue(show_time),'color',[0 0.6 1]); 
xlabel('Time (ms)')
ylabel('488nm AOTF')
linkaxes([ax1 ax2 ax3],'x')
axis tight;




