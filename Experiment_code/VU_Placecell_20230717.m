% Experiment code
% z-stack
%cd 'C:\Users\Labmember\Data\ByungHun'
addpath(genpath('C:\Users\Labmember\Data\ByungHun\GitHub_code'))
Stage=xx.GetDevice("Scientifica_Stage_Controller");
cam=xx.GetDevice('Cam_Controller');
Shutters=xx.GetDevice('Shutter_Device');
dmd = xx.GetDevice('DMD_Device'); %dmd(1,1) is the ALP for the Behavior
Saveto=xx.datafolder;
   cam1 = cam(1,1);
   %cam2 = cam(1,2);
%%
radius=14;
addpath('C:\Users\Labmember\Data\ByungHun')
fpath='C:\Users\Labmember\Data\ByungHun\Optopatch\20230717\Snaps';
files = dir(fullfile(fpath, '*.mat'));
files = files(~[files.isdir]);
[~, idx] = max([files.datenum]);
recent_file = files(idx).name;
load(fullfile(fpath,recent_file));
ref_im=snap.img; sz_ref=size(ref_im);

centers=Cell_segment_circle_10x_VU(ref_im,0.92);
centers=cell_detection_manual(ref_im,centers,[]);
figure;
imshow2(ref_im,[]);
hold all
plot(centers(:,1),centers(:,2),'ro')

dmd_now=xx.GetDevice('DMD_Device');
mask_im = zeros(sz_ref);
for i = 1:size(centers, 1)
    [columnsInImage, rowsInImage] = meshgrid(1:sz_ref(2), 1:sz_ref(1));
    circlePixels = (rowsInImage - centers(i,2)).^2 ...
        + (columnsInImage - centers(i,1)).^2 <= radius.^2;
    mask_im(circlePixels) = 1; % change pixel intensity to 255 (white) inside the circle
end

%Orange
transform=dmd_now(1,1).tform;
%transform=snap.tform;
Rfixed = imref2d([dmd_now(1,1).Dimensions(2) dmd_now(1,1).Dimensions(1)]);

newSize = [2304 2304];
rmask = imresize(mask_im, newSize);
% apply transformation from image space to patterning device space
tmask=imwarp(rmask,transform,'OutputView',Rfixed);

dmd_now(1,1).Target=tmask;
dmd_now(1,1).Write_Static();

if  Shutters(1, 1).State
    else
    Shutters(1, 1).State=1;
end
pause(cam1.exposuretime)
CurrentIm = cam1.Snap();
Shutters(1, 1).State=0;

figure(2); clf;
ax1=nexttile([1 1]);
imshow2(ref_im,[]); hold all;
plot(centers(:,1),centers(:,2),'ro')
title(['Before masking'])
ax2=nexttile([1 1]);
imshow2(CurrentIm,[]); hold all
plot(centers(:,1),centers(:,2),'ro')
title(['After masking'])
linkaxes([ax1 ax2],'xy')
tt=clock; t = datetime(tt);
t.Format='ddHHmm';
saveas(gca,[Saveto '\' char(t) 'Masking.fig'])
saveas(gca,[Saveto '\' char(t) 'Masking.png'])

%% Clicky
ROI_sz=10;
fpath=uigetdir;

load([fpath '/output_data.mat'])
sz=double(Device_Data{1, 3}.ROI([2 4]));
%mov=double(readBinMov([fpath '/frames1.bin'],sz(2),sz(1)));
DAQ_rate=Device_Data{1, 2}.buffered_tasks(1, 1).rate;
CamDAQ_rate=Device_Data{1, 2}.Counter_Inputs.rate;
CamTrig=Device_Data{1, 2}.Counter_Inputs.data;
CamTrig=find(CamTrig(2:end)-CamTrig(1:end-1)>0);
Frm_rate=(CamTrig(2)-CamTrig(1))/CamDAQ_rate;
mov=double(readBinMov_times([fpath '/frames1.bin'],sz(2),sz(1),[1:10]));
%%
figure(3); clf;
imshow2(mean(mov,3),[])
g=1; coord=[];
while g
[x y]=ginput(1);
if isempty(x)
    g=0;
end
coord=[coord; [x y]];
end
close(figure(3));
ROI=round([coord(:,1)-ROI_sz coord(:,1)+ROI_sz coord(:,2)-ROI_sz coord(:,2)+ROI_sz]);
mov_trace=[];
if length(CamTrig)>1e5
   seg=[1:1e5:length(CamTrig)];
   seg=[seg length(CamTrig)];
   t_seg=[[1 seg(2:end-1)]' [seg(2:end-1)-1 seg(end)]'];
for n=1:size(coord,1)   
tmp=[];     
for t=1:size(t_seg,1)
%mov_seg=double(readBinMov_times_ROI([fpath '/frames1.bin'],sz(1),sz(2),[1:10],ROI));     
mov_seg=double(readBinMov_times_ROI([fpath '/frames1.bin'],sz(1),sz(2),[t_seg(t,1):t_seg(t,2)],ROI(n,:)));    
tmp=[tmp; squeeze(mean(mov_seg,[1 2]))];
end
mov_trace(n,:)=tmp;
end
else
for n=1:size(coord,1)
mov_seg=double(readBinMov_times_ROI([fpath '/frames1.bin'],sz(2),sz(1),[1:length(CamTrig)],ROI));
mov_trace(n,:)=squeeze(mean(mov_seg,[1 2]));
end
end

figure;
plot(-mov_trace')
%% Load Virmen data
[VirFile]=uigetfile_n_dir();
clear Dat
for i=1:length(VirFile)
    fid = fopen(VirFile{i});
    Dat{i} = fread(fid,[12 inf],'double');
end
%%
mperVR=6.44*1e-3;

Lap_FR=[]; Lap_V=[];  
Reward_lap=[];
rewardPos=0.8; LapDist=115;
for i=1:length(VirFile)
    VRdata=Dat{i}(:,find(Dat{i}(10,:)));
    %VRdata(7,:)=zscore(abs(VRdata(7,:)-movmedian(VRdata(7,:),20)));
    %HF=abs(VRdata(7,:)-movmedian(VRdata(7,:),20));
    HF=abs(VRdata(9,:)-movmedian(VRdata(9,:),500));
    VRdata(9,:)=(HF-mean(HF))/std(HF);
    [pks lick]=findpeaks(VRdata(9,:));
    Artifact=abs(VRdata(5,:)-rewardPos*LapDist)<0.025*LapDist;
    bwArtifact=bwlabel(Artifact);
    ArtiArg=[];

    for b=1:max(bwArtifact)
        region=find(bwArtifact==b);
        pks_in_region=find(ismember(lick,region));
        if ~isempty(pks_in_region)
            [m, a]=max(pks(pks_in_region));
            if m>5
                ArtiArg=[ArtiArg lick(pks_in_region(a))-4:lick(pks_in_region(a))];
            end
        end
    end
    Lick_trace=zeros(1,size(VRdata,2));
    Lick_trace(lick(pks>2))=1;
    Lick_trace(ArtiArg)=0;
    % plot(rescale(VRdata(7,:)))
    % hold all
    % plot(rescale(VRdata(4,:)))
    t_VR = datetime(datetime(VRdata(1,:),'ConvertFrom','datenum'), 'InputFormat', 'yyyy-MM-dd HH:mm:ss.SSS');
    t_VR= t_VR-t_VR(1);
    t_VR= milliseconds(t_VR)/1000;

    %[Lap_FR Lap_V]=get_LapFR_VU(double(VRdata(7,:)>0.1),100,VRdata([1 1 2 3 4 5 6 8],:),0.005,1);
    mov_trace_hi=-(mov_trace-movmedian(mov_trace,300,2));
    mov_trace_hi=mov_trace_hi./get_threshold(mov_trace_hi,1);
    sp=find_spike_bh(mov_trace_hi,4,3);
    %mov_trace_hi=-movmean((mov_trace-movmedian(mov_trace,500,2)),20,2);
    
    t_DAQ=[CamTrig/CamDAQ_rate];
    clear FR
    for n=1:size(mov_trace_hi)
    %[FR(:,:,n) V]=get_LapFR_VU(double(mov_trace_hi(n,:)),t_DAQ,30,VRdata,0.005,1,115);
    [FR(:,:,n) V]=get_LapFR_VU_DAQ(double(movsum(sp(n,:),1000)),t_DAQ,30,VRdata,0,1,115);
    end
    [LickFR V]=get_LapFR_VU(double(Lick_trace),t_VR,30,VRdata,0,1,115);
Lap_FR=[Lap_FR; FR];
Lap_V=[Lap_V; V*mperVR;];

end
%%
noi=1;
figure(3); clf;
ax1=[];
ax1=[ax1 nexttile([1 1])];
imagesc(Lap_FR(:,:,noi)); hold all;
title('Firing rate')
ax1=[ax1 nexttile([1 1])];
plot(mean(Lap_FR(:,:,noi),1,'omitnan'))

ax1=[ax1 nexttile([1 1])];
imagesc(LickFR); hold all;
title('Lick rate')
ax1=[ax1 nexttile([1 1])];
plot(mean(LickFR,1,'omitnan'))

ax1=[ax1 nexttile([1 1])];
imagesc(Lap_V,[0 0.01]); hold all;
colormap('turbo')
title('Speed')
colorbar

ax1=[ax1 nexttile([1 1])];
plot(mean(Lap_V,1,'omitnan'))

linkaxes(ax1,'x')
axis tight

nexttile([1 1])
t_DAQ_scaled=[t_VR(1): (t_VR(end)-t_VR(1))/(size(mov_trace,2)-1):t_VR(end)];
plot(t_DAQ_scaled,mov_trace_hi(noi,:),'k',t_DAQ_scaled(find(sp(noi,:))),mov_trace_hi(noi,find(sp(noi,:))),'r.')
hold all
plot(t_VR,VRdata(12,:)*15)
plot(t_VR,rescale(VRdata(5,:))*15)

filename=char(VirFile);
saveas(gca,[filename(1:end-5) '.fig'])


