clear

ROI_sz=7;
fpath=uigetdir;

fileList = dir(fullfile(fpath, '*.data'));

    if isempty(fileList)
disp('No Virmen data file in the folder')
[VirFile]=uigetfile_n_dir();
    else
       VirFile=fullfile(fpath,fileList.name);
    end

clear Dat        
    fid = fopen(VirFile);
    Dat = fread(fid,[12 inf],'double');



load([fpath '/output_data.mat'])
sz=double(Device_Data{1, 3}.ROI([2 4]));

DAQ_rate=Device_Data{1, 2}.buffered_tasks(1, 1).rate;
CamDAQ_rate=Device_Data{1, 2}.Counter_Inputs.rate;
CamTrig=Device_Data{1, 2}.Counter_Inputs.data;
CamTrig=find(CamTrig(2:end)-CamTrig(1:end-1)>0);
Frm_rate=(CamTrig(2)-CamTrig(1))/CamDAQ_rate;
mov=double(readBinMov_times([fpath '/frames1.bin'],sz(2),sz(1),[1:10]));
%%
figure(3); clf;
Rfixed = imref2d(size(Device_Data{1, 6}.refimage.img));
inverseTform = invert(Device_Data{1, 6}.tform);
revertedImage = imwarp(Device_Data{1, 6}.Target, inverseTform,'OutputView',Rfixed);

blueDMDimg=imcrop(revertedImage,Device_Data{1, 3}.ROI([1 3 2 4]));
bluePatt = bwboundaries(blueDMDimg);
imshow2(mean(mov,3),[])
hold all
plot(bluePatt{1}(:,2),bluePatt{1}(:,1),'color',[0 0.5 1])

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

%%
mperVR=6.44*1e-3;

Lap_FR=[]; Lap_V=[];  
Reward_lap=[];
rewardPos=0.8; LapDist=115;

    VRdata=Dat(:,find(Dat(10,:)));
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
    t_DAQ_scaled=rescale(t_DAQ)*t_VR(end);
    clear FR
    for n=1:size(mov_trace_hi)
    %[FR(:,:,n) V]=get_LapFR_VU(double(mov_trace_hi(n,:)),t_DAQ,30,VRdata,0.005,1,115);
    [FR(:,:,n) V]=get_LapFR_VU(double(movsum(sp(n,:),1000)),t_DAQ_scaled,30,VRdata,0.005,1,115);
    end
    [LickFR V]=get_LapFR_VU(double(Lick_trace),t_VR,30,VRdata,0,1,115);
Lap_FR=[Lap_FR; FR];
Lap_V=[Lap_V; V*mperVR;];

%%

figure(3); clf;
if size(coord,1) >5
tiledlayout(ceil((size(coord,1)+1)/2),6)
else
tiledlayout(size(coord,1)+1,3)
end
ax1=[];
ax2=[];

for noi=1:size(coord,1)
ax1=[ax1 nexttile([1 1])];
imagesc(Lap_FR(:,:,noi)); hold all;
title(['Firing rate #' num2str(noi)])
ax1=[ax1 nexttile([1 1])];
plot(mean(Lap_FR(:,:,noi),1,'omitnan'))

ax2=[ax2 nexttile([1 1])];
t_DAQ_scaled=[t_VR(1): (t_VR(end)-t_VR(1))/(size(mov_trace,2)-1):t_VR(end)];
plot(t_DAQ_scaled,mov_trace_hi(noi,:),'k',t_DAQ_scaled(find(sp(noi,:))),mov_trace_hi(noi,find(sp(noi,:))),'r.')
hold all
plot(t_VR,VRdata(12,:)*15)
plot(t_VR,rescale(VRdata(5,:))*15)

end

ax1=[ax1 nexttile([1 1])];
imagesc(LickFR); hold all;
title('Lick rate')

ax1=[ax1 nexttile([1 1])];
imagesc(Lap_V,[0 0.01]); hold all;
colormap('turbo')
title('Speed')
colorbar

nexttile([1 1]);
imshow2(mean(mov,3),[])
hold all
plot(coord(:,1),coord(:,2),'ro')
text(coord(:,1)+2,coord(:,2),num2str([1:size(coord,1)]'),'color','r','FontSize',10)

linkaxes(ax1,'x')
axis tight
linkaxes(ax2,'xy')
axis tight

saveas(gca,[fullfile(fpath, 'Brief.fig')])

