%% Load data

clear
clc;
sourcePath='/Volumes/BHL_WD18TB/20230917_BHLm077_78_vrPrism';
cd(sourcePath)
[fpath] = uigetfile_n_dir(); %only Treadmill data
%% Clicky

ROI_sz=10;

load([fpath{1} '/output_data.mat'])
sz=double(Device_Data{1, 3}.ROI([2 4]));
%mov=double(readBinMov([fpath '/frames1.bin'],sz(2),sz(1)));
DAQ_rate=Device_Data{1, 2}.buffered_tasks(1, 1).rate;
CamDAQ_rate=Device_Data{1, 2}.Counter_Inputs.rate;
CamTrig=Device_Data{1, 2}.Counter_Inputs.data;
CamTrig=find(CamTrig(2:end)-CamTrig(1:end-1)>0);
Frm_rate=(CamTrig(2)-CamTrig(1))/CamDAQ_rate;
mov=double(readBinMov_times([fpath{1} '/frames1.bin'],sz(2),sz(1),[1:10]));

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
cmap=distinguishable_colors(size(ROI,1));
figure(3); clf;
imshow2(mean(mov,3),[])
hold all
for n=1:size(ROI,1)
    ROI_box=[[ROI(n,1) ROI(n,3)]; [ROI(n,1) ROI(n,4)]; [ROI(n,2) ROI(n,4)]; [ROI(n,2) ROI(n,3)]; [ROI(n,1) ROI(n,3)];];
    plot(ROI_box(:,1),ROI_box(:,2),'color',cmap(n,:))
    text(mean(ROI_box(:,1)),mean(ROI_box(:,2)),num2str(n),'color',cmap(n,:))
end

%% Load Virmen data
[VirFile]=uigetfile_n_dir();
clear Dat
for i=1:length(VirFile)
    fid =   fopen(VirFile{i});
    VRdata = fread(fid,[12 inf],'double');
end

WorldITrack=find(VRdata(2,:)==1); %world 1
VRdata(5,WorldITrack)=(VRdata(5,WorldITrack)+6)*115/121;
%%
mov_trace=[];
if length(CamTrig)>5e4
   seg=[1:5e4:length(CamTrig)];
   seg=[seg length(CamTrig)];
   t_seg=[[1 seg(2:end-1)]' [seg(2:end-1)-1 seg(end)]'];
for n=1:size(coord,1)   
tmp=[];     
for t=1:size(t_seg,1)
%mov_seg=double(readBinMov_times_ROI([fpath '/frames1.bin'],sz(1),sz(2),[1:10],ROI));     
mov_seg=double(readBinMov_times_ROI([fpath{1} '/frames1.bin'],sz(1),sz(2),[t_seg(t,1):t_seg(t,2)],ROI(n,:)));    
tmp=[tmp; squeeze(mean(mov_seg,[1 2]))];
end
mov_trace(n,:)=tmp;
end
else
for n=1:size(coord,1)
mov_seg=double(readBinMov_times_ROI([fpath{1} '/frames1.bin'],sz(2),sz(1),[1:length(CamTrig)],ROI));
mov_trace(n,:)=squeeze(mean(mov_seg,[1 2]));
end
end

figure;
for n=1:size(coord,1)
plot(mov_trace(n,:),'color',cmap(n,:))
hold all
end
%%
clear Lap_FR Lap_V LickFR

    VRdata=VRdata(:,VRdata(10,:)>0);
    DAQ_rate=Device_Data{1, 2}.Counter_Inputs.rate;
    CamCounter=Device_Data{1, 2}.Counter_Inputs(1, 1).data;
    CamTrigger=find(CamCounter(2:end)-CamCounter(1:end-1));
    t_DAQ=CamTrigger/DAQ_rate;
    t_VR = datetime(datetime(VRdata(1,:),'ConvertFrom','datenum'), 'InputFormat', 'yyyy-MM-dd HH:mm:ss.SSS');
    t_VR= t_VR-t_VR(1);
    t_VR= milliseconds(t_VR)/1000;
    t_VR= t_VR*t_DAQ(end)/t_VR(end);
    VRdata(1,:)=t_VR;
    [Virmen_data_int vel_trace]=virmen_interpolate(VRdata,115,t_DAQ);
    Virmen_data_int(end+1,:)=vel_trace;

    HF=abs(Virmen_data_int(9,:)-movmedian(Virmen_data_int(9,:),500));
    Virmen_data_int(9,:)=(HF-mean(HF))/std(HF);

    Lick_trace=zeros(1,size(Virmen_data_int,2));
    Lick_trace(Virmen_data_int>2)=1;
     
    mov_trace_hi=-(mov_trace-movmedian(mov_trace,300,2));
    mov_trace_hi=mov_trace_hi./get_threshold(mov_trace_hi,1);
    sp=find_spike_bh(mov_trace_hi,4,3);

    for n=1:size(mov_trace_hi)
    [Lap_FR(:,:,n) Lap_V]=PlaceTrigger_average(sp(n,:),50,Virmen_data_int,0.003,115);
    end

    [LickFR]=PlaceTrigger_average(Virmen_data_int(9,:),50,Virmen_data_int,0.003,115);

%%

save_to='/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Updates/2023/20230918_Movs_Figs/';

figure(3); clf;
Nneuron=size(mov_trace_hi,1);
tiledlayout(3,4)
ax1=[];
for n=2%1:Nneuron
ax1=[ax1 nexttile([1 1])];
imagesc(Lap_FR(:,:,n)); hold all;
colormap('turbo')
title(['Firing rate of ROI#' num2str(n)])
set(gca,'FontSize',18)
end

ax1=[ax1 nexttile([1 1])];
imagesc(LickFR); hold all;
title('Lick rate')
set(gca,'FontSize',18)

ax1=[ax1 nexttile([1 1])];
imagesc(Lap_V,[0 0.01]); hold all;
colormap('turbo')
title('Speed')
colorbar
set(gca,'FontSize',18)

linkaxes(ax1,'x')
axis tight

nexttile([1 1])
imshow2(mean(mov,3),[])
hold all
for n=1:size(ROI,1)
    ROI_box=[[ROI(n,1) ROI(n,3)]; [ROI(n,1) ROI(n,4)]; [ROI(n,2) ROI(n,4)]; [ROI(n,2) ROI(n,3)]; [ROI(n,1) ROI(n,3)];];
    plot(ROI_box(:,1),ROI_box(:,2),'color',cmap(n,:))
    text(mean(ROI_box(:,1)),mean(ROI_box(:,2)),num2str(n),'color',cmap(n,:),'FontSize',18)
end


nexttile([2 4])
t_DAQ_scaled=[t_VR(1): (t_VR(end)-t_VR(1))/(size(mov_trace,2)-1):t_VR(end)];
for n=1:Nneuron
tr=rescale(mov_trace_hi(n,:));
plot(t_DAQ_scaled,tr+n,'color',cmap(n,:))
hold all
plot(t_DAQ_scaled(find(sp(n,:))),tr(find(sp(n,:)))+n,'r.')
end
plot(t_VR,VRdata(12,:)*Nneuron)
plot(t_VR,rescale(VRdata(5,:))*Nneuron)
ylabel('High Pass traces'); xlabel('Frames'); axis tight
title(fpath{1},'Interpreter','none')
set(gca,'FontSize',18)

set(gcf,'color','w');
set(gcf,'position',[100 1000 2000 1500]);

saveas(gcf,fullfile(save_to,'VR_Prism_BHL078Trace.fig'))
saveas(gcf,fullfile(save_to,'VR_Prism_BHL078Trace.png'))


