% Analysis on AAV expression sample and plot, in house YQ201
% 2022/11/09, Byung Hun Lee

clear
[fpath] = uigetfile_n_dir;
%%

DAQ_rate=0.00001;
for i=1:length(fpath)
load([fpath{i} '/output_data.mat'])
Blue=Device_Data{1, 2}.buffered_tasks(1, 1).channels(1, 2).data;
Airp=Device_Data{1, 2}.buffered_tasks(1, 1).channels(1, 3).data;

a=find(Device_Data{1, 2}.Counter_Inputs.data(2:end)-Device_Data{1, 2}.Counter_Inputs.data(1:end-1));
frm_rate=(a(2)-a(1))*DAQ_rate;
sz=Device_Data{1, 4}.ROI([2 4]);

mov=double(readBinMov([fpath{i} '/frames1.bin'],sz(2),sz(1)));
t=[1:size(mov,3)-1]*frm_rate; t(t>length(Blue)*DAQ_rate)=[];

[roi int]=clicky(mov); int=int-movmedian(int,500,1);
close all
figure; ax1=[];
tiledlayout(5,1)
nexttile([1 1])
imshow2(mean(mov,3),[]); hold all
for r=1:length(roi); plot(roi{r}(:,1),roi{r}(:,2),'linewidth',2); hold all; end
ax1=[ax1 nexttile([3 1])]; 
plot(t,-int(1:length(t),:)+[1:length(roi)]*10); set(gca,'FontSize',18)
title(fpath_mac2window([fpath{i}]))
ylim([0 (length(roi)+1)*10])
ax1=[ax1 nexttile([1 1])];
plot([1:length(Blue)]*0.00001,Blue,'linewidth',2); set(gca,'FontSize',18);
hold all
plot([1:length(Airp)]*0.00001,Airp,'linewidth',2); set(gca,'FontSize',18);
linkaxes(ax1,'x')
axis tight
set(gcf,'color','white')
saveas(gcf,[fpath{i} '/clicky_plot.fig'])
saveas(gcf,[fpath{i} '/clicky_plot.png'])
end
