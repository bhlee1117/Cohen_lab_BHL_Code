 clear
Path_n='H:\Image_data\Mouse_behavior\20180607-08\';
File_n='MVI_7382.MOV';
File_R='MVI_7384.MOV';

data_exp=importdata([Path_n File_n '.mat']);
data_expR=importdata([Path_n File_R '.mat']);
thres=1300;
control_time=[1 180];
shock_time=[180 240];
ret_time=[1 180];
%%

for i=1:size(data_exp,3)
fz_e(i)=sum(sum(data_exp(:,:,i)));

if fz_e(i)>thres
fz_bino_e(i)=0;
else
fz_bino_e(i)=1;   
end
end
% figure(1)
% plot(fz_bino_e,'o')
% title('Freezing Session')

figure(1)
subplot(3,1,1)
plot(fz_e(1:240),'color','k')
ylim([0 8000])
hold all
 line([shock_time(1,1) shock_time(1,1)],[0 max(fz_e(1:180))],'color','r')
 line([0 shock_time(1,2)],[thres thres],'color','r','LineStyle','--')
ylabel('Movement (Pixel)')
xlabel('Time (sec)')
subplot(3,1,2:3)
bar([1 2],[mean(fz_bino_e(control_time(1,1):control_time(1,2))) mean(fz_bino_e(shock_time(1,1):shock_time(1,2)))]*100,'FaceColor',[0.6 0.6 0.6],'BarWidth',0.8)
set(gca,'FontSize',15,'LineWidth',1,'XTick',[1 2],'XTickLabel',...
    {'Before Shock','After Shock'});

% bar([1],[mean(fz_bino_e(control_time(1,1):control_time(1,2))) ]*100,'FaceColor',[0.6 0.6 0.6],'BarWidth',0.8)
% set(gca,'FontSize',15,'LineWidth',1,'XTick',[1 2],'XTickLabel',...
%     {'Retrieval in context'});
xlim([0.5 2.5])
ylabel('Freezing (%)');


   
        
        
     for i=1:size(data_expR,3)
fz_eR(i)=sum(sum(data_expR(:,:,i)));

if fz_eR(i)>thres
fz_bino_eR(i)=0;
else
fz_bino_eR(i)=1;   
end
end
% figure(1)
% plot(fz_bino_e,'o')
% title('Freezing Session')

figure(2)
subplot(3,1,1)
plot(fz_eR(1:180),'color','k')
ylim([0 8000])
hold all
 line([0 ret_time(1,2)],[thres thres],'color','r','LineStyle','--')
ylabel('Movement (Pixel)')
xlabel('Time (sec)')
subplot(3,1,2:3)
bar([1 2],[mean(fz_bino_eR(ret_time(1,1):ret_time(1,2))) 0]*100,'FaceColor',[0.6 0.6 0.6],'BarWidth',0.8)
set(gca,'FontSize',15,'LineWidth',1,'XTick',[1 2],'XTickLabel',...
    {'Retrieval',''});

% bar([1],[mean(fz_bino_e(control_time(1,1):control_time(1,2))) ]*100,'FaceColor',[0.6 0.6 0.6],'BarWidth',0.8)
% set(gca,'FontSize',15,'LineWidth',1,'XTick',[1 2],'XTickLabel',...
%     {'Retrieval in context'});
xlim([0.5 2.5])
ylabel('Freezing (%)');
   

%%