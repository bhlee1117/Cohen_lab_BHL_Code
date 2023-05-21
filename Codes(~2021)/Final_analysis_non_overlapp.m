%% Load the Final Mouse data
clear
[fnm pth]= uigetfile('.mat','Choose the final data','multiselect','on');  % Cropped image folder

%%
if ischar(fnm)
fnm={fnm};
end

    load([pth char(fnm{1,1})])
    g=1;
   
   TXN=[sum(Final.Data{1}(:,4)==2) sum(Final.Data{2}(:,4)==2) sum(Final.Data{3}(:,4)==2)];
   Freeze=[Final.Freeze.Baseline Final.Freeze.CFC Final.Freeze.Retrieval];
   for i=1:3
       TXN_rate(1,i)=TXN(1,i)/size(Final.Data{i},1); % Percentage
   end
   
   figure(3*i-2)
   subplot(2,2,1) % Mouse Track
   plot(Final.Distance.track{1,1}(1:540,1),Final.Distance.track{1,1}(1:540,2),'g-')
   hold all
   plot(Final.Distance.track{1,2}(1:540,1),Final.Distance.track{1,2}(1:540,2),'r-')
   axis equal
   legend({'Baseline','Retrieval'},'Location','northwest','Orientation','horizontal','FontSize',10)
   
   subplot(2,2,2) % Freezing Rate
   bar([1:1:3],Freeze*100,'FaceColor',[0.6 0.6 0.6],'BarWidth',0.8)
set(gca,'FontSize',8,'LineWidth',1,'XTick',[1:1:3],'XTickLabel',...
    {'Before Shock','After Shock','Retrieval'});
   ylabel('Freezing (%)')
   xlim([0.5 3.5])
   
   subplot(2,2,3) % TXN number
     bar([1:1:3],TXN,'FaceColor',[0.6 0.6 0.6],'BarWidth',0.8)
set(gca,'FontSize',10,'LineWidth',1,'XTick',[1:1:3],'XTickLabel',...
    {'HC','CFC','Ret.'})
  ylabel('# of TXN Cells')
  xlim([0.5 3.5])
  
   subplot(2,2,4) % TXM rate
   bar([1:1:3],TXN_rate,'FaceColor',[0.6 0.6 0.6],'BarWidth',0.8)
set(gca,'FontSize',10,'LineWidth',1,'XTick',[1:1:3],'XTickLabel',...
    {'HC','CFC','Ret.'})
ylabel('Fraction of TXN Cells (%)');
xlim([0.5 3.5])

%%
% rearr=rearrange_actmap(cred_data,4,[2 1 3]);
% subplot(1,3,1)
% imagesc(rearr(:,4:6))
% set(gca,'FontSize',10,'LineWidth',1,'XTick',[1:1:3],'XTickLabel',...
%     {'HC','CFC','Retrieval'});
% ylabel('# of TXN Cells')
% 
% hold all

rearr=rearrange_actmap(cred_data,5,[2 1 3]);
subplot(1,2,1)
imagesc(rearr(:,4:6))
set(gca,'FontSize',10,'LineWidth',1,'XTick',[1:1:3],'XTickLabel',...
    {'HC','CFC','Retrieval'});
ylabel('# of TXN Cells')
colormap('copper')
rearr=rearrange_actmap(cred_data,6,[2 1 3]);
subplot(1,2,2)
imagesc(rearr(:,4:6))
set(gca,'FontSize',10,'LineWidth',1,'XTick',[1:1:3],'XTickLabel',...
    {'HC','CFC','Retrieval'});
ylabel('# of TXN Cells')

