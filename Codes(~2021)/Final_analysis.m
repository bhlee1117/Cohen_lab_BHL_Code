%% Load the Final Mouse data
clear
[fnm pth]= uigetfile('.mat','Choose the final data','multiselect','on');  % Cropped image folder
ref=1;
%%
if ischar(fnm)
fnm={fnm};
end

for i=1:length(fnm)
    load([pth char(fnm{1,i})])
    session=size(Final.Data,2)-3;
ss_name{1,1}=char('HC');
ss_name{1,2}=char('CFC');
ss_name{1,3}=char('Ret.');
if session>3
    for j=1:session-3
        ss_name{1,j+3}=char(['Ret.' num2str(j+1)]);
    end
end

    im=imread([pth Final.Mouse '\MAX.tif']);
    g=1;
   for j=1:size(Final.Data,1)
       if isempty(find(Final.Data(j,4:end)==3))
           cred_data(g,1:session+4)=[Final.Data(j,1:end) j];
           g=g+1;
       end
   end
   TXN=[];
   
   for j=1:session
   TXN=[TXN sum(cred_data(:,j+3)==2)];
   end
   %Freeze=[Final.Freeze.Baseline Final.Freeze.CFC Final.Freeze.Retrieval];
   Freeze=[Final.Freeze.Baseline Final.Freeze.CFC Final.Freeze.Retrieval Final.Freeze.Retrieval2 Final.Freeze.Retrieval3 Final.Freeze.Retrieval4];% Final.Freeze.Retrieval5 Final.Freeze.Retrieval6];
   TXN_rate=TXN/size(cred_data,1)*100; % Percentage
   overlap3=[sum(cred_data(:,4)==2 & cred_data(:,5)==2 & cred_data(:,6)==2),....
             sum(cred_data(:,5)==2 & cred_data(:,6)==2 & cred_data(:,7)==2),...
             sum(cred_data(:,5)==2 & cred_data(:,6)==2 & cred_data(:,7)==2 & cred_data(:,8)==2),...
             sum(cred_data(:,5)==2 & cred_data(:,6)==2 & cred_data(:,7)==2 & cred_data(:,8)==2 & cred_data(:,9)==2)];
             %sum(cred_data(:,5)==2 & cred_data(:,6)==2 & cred_data(:,8)==2),...
             %sum(cred_data(:,5)==2 & cred_data(:,7)==2 & cred_data(:,8)==2),...
             %sum(cred_data(:,6)==2 & cred_data(:,7)==2 & cred_data(:,8)==2),...
             %sum(cred_data(:,4)==2 & cred_data(:,5)==2 & cred_data(:,6)==2 & cred_data(:,7)==2),...
             %sum(cred_data(:,4)==2 & cred_data(:,5)==2 & cred_data(:,6)==2 & cred_data(:,8)==2),...
             %sum(cred_data(:,4)==2 & cred_data(:,5)==2 & cred_data(:,7)==2 & cred_data(:,8)==2),...
             %sum(cred_data(:,4)==2 & cred_data(:,6)==2 & cred_data(:,7)==2 & cred_data(:,8)==2),...
             
%    overlap=[sum(cred_data(:,4)==2 & cred_data(:,5)==2) sum(cred_data(:,4)==2 & cred_data(:,6)==2) sum(cred_data(:,6)==2 & cred_data(:,5)==2) sum(cred_data(:,4)==2 & cred_data(:,5)==2 & cred_data(:,6)==2)]; %HC&CFC HC&Ret CFC&Ret HC&CFC&Ret
%    overlap_rate=[overlap(1,1)/TXN(1,1) overlap(1,1)/TXN(1,2) overlap(1,2)/TXN(1,1) overlap(1,2)/TXN(1,3) overlap(1,3)/TXN(1,2) overlap(1,3)/TXN(1,3)]; %P(HC&CFC|HC) P(HC&CFC|CFC) P(HC&Ret|HC) P(HC&Ret|Ret) P(CFC&Ret|CFC) P(CFC&Ret|Ret)
%
for k=1:session
    for j=1:session
        ov_mat(k,j)=sum(cred_data(:,k+3)==2 & cred_data(:,j+3)==2);
        if k>=j
        cond_mat(k,j)=ov_mat(k,j)/sum(cred_data(:,k+3)==2);
        else
        cond_mat(k,j)=ov_mat(k,j)/sum(cred_data(:,k+3)==2);
        end
    end
end
   figure(3*i-2)
   subplot(2,2,1) % Mouse Track
   for j=1:session-1
   plot(Final.Distance.track{1,j}(1:540,1),Final.Distance.track{1,j}(1:540,2))
   hold all
   end
   axis equal
   legend({'Baseline','Retrieval','2nd Ret.','3rd Ret.','Rem'},'Location','northwest','Orientation','horizontal','FontSize',10)
   fz_name=ss_name;
   fz_name{1,1}=char('Before Shock');
   fz_name{1,2}=char('After Shock');
   subplot(2,2,2) % Freezing Rate
   bar(categorical(fz_name),Freeze(1,1:session)*100,'FaceColor',[0.6 0.6 0.6],'BarWidth',0.8)
set(gca,'FontSize',8,'LineWidth',1);
   ylabel('Freezing (%)')
 %  xlim([0.5 session+0.5])
   
   subplot(2,2,3) % TXN number
     bar(categorical(ss_name),TXN,'FaceColor',[0.6 0.6 0.6],'BarWidth',0.8)
set(gca,'FontSize',10,'LineWidth',1)
  ylabel('# of TXN Cells')
 % xlim([0.5 session+0.5])
  
   subplot(2,2,4) % TXM rate
   bar(categorical(ss_name),TXN_rate,'FaceColor',[0.6 0.6 0.6],'BarWidth',0.8)
set(gca,'FontSize',10,'LineWidth',1)
ylabel('Fraction of TXN Cells (%)');
%xlim([0.5 session+0.5])

   figure(3*i-1)
subplot(2,2,1)   
imagesc(ov_mat)
set(gca,'FontSize',10,'LineWidth',1,'XTick',[1:1:session],'YTick',[1:1:session],'YTickLabel',ss_name,'XTickLabel',...
    ss_name)
axis equal
xlim([0.5 session+0.5])
ylim([0.5 session+0.5])
for i=1:session
    for j=1:session
   text(i-0.2,j,num2str(ov_mat(i,j)))
    end
end

subplot(2,2,3)   
imagesc(cond_mat)
set(gca,'FontSize',10,'LineWidth',1,'XTick',[1:1:session],'YTickLabel',ss_name,'XTickLabel',...
    ss_name,'YTick',[1:1:session])
axis equal 
for i=1:session
    for j=1:session
   text(i-0.2,j,num2str(cond_mat(j,i),2))
    end
end

xlim([0.5 session+0.5])
ylim([0.5 session+0.5])
subplot(1,2,2)
clear act_map
% act_map{1}=rearrange_actmap(Final.Data(:,4:end),ref,[2 1 3]);
% l{1}=act_map{1}(:,ref)==2;
% for j=1:session-ref
%     if sum(l{j})>0
% l{j+1}=act_map{j}(1:max(find(l{j}==1)),j+ref-1)==2;
% act_map{j+1}=rearrange_actmap(act_map{j}(1:max(find(l{j}==1)),:),j+ref,[2 1 3]);
% act_map{j+1}(max(find(l{j}==1)):size(act_map{j},1),:)=act_map{j}(max(find(l{j}==1)):size(act_map{j},1),:);
%     else
%         act_map{j}=act_map{j-1};
%         l{j}=l{j-1};
%     end
% end
ref=2;
act_map{1}=rearrange_actmap(Final.Data(:,4:end),ref,[2 1 3]);
l{1}=act_map{1}(:,ref)==2;
m_row(1)=max(find(l{1}==1));
for j=1:3
    if sum(l{j})>0
        m_row(j)=max(find(l{j}==1));
act_map{j+1}=rearrange_actmap(act_map{j}(1:m_row(j),:),j+ref,[2 1 3]);
act_map{j+1}(m_row(j)+1:size(act_map{j},1),:)=rearrange_actmap(act_map{j}(m_row(j)+1:size(act_map{j},1),:),j+ref,[2 1 3]);
        l{j+1}=act_map{j}(:,ref)==2;
        
    else
        act_map{j}=act_map{j-1};
        l{j}=l{j-1};
    end
end
imagesc(act_map{end})
set(gca,'FontSize',10,'LineWidth',1,'XTick',[1:1:session],'XTickLabel',...
    ss_name)

% bar([1:1:6],overlap_rate*100,'FaceColor',[0.6 0.6 0.6],'BarWidth',0.8)
% set(gca,'FontSize',10,'LineWidth',1,'XTick',[1:1:6],'XTickLabel',...
%     {'P(HC&CFC|HC)','P(HC&CFC|CFC)','P(HC&Ret|HC)','P(HC&Ret|Ret)','P(CFC&Ret|CFC)','P(CFC&Ret|Ret)'})
%    ylabel('Conditional Prob. (%)')
%    xlim([0.5 6.5])
   
% subplot(2,2,3)
% bar([1:1:4],overlap,'FaceColor',[0.6 0.6 0.6],'BarWidth',0.8)
% set(gca,'FontSize',10,'LineWidth',1,'XTick',[1:1:4],'XTickLabel',...
%     {'HC & CFC','HC & Ret.','CFC & Ret.','HC & CFC & Ret.'});
% ylabel('# of TXN Cells')
% xlim([0.5 4.5])
% 
% subplot(2,2,4)
% venn_data=[TXN overlap];
% venn(venn_data, 'ErrMinMode', 'ChowRodgers')
%  legend({'HC','CFC','Ret.'},'Location','northwest','FontSize',10)
 
%    figure(3*i)
%    imagesc(im)
%    colormap('gray')
%    hold all
%    for j=1:session-1
%        clear tmp
%        color=zeros(1,3);
%        color(1,j)=1;
%        mask=cred_data(:,j+4)==2;
%        tmp=cred_data(:,1:3);
%        tmp(~mask,1:3)=NaN;
%        plot3(tmp(:,1),tmp(:,2),tmp(:,3),'color',color,'marker','o','linestyle','none','markersize',3*j-2,'linewidth',6-j)
%    end
%    legend({'CFC','Ret.','2nd Ret.','3rd Ret.','4th Ret.'},'Location','northwest','FontSize',10)
end


