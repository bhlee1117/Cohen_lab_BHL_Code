%% Load mouse data

clear
[fnm pth]=uigetfile('.mat','Select the analyzed mouse data','Multiselect','on');

%%
if ischar(fnm)
    fnm={fnm};
end
lg=[];
for i=1:length(fnm)
    mousedat{i}=importdata([pth fnm{1,i}]);
    g=1;
    if iscell(mousedat{1,i}.Data)
    for j=1:size(mousedat{1,i}.Data,2)
        fracTXN(i,j)=sum(mousedat{1,i}.Data{1,j}(:,4)==2)/sum(mousedat{1,i}.Data{1,j}(:,4)~=3); %columns for HC Tr Ret... Rows for each mouse
        numTXN(i,j)=sum(mousedat{1,i}.Data{1,j}(:,4)==2);
    end
    else
        for j=4:size(mousedat{1,i}.Data,2)
        fracTXN(i,j-3)=sum(mousedat{1,i}.Data(:,j)==2)/sum(mousedat{1,i}.Data(:,j)~=3);
        numTXN(i,j-3)=sum(mousedat{1,i}.Data(:,j)==2);
        for k=1:size(mousedat{1,i}.Data,2)-j
        overlapTXN(i,g)=sum(mousedat{1,i}.Data(:,j)==2 & mousedat{1,i}.Data(:,j+k) ==2);
        g=g+1;
        end
        end
        cond_prob(i,:)=[overlapTXN(i,1)/numTXN(i,2) overlapTXN(i,2)/numTXN(i,3) overlapTXN(i,3)/numTXN(i,4) overlapTXN(i,4)/numTXN(i,3) overlapTXN(i,5)/numTXN(i,2) overlapTXN(i,6)/numTXN(i,3)];% HC&CFC/CFC, HC&Ret/Ret, HC&Rem/Rem, CFC&Ret/Ret, CFC&Rem/CFC, Ret&Rem/Ret
    end
    f=table2array(struct2table(mousedat{1,i}.Freeze));
 Freeze(i,1:size(f,2))=f;
 lg=[lg; mousedat{1,i}.Mouse];
 
end
fracTXN(fracTXN==0)=NaN;
Freeze(Freeze==0)=NaN;
cond_prob(cond_prob==0)=NaN;
%% 'Fraction of TXN Cell (%)'

bar([1:1:6],mean(fracTXN,'omitnan')*100,'FaceColor',[0.6 0.6 0.6],'BarWidth',0.8)
set(gca,'FontSize',8,'LineWidth',1,'XTick',[1:1:5],'XTickLabel',...
    {'HC','CFC','1^s^t Ret.','2^n^d Ret.','3^r^d Ret.'});
   ylabel('Fraction of TXN Cell (%)')
   xlim([0.5 5.5])
   hold all
   errorbar([1:1:6],mean(fracTXN,'omitnan')*100,std(fracTXN,'omitnan')/sqrt(size(fracTXN,1))*100,'k','linestyle','none')
   
% for i=1:size(fracTXN,1)
%    clear color
%    color=zeros(1,3);
%     color(1,i)=1;
%     p{i}=plot([1:1:sum(fracTXN(i,:)~=0)],[fracTXN(i,1:sum(fracTXN(i,:)~=0))]*100,'color',color,'marker','.','linestyle','none','markersize',30);
%     hold all
% end
% legend([p{1} p{2} p{3}],lg)
   
%% 'Freezing (%)'

bar([1:1:4],mean([Freeze(:,1) Freeze(:,3:5)],'omitnan')*100,'FaceColor',[0.6 0.6 0.6],'BarWidth',0.8)
set(gca,'FontSize',8,'LineWidth',1,'XTick',[1:1:5],'XTickLabel',...
    {'HC','1^s^t Ret.','2^n^d Ret.','3^r^d Ret.'});
   ylabel('Freezing (%)')
   xlim([0.5 4.5])
   hold all
   errorbar([1:1:4],mean([Freeze(:,1) Freeze(:,3:5)],'omitnan')*100,std([Freeze(:,1) Freeze(:,3:5)],'omitnan')/sqrt(size(fracTXN,1))*100,'k','linestyle','none')
   
% for i=1:size(Freeze,1)
%    clear color
%    color=zeros(1,3);
%     color(1,i)=1;
%     p{i}=plot([1:1:sum(Freeze(i,:)~=0)],[Freeze(i,1:sum(Freeze(i,:)~=0))]*100,'color',color,'marker','.','linestyle','none','markersize',30);
%     hold all
% end
% legend([p{1} p{2} p{3}],lg)
%% 'Conditional Probability (%)'

bar([1:1:6],mean(cond_prob,'omitnan')*100,'FaceColor',[0.6 0.6 0.6],'BarWidth',0.8)
set(gca,'FontSize',8,'LineWidth',1,'XTick',[1:1:6],'XTickLabel',...
    {'HC&CFC/CFC', 'HC&Ret/Ret', 'HC&Rem/Rem', 'CFC&Ret/Ret', 'CFC&Rem/CFC', 'Ret&Rem/Ret'});
   ylabel('Conditional Probability (%)')
   xlim([0.5 6.5])
   hold all
   errorbar([1:1:6],mean(cond_prob,'omitnan')*100,std(cond_prob,'omitnan')*100,'k')
   
for i=1:size(cond_prob,1)
   clear color
   color=zeros(1,3);
    color(1,i)=1;
    p{i}=plot([1:1:sum(cond_prob(i,:)~=0)],[cond_prob(i,1:sum(cond_prob(i,:)~=0))]*100,'color',color,'marker','.','linestyle','none','markersize',30);
    hold all
end
legend([p{2} p{3}],lg)
%%
for i=1:size(fracTXN,1)
p1=plot(100*Freeze(i,3),fracTXN(i,3),'r.','linestyle','none','markersize',20);
text(100*Freeze(i,3)+3,fracTXN(i,3),mousedat{1,i}.Mouse,'color','r')
hold all
end
for i=1:size(fracTXN,1)
p2=plot(100*Freeze(i,2),fracTXN(i,2),'b.','linestyle','none','markersize',20);
text(100*Freeze(i,2)+3,fracTXN(i,2),mousedat{1,i}.Mouse,'color','b')
hold all
end
for i=1:size(fracTXN,1)
p3=plot(100*Freeze(i,1),fracTXN(i,1),'g.','linestyle','none','markersize',20);
text(100*Freeze(i,1)+3,fracTXN(i,1),mousedat{1,i}.Mouse,'color','g')
hold all
end
for i=1:size(fracTXN,1)
p4=plot(100*Freeze(i,4),fracTXN(i,4),'k.','linestyle','none','markersize',20);
text(100*Freeze(i,4)+3,fracTXN(i,4),mousedat{1,i}.Mouse,'color','k')
hold all
end
xlim([0 100])
ylim([0 0.5])
xlabel('Freezing (%)')
ylabel('Fraction of TXN Cell')
legend([p1 p2 p3 p4],'Retrieval','CFC','HC','Remote')

%%
for i=1:size(fracTXN,1)
    if ~iscell(mousedat{1,i}.Data)
    pie_x(i,:)=[sum(mousedat{1,i}.Data(:,4)==2 & mousedat{1,i}.Data(:,5)~=2 & mousedat{1,i}.Data(:,6)~=2 & mousedat{1,i}.Data(:,7)~=2),...
                sum(mousedat{1,i}.Data(:,4)~=2 & mousedat{1,i}.Data(:,5)==2 & mousedat{1,i}.Data(:,6)~=2 & mousedat{1,i}.Data(:,7)~=2),...
                sum(mousedat{1,i}.Data(:,4)~=2 & mousedat{1,i}.Data(:,5)~=2 & mousedat{1,i}.Data(:,6)==2 & mousedat{1,i}.Data(:,7)~=2),...
                sum(mousedat{1,i}.Data(:,4)~=2 & mousedat{1,i}.Data(:,5)~=2 & mousedat{1,i}.Data(:,6)~=2 & mousedat{1,i}.Data(:,7)==2),...
                sum(mousedat{1,i}.Data(:,4)==2 & mousedat{1,i}.Data(:,5)==2 & mousedat{1,i}.Data(:,6)~=2 & mousedat{1,i}.Data(:,7)~=2),...
                sum(mousedat{1,i}.Data(:,4)==2 & mousedat{1,i}.Data(:,5)~=2 & mousedat{1,i}.Data(:,6)==2 & mousedat{1,i}.Data(:,7)~=2),...
                sum(mousedat{1,i}.Data(:,4)==2 & mousedat{1,i}.Data(:,5)~=2 & mousedat{1,i}.Data(:,6)~=2 & mousedat{1,i}.Data(:,7)==2),...
                sum(mousedat{1,i}.Data(:,4)~=2 & mousedat{1,i}.Data(:,5)==2 & mousedat{1,i}.Data(:,6)==2 & mousedat{1,i}.Data(:,7)~=2),...
                sum(mousedat{1,i}.Data(:,4)~=2 & mousedat{1,i}.Data(:,5)==2 & mousedat{1,i}.Data(:,6)~=2 & mousedat{1,i}.Data(:,7)==2),...
                sum(mousedat{1,i}.Data(:,4)~=2 & mousedat{1,i}.Data(:,5)~=2 & mousedat{1,i}.Data(:,6)==2 & mousedat{1,i}.Data(:,7)==2),...
                sum(mousedat{1,i}.Data(:,4)==2 & mousedat{1,i}.Data(:,5)==2 & mousedat{1,i}.Data(:,6)==2 & mousedat{1,i}.Data(:,7)~=2),...
                sum(mousedat{1,i}.Data(:,4)==2 & mousedat{1,i}.Data(:,5)==2 & mousedat{1,i}.Data(:,6)~=2 & mousedat{1,i}.Data(:,7)==2),...
                sum(mousedat{1,i}.Data(:,4)~=2 & mousedat{1,i}.Data(:,5)==2 & mousedat{1,i}.Data(:,6)==2 & mousedat{1,i}.Data(:,7)==2),...
                sum(mousedat{1,i}.Data(:,4)==2 & mousedat{1,i}.Data(:,5)~=2 & mousedat{1,i}.Data(:,6)==2 & mousedat{1,i}.Data(:,7)==2),...
                sum(mousedat{1,i}.Data(:,4)==2 & mousedat{1,i}.Data(:,5)==2 & mousedat{1,i}.Data(:,6)==2 & mousedat{1,i}.Data(:,7)==2),...
                sum(mousedat{1,i}.Data(:,4)~=2 & mousedat{1,i}.Data(:,5)~=2 & mousedat{1,i}.Data(:,6)~=2 & mousedat{1,i}.Data(:,7)~=2)
                ];%HC CFC Ret Rem only,HC&CFC HC&Ret HC&Rem CFC&Ret CFC&Rem Ret&Rem HC&CFC&Ret HC&CFC&Rem CFC&Ret&Rem HC&Ret&Rem all
    figure(i)
pie(pie_x(i,:),{'HC', 'CFC', 'Ret', 'Rem','HC&CFC', 'HC&Rt', 'HC&Rm', 'CFC&Rt', 'CFC&Rm', 'Rt&Rm', '','','','','all','None'})
lg_con={'HC', 'CFC', 'Ret', 'Rem','HC&CFC', 'HC&Rt', 'HC&Rm', 'CFC&Rt', 'CFC&Rm', 'Rt&Rm', 'H&C&Rt','H&C&Rm','C&Rt&Rm','H&Rt&Rm','all','None'};
g=1;
for ii=1:size(lg_con,2)
    if pie_x(i,ii)~=0
    lg_con2(1,g)=lg_con(1,ii);
    g=g+1;
    end
end
legend(lg_con2)
    title(mousedat{1,i}.Mouse)
    end
    
end

%%
clear chance_overlap
        g=1;
        gg=1;
for i=1:size(mousedat,2)
    if mousedat{1,i}.overlap==1
        chance_overlap{1}(gg,1:2)=[numTXN(i,2)/size(mousedat{1,i}.Data,1)*numTXN(i,3)/size(mousedat{1,i}.Data,1) sum(mousedat{1,i}.Data(:,5)==2&mousedat{1,i}.Data(:,6)==2)/size(mousedat{1,i}.Data,1) ];
        chance_overlap_mouseinfo{1}{gg,1}=mousedat{1,i}.Mouse;
        if ~isempty(findstr(cell2mat(mousedat{1,i}.behavior),'Ret2'))
            chance_overlap_mouseinfo{2}{g,1}=mousedat{1,i}.Mouse;  
        chance_overlap{2}(g,1:2)=[numTXN(i,2)/size(mousedat{1,i}.Data,1)*numTXN(i,4)/size(mousedat{1,i}.Data,1)...      %CFC&Ret2
                              sum(mousedat{1,i}.Data(:,5)==2&mousedat{1,i}.Data(:,7)==2)/size(mousedat{1,i}.Data,1) ];  
        chance_overlap{3}(g,1:2)=[numTXN(i,2)/size(mousedat{1,i}.Data,1)*numTXN(i,5)/size(mousedat{1,i}.Data,1)...      %CFC&Ret3
                               sum(mousedat{1,i}.Data(:,5)==2&mousedat{1,i}.Data(:,8)==2)/size(mousedat{1,i}.Data,1) ];
        chance_overlap{4}(g,1:2)=[numTXN(i,3)/size(mousedat{1,i}.Data,1)*numTXN(i,4)/size(mousedat{1,i}.Data,1)...      %Ret&Ret2
                               sum(mousedat{1,i}.Data(:,6)==2&mousedat{1,i}.Data(:,7)==2)/size(mousedat{1,i}.Data,1) ];
        chance_overlap{5}(g,1:2)=[numTXN(i,3)/size(mousedat{1,i}.Data,1)*numTXN(i,5)/size(mousedat{1,i}.Data,1)...      %Ret&Ret3
                               sum(mousedat{1,i}.Data(:,6)==2&mousedat{1,i}.Data(:,8)==2)/size(mousedat{1,i}.Data,1) ];
        chance_overlap{6}(g,1:2)=[numTXN(i,4)/size(mousedat{1,i}.Data,1)*numTXN(i,5)/size(mousedat{1,i}.Data,1)...      %Ret2&Ret3
                               sum(mousedat{1,i}.Data(:,7)==2&mousedat{1,i}.Data(:,8)==2)/size(mousedat{1,i}.Data,1) ];
        chance_over_tri{1}(g,1:2)=[numTXN(i,4)/size(mousedat{1,i}.Data,1)*numTXN(i,5)/size(mousedat{1,i}.Data,1)...      %Ret2&Ret3
                               sum(mousedat{1,i}.Data(:,7)==2&mousedat{1,i}.Data(:,8)==2)/size(mousedat{1,i}.Data,1) ];
        chance_over_tri{2}(g,1:2)
        chance_over_tri{3}(g,1:2)
        chance_over_tri{4}(g,1:2)
                           g=g+1;
        end
    gg=gg+1;    
    end
end

