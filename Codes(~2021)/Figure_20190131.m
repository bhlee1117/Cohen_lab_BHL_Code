%% Load mouse data

clear
[fnm pth]=uigetfile('.mat','Select the analyzed mouse data','Multiselect','on');
% Group 1 = HC CFC Ret 
% Group 2 = HC CFC Ret Rem
% Group 3 = HC CFC Ret Ret2 Ret3 Rem
% Group 4 = HC CFC RetB Ret2A
% Group 5 = RSC
%%
if ischar(fnm)
    fnm={fnm};
end
lg=[];
for i=1:length(fnm)
    mousedat{i}=importdata([pth fnm{1,i}]);
    g=1;
    if size((mousedat{1,i}.Data),1)>1
    for j=1:size(mousedat{1,i}.Data,1)
        fracTXN{1}(i,j)=sum(mousedat{1,i}.Data{j,1}(:,4)==2)/sum(mousedat{1,i}.Data{j,1}(:,4)~=3); %columns for HC Tr Ret... Rows for each mouse
        fracTXN{2}(i,j)=sum(mousedat{1,i}.Data{j,2}(:,4)==2)/sum(mousedat{1,i}.Data{j,2}(:,4)~=3);
    end
    else
        for j=4:size(mousedat{1,i}.Data{1,1},2)
        fracTXN{1}(i,j-3)=sum(mousedat{1,i}.Data{1,1}(:,j)==2)/sum(mousedat{1,i}.Data{1,1}(:,j)~=3);
        fracTXN{2}(i,j-3)=sum(mousedat{1,i}.Data{1,2}(:,j)==2)/sum(mousedat{1,i}.Data{1,2}(:,j)~=3);
        
%         for l=1:2
%         for k=1:size(mousedat{1,i}.Data{1,l},2)-j
%         overlapTXN{l}(i,g)=sum(mousedat{1,i}.Data{1,l}(:,j)==2 & mousedat{1,i}.Data{1,l}(:,j+k) ==2);
%         g=g+1;
%         end
%         end
        end
        %cond_prob(i,:)=[overlapTXN(i,1)/numTXN(i,2) overlapTXN(i,2)/numTXN(i,3) overlapTXN(i,3)/numTXN(i,4) overlapTXN(i,4)/numTXN(i,3) overlapTXN(i,5)/numTXN(i,2) overlapTXN(i,6)/numTXN(i,3)];% HC&CFC/CFC, HC&Ret/Ret, HC&Rem/Rem, CFC&Ret/Ret, CFC&Rem/CFC, Ret&Rem/Ret
    end
    f=table2array(struct2table(mousedat{1,i}.Freeze));
 Freeze(i,1:size(f,2))=f;
 %lg=[lg; cell(mousedat{1,i}.Mouse)];
 lg{i}=mousedat{1,i}.Mouse;
end
fracTXN{1,1}(fracTXN{1,1}==0)=NaN;
fracTXN{1,2}(fracTXN{1,2}==0)=NaN;
fracTXN{1,1}(isnan(fracTXN{1,1}(:,1)),1)=0;
fracTXN{1,2}(isnan(fracTXN{1,1}(:,1)),1)=0;

Freeze(Freeze==0)=NaN;
if size(fracTXN{1,1},2)<9
    fracTXN{1,1}(:,size(fracTXN{1,1},2)+1:9)=NaN;
    fracTXN{1,2}(:,size(fracTXN{1,2},2)+1:9)=NaN;
    Freeze(:,size(Freeze,2)+1:9)=NaN;
else
for i=1:size(mousedat,2)
    if ~isnan(fracTXN{1,1}(i,9))
        a=isnan(Freeze(i,1:8));
        Freeze(i,9)=Freeze(i,max(find(a==0)));
Freeze(i,max(find(a==0)))=NaN;
    else
        Freeze(i,9)=NaN;
    end
end
end
Freeze(isnan(Freeze(:,1)),1)=0;

%cond_prob(cond_prob==0)=NaN;
%% 'Fraction of TXN Cell (%)'
for l=1:2
    figure(l)

   p=[];
for i=1:size(fracTXN{1,l},1)
   clear color
%    color=zeros(1,3);
%     color(1,i)=1;
    ptm=plot([1:1:6],[fracTXN{1,l}(:,1:5) fracTXN{1,l}(:,9)]*100,'marker','.','markersize',20,'color',[0.6 0.6 0.6]); %'color',color,
    p=[p ptm];
    hold all
end
%legend(p,lg{1,1:size(fracTXN{1,l},1)})
plot([1:1:6],mean([fracTXN{1,l}(:,1:5) fracTXN{1,l}(:,9)],'omitnan')*100,'r','linewidth',4)%'FaceColor',[0.6 0.6 0.6],'BarWidth',0.8)
set(gca,'FontSize',8,'LineWidth',1,'XTick',[1:1:6],'XTickLabel',...
    {'HC','CFC','1^s^t Ret.','2^n^d Ret.','3^r^d Ret.','Remote'});
   ylabel('Fraction of TXN Cell (%)')
   xlim([0.5 6.5])
   hold all
   errorbar([1:1:6],mean([fracTXN{1,l}(:,1:5) fracTXN{1,l}(:,9)],1,'omitnan')*100,std([fracTXN{1,l}(:,1:5) fracTXN{1,l}(:,9)],0,1,'omitnan')/sqrt(size(fracTXN{1,l}(:,1:3),1))*100,'r','linewidth',4)
end
%% 'Freezing (%)'

bar([1:1:5],mean([Freeze(:,1) Freeze(:,3:5) Freeze(:,9)],1,'omitnan')*100,'FaceColor',[0.6 0.6 0.6],'BarWidth',0.8)
set(gca,'FontSize',8,'LineWidth',1,'XTick',[1:1:5],'XTickLabel',...
    {'HC','1^s^t Ret.','2^n^d Ret.','3^r^d Ret.','Remote'});
   ylabel('Freezing (%)')
   xlim([0.5 5.5])
   ylim([0 80])
   hold all
   errorbar([1:1:5],mean([Freeze(:,1) Freeze(:,3:5) Freeze(:,9)],1,'omitnan')*100,std([Freeze(:,1) Freeze(:,3:5) Freeze(:,9)],0,1,'omitnan')/sqrt(size(fracTXN{1,1},1))*100,'k','linestyle','none')
   p=[];
for i=1:size(Freeze,1)
   clear color
   color=zeros(1,3);
    color(1,i)=1;
    ptm=plot([1:1:5],[Freeze(i,1) Freeze(i,3:5) Freeze(i,9)]*100,'marker','.','linestyle','none','markersize',30);%,'color',color
    p=[p ptm];
    hold all
end
legend(p,lg{1,1:size(fracTXN{1,l},1)})
%% 'Conditional Probability (%)'
clear omitna_data overlap over_name overlapTXN
% HC=1, CFC=2; Ret=3; Ret2=4; Ret3=5; Rem=9;
sess_name={'HC','CFC','Ret','Ret2','Ret3','Ret4','Ret5','Ret6','Rem.'};
overlap=[1 2;1 3;2 3;3 4;3 5;4 5;1 9;2 9;3 9];
for i=1:size(overlap,1)
    over_name{i,1}=[char(sess_name(overlap(i,1))) ' ¡û ' char(sess_name(overlap(i,2)))];
end
for l=1:2  % Make NA cell omited data
for i=1:size(mousedat,2)
       g=1;
    if size((mousedat{1,i}.Data),1)>1
    else
 
    for j=1:size(mousedat{1,i}.Data{1,l},1)
        if sum(mousedat{1,i}.Data{1,l}(j,:)==3)==0
    omitna_data{i,l}(g,:)=mousedat{1,i}.Data{1,l}(j,:);
    g=g+1;
        end
    end
    end
end
end
for l=1:2
for i=1:size(omitna_data,1)
    if ~isempty(omitna_data{i,l})
        for j=1:size(overlap,1)
            if size(omitna_data{i,l},2)>=overlap(j,1)+3 && size(omitna_data{i,l},2)>=overlap(j,2)+3 ...
            && sum(omitna_data{i,l}(:,overlap(j,1)+3)==0)==0 && sum(omitna_data{i,l}(:,overlap(j,2)+3)==0)==0
            overlapTXN{l}(i,j)=sum(omitna_data{i,l}(:,overlap(j,1)+3)==2 & omitna_data{i,l}(:,overlap(j,2)+3)==2)/size(omitna_data{i,l},1);
            else
                overlapTXN{l}(i,j)=NaN;
            end
        end
    else
        overlapTXN{l}(i,1:size(overlap,1))=NaN;
    end
end
end
    
% plot
clear p
for i=1:2
    figure(i)
    for j=1:size(overlapTXN{1,i},1)
        p{j}=plot([1:1:size(overlapTXN{1,i},2)],overlapTXN{1,i}(j,:),'.','color',[0.7 0.7 0.7],'markersize',20);
        hold all
    end
    
    bar([1:1:size(overlap,1)],mean(overlapTXN{1,i},1,'omitnan')*100,'FaceColor',[0.7 0.7 0.7],'BarWidth',0.8)
set(gca,'FontSize',8,'LineWidth',1,'XTick',[1:1:size(overlap,1)],'XTickLabel',...
    over_name);
xtickangle(45);

     errorbar([1:1:size(overlap,1)],mean(overlapTXN{1,i},1,'omitnan')*100,std(overlapTXN{1,i},0,1,'omitnan')/sqrt(size(overlapTXN{1,i},1))*100,'r','linewidth',4,'linestyle','none')
end
%%
clear overlap over_name overlapTXN overlapTXN_freez
sess_name={'HC','CFC','Ret','Ret2','Ret3','Ret4','Ret5','Ret6','Rem.'};
overlap=[2 3;2 4;2 5;2 9;3 9];
for i=1:size(overlap,1)
    over_name{i,1}=[char(sess_name(overlap(i,1))) ' ¡û ' char(sess_name(overlap(i,2)))];
end
for l=1:2
for i=1:size(omitna_data,1)
    if ~isempty(omitna_data{i,l})
        for j=1:size(overlap,1)
            if size(omitna_data{i,l},2)>=overlap(j,1)+3 && size(omitna_data{i,l},2)>=overlap(j,2)+3 ...
            && sum(omitna_data{i,l}(:,overlap(j,1)+3)==0)==0 && sum(omitna_data{i,l}(:,overlap(j,2)+3)==0)==0
            overlapTXN_freez{l}(i,j)=sum(omitna_data{i,l}(:,overlap(j,1)+3)==2 & omitna_data{i,l}(:,overlap(j,2)+3)==2)/size(omitna_data{i,l},1)*Freeze(i,overlap(j,2));
            else
                overlapTXN_freez{l}(i,j)=NaN;
            end
        end
    else
        overlapTXN_freez{l}(i,1:size(overlap,1))=NaN;
    end
end
end

% plot

for i=1:2
    figure(i)
    for j=1:size(overlapTXN_freez{1,i},1)
        p{j}=plot([1:1:size(overlapTXN_freez{1,i},2)],overlapTXN_freez{1,i}(j,:),'.','color',[0.7 0.7 0.7],'markersize',20);
        hold all
    end
    
    bar([1:1:size(overlap,1)],mean(overlapTXN_freez{1,i},'omitnan')*100,'FaceColor',[0.7 0.7 0.7],'BarWidth',0.8)
set(gca,'FontSize',8,'LineWidth',1,'XTick',[1:1:size(overlap,1)],'XTickLabel',...
    over_name);
xtickangle(45);

     errorbar([1:1:size(overlap,1)],mean(overlapTXN_freez{1,i},'omitnan')*100,std(overlapTXN_freez{1,i},'omitnan')/sqrt(size(overlapTXN_freez{1,i},1))*100,'r','linewidth',4,'linestyle','none')
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

