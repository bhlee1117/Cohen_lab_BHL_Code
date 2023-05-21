clear 

[fnm pth]=uigetfile('*.mat','Multiselect','on');

%%
for i=1:length(fnm)
    color=rand(1,3);
    load([pth fnm{1,i}])

end

%%

for i=1:length(fnm)
    color=rand(1,3);
    load([pth fnm{1,i}])
    grid{1}=randn(size(single,1),1)*0.1+1;
    grid{2}=randn(size(TXN,1),1)*0.1+2;
I{i,1}=intensity(single,1);
I{i,2}=intensity(TXN,1);
M(i,1:3)=[mean(intensity(single,1)) median(intensity(single,1)) std(intensity(single,1))];
N_TXN{i,1}=intensity(TXN,1)/mean(intensity(single,1));
N_TXN{i,2}=intensity(TXN,1)/median(intensity(single,1));
N_TXN{i,1}(N_TXN{i,1}>500)=NaN;
N_TXN{i,2}(N_TXN{i,2}>500)=NaN;
plot([grid{1};grid{2}],[I{i,1};I{i,2}],'color',color,'marker','.','linestyle','none')
hold all

end

set(gca,'YScale','log','XTick',[1 2],'XTickLabel',{'Single','TXN'})

%%
g=1;
for j=1:2
for i=1:2
 rep=cell2mat(N_TXN([j:2:22+j],i));
 size(rep)
grid{2*(j-1)+i}=randn(sum(~isnan(rep)),1)*0.1+2*(j-1)+i;
plot(grid{2*(j-1)+i},rep(find(~isnan(rep))),'k.')
hold all
end
end
errorbar([1:4],[mean(cell2mat(N_TXN([1:2:23],1)),'omitnan') mean(cell2mat(N_TXN([1:2:23],2)),'omitnan') mean(cell2mat(N_TXN([2:2:24],1)),'omitnan') mean(cell2mat(N_TXN([2:2:24],2)),'omitnan')],...
[std(cell2mat(N_TXN([1:2:23],1)),0,1,'omitnan') std(cell2mat(N_TXN([1:2:23],2)),0,1,'omitnan') std(cell2mat(N_TXN([2:2:24],1)),0,1,'omitnan') std(cell2mat(N_TXN([2:2:24],2)),0,1,'omitnan')],'r','linestyle','none')
% for i=1:size(N_TXN)
%     grid{1}=randn(size(N_TXN{i,2},1),1)*0.1+1;
% %     plot(grid{1},N_TXN{i,2},'k.')
% %     hold all
% %     errorbar(i,mean(N_TXN{i,2}),std(N_TXN{i,2}),'r')
%     
% end
set(gca,'XTick',[1 2 3 4],'XTickLabel',{'CDS (/mean)','CDS (/median)','PP7 (/mean)','PP7 (/median)'})
ylim([0 200])
mRNA_per_TXN=cell2mat(N_TXN);
xtickangle(45)
% 
% histogram(mRNA_per_TXN(:,1),[0:3:180])
% hold all
% histogram(mRNA_per_TXN(:,2),[0:3:180])

%%
