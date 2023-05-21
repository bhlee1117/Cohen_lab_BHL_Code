function hh=plot_percent_TXN(FR,bin_number,cat,t_mouse)
switch cat
    case 1 
        xlab='Mean \DeltaF/F';
    case 2 
        xlab='\Sigma peak';
    case 3 
        xlab='\Sigma transient';
    case 4 
        xlab='Ca^2^+ event rate (Hz)';
    case 5 
        xlab='Ca^2^+ event rate (Hz)';    
end

if t_mouse
max_rate=max(cell2mat(FR(:)))*1.1;
for m=1:size(FR,1)
for i=1:2
hh{m,i}=histcounts(FR{m,i},[0:max_rate/(bin_number):max_rate]);
end
datum2(m,:)=(hh{m,2})./(hh{m,1}+hh{m,2});
end
figure1 = figure('InvertHardcopy','off','PaperUnits','centimeters',...
    'Color',[1 1 1],...
    'Renderer','painters','position',[100 100 300 300]);
lineProps.col{1}=[1 0 0];
mseb([0:max_rate/(bin_number):max_rate/(bin_number)*(size(hh{1,1},2)-1)],...
         mean(datum2,1,'omitnan'),std(datum2,0,1,'omitnan')/sqrt(size(datum2,1)),lineProps,1);
xlim([0 max_rate])
ylim([0 1.1])
set(gca,'LineWidth',2,'FontSize',13,...
    'FontName','arial rounded mt bold')
xlabel(xlab,'LineWidth',2,'FontSize',13,...
    'FontName','arial rounded mt bold')
ylabel('Fraction of TXN neurons','LineWidth',2,'FontSize',13,...
    'FontName','arial rounded mt bold')
else
    figure1 = figure('InvertHardcopy','off','PaperUnits','centimeters',...
    'Color',[1 1 1],...
    'Renderer','painters','position',[100 100 300 300]);
m=max(cell2mat(FR'));
for i=1:2
hh{i}=histcounts(FR{i},[0:m/(bin_number):m]);
end
plot([0:m/(bin_number):m/(bin_number)*(size(hh{1},2)-1)],hh{2}./(hh{1}+hh{2}),'LineWidth',2)
xlim([0 m])
ylim([0 1.1])
set(gca,'LineWidth',2,'FontSize',13,...
    'FontName','arial rounded mt bold')
xlabel(xlab,'LineWidth',2,'FontSize',13,...
    'FontName','arial rounded mt bold')
ylabel('Fraction of TXN neurons','LineWidth',2,'FontSize',13,...
    'FontName','arial rounded mt bold')
end
end