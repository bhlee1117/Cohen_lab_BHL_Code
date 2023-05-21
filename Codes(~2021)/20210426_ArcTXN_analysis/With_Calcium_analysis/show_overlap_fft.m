function show_overlap_fft(Full_result,mouse,method,comb,cmap)
for m=mouse %mouse

[a aa]=cellfun(@size,Full_result{m}.Calcium_identified);
overlap_list=find(min(aa,[],2)>0); %Detected in All days
datum=Full_result{m}.Arc_class(Full_result{m}.list_identified(overlap_list,1),4:end);
g=1;
for i=overlap_list'
    for j=1:3 %day
    Cal{g,j}=Full_result{1, 1}.spike_identified{i,j};  
    end
    g=g+1;
end
figure1 = figure('InvertHardcopy','off','PaperUnits','centimeters',...
    'Color',[1 1 1],...
    'Renderer','painters','position',[100 100 2100 900]);

day_ti={'Day 1 (Ctx A)','Day 2 (Ctx A)','Day 3 (Ctx B)'};
for day=1:3
    subplot(1,3,day)
    hold all
    title(day_ti{1,day})
    xlabel('Frequency (Hz)')
    ylabel('Cell ID')
    g=1;
for i=1:size(comb,1) %conditions
    comb_list{i}=find(datum(:,2)==comb(i,1) & datum(:,4)==comb(i,2) & datum(:,6)==comb(i,3)) ;
switch method
    case 'line'
    for j=1:size(comb_list{i},1) %cell
        [time_vec FT]=run_fft([0:1/30:1/30*(size(Cal{comb_list{i}(j,1),day},2)-1)],full(Cal{comb_list{i}(j,1),day}));
        plot(time_vec,FT+g*20,'color',cmap(i,:),'linewidth',0.5)
        hold all
        g=g+1;
    end
    case 'image'
    for j=1:size(comb_list{i},1)
        c_df_im(g,1:size(datum{m}.cal{comb_list{i}(j,1),1},2))=datum{m}.cal{comb_list{i}(j,1),1};
        hold all
        g=g+1;
    end
        
end
end
ylim([0 800])
xlim([0 15])
% set(gca,'FontSize',8,'LineWidth',1,'YTick',[0.005:0.025:g*0.005],'YTickLabel',...
%  [1:5:g],'FontName','arial rounded mt bold','FontSize',13,'LineWidth',2);

end
switch method
    case 'image'
imagesc(c_df_im)
axis tight off
[sz1 sz2]=cellfun(@size,comb_list);
cline([0 0 0 0 0],[1 cumsum(sz1)],zeros(5,1),[1 2 3 4],cmap,10);
    case 'line'

end
end
end
function [f P1]=run_fft(time_vector,signal)
 
L = size(time_vector,2);             % Length of signal

Y=fft(signal);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = 1/time_vector(1,2)*(0:(L/2))/L;
end