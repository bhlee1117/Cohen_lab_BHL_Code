function show_overlap_spike(Full_result,mouse,method,comb,cmap)
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

%cmap=[distinguishable_colors(7)];

day_ti={'Day 1 (Ctx A)','Day 2 (Ctx A)','Day 3 (Ctx B)'};
for day=1:3
    subplot(1,3,day)
    hold all
    title(day_ti{1,day})
    xlabel('Time (sec)')
    ylabel('Cell ID')
    g=1;
    
    
for i=1:size(comb,1)
    comb_list{i}=find(datum(:,2)==comb(i,1) & datum(:,4)==comb(i,2) & datum(:,6)==comb(i,3)) ;

switch method
    case 'line'
    for j=1:size(comb_list{i},1)

        plot([1/30:1/30:1/30*size(Cal{comb_list{i}(j,1),day},2)],Cal{comb_list{i}(j,1),day}+3000*g,'color',cmap(i,:),'linewidth',0.5)
        hold all
        g=g+1;
    end
    case 'image'
    for j=1:size(comb_list{i},1)
        c_df_im(g,1:size(Cal{comb_list{i}(j,1),day},2))=Cal{comb_list{i}(j,1),day};
        hold all
        g=g+1;
    end
        
end
end
    if size(comb,1)==1
       comb_list{2}=setdiff(1:size(datum,1),comb_list{1})';
    for j=1:size(comb_list{2},1)
        plot([1/30:1/30:1/30*size(Cal{comb_list{2}(j,1),day},2)],Cal{comb_list{2}(j,1),day}+g*1000,'color',cmap(2,:),'linewidth',0.5)
        hold all
        g=g+1;
    end
    end
switch method
    case 'image'
imagesc(c_df_im)
axis tight off
[sz1 sz2]=cellfun(@size,comb_list);
cline(zeros(size(comb,1)+1,1),[1 cumsum(sz1)],zeros(size(comb,1)+1,1),[1:1:size(comb,1)],cmap,10);
    case 'line'

ylim([0.5 3000*(g+1)])
xlim([0 1/30*size(Cal{comb_list{1}(1,1),day},2)])
end

end

end
end