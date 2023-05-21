function output_data=show_overlap_cdf(Full_result,mouse,comb,cmap)
scale=20;
for m=mouse %mouse

[a aa]=cellfun(@size,Full_result{m}.Calcium_identified);
overlap_list=find(min(aa,[],2)>0); %Detected in All days
overlap_Arc_ref=Full_result{m}.list_identified(overlap_list,1);
datum{m}.Arc=Full_result{m}.Arc_class(overlap_Arc_ref,4:end); % Arc TXN matrix
g=1;
for i=overlap_list' % overlap cells
    for j=1:3 %day
    datum{m}.Cal{g,j}=Full_result{1, 1}.Calcium_identified{i,j}; 
    datum{m}.cal_sigma{g,j}=Full_result{1, 1}.Calcium{1,j}.cal_sigma(Full_result{1, 1}.list_identified(i,j+1),:);
    ini_fin=Full_result{1, 1}.Calcium{1,j}.ini_fin{1,Full_result{1, 1}.list_identified(i,j+1)};
    datum{m}.Peaks{g,j}=peak_finding(full(datum{m}.cal_sigma{g,j}),datum{m}.Cal{g,j},ini_fin);
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
    xlabel('Time (sec)')
    ylabel('Cell ID')
    g=1;
    
for i=1:size(comb,1) % Condition
    comb_list{i}=find(datum{m}.Arc(:,2)==comb(i,1) & datum{m}.Arc(:,4)==comb(i,2) & datum{m}.Arc(:,6)==comb(i,3)) ;

    for j=1:size(comb_list{i},1) % Neuron

        plot([1/30:1/30:1/30*size(datum{m}.cal_sigma{comb_list{i}(j,1),day},2)],datum{m}.cal_sigma{comb_list{i}(j,1),day}+g*scale,'color',cmap(i,:),'linewidth',0.5)
        hold all
        if ~isempty(datum{m}.Peaks{comb_list{i}(j,1),day})
        plot(datum{m}.Peaks{comb_list{i}(j,1),day}(:,3),datum{m}.Peaks{comb_list{i}(j,1),day}(:,1)+g*scale,'color',[1 0 0],'linestyle','none','linewidth',0.5,'marker','o','markersize',3)
        
        
        output_data{m}{g,i}(:,day)=full([sum(datum{m}.Peaks{comb_list{i}(j,1),day}(:,1)) sum(datum{m}.Peaks{comb_list{i}(j,1),day}(:,2))...
                              mean(datum{m}.Peaks{comb_list{i}(j,1),day}(:,1)) mean(datum{m}.Peaks{comb_list{i}(j,1),day}(:,2))...
                               size(datum{m}.Peaks{comb_list{i}(j,1),day},1)]');
                      %  Sum of peaks (sigma) ;  Sum of peaks (dF/F) ;  
                      %  Mean of peaks (dF/F) ;  Mean of peaks (dF/F) ;  number of peaks
        end
        g=g+1;
    end
        
end
  %Find the row that satisfy each condition
    if size(comb,1)==1 % Just one condition
       comb_list{2}=setdiff(1:size(datum,1),comb_list{1})';
    for j=1:size(comb_list{2},1)
        plot([1/30:1/30:1/30*size(Cal{comb_list{2}(j,1),day},2)],Cal{comb_list{2}(j,1),day}+g,'color',cmap(2,:),'linewidth',0.5)
        hold all
        g=g+1;
    end
    end
ylim([0.5 scale*(g+1)])
xlim([0 1/30*size(datum{m}.Cal{comb_list{1}(1,1),day},2)])
end
end
end