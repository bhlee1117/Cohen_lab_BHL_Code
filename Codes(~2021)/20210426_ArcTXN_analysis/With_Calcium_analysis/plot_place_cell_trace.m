function plot_place_cell_trace(Full_result,mouse,neuron,min_track)

for m=mouse %mouse
figure1 = figure('InvertHardcopy','off','PaperUnits','centimeters',...
    'Color',[1 1 1],...
    'Renderer','painters','position',[100 100 2500 400]);
[a aa]=cellfun(@size,Full_result{m}.Calcium_identified);
overlap_list=find(min(aa,[],2)>0); %Detected in All days
datum=Full_result{m}.Arc_class(Full_result{m}.list_identified(overlap_list,1),4:end);
g=1;
for j=1:3%day
      subplot(1,3,j)
    for i=overlap_list(neuron,1) 
    Cal{g,j}=Full_result{1, 1}.Calcium_identified{i,j};  
    step=Full_result{1, 1}.VR{1, j}(2:end,1)-Full_result{1, 1}.VR{1, j}(1:end-1,1);
    l=find(abs(step)>abs(min_track));
    for lap=1:size(l,1)-1
    P{lap,j}=[Full_result{1, 1}.VR{1, j}(l(lap,1)+1:l(lap+1,1),1) Full_result{1, 1}.Calcium_identified{i,j}(1,l(lap,1)+1:l(lap+1,1))'];
    plot(P{lap,j}(:,1),P{lap,j}(:,2))
        hold all
    end
    end
end
end