function Full_result=combine_var(pthfile)
tic;
for i=1:length(pthfile)
    load([pthfile{1,i}])
end
toc;disp([pthfile{1,2} ' is loaded'])
Full_result.Calcium=ref_result;
for p=1:3
    figure1 = figure('InvertHardcopy','off','PaperUnits','centimeters',...
    'Color',[1 1 1],...
    'Renderer','painters','position',[100 100 1500 1000]); 
    for c1=10:40
plot([1:size(Full_result.Calcium{p}.C_df,2)],Full_result.Calcium{p}.C_df(c1,:)+c1)
hold all    
    end
frame_end=round(ginput(1)); frame_end=frame_end(1,1);
close all
Full_result.Calcium{p}.C_df=[];
Full_result.Calcium{p}.C_or=full(Full_result.Calcium{p}.C_or(:,1:frame_end));
Full_result.Calcium{p}.cal_transient=Full_result.Calcium{p}.cal_transient(:,1:frame_end);
Full_result.Calcium{p}.cal_sigma=Full_result.Calcium{p}.cal_sigma(:,1:frame_end);
[BaseX{p} Full_result.Calcium{p}]=baseline_correction(Full_result.Calcium{p});
save('20temp.mat','BaseX','Full_result')
close all    
    
[Full_result.Calcium{p}.cal_sigma Full_result.Calcium{p}.cal_transient Full_result.Calcium{p}.ini_fin]=...
                 cal_transient_detection(Full_result.Calcium{p},3,0.5,baseX{p},1);
[~,~,Full_result.Calcium{p}.off_focus]=...
                 cal_transient_detection(Full_result.Calcium{p},3,0.5,baseX{p},-1); % detect off focus             

%show_transient(Full_result.Calcium{1, p},[10:40],1,0)

end


%Full_result.options=result{1}.options;
Full_result.regist=regist;
lim=[3000 3000 3150];
for i=1:size(position,2)
int_position{i}=[[1/30:1/30:position{1,i}(end,1)]' interp1q(position{1,i}(:,1),position{1,i}(:,2),[1/30:1/30:position{1,i}(end,1)]')];
step=int_position{i}(2:end,2)-int_position{i}(1:end-1,2);
[odd_list n]=bwlabel(step>-lim(1,i)-50 & step<-100);
for k=1:n
 if size(find(odd_list==k),1)==2
     l=find(odd_list==k);
  int_position{i}(l(2,1),2)=mean([int_position{i}(l(2,1)-1,2) int_position{i}(l(2,1)+1,2)+lim(1,i)])-lim(1,i)*(mean([int_position{i}(l(2,1)-1,2) int_position{i}(l(2,1)+1,2)+lim(1,i)])>(lim(1,i)-2000));
 end
end
end
Full_result.VR=int_position;
Full_result.Calcium_identified=identified_calcium;
Full_result.list_identified=identified_list;
Full_result.Arc_class=Data(:,[1 2 3 5 4 7 6 9 8]);  % arrange 1B 1A 2B 2A 3B 3A
Full_result.cell_list=cell_list;
Full_result.matched_list=match_CS;
Full_result.SP_identified=identified_SpC;
Full_result.spike_identified=identified_spike;
end