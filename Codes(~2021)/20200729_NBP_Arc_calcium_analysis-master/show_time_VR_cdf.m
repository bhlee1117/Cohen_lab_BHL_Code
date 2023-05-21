function show_time_VR_cdf(Full_result,m,d,cell,range)
figure
cmap=jet(100);
C=full(Full_result{m}.Calcium{1,d}.C_df(cell,:));
C=(C-min(C))/(max(C)-min(C))+0.01;
cline([range]*1/30,Full_result{m}.VR{d}(range,2),zeros(size(range,2),1),round(C*99),cmap,2);

end