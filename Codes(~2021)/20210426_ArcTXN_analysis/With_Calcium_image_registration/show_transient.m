function show_transient(ref_result,range,timescale,sw)
figure1 = figure('InvertHardcopy','off','PaperUnits','centimeters',...
    'Color',[1 1 1],...
    'Renderer','painters','position',[100 100 500 20*size(range,2)]);

g=1;

if sw==1
    scale=20;
for i=range
l=bwlabel(ref_result.cal_transient(i,:));
plot([timescale:timescale:timescale*ref_result.ini_fin{i}(1,1)],ref_result.cal_sigma(i,1:ref_result.ini_fin{i}(1,1))+g*scale,'k')
hold all
for j=1:size(ref_result.ini_fin{i},1)-1
plot([timescale*(ref_result.ini_fin{i}(j,2)+1):timescale:timescale*(ref_result.ini_fin{i}(j+1,1)-1)],ref_result.cal_sigma(i,ref_result.ini_fin{i}(j,2)+1:ref_result.ini_fin{i}(j+1,1)-1)+g*scale,'k')
end
plot([timescale*(ref_result.ini_fin{i}(end,2)+1):timescale:timescale*(size(ref_result.cal_transient(i,:),2))],ref_result.cal_sigma(i,ref_result.ini_fin{i}(end,2)+1:end)+g*scale,'k')
for j=1:size(ref_result.ini_fin{i},1)
if ref_result.ini_fin{i}(j,1)-1>0
plot([timescale*(ref_result.ini_fin{i}(j,1)-1):timescale:timescale*ref_result.ini_fin{i}(j,2)],ref_result.cal_sigma(i,ref_result.ini_fin{i}(j,1)-1:ref_result.ini_fin{i}(j,2))+g*scale,'r')
else
plot([timescale*(ref_result.ini_fin{i}(j,1)):timescale:timescale*ref_result.ini_fin{i}(j,2)],ref_result.cal_sigma(i,ref_result.ini_fin{i}(j,1):ref_result.ini_fin{i}(j,2))+g*scale,'r')
end
hold all
end
g=g+1;
end
xlabel('Time (sec)','FontSize',13,'FontName','arial rounded mt bold')
ylabel('\sigma','FontSize',13,'FontName','arial rounded mt bold')    
    
else
    scale=2;
for i=range
l=bwlabel(ref_result.cal_transient(i,:));
plot([timescale:timescale:timescale*ref_result.ini_fin{i}(1,1)],ref_result.C_df(i,1:ref_result.ini_fin{i}(1,1))+g*scale,'k')
hold all
for j=1:size(ref_result.ini_fin{i},1)-1
plot([timescale*(ref_result.ini_fin{i}(j,2)+1):timescale:timescale*(ref_result.ini_fin{i}(j+1,1)-1)],ref_result.C_df(i,ref_result.ini_fin{i}(j,2)+1:ref_result.ini_fin{i}(j+1,1)-1)+g*scale,'k')
end
plot([timescale*(ref_result.ini_fin{i}(end,2)+1):timescale:timescale*(size(ref_result.cal_transient(i,:),2))],ref_result.C_df(i,ref_result.ini_fin{i}(end,2)+1:end)+g*scale,'k')
for j=1:size(ref_result.ini_fin{i},1)
if ref_result.ini_fin{i}(j,1)-1>0
plot([timescale*(ref_result.ini_fin{i}(j,1)-1):timescale:timescale*ref_result.ini_fin{i}(j,2)],ref_result.C_df(i,ref_result.ini_fin{i}(j,1)-1:ref_result.ini_fin{i}(j,2))+g*scale,'r')
else
plot([timescale*(ref_result.ini_fin{i}(j,1)):timescale:timescale*ref_result.ini_fin{i}(j,2)],ref_result.C_df(i,ref_result.ini_fin{i}(j,1):ref_result.ini_fin{i}(j,2))+g*scale,'r')
end
hold all
end
g=g+1;
end
xlabel('Time (sec)','FontSize',13,'FontName','arial rounded mt bold')
ylabel('\sigma','FontSize',13,'FontName','arial rounded mt bold')
end
end