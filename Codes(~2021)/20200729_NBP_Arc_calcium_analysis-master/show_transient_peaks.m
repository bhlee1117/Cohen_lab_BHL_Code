function show_transient_peaks(dat,mouse,day,range,timescale,cmap)
figure1 = figure('InvertHardcopy','off','PaperUnits','centimeters',...
    'Color',[1 1 1],...
    'Renderer','painters','position',[100 100 500 20*size(range,2)]);

g=0;
scale=1.5;

for i=range    
rmv=find(sum(dat{mouse}.Cal{i,6+day}>size(dat{mouse}.Cal{i,day},2),2));
dat{mouse}.Cal{i,6+day}(rmv,:)=[];
    if size(dat{mouse}.Cal{i,day},2)>1
plot([timescale:timescale:timescale*dat{mouse}.Cal{i,6+day}(1,1)],dat{mouse}.Cal{i,day}(1,1:dat{mouse}.Cal{i,6+day}(1,1))+g*scale,'color',cmap)
hold all
for j=1:size(dat{mouse}.Cal{i,6+day},1)-1
plot([timescale*(dat{mouse}.Cal{i,6+day}(j,2)):timescale:timescale*(dat{mouse}.Cal{i,6+day}(j+1,1)-1)],dat{mouse}.Cal{i,day}(1,dat{mouse}.Cal{i,6+day}(j,2):dat{mouse}.Cal{i,6+day}(j+1,1)-1)+g*scale,'color',cmap)
end
plot([timescale*(dat{mouse}.Cal{i,6+day}(end,2)+1):timescale:timescale*(size(dat{mouse}.Cal{i,day},2))],dat{mouse}.Cal{i,day}(1,dat{mouse}.Cal{i,6+day}(end,2)+1:end)+g*scale,'color',cmap)
for j=1:size(dat{mouse}.Cal{i,6+day},1)
if dat{mouse}.Cal{i,6+day}(j,1)-1>0
plot([timescale*(dat{mouse}.Cal{i,6+day}(j,1)-1):timescale:timescale*dat{mouse}.Cal{i,6+day}(j,2)],dat{mouse}.Cal{i,day}(1,dat{mouse}.Cal{i,6+day}(j,1)-1:dat{mouse}.Cal{i,6+day}(j,2))+g*scale,'color',cmap)
else
plot([timescale*(dat{mouse}.Cal{i,6+day}(j,1)):timescale:timescale*dat{mouse}.Cal{i,6+day}(j,2)],dat{mouse}.Cal{i,day}(1,dat{mouse}.Cal{i,6+day}(j,1):dat{mouse}.Cal{i,6+day}(j,2))+g*scale,'color',cmap)
end
hold all
 end
%line([0 timescale*size(dat{mouse}.Cal{i,day},2)],[g*scale g*scale],'color',[1 0.4 0])
% if ~isempty(dat{mouse}.Peak{i,day})
% plot(dat{mouse}.Peak{i,day}(:,3),dat{mouse}.Peak{i,day}(:,2)+g*scale*1.01,'g.') 
% end
% plot([1:size(dat{mouse}.S{i,day},2)]*timescale,dat{mouse}.S{i,3+day}/10+g*scale,'r') 
S=dat{mouse}.S{i,day}(find(dat{mouse}.S{i,day}>2.5)); SL=find(dat{mouse}.S{i,day}>2.5);
%plot([SL]*timescale,S/15+g*scale,'co') 
for of=1:size(dat{mouse}.Cal{i,9+day},1)
    if sum(dat{mouse}.Cal{i,9+day}(of,:)==0)==0
plot([dat{mouse}.Cal{i,9+day}(of,1):dat{mouse}.Cal{i,9+day}(of,2)]*1/30,...
      dat{mouse}.Cal{i,day}(1,dat{mouse}.Cal{i,9+day}(of,1):dat{mouse}.Cal{i,9+day}(of,2))+g*scale,'color',cmap);
 hold all
    end
end
g=g+1;
    end
end
xlabel('Time (sec)','FontSize',13,'FontName','arial rounded mt bold')
ylabel('\DeltaF/F','FontSize',13,'FontName','arial rounded mt bold')    
end