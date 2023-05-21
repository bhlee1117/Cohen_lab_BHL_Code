function [BaseX ref_Cal]=baseline_correction(Cal)
ref_Cal=Cal;
   
for i=1:size(Cal.C_or,1)
    figure1 = figure('InvertHardcopy','off','PaperUnits','centimeters',...
    'Color',[1 1 1],...
    'Renderer','painters','position',[100 100 1500 300]); 
plot([1:1:size(Cal.C_or(i,:),2)],Cal.C_or(i,:))
[x(i,:) y]=ginput(2);
line([x(i,1) x(i,1)],[0 max(Cal.C_or(i,:))],'color','r')
line([x(i,2) x(i,2)],[0 max(Cal.C_or(i,:))],'color','r')
[x2(i,:) y]=ginput(2);
m=mean(full(Cal.C_or(i,round(x(i,1)):round(x(i,2)))));
ref_Cal.C_df(i,:)=(Cal.C_or(i,:)-m)/m;
figure(2)
clf('reset')
plot(ref_Cal.C_df(i,:))
m2=mean(ref_Cal.C_df(i,round(x2(i,1)):round(x2(i,2))));
aa=interp1([mean(x(i,:)) mean(x2(i,:))],[0 m2],[1:1:size(Cal.C_or(i,:),2)],'linear','extrap');
ref_Cal.C_df(i,:)=ref_Cal.C_df(i,:)-aa;
hold all
plot(ref_Cal.C_df(i,:))
plot(aa,'c')
plot([1 size(ref_Cal.C_df,2)],[0 0],'r')
close(figure(1))
end
BaseX=[x x2];
close all
end