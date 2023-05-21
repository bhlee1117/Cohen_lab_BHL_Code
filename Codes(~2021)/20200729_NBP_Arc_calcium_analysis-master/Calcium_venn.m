%function
g=1;
cd=[1 1 1;2 1 1;1 2 1;1 1 2;2 2 1;2 1 2;1 2 2;2 2 2];
clear Frac_Arc
 for Cal_threshold=[0.25]
 figure1 = figure('InvertHardcopy','off','PaperUnits','centimeters',...
     'Color',[1 1 1],'Renderer','painters','position',[100 100 300 300]);
clf('reset')   
for m=1:6
for d=1:3
BB{m}(:,d)=IBI.Brst{m,d}(:,1); end
Arc=dat{m}.Arc_post_ref(:,2:end);
BB_filt=BB{m}(find(sum(isnan(BB{m}),2)==0),:);
Arc_filt{m}=Arc(find(sum(isnan(BB{m}),2)==0),:);
B_class{m}=double(BB_filt>Cal_threshold)+1;
end
Venn_data=plot_venn_Arcal(B_class,1,[2]);
%text(0,0.7,['Threshold :',num2str(Cal_threshold),'Hz'],'Fontname','Arial','Fontsize',10,'HorizontalAlignment','center')
%F=getframe;
%imwrite(uint8(imresize(F.cdata,[255 255])),['\\Neurobiophysics\Byunghun_Lee\Reports\20210209_자료\Brst\' num2str(g) '.jpeg']);
g=g+1;
%close all
for m=1:6
for c=1:size(cd,1)
    clear list
for d=1:3
list(:,d)=B_class{m}(:,d)==cd(c,d);
end
Frac_Arc(:,c,m)=sum(Arc_filt{m}(find(sum(list,2)==3),:)==2,1)';
n(m,c)=sum(sum(list,2)==3);
end
end
 end
M=(sum(Frac_Arc,3)./sum(n,1))'
%end

%%

%function
g=1;
clear Frac_Arc
 for Cal_threshold=2%[0.2:0.1:4]
clf('reset')   
for m=1:6
SP{m}=cellfun(@decon_sum,dat{m}.S(:,4:6));
SP_filt=SP{m}(find(sum(isnan(SP{m}),2)==0),:);
S_class{m}=double(SP_filt>Cal_threshold)+1;
end

Venn_data=plot_venn_Arcal(S_class,1,[2]);
%text(0,0.7,['Threshold :',num2str(Cal_threshold),'Hz'],'Fontname','Arial','Fontsize',10,'HorizontalAlignment','center')
F=getframe;
%imwrite(uint8(imresize(F.cdata,[255 255])),['\\Neurobiophysics\Byunghun_Lee\Reports\20210209_자료\Calcium\' num2str(g) '.jpeg']);
g=g+1;
%close all

for m=1:6
for c=1:size(cd,1)
    clear list
for d=1:3
list(:,d)=S_class{m}(:,d)==cd(c,d);
end
Frac_Arc(:,c,m)=sum(Arc_filt{m}(find(sum(list,2)==3),:)==2,1)';
n(m,c)=sum(sum(list,2)==3);
end
end

 end
M=(sum(Frac_Arc,3)./sum(n,1))'
%end
%%
g=1;
clear Frac_Arc
 for Cal_threshold=0.03%[0.01:0.01:0.25]
clf('reset')   
for m=1:6
Theta{m}=cell2mat(IBI.theta_Brst(m,:));
Theta_filt=Theta{m}(find(sum(isnan(Theta{m}),2)==0),:);
S_class{m}=double(Theta_filt>Cal_threshold)+1;
end

Venn_data=plot_venn_Arcal(S_class,1,[2]);
%text(0,0.7,['Threshold :',num2str(Cal_threshold),'Hz'],'Fontname','Arial','Fontsize',10,'HorizontalAlignment','center')
F=getframe;
%imwrite(uint8(imresize(F.cdata,[255 255])),['\\Neurobiophysics\Byunghun_Lee\Reports\20210209_자료\Theta\' num2str(g) '.jpeg']);
g=g+1;
%close all

for m=1:6
for c=1:size(cd,1)
    clear list
for d=1:3
list(:,d)=S_class{m}(:,d)==cd(c,d);
end
Frac_Arc(:,c,m)=sum(Arc_filt{m}(find(sum(list,2)==3),:)==2,1)';
n(m,c)=sum(sum(list,2)==3);
end
end

 end
M=(sum(Frac_Arc,3)./sum(n,1))'