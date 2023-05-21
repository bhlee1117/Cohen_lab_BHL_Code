function [norm_d polar_d]=plot_3d_FR(dat,cat,post_ref,doi,bin)


% groups={[1 1 2 1 3 1],[1 2 2 1 3 1],[1 2 2 2 3 1],[1 1 2 2 3 1],[1 1 2 2 3 2],[1 1 2 1 3 2]...
%     ,[1 2 2 1 3 2],[1 2 2 2 3 2]};
 groups={[1 1 2 1 3 1],[ 1 1 2 2 3 1],[1 1 2 2 3 2],[1 1 2 1 3 2]};
%tickss={'---','+--','++-','-+-','-++','--+','+-+','+++'};
tickss={'---','-+-','-++','--+'};
for m=1:size(dat,2) %mouse
    switch cat
        case 1 %Mean dFF
            ylab='Mean \DeltaF/F';
            for c1=1:size(dat{m}.Cal,1)
                for d1=1:3
                    if size(dat{m}.Cal{c1,d1},2)<2
                        SP{m}(c1,d1)=NaN;
                    else
                        SP{m}(c1,d1)=mean(dat{m}.Cal{c1,3+d1},'omitnan');
                    end
                end
            end
        case 2 % Sum peak
            ylab='\Sigma Peaks';
            SP{m}=cellfun(@extrack_sigma_peak,dat{m}.Peak);
        case 3 % Sum Transient
            ylab='\Sigma Transient (\DeltaF/F)';
            for c1=1:size(dat{m}.Cal,1)
                for d1=1:3
                    if size(dat{m}.Cal{c1,d1},2)<2
                        SP{m}(c1,d1)=NaN;
                    else
                        tr_area=0;
                        dat{m}.Cal{c1,6+d1}(find(dat{m}.Cal{c1,6+d1}(:,1)>size(dat{m}.Cal{c1,3+d1},2)),:)=[];
                        for tr=1:size(dat{m}.Cal{c1,6+d1},1)
                            if dat{m}.Cal{c1,6+d1}(tr,2)>size(dat{m}.Cal{c1,3+d1},2)
                                dat{m}.Cal{c1,6+d1}(tr,2)=size(dat{m}.Cal{c1,3+d1},2);
                            end
                            if sum(dat{m}.Cal{c1,6+d1}(tr,:)==0)==0
                                tr_area=tr_area+sum(dat{m}.Cal{c1,3+d1}(dat{m}.Cal{c1,6+d1}(tr,1):dat{m}.Cal{c1,6+d1}(tr,2)),'omitnan');
                            end
                        end
                        SP{m}(c1,d1)=tr_area;
                    end
                end
            end
        case 4 % Number of peaks
            ylab='Ca^2^+ event rate (Hz)';
            [cellnumber length]=cellfun(@size,dat{m}.Cal(:,1:3));
            length=max(length,1);
            [SP{m} rw]=cellfun(@size,dat{m}.Peak);
            SP{m}(rw==1)=NaN;
            SP{m}=SP{m}./length*30;
        case 5 % Number of big transients
             ylab='Ca^2^+ event rate (Hz)';
             %SP{m}=cellfun(@decon_sum,dat{m}.S(:,1:3));
             SP{m}=cellfun(@decon_sum,dat{m}.S(:,4:6));
             [cw rw]=cellfun(@size,dat{m}.Peak);
             SP{m}(rw==1)=NaN;
    end
    if post_ref
        Arc_class{m}=[dat{m}.Arc_post_ref(:,2:end)];
    else
        Arc_class{m}=[dat{m}.Arc(:,2:end)];
    end
    for g=1:size(groups,2)
        for cd=1:size(groups{g},2)/2
            list{m,g}(:,cd)=Arc_class{m}(:,groups{g}(1,2*cd-1))==groups{g}(1,2*cd);
        end
        list{m,g}=sum(list{m,g},2)==size(groups{g},2)/2; % cells that match the conditions
        FR{m,g}=[SP{m}(find(list{m,g}),:)];
        % Matrix FR : cell array consist of n X 3 matrix, each column represents
        % days.
        % Quantity of interest is arranged in mouse and groups.
        
    end
end
%p_list=[1 2;1 3;2 3;1 4;2 4;3 4];

%%%%% Data arrange done %%%%%
%%%%% Significance check & plotting start %%%%%
%cmap=abs(jet(8)-0.2)+0.2;
%cmap=[0 0 0; jet(size(groups,2)-1)];
cmap=[0 0 0; distinguishable_colors(size(groups,2)-1)];
p=[];
figure1 = figure('InvertHardcopy','off','PaperUnits','centimeters',...
    'Color',[1 1 1],'Renderer','painters','position',[100 100 300 300]);
for g=1:size(FR,2)
    d=cell2mat(FR(:,g));
    norm_d{g}=d./sqrt(sum(d.^2,2));
    plot3(norm_d{g}(:,1),norm_d{g}(:,2),norm_d{g}(:,3),'marker','o','linestyle','none','MarkerFaceColor',cmap(g,:),...
         'MarkerEdgeColor','none','markersize',10);
    %pm.FaceColorData=[cmap(g,:) 50];
    hold all
    md=[0 0 0;mean(norm_d{g},1,'omitnan')];
    p=[p plot3(md(:,1),md(:,2),md(:,3),'linewidth',2,'color',cmap(g,:))];
end
xlabel('Day 1 (A1)','LineWidth',2,'FontSize',13,...
    'FontName','arial rounded mt bold')
ylabel('Day 2 (A2)','LineWidth',2,'FontSize',13,...
    'FontName','arial rounded mt bold')
zlabel('Day 3 (B)','LineWidth',2,'FontSize',13,...
    'FontName','arial rounded mt bold')
axis equal
lgd=legend(p,tickss,'Location','North');
lgd.NumColumns = 4;
grid on

%%
figure1 = figure('InvertHardcopy','off','PaperUnits','centimeters',...
    'Color',[1 1 1],'Renderer','painters','position',[100 100 300 300]);
%polar plots

for g=1:size(FR,2)
    d=cell2mat(FR(:,g));
    if size(doi,2)==3 % if comparing three days
    norm_d{g}=d(:,doi)./sqrt(sum(d(:,doi).^2,2));    
    im_d=norm_d{g}(:,1)+norm_d{g}(:,2)*exp(i*pi*2/3)+norm_d{g}(:,3)*exp(i*pi*4/3);
     polar_d{g}=[angle(im_d) abs(im_d)]; polar_d{g}(find(polar_d{g}(:,1)<0),1)=polar_d{g}(find(polar_d{g}(:,1)<0),1)+2*pi;
    bin_theta=floor(polar_d{g}(:,1)/bin); polar_weightd=[];
    %polarplot(polar_d(:,1),polar_d(:,2),'marker','.','markersize',12,'color',cmap(g,:),'linestyle','none')
    for b=unique(bin_theta(find(~isnan(bin_theta)),1))'
    polar_weightd=[polar_weightd; [b sum(polar_d{g}(find(bin_theta==b),2)) size(find(bin_theta==b),1)]];
    end
    polar_weightd(:,2:3)=polar_weightd(:,2:3)/sum(~isnan(norm_d{g}(:,1)));
    polarplot([polar_weightd(:,1);polar_weightd(1,1)]*bin+bin/2,...
              [polar_weightd(:,2);polar_weightd(1,2)],'color',cmap(g,:),'linewidth',1.5)
    hold all
    else %if comparing just two days
   norm_d{g}=d(:,doi)./sqrt(sum(d(:,doi),2));
   %plot(norm_d{g}(:,1),norm_d{g}(:,2),'marker','.','linestyle','none','color',cmap(g,:),'markersize',15)
   S=std(norm_d{g},0,1,'omitnan')./sqrt(sum(~isnan(norm_d{g}(:,1))));
   M=mean(norm_d{g},1,'omitnan');
   errorbar(M(1,1),M(1,2),S(1,2),S(1,2),S(1,1),S(1,1),...
                  'marker','o','linestyle','none','color',cmap(g,:),'markersize',10,'linewidth',2)
   hold all
    end
end

if size(doi,2)==3
set(gca,'FontSize',8,'LineWidth',1,'ThetaTick',[0 120 240],'ThetaTickLabel',...
    {'A1','A2','B'},'RGrid','off','RTick',[0:0.2],'FontName','arial rounded mt bold','FontSize',13,'LineWidth',2);
else
    d_name={'A1','A2','B'};
% set(gca,'FontSize',8,'LineWidth',1,'ThetaTick',[0 90],'ThetaTickLabel',...
%     d_name(doi),'RGrid','off','FontName','arial rounded mt bold','FontSize',13,'LineWidth',2);    
 set(gca,'FontName','arial rounded mt bold','FontSize',13,'LineWidth',2); 
 axis equal tight
%  zoom(0.5)
 T = axis;
 line([0 2],[0 2],'color',[0.6 0.6 0.6],'linestyle','--','linewidth',2)
 axis(T);
 xlabel(d_name(doi(1)),'FontName','arial rounded mt bold','FontSize',13)
 ylabel(d_name(doi(2)),'FontName','arial rounded mt bold','FontSize',13)
 polar_d=[];
end
 zoom(0.5)
lgd=legend(tickss,'Location','North');
lgd.NumColumns = 2;
%end
% function length=timelength(tr)
% length=sum(~isnan(tr));
% end
% function newpeak=rmvnanPeak(Peak,isnanPeriod)
% newpeak=Peak;
% for c=1:size(Peak,1)
%     for d1=1:size(Peak,2)
%         if size(Peak{c,d1},2)>1
%         newpeak{c,d1}(find(isnanPeriod{c,d1}(round(Peak{c,d1}(:,3)*30))),:)=[];
%         end
%     end
% end
% end


end