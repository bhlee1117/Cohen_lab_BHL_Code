%function plot_population_vector(dat,groups,xtick,y_lim,sw,cat)
cond=[1 2;1 3;2 3];
clear list 
ylab='Mean firing rate correlation';
for g=1:size(groups,2)
    for cd=1:size(cond,1)
FR_C{g,cd}=[];end 
end

for m=1:size(dat,2) %mouse
    for c1=1:size(dat{m}.Cal,1)
        for d1=1:3
            if size(dat{m}.Cal{c1,d1},2)<2
            firing_rate{m}(c1,d1)=NaN;
            else
            switch cat
                case 1 % Mean dFF
            firing_rate{m}(c1,d1)=mean(dat{m}.Cal{c1,3+d1},'omitnan');
                case 2 % Sigma Peak
                    if isempty(dat{m}.Peak{c1,d1})
            firing_rate{m}(c1,d1)=0;
                    else
            firing_rate{m}(c1,d1)=sum(dat{m}.Peak{c1,d1}(:,2),'omitnan');
                    end
                case 3 % Sigma Tr
                    tr_area=0;
                    for tr=1:size(dat{m}.Cal{c1,6+d1},1)
                        if sum(dat{m}.Cal{c1,6+d1}(tr,:)==0)==0
                            if dat{m}.Cal{c1,6+d1}(tr,1)>size(dat{m}.Cal{c1,3+d1},2)
                            else
                                if dat{m}.Cal{c1,6+d1}(tr,2)>size(dat{m}.Cal{c1,3+d1},2)
                                   dat{m}.Cal{c1,6+d1}(tr,2)=size(dat{m}.Cal{c1,3+d1},2);
                                end
                         tr_area=tr_area+sum(dat{m}.Cal{c1,3+d1}(dat{m}.Cal{c1,6+d1}(tr,1):dat{m}.Cal{c1,6+d1}(tr,2)));
                            end
                        end
                    end
                    firing_rate{m}(c1,d1)=tr_area;
                case 4 % Number of Peak
            firing_rate{m}(c1,d1)=size(dat{m}.Peak{c1,d1},1);
            end
            end
        end
    end
    Arc_class{m}=[dat{m}.Arc(:,2:end)];
    for g=1:size(groups,2)
        for cd=1:size(groups{g},2)/2
            list{m,g}(:,cd)=Arc_class{m}(:,groups{g}(1,2*cd-1))==groups{g}(1,2*cd);
        end
        list{m,g}=find(sum(list{m,g},2)==size(groups{g},2)/2); % cells that match the conditions
        for cd=1:size(cond,1)
      tmp_fr=firing_rate{m}(list{m,g},cond(cd,:));      
      x=find(sum(isnan(tmp_fr),2)~=0);
      tmp_fr(x,:)=[];
      pop_cor{g,1}(m,cd)=corr(tmp_fr(:,1),tmp_fr(:,2),'Type','Pearson');
      %pop_cor{g,1}(m,cd)=abs(tmp_fr(:,1)-tmp_fr(:,2))./abs(tmp_fr(:,1)+tmp_fr(:,2));
      FR_list{m}{g,cd}=tmp_fr;
      FR_C{g,cd}=[FR_C{g,cd}; tmp_fr];
        end
    end
end

%%%%% Data arrange done %%%%%
%%%%% Significance check & plotting start %%%%%
cmap=distinguishable_colors(6);
for g=1:size(groups,2)
 figure1 = figure('InvertHardcopy','off','PaperUnits','centimeters',...
            'Color',[1 1 1],...
            'Renderer','painters','position',[100 100 600 300]);
        for  cd=1:size(cond,1)
            for m=1:6
        subplot(1,size(cond,1),cd)
        plot(FR_list{m}{g,cd}(:,1),FR_list{m}{g,cd}(:,2),'.','color',cmap(m,:))
        hold all
        axis equal tight
        xlim([0 3000])
        ylim([0 3000])
            end
        end
end

figure1 = figure('InvertHardcopy','off','PaperUnits','centimeters',...
            'Color',[1 1 1],...
            'Renderer','painters','position',[100 100 600 300]);
for m=1:size(FR_list,2)
    plot([1 2],pop_cor{1}(m,[1 2]),'o','color',cmap(m,:),'linestyle','-')
    hold all
    plot([3 4],pop_cor{2}(m,[1 2]),'o','color',cmap(m,:),'linestyle','-')
end

%%
subplot(1,2,1)
plot(FR_C{1,1}(:,1),FR_C{1,1}(:,2),'.')
hold all
plot(FR_C{2,1}(:,1),FR_C{2,1}(:,2),'.')
axis equal tight
subplot(1,2,2)
plot(FR_C{1,3}(:,1),FR_C{1,3}(:,2),'.')
hold all
plot(FR_C{2,3}(:,1),FR_C{2,3}(:,2),'.')
axis equal tight