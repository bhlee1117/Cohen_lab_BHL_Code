function [p_value datum L]=plot_inform_theta(dat,IBI,post_ref,y_lim,day_sw)
conds=[1 2 2;2 2 2];
%conds=[2 1;2 2];
xtick={'A1 Arc^-','A1 Arc^+','-+-','--+','++-','+-+','-++','+++'};
%xtick={'--','+-','++','-+','++-','+-+','-++','+++'};
%xtick={'--','+-','++','-+','++-','+-+','-++','+++'};
%doi=[1:size(conds,2)];
doi=[1 2 3];
for m=1:size(dat,2)
    [r cc]=cellfun(@size,dat{m}.Cal(:,1:3));
    avail_list=find(sum(cc==1,2)==0);
    P=dat{m}.Peak(avail_list,:);
    S=dat{m}.S(avail_list,4:6);
    tmp_theta=cell2mat(IBI.theta_Brst(m,:));
    tmp_theta=cell2mat(cellfun(@brst_sum,IBI.IBI(m,:),'UniformOutput',false));
    [deconSum]=tmp_theta(avail_list,:)+1;
    nanS=cellfun(@isnan,S,'UniformOutput',false);
    [a sizeS]=cellfun(@size,nanS); nanSumS=cellfun(@sum,nanS); t_list=(sizeS-nanSumS)/30;
    %deconSum=deconMean.*t_list;
    if day_sw
        lambda{m}=[sum(deconSum(:,1:2),2)./sum(t_list(:,1:2),2) deconSum(:,3)./t_list(:,3) sum(deconSum,2)./sum(t_list,2)]; % lambda_A lambda_B lambda_all
        %         M_inform{m,1}=lambda{m}(:,1).*log(lambda{m}(:,1)./lambda{m}(:,3))/log(2).*sum(t_list(:,1:2),2)./sum(t_list,2)+...
        %                                     lambda{m}(:,2).*log(lambda{m}(:,2)./lambda{m}(:,3))/log(2).*sum(t_list(1,3),2)./sum(t_list,2);
        %Corrected
             M_inform{m,1}=lambda{m}(:,1)./lambda{m}(:,3).*log(lambda{m}(:,1)./lambda{m}(:,3))/log(2).*sum(t_list(:,1:2),2)./sum(t_list,2)+...
                                        lambda{m}(:,2)./lambda{m}(:,3).*log(lambda{m}(:,2)./lambda{m}(:,3))/log(2).*sum(t_list(:,3),2)./sum(t_list,2);
% Raw
%         M_inform{m,1}=lambda{m}(:,1).*log(lambda{m}(:,1)./lambda{m}(:,3))/log(2).*sum(t_list(:,1:2),2)./sum(t_list,2)+...
%             lambda{m}(:,2).*log(lambda{m}(:,2)./lambda{m}(:,3))/log(2).*sum(t_list(:,3),2)./sum(t_list,2);
    else
        %lambda{m}=[sum(deconMean(:,2),2)/sum(t_list(1,2)) deconMean(:,3)/t_list(1,3) sum(deconMean(:,2:3),2)/sum(t_list(1,2:3))]; % lambda_A lambda_B lambda_all
        lambda{m}=[deconSum(:,2)./t_list(:,2) deconSum(:,3)./t_list(:,3) sum(deconSum(:,2:3),2)./sum(t_list(:,2:3),2)]; % lambda_A lambda_B lambda_all
        %     M_inform{m,1}=lambda{m}(:,1).*log(lambda{m}(:,1)./lambda{m}(:,3))/log(2)*sum(t_list(1,2))/sum(t_list(2:3))+...
        %        lambda{m}(:,2).*log(lambda{m}(:,2)./lambda{m}(:,3))/log(2)*sum(t_list(1,3))/sum(t_list(2:3));
        %         M_inform{m,1}=lambda{m}(:,1)./lambda{m}(:,3).*log(lambda{m}(:,1)./lambda{m}(:,3))/log(2).*t_list(:,2)./sum(t_list,2)+...
        %             lambda{m}(:,2)./lambda{m}(:,3).*log(lambda{m}(:,2)./lambda{m}(:,3))/log(2).*t_list(:,3)./sum(t_list,2);
        % Corrected
%                 M_inform{m,1}=lambda{m}(:,1)./lambda{m}(:,3).*log(lambda{m}(:,1)./lambda{m}(:,3))/log(2).*t_list(:,2)./sum(t_list(2:3),2)+...
%                     lambda{m}(:,2)./lambda{m}(:,3).*log(lambda{m}(:,2)./lambda{m}(:,3))/log(2).*t_list(:,3)./sum(t_list(2:3),2);
        % Raw
        M_inform{m,1}=lambda{m}(:,1).*log(lambda{m}(:,1)./lambda{m}(:,3))/log(2).*t_list(:,2)./sum(t_list(2:3),2)+...
            lambda{m}(:,2).*log(lambda{m}(:,2)./lambda{m}(:,3))/log(2).*t_list(:,3)./sum(t_list(2:3),2);
        
    end
    % M_inform{m,1}=lambda{m}(:,1)./lambda{m}(:,3).*log(lambda{m}(:,1)./lambda{m}(:,3))/log(2)*sum(t_list(1,1:2))/sum(t_list)+...
    %         lambda{m}(:,2)./lambda{m}(:,3).*log(lambda{m}(:,2)./lambda{m}(:,3))/log(2)*sum(t_list(1,3))/sum(t_list);
    
    if post_ref
        Arc_class=dat{m}.Arc_post_ref(avail_list,2:end);
        Arc_list{m,1}=zeros(size(avail_list,1),1);
    else
        Arc_class=dat{m}.Arc(avail_list,2:end);
        Arc_list{m,1}=zeros(size(avail_list,1),1);
    end
    for cd=1:size(conds,1)
        list=[];
        for d1=1:size(doi,2)
            list(:,d1)=Arc_class(:,doi(d1))==conds(cd,d1);
        end
        Arc_list{m,1}(find(sum(list,2)==size(conds,2)),1)=cd;
    end
end
Arc_list_pool=cell2mat(Arc_list);  M_inform_pool=cell2mat(M_inform); lambda_pool=cell2mat(lambda');

%plot
%%
figure1 = figure('InvertHardcopy','off','PaperUnits','centimeters',...
    'Color',[1 1 1],'Renderer','painters','position',[100 100 300 300]);

for cd=1:size(conds,1)
    datum{cd}=M_inform_pool(find(Arc_list_pool==cd),1);
    L{cd}=lambda_pool(find(Arc_list_pool==cd),:);
    %datum{cd}(isnan(datum{cd}))=0;
    plot(cd/2+randn(size(datum{cd},1),1)*0.03,datum{cd},'.','color',[0.6 0.6 0.6],'markersize',8)
    hold all
end
for i=1:size(conds,1)
    for j=1:size(conds,1)
        %[p_value(i,j) a ]=ranksum(datum{i},datum{j});
        %[a p_value(i,j)]=ttest2(datum{i},datum{j});
        %[a p_value(i,j)]=kstest2(datum{i},datum{j});
        [p_value(i,j) a]=ranksum(datum{i},datum{j});
    end
end
p_value=triu(p_value,1); p_value(p_value==0)=NaN;
[N N2]=cellfun(@size,datum);
M=cellfun(@mean_omitnan,datum); S=cellfun(@std_omitnan,datum)./sqrt(N);
errorbar([1:size(conds,1)]/2,M,S,'linestyle','none','LineWidth',2,'color','r','Capsize',10,'marker','+','markersize',8)
set(gca,'FontSize',8,'LineWidth',1,'XTick',[1:size(conds,1)]/2,'XTickLabel',...
    xtick,'FontName','arial rounded mt bold','FontSize',13,'LineWidth',2);
ylabel('Mutual information (bit/spike)','LineWidth',2,'FontSize',13,...
    'FontName','arial rounded mt bold')
xlim([0.3 size(conds,1)/2+0.3])
ylim([0 y_lim])
g=1;
p_value
for i=find(p_value<0.05)'
    if p_value(i)<0.05
        star='*';
        if p_value(i)<0.01
            star='**';
            if p_value(i)<0.001
                star='***';
            end
        end
        [r1 c1]=ind2sub(size(conds,1),i);
        line([r1 c1]/2,[y_lim*0.77+y_lim*0.04*g y_lim*0.77+y_lim*0.04*g],'color','k','linewidth',2);
        line([r1 r1]/2,[y_lim*0.70+y_lim*0.04*g y_lim*0.77+y_lim*0.04*g],'color','k','linewidth',2);
        line([c1 c1]/2,[y_lim*0.70+y_lim*0.04*g y_lim*0.77+y_lim*0.04*g],'color','k','linewidth',2);
        text(mean([r1 c1]/2),y_lim*0.79+y_lim*0.04*g,star,'FontSize',13,'FontName','arial rounded mt bold',...
            'HorizontalAlignment', 'center')
    end
    g=g+1;
end

end
function B=brst_sum(IBI)
    B=sum(IBI(:,3:100),2);
end
function M=mean_omitnan(D)
M=mean(D,'omitnan');
end
function M=std_omitnan(D)
M=std(D,0,1,'omitnan');
end
