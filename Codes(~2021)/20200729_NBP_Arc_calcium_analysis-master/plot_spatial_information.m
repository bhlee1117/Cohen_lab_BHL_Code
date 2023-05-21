%% This code is written for calculating spatial information, according to 'An Information-Theoretic Approach to Deciphering the Hippocampal Code'

% MODIFICATION HISTORY :
%           Written by Byung Hun Lee, Deptartment of Physics and Astronomy,
%           Seoul National University, 2020/11/30
function [FR lambda_p]=plot_spatial_information(dat,bin,Full_result,groups,post_ref,PC_day,t_mouse,sw,y_lim,xtick)
SI_pool=[];
for m=1:size(dat,2)
    min_track=-2000;
    [r frm_max]=cellfun(@size,dat{m}.Cal(:,1:3));
    frm_max=max(frm_max);
    P=dat{m}.S;
    clear P_tr
    for d1=1:3
        clear p_bin 
        place_cdf=[Full_result{m}.VR{d1}(1:end,2)'];
        place_cdf(1,1)=min_track;
        
        for c=1:size(P,1)
            if size(P{c,d1},2)>1
                P_tr{c,d1}(:,1)=place_cdf;
                P_tr{c,d1}(:,2)=P{c,3+d1};
                %P_tr{c,d1}(round(P{c,d1}(:,3)*30),2)=1;%P{c,d1}(:,2);
                %cal=dat{m}.Cal{c,3+d1};
                %P_tr{c,d1}(1:size(dat{m}.run{d1},1),2)=cal(1:size(dat{m}.run{d1},1));
                P_tr{c,d1}(find(dat{m}.run{d1}==0),:)=NaN;
                bin_place=round((P_tr{c,d1}(:,1)-min_track)/bin);
                lambda_all=sum(P_tr{c,d1}(:,2),'omitnan')/(sum(~isnan(P_tr{c,d1}(:,1)))/30);
                for p=1:150%max(bin_place)
                    lambda_p{m}{c,d1}(p,1)=sum(P_tr{c,d1}(find(bin_place==p),2),'omitnan')/(sum(bin_place==p,'omitnan')/30);
                    p_bin(p,1)=sum(bin_place==p,'omitnan')/sum(~isnan(P_tr{c,d1}(:,1)));
                end
                lambda_p{m}{c,d1}(lambda_p{m}{c,d1}<0)=NaN;
                SI{m}(c,d1)=sum(lambda_p{m}{c,d1}./abs(lambda_all).*log(lambda_p{m}{c,d1}./abs(lambda_all))/log(2).*p_bin,'omitnan');
            else
                SI{m}(c,d1)=NaN;
            end
        end
    end
    
    if post_ref
        Arc_class{m}=[dat{m}.Arc_post_ref(:,2:end)];
    else
        Arc_class{m}=[dat{m}.Arc(:,2:end)];
    end
    for g=1:size(groups,2)
        for cd=1:size(groups{g},2)/2
            list{m,g}(:,cd)=double(Arc_class{m}(:,groups{g}(1,2*cd-1))==groups{g}(1,2*cd));
        end
        if isempty(PC_day)
            %list{m,g}(:,end+1)=cell2mat(dat{m}.isPC(:,PC_day))==1;
            list{m,g}=sum(list{m,g},2)==size(groups{g},2)/2;
        else
            list{m,g}(:,end+1)=cell2mat(dat{m}.isPC(:,PC_day))==1;
            list{m,g}=sum(list{m,g},2)==size(groups{g},2)/2+1; % cells that match the conditions
        end
        FR{m,g}=[SI{m}(find(list{m,g}),:)];
        SI_pool=[SI_pool;SI{m}];
        
    end
end

%%
p_list=[1 2;1 3;2 3;1 4;2 4;3 4];
if sw  % arrange  day1, day2, day3
    p_list=p_list(1:3,:);
    xtick={'Day 1','Day 2','Day 3'};
    for g=1:size(groups,2)
        figure1 = figure('InvertHardcopy','off','PaperUnits','centimeters',...
            'Color',[1 1 1],...
            'Renderer','painters','position',[100 100 300 300]);
        if t_mouse % mouse by mouse
            for mm=1:size(FR,1)
                d(mm,:)=mean(FR{mm,g},'omitnan');
            end
        else % pooled
            d=cell2mat(FR(:,g));
        end
        for pp=1:size(p_list,1)
            [a p_value(pp)]=ttest(d(:,p_list(pp,1)),d(:,p_list(pp,2))); end
        
        M=mean(d,1,'omitnan');   n=sum(~isnan(d),1);     S=std(d,0,1,'omitnan')./sqrt(n);      G=[1 2 3];
        for i=1:size(d,1)
            plot(G,d(i,:),'color',[0.7 0.7 0.7],'marker','.')
            hold all
        end
        errorbar(G,M,S,'LineWidth',2,'linestyle','-','color','k','Capsize',10)
        hold all
        set(gca,'FontSize',8,'LineWidth',1,'XTick',G,'XTickLabel',xtick,'FontName','arial rounded mt bold','FontSize',13,'LineWidth',2);
        ylabel('Spatial information','LineWidth',2,'FontSize',13,'FontName','arial rounded mt bold')
        xlim([0.3 3.7])
        ylim([0 y_lim])
        
        for i=1:size(p_value,2)
            if p_value(i)<0.05
                star='*';
                if p_value(i)<0.01
                    star='**';
                    if p_value(i)<0.001
                        star='***';
                    end
                end
                
                line(p_list(i,:),[y_lim*0.8+y_lim*0.05*i y_lim*0.8+y_lim*0.05*i],'color','k','linewidth',2);
                line([p_list(i,1) p_list(i,1)],[y_lim*0.7+y_lim*0.05*i y_lim*0.8+y_lim*0.05*i],'color','k','linewidth',2);
                line([p_list(i,2) p_list(i,2)],[y_lim*0.7+y_lim*0.05*i y_lim*0.8+y_lim*0.05*i],'color','k','linewidth',2);
                text(mean(p_list(i,:)),y_lim*0.82+y_lim*0.05*i,star,'FontSize',13,'FontName','arial rounded mt bold',...
                    'HorizontalAlignment', 'center')
            end
        end  % Calculate p-value
    end
    
else  % Day -> groups
    
    if size(groups,2)==2
        p_list=p_list(1,:); end
    if size(groups,2)==3
        p_list=p_list(1:3,:); end
    if size(groups,2)==4
        p_list=p_list(1:6,:); end
    
    figure1 = figure('InvertHardcopy','off','PaperUnits','centimeters',...
        'Color',[1 1 1],...
        'Renderer','painters','position',[100 100 350 250]);
    G=[0.8 1.2 1.8 2.2 2.8 3.2];
    for day=1:3
        if t_mouse
            for g=1:size(groups,2)
                for kk=1:size(FR,1)
                    tmp(kk,:)=mean(FR{kk,g},1,'omitnan'); end
                d{day,g}=tmp(:,day);
            end
        else
            for g=1:size(groups,2)
                tmp=cell2mat(FR(:,g));
                d{day,g}=tmp(:,day);
            end
        end
        M=[]; S=[];
        cmap=[0.1 0.1 0.1;1 0 0];
%         for g=1:size(groups,2)
%             M=[M mean(d{day,g},1,'omitnan')];
%             S=[S  std(d{day,g},0,1,'omitnan')./sqrt(sum(~isnan(d{day,g})))]; end
%         
%         b=bar(G(2*day-1:2*day),M,'Barwidth',0.7,'LineWidth',2);
%         hold all
%         b.FaceColor='flat';
%         for j=1:2
%             b.CData(j,:)= cmap(j,:);
%         end
plot_errorbar3(G(2*day-1:2*day),d(day,:),'ranksum',5,'Spatial information',{'A1-Arc^-','A1-Arc^+'},[0.1 0.1 0.1; 1 0 0]);
        
        if t_mouse
            for pp=1:size(p_list,1)
                [a p_value(pp,1)]=ttest(d{day,p_list(pp,1)},d{day,p_list(pp,2)});
            end
            for mice=1:size(d{day,g},1)
                tmp2=cell2mat(d(day,:));
                %plot([1:size(groups,2)]+randn(1,1)*0,tmp2(mice,:),'marker','o','color',[0.4 0.4 0.4],'linewidth',1)
                plot(G([2*day-1:2*day]),tmp2(mice,:),'marker','o','color',[0.7 0.7 0.7],'linewidth',1)
                hold all
            end
        else
            
            for pp=1:size(p_list,1)
                [a p_value(pp,1)]=kstest2(d{day,p_list(pp,1)},d{day,p_list(pp,2)});
                [a p_value(pp,2)]=ttest2(d{day,p_list(pp,1)},d{day,p_list(pp,2)});
                [p_value(pp,3) a]=ranksum(d{day,p_list(pp,1)},d{day,p_list(pp,2)});
            end
            for g=1:size(groups,2)
                plot(G(2*day-2+g)+randn(size(d{day,g},1),1)*0.02,d{day,g},'marker','.','color',[0.7 0.7 0.7],'linestyle','none')
                hold all
            end
        end
        p_value
%         errorbar(G(2*day-1:2*day),M,S,...
%             'LineWidth',2,'linestyle','-','color','k','Capsize',10)
        
        set(gca,'FontSize',8,'LineWidth',1,'XTick',G,'XTickLabel',...
            xtick,'FontName','arial rounded mt bold','FontSize',13,'LineWidth',2);
        ylabel('Spatial information','LineWidth',2,'FontSize',13,...
            'FontName','arial rounded mt bold')
        %xlim([0.3 size(groups,2)+0.7])
        xlim([0.5 G(end)+0.3])
        ylim([0 y_lim])
        if ~t_mouse
            ref=3;
        else
            ref=1; end
        for i=1:size(p_list,1)
            if p_value(i,ref)<0.05
                star='*';
                if p_value(i,ref)<0.01
                    star='**';
                    if p_value(i,ref)<0.001
                        star='***';
                    end
                end
                line(day-1+G(p_list(i,:)),[y_lim*0.77+y_lim*0.04*i y_lim*0.77+y_lim*0.04*i],'color','k','linewidth',2);
                line(day-1+G([p_list(i,1) p_list(i,1)]),[y_lim*0.70+y_lim*0.04*i y_lim*0.77+y_lim*0.04*i],'color','k','linewidth',2);
                line(day-1+G([p_list(i,2) p_list(i,2)]),[y_lim*0.70+y_lim*0.04*i y_lim*0.77+y_lim*0.04*i],'color','k','linewidth',2);
                text(day-1+mean(G(p_list(i,:))),y_lim*0.79+y_lim*0.04*i,star,'FontSize',13,'FontName','arial rounded mt bold',...
                    'HorizontalAlignment', 'center')
            end
        end
        xtickangle(45)
    end
end

end