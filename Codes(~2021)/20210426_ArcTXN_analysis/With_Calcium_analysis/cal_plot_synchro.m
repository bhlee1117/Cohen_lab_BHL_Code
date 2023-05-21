%function cal_plot_synchro(dat,groups,sprt)
bin=1;
clear list
for m=1:size(dat,2)
    Arc_class{m}=[dat{m}.Arc(:,2:end)];
    for g=1:size(groups,2)
        for d1=1:3
            list{m}{g,d1}=find(Arc_class{m}(:,d1)==[groups(1,g)]);
        end
    end
    for d1=1:3
        binned{m,d1}=[]; group_div{m,d1}=[]; cal_source{m,d1}=[];
        [r c]=cellfun(@size,dat{m}.Cal(:,d1));  max_bin=ceil(max(c)*1/30/bin); % Calculate the length of trace.
        Peaks=dat{m}.Peak(cell2mat(list{m}(:,d1))',d1); % Cells consist of Peak arrays
        Cal_trace=dat{m}.Cal(cell2mat(list{m}(:,d1))',d1);
        [gs cs]=cellfun(@size,list{m}(:,d1));
        gs=cumsum(gs);
        for p=1:size(Peaks,1) %Cell
            if sum(sum(isnan(Peaks{p,1})))==0
                bin_tmp=zeros(1,max_bin);
                if isempty(Peaks{p,1}) % no Peaks
                    binned{m,d1}=[binned{m,d1}; bin_tmp];
                else
                    bin_tmp(1,ceil(Peaks{p,1}(:,3)/bin))=1;
                    binned{m,d1}=[binned{m,d1}; bin_tmp];
                    cal_source{m,d1}=[cal_source{m,d1}; Cal_trace{p,1}];
                end
            end
            if ~isempty(find(gs==p)) % Where the each conditions end.
                group_div{m,d1}=[group_div{m,d1}; size(binned{m,d1},1)];
            end
        end
        %Calculate correlation matrix
        for i=1:size(binned{m,d1},1) %Cell 
            for j=1:size(binned{m,d1},1) %Cell
                if isempty(sprt)
                corr_mat{m,d1}(i,j)=corr(binned{m,d1}(i,:)',binned{m,d1}(j,:)','Type','Pearson');
                else
                    for t=1:ceil(max_bin/bin/sprt)
                        if t~=ceil(max_bin/bin/sprt)
                corr_mat{m,d1}(i,j,t)=corr(binned{m,d1}(i,sprt*(t-1)+1:sprt*t)',binned{m,d1}(j,sprt*(t-1)+1:sprt*t)','Type','Pearson');    
                        else
                corr_mat{m,d1}(i,j,t)=corr(binned{m,d1}(i,sprt*(t-1)+1:end)',binned{m,d1}(j,sprt*(t-1)+1:end)','Type','Pearson');                
                        end
                    end
                end
            end
        end
    end
end
%% Show example adjacent matrix
m=6; d=1;
for t=1:size(corr_mat{m,d1},3)
    ff=figure(1);
imagesc(corr_mat{m,d}(:,:,t),[0 1])
axis equal tight
colormap('jet')
colorbar
hold all
for i=1:size(group_div{2,1},1)-1
line([group_div{m,d}(i,1) group_div{m,d}(i,1)],[0 size(corr_mat{m,d},1)],'color','w','linewidth',1)
line([0 size(corr_mat{m,d},1)],[group_div{m,d}(i,1) group_div{m,d}(i,1)],'color','w','linewidth',1)
end
F=getframe(ff);
imwrite(F.cdata,'corr_mat.tif','Writemode','Append')
end

%%
clear corr_mat_pooled
for m=1:size(corr_mat,1)
    for d1=1:3
        z=zeros(size(corr_mat{m,d1},1),size(corr_mat{m,d1},2))+1;
        z=triu(z,1);
        cor_upmat=triu(corr_mat{m,d1},1);
        cor_upmat(z==0)=NaN;
        tmp=[1;group_div{m,d1}];
        ind=find((tmp(2:end)-tmp(1:end-1))~=0);
        divide=[[1;group_div{m,d1}(1:end-1,1)+1] [group_div{m,d1}]];
        for divx=ind'  %conditions
            for divy=ind'
        frag_mat=cor_upmat(divide(divx,1):divide(divx,2),divide(divy,1):divide(divy,2));
        corr_mat_pooled{divx,divy}{m,d1}=frag_mat(:);
            end
        end
    end
end

%%
cmap=distinguishable_colors(size(corr_mat_pooled,1)*size(corr_mat_pooled,2));
g=1;
for i=1:size(corr_mat_pooled,1)
    for j=1:size(corr_mat_pooled,2)
M(i,j)=mean(cell2mat(corr_mat_pooled{i,j}(:)),'omitnan');
X=randn(size(cell2mat(corr_mat_pooled{i,j}(:)),1),1)*0.05+i;
Y=randn(size(cell2mat(corr_mat_pooled{i,j}(:)),1),1)*0.05+j;
plot3(X,Y,cell2mat(corr_mat_pooled{i,j}(:)),'color',cmap(g,:),'marker','.','linestyle','none')
hold all
g=g+1;
    end
end
imagesc(imrotate(flipud(M),270),[0 0.04])
colormap('jet')
axis equal tight
grid on
set(gca,'XTick',[1 2 3 4],'XTickLabel',{'~B & ~A','~B & A','B & ~A','B & A'},...
        'YTick',[1 2 3 4],'YTickLabel',{'~B & ~A','~B & A','B & ~A','B & A'})


%end