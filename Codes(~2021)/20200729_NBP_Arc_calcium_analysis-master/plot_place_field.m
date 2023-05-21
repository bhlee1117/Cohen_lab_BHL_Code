function PV_corr=plot_place_field(dat,m,ref_day,sw,image_on)
if image_on
figure1 = figure('InvertHardcopy','off','PaperUnits','centimeters',...
    'Color',[1 1 1],...
    'Renderer','painters','position',[100 100 500 800]);
end
cond={[1 3],[2 4]};
if sw % TXN condition
    for g=1:3
        im{g}=[];  
    end
    for i=1:size(cond,2) % group
        for j=1:size(cond{i},2) % conditions for group
            Arc_class{i}(:,j)=dat{m}.Arc(:,ref_day+1)==cond{i}(1,j);
        end
        Arc_class{i}=sum(Arc_class{i},2)>0;
        list{i}=find(cell2mat(dat{m}.isPC(:,ref_day))==1 & ~isnan(cell2mat(dat{m}.isPC(:,ref_day+1))) ...
            & Arc_class{i});
%            list{i}=find(cell2mat(dat{m}.isPC(:,ref_day))==1 & ~isnan(cell2mat(dat{m}.isPC(:,ref_day+1))) & ~isnan(cell2mat(dat{m}.isPC(:,ref_day+2))) ...
%                       & Arc_class{i});
        if ~isempty(list{i})
        for c=1:size(list{i}) % cells in the group
            C=dat{m}.place_FM{list{i}(c),ref_day};
            %C(dat{m}.place_field{list{i}(c),ref_day}==0)=0;
            [max_fr max_arg{i}(c)]=max(C);
            %max_arg{i}(c) = sum(C'.*[1:size(C,1)])/sum(C);
        end
        [B order]=sort(max_arg{i},'ascend');
        s=setdiff([1 2],ref_day);
        g=1; % Day
        for k=[ref_day s]
            im{g}=[im{g}; cell2mat(dat{m}.place_FM(list{i}(order),k)')'];
            g=g+1;
        end
        end
    end  % group end
    g=1;
    for k=[ref_day s]
        
        %norm_im{g}=im{g};%./(max(im{g},[],2));%-min(im{g},[],2));
        for ccc=1:size(im{g},1)
            im{g}(ccc,:)=imgaussfilt(im{g}(ccc,:),2.7);
        end
        norm_im{g}=(im{g}-min(im{g},[],2))./(max(im{g},[],2)-min(im{g},[],2));
        if image_on
            subplot(6,1,2*g-1) 
        imagesc(norm_im{g}(1:length(list{1}),:),[0.2 1])
        colormap('jet')
        axis tight off equal
        subplot(6,1,2*g) 
        imagesc(norm_im{g}(length(list{1})+1:end,:),[0.2 1])
        colormap('jet')
        axis equal tight off
        g=g+1;
        end
    end
    
    for i=1:size(norm_im{1},1) % calculate PV correlation between Arc-on vs Arc-off
        minimum_bin=min([size(norm_im{1},2) size(norm_im{2},2) ]);
        %d1=imgaussfilt(im{1}(i,1:minimum_bin)',0.5); d2=imgaussfilt(im{2}(i,1:minimum_bin)',0.5);
        d1=(norm_im{1}(i,1:minimum_bin)'); d2=(norm_im{2}(i,1:minimum_bin)');
        PV(i,1)=corr(d1,d2,'Type','Pearson');
    end
    PV_corr{1}=PV(1:size(list{1},1),:); PV_corr{2}=PV(size(list{1},1)+1:end,:);
    
    
else  % Non TXN condition
    list=find(cell2mat(dat{m}.isPC(:,ref_day))==1 & sum(isnan(cell2mat(dat{m}.isPC)),2)==0);
    for i=1:size(list)
        C=dat{m}.place_FM{list(i),ref_day};
        %C(dat{m}.place_field{list(i),ref_day}==0)=0;
        [max_fr max_arg(i)]=max(C);
        %max_arg(i) = sum(C'.*[1:size(C,1)])/sum(C);
    end
    [B order]=sort(max_arg,'ascend');
    s=setdiff([1 2 3],ref_day);
    
    g=1;
    for i=[ref_day s]
        im{g}=cell2mat(dat{m}.place_FM(list(order),i)')';
        norm_im{g}=(im{g}-min(im{g},[],2))./(max(im{g},[],2)-min(im{g},[],2));
        
        if image_on
        subplot(3,1,g)
        imagesc(norm_im{g},[0.2 1])
        colormap('jet')
        axis tight off 
        end
        g=g+1;
    end
    for i=1:size(norm_im{1},1)
        minimum_bin=min([size(norm_im{1},2) size(norm_im{2},2) size(norm_im{3},2)]);
        PV_corr(i,1)=corr(norm_im{1}(i,1:minimum_bin)',norm_im{2}(i,1:minimum_bin)','Type','Pearson');
        PV_corr(i,2)=corr(norm_im{1}(i,1:minimum_bin)',norm_im{3}(i,1:minimum_bin)','Type','Pearson');
        PV_corr(i,3)=corr(norm_im{3}(i,1:minimum_bin)',norm_im{2}(i,1:minimum_bin)','Type','Pearson');
    end
    
end

end
