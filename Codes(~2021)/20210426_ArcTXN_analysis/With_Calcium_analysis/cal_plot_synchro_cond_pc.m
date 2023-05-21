function [corr_mat corr_mat_pooled group_div p_val]=cal_plot_synchro_cond_pc(dat,groups,sprt)
bin=1;
clear list corr_mat
for m=1:size(dat,2)
    PC_list{m}=cell2mat(dat{m}.isPC);
    for g=1:size(groups,2)
        for cd=1:size(groups{g},2)/2
            list{m}{g,1}(:,cd)=PC_list{m}(:,groups{g}(1,2*cd-1))==groups{g}(1,2*cd);
        end
        list{m}{g,1}=find(sum(list{m}{g,1},2)==size(groups{g},2)/2); % cells that match the conditions
    end
    list{m}{g+1,1}=setdiff([1:size(PC_list{m},1)],cell2mat(list{m}))';
    for d1=1:3
        binned{m,d1}=[]; group_div{m,d1}=[]; cal_source{m,d1}=[]; pc{m,d1}=[];
        [r c]=cellfun(@size,dat{m}.Cal(:,3+d1));  max_bin=ceil(max(c)*1/30/bin); % Calculate the length of trace.
        %Peaks=dat{m}.Peak(cell2mat(list{m})',d1); % Cells consist of Peak arrays
        Peaks=cellfun(@spike2binnedarray, dat{m}.S(cell2mat(list{m})',3+d1),'UniformOutput',false); % Cells consist of Peak arrays
        Cal_trace=dat{m}.Cal(cell2mat(list{m})',3+d1); 
        [gs cs]=cellfun(@size,list{m});
        gs=cumsum(gs);
        for p=1:size(Peaks,1) %Cell
            if size(Peaks{p,1},2)>1%sum(sum(isnan(Peaks{p,1})))==0
               % bin_tmp=zeros(1,max_bin);
                if isempty(Peaks{p,1}) % no Peaks
                   %nan_time=find(isnan(Cal_trace{p,1}))/30;
                    %bin_tmp(1,ceil(nan_time/bin))=NaN; %out-of focus & grooming convert to NaN
                    
                    binned{m,d1}=[binned{m,d1}; Peaks{p,1}];
                    
                else
%                     bin_tmp(1,ceil(Peaks{p,1}(:,2)/bin))=1;
                 % nan_time=find(isnan(Cal_trace{p,1}))/30;
                 % bin_tmp(1,ceil(nan_time/bin))=NaN;
                  binned{m,d1}=[binned{m,d1}; Peaks{p,1}];
                    
                    %                    cal_source{m,d1}=[cal_source{m,d1}; Cal_trace{p,1}];
                end
            end
            if ~isempty(find(gs==p)) % Where the each conditions end.
                group_div{m,d1}=[group_div{m,d1}; size(binned{m,d1},1)];
            end
        end
        for i=1:size(binned{m,d1},1) %Cell
            for j=1:size(binned{m,d1},1) %Cell
                if isempty(sprt)
              [corr_mat{m,d1}(i,j) p_val{m,d1}(i,j)]=corr(binned{m,d1}(i,:)',binned{m,d1}(j,:)'...
                                                   ,'Type','Pearson','Rows','pairwise');  % not include NaN value
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

clear corr_mat_pooled cor_upmat
for m=1:size(corr_mat,1)
    for d1=1:3
        clear cor_upmat
        z=zeros(size(corr_mat{m,d1},1),size(corr_mat{m,d1},2))+1;
        z=triu(z,1);
        for t=1:size(corr_mat{m,d1},3)
            tp=triu(corr_mat{m,d1}(:,:,t),1);
            tp(z==0)=NaN;
            cor_upmat(:,:,t)=tp;
        end
        tmp=[1;group_div{m,d1}];
        ind=find((tmp(2:end)-tmp(1:end-1))~=0);
        divide=[[1;group_div{m,d1}(1:end-1,1)+1] [group_div{m,d1}]];
        for divx=ind'  %conditions
            for divy=ind'
                if isempty(sprt)
                    frag_mat=cor_upmat(divide(divx,1):divide(divx,2),divide(divy,1):divide(divy,2));
                    corr_mat_pooled{divx,divy}{m,d1}=frag_mat(:);
                else
                    frag_mat=cor_upmat(divide(divx,1):divide(divx,2),divide(divy,1):divide(divy,2),:);
                    corr_mat_pooled{divx,divy}{m,d1}=frag_mat;
                end
            end
        end
    end
end
end
function P=spike2array(S)
if size(S,2)>1
sl=find(S>0);
P=[S(sl)' sl'/30];
else
    P=NaN;
end
end