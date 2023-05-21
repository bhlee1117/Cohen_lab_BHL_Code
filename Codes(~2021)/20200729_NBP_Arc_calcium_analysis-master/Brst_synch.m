groups={[1 1 2 1],[1 2 2 1],[1 1 2 2],[1 2 2 2]};

for m=1:6
    for d1=1:3
        for c1=1:size(dat{m}.S,1)
            if size(dat{m}.S{c1,d1},2)>1
                moS=movsum(dat{m}.S{c1,3+d1},2);
                [val loc]=findpeaks(moS);
                S=zeros(1,size(dat{m}.S{c1,3+d1},2)); S(loc)=val;
                S_th{m}{c1,d1}=double(S>2.5); S_th{m}{c1,d1}(isnan(dat{m}.S{c1,d1}))=NaN;
                MovS{m}{c1,d1}=movsum(S_th{m}{c1,d1},3,'omitnan')>0;
            else
                S_th{m}{c1,d1}=[];
                MovS{m}{c1,d1}=[];
            end
        end
    end
    
    Arc_class{m}=dat{m}.Arc_post_ref(:,2:end);
    
    for g=1:size(groups,2)
        for cd=1:size(groups{g},2)/2
            list{m,g}(:,cd)=Arc_class{m}(:,groups{g}(1,2*cd-1))==groups{g}(1,2*cd);
        end
        list{m,g}=double(sum(list{m,g},2)==size(groups{g},2)/2); % cells that match the conditions
    end
    list{m,size(groups,2)+1}=zeros(size(dat{m}.S,1),1);
    list{m,size(groups,2)+1}(setdiff([1:size(dat{m}.Arc,1)],find(sum(cell2mat(list(m,:)),2)>0)),1)=1;
end

%%

for m=1:6
    for d=1:3
        [sz1 sz2]=cellfun(@size, MovS{m}(:,d));
        avail=find(sz2>0);
        %sync_N{m,d}=sum(cell2mat(MovS{m}(avail,d)),'omitnan');
        for g=1:size(list,2)
            l=(sz2>0) & list{m,g}; % calcium available & satisfy TXN conditions
            sync_N{m,d}(g,:)=sum(cell2mat(MovS{m}(l,d)),'omitnan');
        end
    end
end

%%
