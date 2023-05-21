function export_to_Gephi(corr_mat,p_val,m,d1,corr_th,group_div,pc,s_index,fnm)

        A=round(corr_mat{m,d1},3);
        A(find(diag(zeros(size(A,1),1)+1)==1))=0;
        A=double(A>corr_th & p_val{m,d1}<0.05);
        A=triu(A,1);
        divide=[[1;group_div{m,d1}(1:end-1,1)+1] [group_div{m,d1}]];
        if size(divide,1)~= size(s_index,2)
            error('Check the size of index vector')
        end
        for i=1:size(divide,1)
         s(divide(i,1):divide(i,2),1)=s_index(i);
        end
        pc=cell2mat(pc{m,d1});
        s=s+pc*max(s_index);
EdgeL=adj2gephilab(fnm,A,s);
end