function [centrality cc_dist estrada_index M]=cal_parms_graphs(corr_mat,group_div,corr_th,p_val)

cmap_group=[0 0 0;1 0 0]; s_index=[1 -1];
for m=1:size(corr_mat,1)
    for d1=1:size(corr_mat,2)
        A=round(corr_mat{m,d1},3);
        A(find(diag(zeros(size(A,1),1)+1)==1))=0;
        A=double(A>corr_th & p_val{m,d1}<0.05);
        
        
        divide=[[1;group_div{m,d1}(1:end-1,1)+1] [group_div{m,d1}]];
        dc_norm=sum(A,2)/sum(A(:))*2;
        s=[];
        for g=1:size(group_div{m,d1},1)
            centrality.dcdist{g}{m,d1}=sum(A(divide(g,1):divide(g,2),:),2);
            %centrality.dcnormdist{g}{m,d1}=sum(A(divide(g,1):divide(g,2),:),2)/sum(A(:))*2;
            centrality.dcnormdist{g}{m,d1}=sum(A(divide(g,1):divide(g,2),:),2)/(size(A,1)-1);
            cc=clusteringcoef(A,'old');
            cc_dist{g}{m,d1}=cc(divide(g,1):divide(g,2))';
            nsc=diag(expm(A))/trace(expm(A)); %Normalized subgraph centrality
            estrada_index{g}{m,d1}=nsc(divide(g,1):divide(g,2));
            if g<3
            s(divide(g,1):divide(g,2),1)=s_index(1,g);
            end
        end
        G=graph(double(A(1:size(s,1),1:size(s,1))));
         [M.modular(m,d1) M.modular_norm(m,d1)]=modularity(G, s');
%         for ran_it=1:50
%             [NM rand_order]=sort(randn(size(s,1),1));
%             [ran_mod(m,d1,ran_it) ran_mod_norm(m,d1,ran_it)]=modularity(G, s(rand_order)');
%         end
%         modular{2}(m,d1)=mean(ran_mod(m,d1,:)); modular_norm{2}(m,d1)=mean(ran_mod_norm(m,d1,:));
    end
end
end