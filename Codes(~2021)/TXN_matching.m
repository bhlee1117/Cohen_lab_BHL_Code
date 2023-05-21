%%
clear
[fnm pth]=uigetfile('*.mat','Multiselect','on');
%%
clear TXN_all
for i=1%:length(fnm)/2
    for j=1:2
    load([pth fnm{1,2*(i-1)+j}])
    TXN_all{i,j}=[datum(TXN,[1:2 4]) intensity(TXN,1)];
    Single_all{i,j}=[datum(single,[1:2 4]) intensity(single,1)];
    total_all{i,j}=[datum(:,[1:2 4]) intensity(:,1)];
    end
end

%%
clear list matlist CDS_PP_TXN CDS_PP_Sig Dst Dst_sig
for i=1:length(fnm)/2
    for j=1:size(TXN_all{i,1},1)
        for k=1:size(total_all{i,2},1)
    Dst{i,1}(j,k)=distance_BH(TXN_all{i,1}(j,1:2),total_all{i,2}(k,1:2));
        end
    end
%       for j=1:size(Single_all{i,1},1)
%         for k=1:size(Single_all{i,2},1)
%     Dst_sig{i,1}(j,k)=distance_BH(Single_all{i,1}(j,1:2),Single_all{i,2}(k,1:2));
%         end
%     end
[minimum arg]=min(Dst{i,1},[],2);
list=find(minimum<240); % 1 um
matlist{i,1}=[list' arg(list,1)'];
% 
% [minimum_sig arg_sig]=min(Dst_sig{i,1});
% list_sig=find(minimum_sig<300); % 1 um
% matlist_sig{i,1}=[list_sig' arg_sig(1,list_sig)'];

CDS_PP_TXN{i,1}=[TXN_all{i,1}(matlist{i,1}(:,2),:) total_all{i,2}(matlist{i,1}(:,1),:)];
%CDS_PP_Sig{i,1}=[Single_all{i,1}(matlist_sig{i,1}(:,2),:) Single_all{i,2}(matlist_sig{i,1}(:,1),:)];
mat_prob(i,1:2)=[size(list,1)/size(TXN_all{i,1},1) size(list,1)/size(total_all{i,2},1)];
end
%%
for i=1:length(fnm)/2
    single_int(i,1:2)=median(CDS_PP_Sig{i,1}(:,[4 8]));
    mRNA_per_TXN{i,1}=[CDS_PP_TXN{i,1}(:,4)/single_int(i,1) CDS_PP_TXN{i,1}(:,8)/single_int(i,2)]; 
end
%% plot
load(['H:\Image_data\KIST_RNAscope\40x\20190209_SNU_RNAscope\RNA_scope_confocal\Result\Coord_included_20200129\Total.mat'])
mat_CDS_PP_TXN=cell2mat(CDS_PP_TXN);
mat_CDS_PP_Sig=cell2mat(CDS_PP_Sig);
mat_mRNA_per_TXN=cell2mat(mRNA_per_TXN);
figure(1)
plot(mat_mRNA_per_TXN(:,1),mat_mRNA_per_TXN(:,2),'.')
 xlim([0 250])
 ylim([0 250])
 ylabel('number of mRNA per TXN (PP7)')
 xlabel('number of mRNA per TXN (CDS)')
set(gca,'XScale','log','YScale','log')
% figure(2)
grid{1}=randn(size(mat_mRNA_per_TXN,1),1)*0.1+1;

grid{2}=randn(size(mat_mRNA_per_TXN,1),1)*0.1+2;
plot([grid{1} grid{2}],[mat_mRNA_per_TXN(:,1) mat_mRNA_per_TXN(:,2)],'k.')
hold all
errorbar([1 2],mean(mat_mRNA_per_TXN,1),std(mat_mRNA_per_TXN,0,1),'color','r','linestyle','none')
set(gca,'XTick',[1 2],'XTicklabel',{'CDS','PP7'},'YScale','linear')
ylabel('number of mRNA per TXN')
ylim([0 250])