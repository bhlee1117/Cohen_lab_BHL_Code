function [P_real P_ind patt]=cal_patt_prob(sp)
%clear P_real patt Prob_each P_ind

%sp=movsum(Result{1}.spike,2,2)>0;
%sp=sp(:,[1:30:end]);

bin_fact_sp=sp.*2.^[1:size(sp,1)]';
bin_fact=sum(bin_fact_sp,1);
[bin_unique, bin_loc]=unique(bin_fact);
FR=sum(sp,2,'omitnan')./sum(~isnan(sp),2);
nbytes = fprintf('processing pattern %0 of %d', length(bin_unique));
for i=1:length(bin_unique)
fprintf(repmat('\b',1,nbytes))
nbytes = fprintf('processing pattern %d of %d', i, length(bin_unique));
P_real(i)=length(find(bin_fact==bin_unique(i)))/mean(sum(~isnan(sp),2));
patt(:,i)=sp(:,bin_loc(i));
Prob_each(:,i)=(2*FR-1).*sp(:,bin_loc(i))+(1-FR);
P_ind(i)=prod((2*FR-1).*sp(:,bin_loc(i))+(1-FR),'all');
end



%% plot
% figure;
% ActN=sum(patt,1);
% cmap=hsv(max(ActN)+1);
% dot_color=cmap(ActN+1,:);
% 
% lines=scatter(P_real,P_ind,10,dot_color,'filled');
% %arrayfun(@(l,c) set(l,'Color',c{:}),lines,num2cell(dot_color,2))
% hold all
% plot([1e-8 1],[1e-8 1],'r')
% set(gca,'xscale','log','yscale','log')
end