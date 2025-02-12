function dff_robust_constant = get_robustdFF(STAmov,ftprnt,F0image)

F0=tovec(F0image);
STAmovVec=tovec(STAmov);

clf;
for n=1:size(ftprnt,3)
    
px=find(tovec(ftprnt(:,:,n)>0));
[~, maxfrm]=max(mean(STAmovVec(px,:),1,'omitnan'));

dF=STAmovVec(:,maxfrm);

F0_weight=F0.*(tovec(ftprnt(:,:,n)));
dF_weight=dF.*(tovec(ftprnt(:,:,n)));

px2=px(find(F0(px)>prctile(F0(px),30)));
px2_weight=px(find(F0_weight(px)>prctile(F0_weight(px),30) & F0_weight(px)<prctile(F0_weight(px),98)));

nexttile([1 1])
%plot(F0(px),dF(px),'.'); hold all
[p]=polyfit(F0(px2), dF(px2), 1);
[p_weight]=polyfit(F0_weight(px2_weight), dF_weight(px2_weight), 1);
% Get fitted values
y_fit = polyval(p, F0(px2));
y_fit_weight = polyval(p_weight, F0_weight(px2_weight));

scatter(F0(px),dF(px),10,'filled'); hold all
plot(F0(px2),y_fit,'r')
SS_res = sum((dF(px2) - y_fit).^2);
SS_tot = sum((dF(px2) - mean(dF(px2))).^2);
R_squared = 1 - (SS_res / SS_tot);

SS_res = sum((dF_weight(px2_weight) - y_fit_weight).^2);
SS_tot = sum((dF_weight(px2_weight) - mean(dF(px2_weight))).^2);
R_squared_weight = 1 - (SS_res / SS_tot);

title(['d(\DeltaF)/dF0 : ' num2str(p_weight(1),2) ', R^2 : ' num2str(R_squared_weight,2)])
Fslope(n)=p(1);
Fslope_weight(n)=p_weight(1);
Rsq(n)=R_squared;

%dff(n)=sum((dF).*(tovec(Result.ftprnt(:,:,n))>0),1,'omitnan')./sum((F0_filter).*(tovec(Result.ftprnt(:,:,n))>0),1,'omitnan');
%dff_weight(n)=sum((dF).*(tovec(ftprnt(:,:,n))),1,'omitnan')./sum(F0_weight,1,'omitnan');
weighted_dFSum=dF'*(tovec(ftprnt(:,:,n)));
dff_robust_constant(n)=(weighted_dFSum)/p(1);
end
