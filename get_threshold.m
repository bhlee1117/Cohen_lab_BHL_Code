function thres=get_threshold(volt_hp,ratio,bin)
if nargin<2
    ratio=2.7;
else if nargin<3
        bin=20;
end
end
%[cum bin]=histcounts(volt_hp-median(volt_hp),'Normalization','cumcount'); %cumulative count
%bin=mean([bin(1:end-1)' bin(2:end)'],2);
%fit_range=round(length(bin)*0.7);
%ft=fittype('c*0.5*(1+erf((x-a)/(b*sqrt(2))))'); %fit funtion; cumulative normal distribution
%2*cum(max(find(bin<0)))
%fitresult=fit(bin(1:fit_range),cum(1:fit_range)',ft,...
%    'lower',[-0.5 -10 0],'upper',[0.5 std(volt_hp)*10 inf],'StartPoint',[0 std(volt_hp) size(volt_hp,2)*0.8])
%fitresult=fit(bin(1:fit_range),cum(1:fit_range)',ft,...
%    'lower',[-0.5 std(volt_hp)*0.5 0],'upper',[0.5 std(volt_hp)*2 0.99*size(volt_hp,2)],'StartPoint',[0 std(volt_hp) size(volt_hp,2)*0.8]);
  
%plot(fitresult,bin,cum)
%s=sort(volt_hp,'descend');
%thres=s(round(size(volt_hp,2)-fitresult.c));
%volt_hp = volt_hp - movmean(volt_hp, 5);
volt_hp=volt_hp-movmedian(volt_hp,bin,2);
for i=1:size(volt_hp,1)    
t=find(volt_hp(i,:)<0);
s=std([volt_hp(i,t) -volt_hp(i,t)]);
thres(i,1)=s*ratio;
end
end

