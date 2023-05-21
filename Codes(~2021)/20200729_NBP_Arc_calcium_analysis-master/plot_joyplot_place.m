function plot_joyplot_place(place_cdf,ll)
for i=1:size(ll,1)
try 
int_cdf(i,:)=interp1(place_cdf(ll(i,1):ll(i,2),1)+2000,place_cdf(ll(i,1):ll(i,2),2),[0:5:3000]);
catch
m=-2001;
tmp=place_cdf(ll(i,1):ll(i,2),1);
for j=1:size(tmp,1)
if tmp(j,1)>m
m=tmp(j,1);
else
m=m+10;
tmp(j,1)=m;
end
end
int_cdf(i,:)=interp1(tmp+2000,place_cdf(ll(i,1):ll(i,2),2),[0:5:3000]);
end
end
%int_cdf=int_cdf-min(int_cdf,[],2);
%%
int_cdf(isnan(int_cdf))=0;
joyPlot(int_cdf',[0:5:3000],0.25,'FaceColor',[1:37])
colormap('jet')
axis tight off
colorbar
end
