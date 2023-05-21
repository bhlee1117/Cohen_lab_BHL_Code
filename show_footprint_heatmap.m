function show_footprint_heatmap(c_ftprnt,dat)
%figure; tiledlayout(2,1)
%a=randn(size(c_ftprnt,3),1); [~,b]=sort(a);
if size(c_ftprnt,3)~= length(dat)
    error('Check the size of two variables')
end
dat_norm=(dat-min(dat))./(max(dat)-min(dat)); %normalize
b=[1:size(c_ftprnt,3)];
ax1=nexttile([1 1]);
binsize=50;
cmap=jet(binsize+1);
colr=cmap(ceil(dat_norm./((max(dat_norm)-min(dat_norm))/binsize))+1,:);
%colr = flip(max(colormap(jet(size(c_ftprnt,3))),0),1); colr(colr>1)=1;
imshow2(squeeze(sum(c_ftprnt.*reshape(colr(b,:),1,1,[],3),3)),[]);
coord=get_coord(c_ftprnt);
text(coord(:,1)',coord(:,2)',num2str([1:size(c_ftprnt,3)]'),'color','w')
colormap('jet')
colorbar('Ticks',[0:0.2:1],...
         'TickLabels',num2str([min(dat):(max(dat)-min(dat))/5:max(dat)]','%2.2f'))
end