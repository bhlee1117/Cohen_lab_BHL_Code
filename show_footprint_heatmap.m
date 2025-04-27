function colored_image=show_footprint_heatmap(c_ftprnt,dat,dat_min,dat_max)
%figure; tiledlayout(2,1)
%a=randn(size(c_ftprnt,3),1); [~,b]=sort(a);
if size(c_ftprnt,3)~= length(dat)
    error('Check the size of two variables')
end

if nargin<3
    dat_min=min(dat(:));
    dat_max=max(dat(:));
end
c_ftprnt=mat2gray(c_ftprnt);
%dat_norm=(dat-min(dat))./(max(dat)-min(dat)); %normalize
b=[1:size(c_ftprnt,3)];
ax1=nexttile([1 1]);
colr=grs2rgb(dat,colormap(turbo),dat_min,dat_max);
%colr = flip(max(colormap(jet(size(c_ftprnt,3))),0),1); colr(colr>1)=1;
colored_image=squeeze(sum(c_ftprnt.*reshape(colr(b,:),1,1,[],3),3));
imshow2(colored_image,[]);
coord=get_coord(c_ftprnt);
%text(coord(:,1)',coord(:,2)',num2str([1:size(c_ftprnt,3)]'),'color','w')
%colormap('jet')
%colorbar('Ticks',[0:0.2:1],...
%         'TickLabels',num2str([min(dat):(max(dat)-min(dat))/5:max(dat)]','%2.2f'))
end