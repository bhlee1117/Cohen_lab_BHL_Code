function show_footprnt_contour(c_ftprnt,mov_mc,colr)

if nargin<3
    colr = max(colormap(turbo(size(c_ftprnt,3))),0); colr(colr>1)=1;
end
%figure; tiledlayout(2,1)
%a=randn(size(c_ftprnt,3),1); [~,b]=sort(a);
c_ftprnt=mat2gray(c_ftprnt);
mov_mc=mat2gray(mov_mc);
b=[1:size(c_ftprnt,3)];
%ax1=nexttile([1 1]);
%ax1=subplot(2,1,1);
imshow2(mov_mc); hold all
for i=1:size(c_ftprnt,3)
ROIcontour = bwboundaries(c_ftprnt(:,:,i));
plot(ROIcontour{1}(:,2),ROIcontour{1}(:,1),'color',colr(i,:),'LineWidth',1.5)
end

coord=get_coord(c_ftprnt);
text(coord(:,1)',coord(:,2)',num2str([1:size(c_ftprnt,3)]'),'color','w')

end