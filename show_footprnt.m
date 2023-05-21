function show_footprnt(c_ftprnt,mov_mc,colr)

if nargin<3
    colr = flip(max(colormap(jet(size(c_ftprnt,3))),0),1); colr(colr>1)=1;
end
%figure; tiledlayout(2,1)
%a=randn(size(c_ftprnt,3),1); [~,b]=sort(a);
b=[1:size(c_ftprnt,3)];
%ax1=nexttile([1 1]);
ax1=subplot(2,1,1);

imshow2(squeeze(sum(c_ftprnt.*reshape(colr(b,:),1,1,[],3),3)),[]);
coord=get_coord(c_ftprnt);
text(coord(:,1)',coord(:,2)',num2str([1:size(c_ftprnt,3)]'),'color','w')

ax2=subplot(2,1,2);
imshow2(imfuse(mean(mov_mc,3),sum(c_ftprnt,3)),[])
linkaxes([ax1 ax2],'xy')
end