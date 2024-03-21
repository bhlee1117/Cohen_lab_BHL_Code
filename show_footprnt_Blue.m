function show_footprnt_Blue(c_ftprnt,mov_mc,Blue)


    colr = max(colormap(turbo(size(c_ftprnt,3))),0); colr(colr>1)=1;

%figure; tiledlayout(2,1)
%a=randn(size(c_ftprnt,3),1); [~,b]=sort(a);
c_ftprnt=mat2gray(c_ftprnt);
mov_mc=mat2gray(mov_mc);
b=[1:size(c_ftprnt,3)];
%ax1=nexttile([1 1]);
ax1=subplot(2,1,1);

imshow2(squeeze(sum(c_ftprnt.*reshape(colr(b,:),1,1,[],3),3)),[]);
coord=get_coord(c_ftprnt);
text(coord(:,1)',coord(:,2)',num2str([1:size(c_ftprnt,3)]'),'color','w')

ax2=subplot(2,1,2);
imshow2(imfuse(mean(mov_mc,3),sum(c_ftprnt,3)),[])
hold all
plot(Blue(:,2),Blue(:,1),'color',[0 0.6 1],'LineWidth',2)
linkaxes([ax1 ax2],'xy')


end