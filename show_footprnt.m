function show_footprnt(c_ftprnt,mov_mc,colr)
coord=get_coord(c_ftprnt>0);
if nargin<3
    colr = max(colormap(turbo(size(c_ftprnt,3))),0); colr(colr>1)=1;
end
%figure; tiledlayout(2,1)
%a=randn(size(c_ftprnt,3),1); [~,b]=sort(a);
c_ftprnt=toimg(rescale2(tovec(c_ftprnt),1),size(c_ftprnt,1),size(c_ftprnt,2));
c_ftprnt=mat2gray(c_ftprnt);
mov_mc=mat2gray(mov_mc);
b=[1:size(c_ftprnt,3)];
ax1=nexttile([1 1]);
%ax1=subplot(2,1,1);

imshow2(squeeze(sum(c_ftprnt.*reshape(colr(b,:),1,1,[],3),3)),[]);
text(coord(:,1)',coord(:,2)',num2str([1:size(c_ftprnt,3)]'),'color','w')

%ax2=subplot(2,1,2);
ax2=nexttile([1 1]);
imshow2(imfuse(mean(mov_mc,3),sum(c_ftprnt,3)),[])
linkaxes([ax1 ax2],'xy')
end