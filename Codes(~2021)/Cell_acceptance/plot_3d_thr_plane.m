function plot_3d_thr_plane(Class_dat,f,thres,mov)
clear im
% for i=1:4
%plot3(Class_dat(:,f(1,1)),Class_dat(:,f(1,2)),Class_dat(:,f(1,3)),'color',[1 0 0],'marker','.','linestyle','none')
plot3(Class_dat(:,f(1,1)),Class_dat(:,f(1,2)),Class_dat(:,f(1,3)),'marker','.','linestyle','none')

% hold all
% end
% for i=1:2
% plot3(Class_dat(888+328*(i-1)+1:888+328*(i),f(1,1)),Class_dat(888+328*(i-1)+1:888+328*(i),f(1,2)),Class_dat(888+328*(i-1)+1:888+328*(i),f(1,3)),'color',[0 0 1],'marker','.','linestyle','none')
% end
patch( [140 0 0 140] , [thres(1,2) thres(1,2) thres(1,2) thres(1,2)], [300 300 0 0], [0.4 0.7 0.4])
alpha(0.3)
patch( zeros(1,4)+thres(1,1) , [45 0 0 45], [300 300 0 0], [0.4 0.7 0.4])
alpha(0.3)
patch( [0 140 140 0] , [0 0 45 45], zeros(1,4)+thres(1,3), [0.4 0.7 0.4])
alpha(0.3)
grid on
xlabel('Number of lipofuscin')
ylabel('Mean intensity of lipofuscin / mean Green')
zlabel('Top 10 % mean intensity / mean Green')
xlim([0 140])
ylim([0 15])
zlim([0 100])
view([1,45])
if mov==1
F = getframe;
k=1;
for az=1:135
    view([az,45])
    F = getframe(gcf);
    for j=1:3
    im(:,:,j,k) = imresize(F.cdata(:,:,j),[334 428]);
    end
    k=k+1;
end
for la=45:-0.5:10
    view([az,la])
    F = getframe(gcf);
    for j=1:3
    im(:,:,j,k) = imresize(F.cdata(:,:,j),[334 428]);
    end
    k=k+1;
end
im=reshape(im,334,428,1,3*(k-1));
end
%imwrite(im,'Animation.gif','DelayTime',0.05,'LoopCount',inf)