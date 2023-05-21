%% Evaluating the intensity profile of TXN site.
clear
load(['E:\BACKUP\대학원\연구실\MY_Projects\In_vivo_imaging\Cell_detection\TXN_detection_deep_learning\Short_Term\Test\1.mat'])
for i=1:12
    im(:,:,i)=imread(['E:\BACKUP\대학원\연구실\MY_Projects\In_vivo_imaging\Cell_detection\TXN_detection_deep_learning\Short_Term\Test\','1.tif']);
end

%%
N=3;
n=3;
fit_x=2*N+3;
fit_y=2*N+3;
for i=1:9
    if mod(i,n)==0
mont(1+(floor(i/n)-1)*(2*N-1):(floor(i/n))*(2*N-1),1+(n-1)*(2*N-1):(n)*(2*N-1))=im(round(data.tr(i,2))-ceil(N/2):round(data.tr(i,2))+ceil(N/2),round(data.tr(i,1))-ceil(N/2):round(data.tr(i,1))+ceil(N/2));
    else
mont(1+floor(i/n)*(2*N-1):(floor(i/n)+1)*(2*N-1),1+(mod(i,n)-1)*(2*N-1):(mod(i,n))*(2*N-1))=im(round(data.tr(i,2))-ceil(N/2):round(data.tr(i,2))+ceil(N/2),round(data.tr(i,1))-ceil(N/2):round(data.tr(i,1))+ceil(N/2));
    end
background=mean(min_rank(reshape(mon(:,:,i),(2*N-1)^2,1),N));
int(i,1)=sum(sum(im(round(data.tr(i,2))-ceil(N/2):round(data.tr(i,2))+ceil(N/2),round(data.tr(i,1))-ceil(N/2):round(data.tr(i,1))+ceil(N/2))))-background*(2*ceil(N/2)+1)^2;

end
imagesc(mont)
axis equal
colormap('gray')