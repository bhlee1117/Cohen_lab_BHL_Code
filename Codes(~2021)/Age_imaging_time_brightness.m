%%
clear 
load(['C:\Users\Administrator\Dropbox\In vivo imaging\Analysis\','Mouse_age_imagingTime.mat'])
%%
for i=9:size(Mouse_number,2)
    [f p] = uigetfile('.tif',Mouse_number{1,i},'Multiselect','on');  % Cropped image folder
    pth{i,1}=p;
    for j=1:3
        clear im
        fnm{i,j}=f{1,j};
        for k=1:81
        im(:,:,k)=imread([pth{i,1} fnm{i,j}],k+81);
        end
    brightness(i,j)=mean(mean(mean(im,3)));
    end
end
%%
%load(['C:\Users\Administrator\Dropbox\In vivo imaging\Analysis\','Mouse_age_imagingTime.mat'])
for i=1:size(Mouse_number,2)
    for j=1:3
    norm_bright(i,j)=brightness(i,j)/laser_power(i,j);
    end
end
figure(1)
bar([1:size(Mouse_number,2)-1],mean(norm_bright(1:15,:),2))
hold all
bar([size(Mouse_number,2)],mean(norm_bright(end,:),2))
set(gca,'FontSize',8,'LineWidth',1,'XTick',[1:1:size(Mouse_number,2)],'XTickLabel',...
    Mouse_number);

ylabel('Normalized brightness (A.U.)')
errorbar([1:size(Mouse_number,2)],mean(norm_bright,2),std(norm_bright,0,2)/sqrt(3),'linestyle','none')
xtickangle(40)
figure(2)
plot(Age_image_dev(1,[1:6 8:10 13:15]),mean(norm_bright([1:6 8:10 13:15],:),2),'r.')
hold all
plot(Age_image_dev(1,[7 11 12]),mean(norm_bright([7 11 12],:),2),'b.')
plot(Age_image_dev(1,[16]),mean(norm_bright([16],:),2),'k.')
xlabel('Image day from the birth')
ylabel('Normalized brightness (A.U.)')