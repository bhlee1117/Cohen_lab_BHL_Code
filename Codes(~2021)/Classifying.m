clear 
[fnm pth]=uigetfile('.tif','Pick cropped image files to analysis');  % Cropped image folder
load('OE9_CFC-1_2018-m12-d11-13h21_TNT.mat')

%%
imginf=imfinfo([pth fnm]);
numstack=size(imginf,1);

for i=1:numstack
    im(:,:,i)=imread([pth fnm],i);
end
im_max=max(im,[],3);
%%
% threshold = [276 469;... % Sigma xy
%             700 1383;... % Simga Z
%             400 1000;... % Amp
%             150 700];    % Amp-bg
threshold = [0 9;... % Sigma xy
            0.500 10000;... % Amp
            150 700];    % Amp-bg
clear cred_track cred incred_track incred
g=1;
gg=1;
n=5;
% for i=1:size(refinementData{1,1},1)
%     if threshold(1,1)<mean(refinementData{1,1}(i,6:7)) && threshold(1,2)>mean(refinementData{1,1}(i,6:7))   %Sigma xy
%         if threshold(2,1)<refinementData{1,1}(i,4) && threshold(2,2)>refinementData{1,1}(i,4)
% %             if threshold(3,1)<refinementData{1,1}(i,4)-refinementData{1,1}(i,5) && threshold(3,2)>refinementData{1,1}(i,)(i,4)
% %                 if threshold(4,1)<refinementData{1,1}(i,)(i,4)-refinementData{1,1}(i,)(i,5) && threshold(4,2)>refinementData{1,1}(i,)(i,4)-refinementData{1,1}(i,)(i,5)
%                     cred(g,1:8)=refinementData{1,1}(i,:);
%                     g=g+1;
%                 end
%             end
%         end

for ptl=1:max(trackingData(:,1))
    ptl_row=find(trackingData(:,1)==ptl);
if size(find(trackingData(:,1)==ptl),1)>n
    if max(trackingData(ptl_row,6))>400  && mean(trackingData(ptl_row(1,1)+round(size(ptl_row,1)/2)-1:ptl_row(1,1)+round(size(ptl_row,1)/2)+1,8))<3 && sum(trackingData(ptl_row,6)<trackingData(ptl_row,7))<3
    if sum(trackingData(ptl_row,8)<1.5) && sum(trackingData(ptl_row,8)>1.2 & trackingData(ptl_row,8)<1.9)
        if mod(round(trackingData(ptl_row(1,1),3)),51)>6 && mod(round(trackingData(ptl_row(1,1),4)),51)>6 && mod(round(trackingData(ptl_row(1,1),3)),51)<45 && mod(round(trackingData(ptl_row(1,1),4)),51)<45
cred_track{g}=trackingData(find(trackingData(:,1)==ptl),1:8);
    cred(g,1:2)=[mean(cred_track{g}(:,3)) mean(cred_track{g}(:,4))];
    g=g+1;
        end
    end
    end
else
    incred_track{gg}=trackingData(find(trackingData(:,1)==ptl),1:8);
    incred(gg,1:2)=[mean(incred_track{gg}(:,3)) mean(incred_track{gg}(:,4))];
    gg=gg+1;
end
end
%         
%     end
% end

figure(1)
colormap('gray')
imagesc(im_max)
hold all
axis equal
plot(cred(:,1),cred(:,2),'ro')

%plot(incred(:,1),incred(:,2),'go')

for i=1:size(cred,1)
    text(cred(i,1),cred(i,2)+5,num2str(i),'color','r')
end
% 
% for i=1:size(incred,1)
%     text(incred(i,1),incred(i,2)+5,num2str(i),'color','g')
% end

%%
n=6;
resized_map=zeros(size(im,1)/n,size(im,2)/n);
for i=1:size(refinementData,1)
    clear row
    row=refinementData{i,1}(:,4)>200 &refinementData{i,1}(:,6)<2 & refinementData{i,1}(:,4)>refinementData{i,1}(:,5) &...
        mod(round(refinementData{i,1}(:,1)),51)>6 & mod(round(refinementData{i,1}(:,2)),51)>6 & mod(round(refinementData{i,1}(:,1)),51)<45 & mod(round(refinementData{i,1}(:,2)),51)<45 ;
 for j=1:size(row,1)
   if row(j,1)
     resized_map(round(refinementData{i,1}(j,2)/n),round(refinementData{i,1}(j,1)/n))=resized_map(round(refinementData{i,1}(j,2)/n),round(refinementData{i,1}(j,1)/n))+1;
   end
 end
end
imwrite(uint16(resizem(resized_map,n)),'resize_map.tif')
% imagesc(resized_map)

% axis equal
% alpha(0.5)
% hold all
% imagesc(im_max)
% colormap('gray')
% axis equal