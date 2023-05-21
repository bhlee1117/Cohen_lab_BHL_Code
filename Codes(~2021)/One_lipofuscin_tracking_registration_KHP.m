clear
[fnm,pth]=uigetfile('*.tif','Select the raw image','Multiselect','on');
[fnmtr,pthtr]=uigetfile('*.mat','Select the tracking data','Multiselect','on');
channel=2;
%%
if ischar(fnm)
    fnm={fnm};
end

if ischar(fnmtr)
    fnmtr={fnmtr};
end
h=waitbar(0,'Drift calculation');

for i=1:length(fnm)
imginf=imfinfo([pth char(fnm{1,i})]);
numstack=size(imginf,1);  
sp_name=split(fnm{1,i},'_');
g=1;
clear lip_cred_tr lip_cred_tr_cov
for k=1:length(fnmtr)
    if ~isempty(findstr(fnmtr{1,k},sp_name{3,1}))
       lip_tr=importdata([pthtr fnmtr{1,k}]);
       for ptl=1:max(lip_tr.trackingData(:,1))
           if size(find(lip_tr.trackingData(:,1)==ptl),1)>3
               clear dif row
%                row=find(lip_tr.trackingData(:,1)==ptl);
%                dif(:,1)=lip_tr.trackingData(row,3)-lip_tr.trackingData(row(1,1),3);
%                dif(:,2)=lip_tr.trackingData(row,4)-lip_tr.trackingData(row(1,1),4);
% %                
% %                lip_cred_tr{1}(lip_tr.trackingData(find(lip_tr.trackingData(:,1)==ptl),2),g)=[NaN ; dif(2:end,1)];
% %                lip_cred_tr{2}(lip_tr.trackingData(find(lip_tr.trackingData(:,1)==ptl),2),g)=[NaN ; dif(2:end,2)];

               lip_cred_tr{1}(lip_tr.trackingData(find(lip_tr.trackingData(:,1)==ptl),2),g)=[NaN ; diff(lip_tr.trackingData(find(lip_tr.trackingData(:,1)==ptl),3))];
               lip_cred_tr{2}(lip_tr.trackingData(find(lip_tr.trackingData(:,1)==ptl),2),g)=[NaN ; diff(lip_tr.trackingData(find(lip_tr.trackingData(:,1)==ptl),4))];
                g=g+1;
           end
       end
       lip_cred_tr_cov=lip_cred_tr;
       lip_cred_tr_cov{1}(lip_cred_tr{1}==0)=NaN;
       lip_cred_tr_cov{2}(lip_cred_tr{2}==0)=NaN;
       
       lip_cred_tr_cov{1}(isnan(lip_cred_tr{1}))=0;
       lip_cred_tr_cov{2}(isnan(lip_cred_tr{2}))=0;
       drift{i,1}(1:numstack/channel,:)=round(cumsum(mean(lip_cred_tr_cov{1}(1:numstack/channel,:),2,'omitnan'),'omitnan'));
       drift{i,2}(1:numstack/channel,:)=round(cumsum(mean(lip_cred_tr_cov{2}(1:numstack/channel,:),2,'omitnan'),'omitnan'));
    end
end
waitbar((i)/length(fnm))
end
close(h)
%calculate drift
%%
h=waitbar(0,'Registration');

for i=1:length(fnm)
    clear im2
    imginf=imfinfo([pth char(fnm{1,i})]);
numstack=size(imginf,1);  
for c=1:channel
im2{c}=zeros(imginf(1).Height+abs(min(drift{i,2}(:,1)))*(1-heaviside(min(drift{i,2}(:,1))))+abs(max(drift{i,2}(:,1))),imginf(1).Width+abs(min(drift{i,1}(:,1)))*(1-heaviside(min(drift{i,1}(:,1))))+abs(max(drift{i,1}(:,1))),numstack/channel);
end
for j=1:numstack/channel   

im=zeros(imginf(1).Height+abs(min(drift{i,2}(:,1)))*(1-heaviside(min(drift{i,2}(:,1))))+abs(max(drift{i,2}(:,1))),imginf(1).Width+abs(min(drift{i,1}(:,1)))*(1-heaviside(min(drift{i,1}(:,1))))+abs(max(drift{i,1}(:,1))),1);
im(abs(max(drift{i,2}(:,1)))-drift{i,2}(j,1)+1:abs(max(drift{i,2}(:,1)))-drift{i,2}(j,1)+imginf(1).Height,...
   abs(max(drift{i,1}(:,1)))-drift{i,1}(j,1)+1:abs(max(drift{i,1}(:,1)))-drift{i,1}(j,1)+imginf(1).Width,1)=imread([pth char(fnm{1,i})],j);

for c=1:channel
im2{c}(abs(max(drift{i,2}(:,1)))-drift{i,2}(j,1)+1:abs(max(drift{i,2}(:,1)))-drift{i,2}(j,1)+imginf(1).Height,...
   abs(max(drift{i,1}(:,1)))-drift{i,1}(j,1)+1:abs(max(drift{i,1}(:,1)))-drift{i,1}(j,1)+imginf(1).Width,j)=imread([pth char(fnm{1,i})],(numstack/channel)*(c-1)+j);
end

   
%     imwrite(uint16(imcrop(im,cro_pos)),[pth 'lip_reg_' char(fnm{1,i})],'Writemode','Append','Compression','none');
%imwrite(uint16(mask(:,:,j)),[pth 'mask_' char(fnm{1,1})],'Writemode','Append','Compression','none');

end   

cro_pos=[abs(max(drift{i,1}(:,1)))-drift{i,2}(1,1)+1 abs(max(drift{i,2}(:,1)))-drift{i,1}(1,1)+1 imginf(1).Height-1 imginf(1).Width-1];
for c=1:channel
for j=1:numstack/channel
%     im2(abs(max(drift{i,2}(:,1)))-drift{i,2}(j,1)+1:abs(max(drift{i,2}(:,1)))-drift{i,2}(j,1)+imginf(1).Height,...
%    abs(max(drift{i,1}(:,1)))-drift{i,1}(j,1)+1:abs(max(drift{i,1}(:,1)))-drift{i,1}(j,1)+imginf(1).Width,2)=imread([pth char(fnm{1,i})],numstack/2+j);
%      imwrite(uint16(imcrop(im2(:,:,2),cro_pos)),[pth 'lip_reg_' char(fnm{1,i})],'Writemode','Append','Compression','none');

imwrite(uint16(imcrop(im2{c}(:,:,j),cro_pos)),[pth 'Reg_' char(fnm{1,i})],'Writemode','Append','Compression','none');
end
end
waitbar((i)/length(fnm))
end
close(h)