%%
clear
data_file=['H:\Image_data\OlympusTPM\Invivo_images_CFC\','OE40.mat'];
image_file=['H:\Image_data\OlympusTPM\Invivo_images_CFC\OE40\Transformed2\crop\','OE40_CFC.tif'];
n=5;
session_to_anal=4;  %Ret3 = 5
%%
clear txn
i=1;
    mousedat{i}=importdata(data_file);
       g=1;
       thr_Data{i}=mousedat{1,i}.Cell.Data;
       tmp=thr_Data{i}(:,4:end);
       tmp(mousedat{1,i}.Cell.thres==0)=4;
       tmp(mousedat{1,i}.Cell.Data(:,4:end)==2)=2;
       thr_Data{i}(:,4:end)=tmp;
       
    if iscell(mousedat{1,i}.Cell.Data) || ~isempty(find(mousedat{1,i}.Cell.Data(:,4:end)~=1 & mousedat{1,i}.Cell.Data(:,4:end)~=2 & mousedat{1,i}.Cell.Data(:,4:end)~=3 & mousedat{1,i}.Cell.Data(:,4:end)~=0))
  
    else
 %omitna_data{i,1}=zeros(size(mousedat{1,i}.Cell.Data,1),12);
    for j=1:size(mousedat{1,i}.Cell.Data,1) %Cell 
        if sum(thr_Data{i}(j,:)==3 | thr_Data{i}(j,:)==4)==0
    omitna_data{i,1}(g,1:size(thr_Data{i},2))=thr_Data{i}(j,:);
    g=g+1;
        end
    end
    end

info=imfinfo(image_file);
numstack=size(imfinfo(image_file),1);
load(data_file)
g=1;
  
   cred_data=omitna_data{1,1};
for i=[1:session_to_anal]
   row=find(cred_data(:,i+3)==2);
   for j=1:size(row,1)
       txn{i}(j,1:3)=cred_data(row(j,1),1:3);
   end
   txn_shift=txn{i};
   txn_shift(:,:,3)=txn_shift+n;
   pt_im{i}=point_image(txn_shift,[info(1).Height info(1).Width size(info,1)+2*n],5);
end

%%
all_shift=cred_data(:,1:3);
all_shift(:,3)=all_shift(:,3)+n;
pt_all=point_image(all_shift,[info(1).Height info(1).Width size(info,1)+2*n],15);

% for i=n+1:n+numstack
%     for j=[1:5]
%     imwrite(uint16(100*pt_im{j}(:,:,i)),['E:\BACKUP\대학원\연구실\Lab meeting\20190107_presentation_발표자료\',num2str(j),'.tif'],'Writemode','Append','Compression','none')
% end
% end
% for j=1:6
% imwrite(uint16(100*max(pt_im{j},[],3)),['E:\BACKUP\대학원\연구실\Lab meeting\20200107_presentation_발표자료\',num2str(j),'.tif'],'Writemode','Append','Compression','none')
% end
% %%
% % for i=n+1:n+numstack
% % 
% %     imwrite(uint16(100*pt_all(:,:,i)),['E:\BACKUP\대학원\연구실\Lab meeting\20190107_presentation_발표자료\','All','.tif'],'Writemode','Append','Compression','none')
% % 
% % end
% imwrite(uint16(100*max(pt_all,[],3)),['E:\BACKUP\대학원\연구실\Lab meeting\20200107_presentation_발표자료\','All','.tif'],'Writemode','Append','Compression','none')
% %%
clear pt_im2 txn 
for i=1:3
   row=find(sum(cred_data(:,5:4+i)==2,2)==i);
   for j=1:size(row,1)
       txn{i}(j,1:3)=cred_data(row(j,1),1:3);
   end
   txn_shift=txn{i};
   txn_shift(:,:,3)=txn_shift+n;
   pt_im2{i}=point_image(txn_shift(:,:,1),[info(1).Height info(1).Width size(info,1)+2*n],15);
    imwrite(uint16(100*max(pt_im2{i},[],3)),['S:\Reports\20201117_LM자료\','CFC&Ret','.tif'],'Writemode','Append','Compression','none')
end
imwrite(uint16(100*max(pt_all,[],3)),['S:\Reports\20201117_LM자료\','CFC&Ret','.tif'],'Writemode','Append','Compression','none')