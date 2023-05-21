clear
[fnm,pth]=uigetfile('*.tif','Select the raw image','Multiselect','on');
[fnmtr,pthtr]=uigetfile('*.mat','Select the tracking data','Multiselect','on');
Cell_size=[51 51 51];
cell_list_fnm=['H:\Image_data\OlympusTPM\20181008-20181012\OE9\Transformed_2\','Cell_list_T_sub_lip_reg_20181008_OE9_HC_.tif.mat'];
%%
if ischar(fnm)
    fnm={fnm};
end

if ischar(fnmtr)
    fnmtr={fnmtr};
end
h=waitbar(0,'zstack registration');

for i=1:length(fnm)
imginf=imfinfo([pth char(fnm{1,i})]);
numstack=size(imginf,1);  
sp_name=split(fnm{1,i},'_');
g=1;
clear lip_cred_tr lip_cred_tr_cov
for k=1:length(fnmtr)
    lip_cred_tr=cell(460,2);
    if ~isempty(findstr(fnmtr{1,k},sp_name{2,1}))
       lip_tr=importdata([pthtr fnmtr{1,k}]);
       for ptl=1:max(lip_tr.trackingData(:,1))
           if size(find(lip_tr.trackingData(:,1)==ptl),1)>10
               clear dif row
%                row=find(lip_tr.trackingData(:,1)==ptl);
%                dif(:,1)=lip_tr.trackingData(row,3)-lip_tr.trackingData(row(1,1),3);
%                dif(:,2)=lip_tr.trackingData(row,4)-lip_tr.trackingData(row(1,1),4);
% %                
% %                lip_cred_tr{1}(lip_tr.trackingData(find(lip_tr.trackingData(:,1)==ptl),2),g)=[NaN ; dif(2:end,1)];
% %                lip_cred_tr{2}(lip_tr.trackingData(find(lip_tr.trackingData(:,1)==ptl),2),g)=[NaN ; dif(2:end,2)];
row=find(lip_tr.trackingData(:,1)==ptl);
               lip_cred_tr{floor(lip_tr.trackingData(row(1,1),3)/Cell_size(1,1))*10+floor(lip_tr.trackingData(row(1,1),4)/Cell_size(1,1))+1,1}(lip_tr.trackingData(row,2),end+1)=[NaN ; diff(lip_tr.trackingData(row,3))];
               lip_cred_tr{floor(lip_tr.trackingData(row(1,1),3)/Cell_size(1,1))*10+floor(lip_tr.trackingData(row(1,1),4)/Cell_size(1,1))+1,2}(lip_tr.trackingData(row,2),end+1)=[NaN ; diff(lip_tr.trackingData(row,4))];
  
           end
       end
       
       lip_cred_tr_cov=lip_cred_tr;
       for j=1:size(lip_cred_tr,1)
       lip_cred_tr_cov{j,1}(lip_cred_tr{j,1}==0)=NaN;
       lip_cred_tr_cov{j,2}(lip_cred_tr{j,2}==0)=NaN;
       
%        lip_cred_tr_cov{1}(isnan(lip_cred_tr{1}))=0;
%        lip_cred_tr_cov{2}(isnan(lip_cred_tr{2}))=0;
       drift{j,1,i}=round(cumsum(mean(lip_cred_tr_cov{j,1},2,'omitnan'),'omitnan'));
       drift{j,2,i}=round(cumsum(mean(lip_cred_tr_cov{j,2},2,'omitnan'),'omitnan'));
       end
    end
end
waitbar((i)/length(fnm))
end
close(h)
%calculate drift
%%
load(cell_list_fnm)
clear map_lip map_filt im_exp crop_im Filt_im 
% map=zeros(10*Cell_size(1,1),ceil(size(Cells.list,1)/10)*Cell_size(1,2),Cell_size(:,3));
map_lip=zeros(10*Cell_size(1,1),ceil(size(Cells.list,1)/10)*Cell_size(1,2),Cell_size(:,3));
%%
for i=1:length(fnm)
    sp_name=split(fnm{1,i},'_');
    load(cell_list_fnm)
for j=1:Cells.zstack
    im(:,:,j)=imcrop(imread([pth char(fnm{1,i})],j),round(Cells.crop));
end

    im_exp=zeros(size(im,1)+2*floor(Cell_size(1,1)/2),size(im,2)+2*floor(Cell_size(1,2)/2),size(im,3)+2*floor(Cell_size(1,3)/2));
    im_exp(floor(Cell_size(1,2)/2)+1:floor(Cell_size(1,2)/2)+size(im,1),floor(Cell_size(1,1)/2)+1:floor(Cell_size(1,1)/2)+size(im,2)...
          ,floor(Cell_size(1,3)/2)+1:floor(Cell_size(1,3)/2)+size(im,3))=im;
    for k=1:size(Cell_size,2)
        Cells.list(:,k)=Cells.list(:,k)+floor(Cell_size(1,k)/2);
    end
    
for j=1:size(Cells.list,1)    
    clear crop_im max_crop_im;

    crop_im=im_exp(round(Cells.list(j,2))-floor(Cell_size(1,2)/2):round(Cells.list(j,2))+floor(Cell_size(1,2)/2),...
               round(Cells.list(j,1))-floor(Cell_size(1,1)/2):round(Cells.list(j,1))+floor(Cell_size(1,1)/2),...
               round(Cells.list(j,3))-floor(Cell_size(1,3)/2):round(Cells.list(j,3))+floor(Cell_size(1,3)/2));
           
           
im2=zeros(Cell_size(1,1)+abs(min(drift{j,2,i}(:,1)))*(1-heaviside(min(drift{j,2,i}(:,1))))+abs(max(drift{j,2,i}(:,1))),Cell_size(1,2)+abs(min(drift{j,1,i}(:,1)))*(1-heaviside(min(drift{j,1,i}(:,1))))+abs(max(drift{j,1,i}(:,1))),1);
for jj=1:size(drift{j,2,i},1)
im2(abs(max(drift{j,2,i}(:,1)))-drift{j,2,i}(jj,1)+1:abs(max(drift{j,2,i}(:,1)))-drift{j,2,i}(jj,1)+size(crop_im,1),...
   abs(max(drift{j,1,i}(:,1)))-drift{j,1,i}(jj,1)+1:abs(max(drift{j,1,i}(:,1)))-drift{j,1,i}(jj,1)+size(crop_im,2))=crop_im(:,:,jj);

cro_pos=[abs(max(drift{j,1,i}(:,1)))-drift{j,2,i}(1,1)+1 abs(max(drift{j,2,i}(:,1)))-drift{j,1,i}(1,1)+1 size(crop_im,1)-1 size(crop_im,2)-1];
crop_im_reg(:,:,jj)=imcrop(im2,cro_pos);
end
           if mod(j,10)==0
map_lip(1+(10-1)*Cell_size(1,2):(10)*Cell_size(1,2),1+floor(j/10)*Cell_size(1,1):(floor(j/10)+1)*Cell_size(1,1),1:size(crop_im,3))=crop_im_reg;               
           else
map_lip(1+(mod(j,10)-1)*Cell_size(1,2):(mod(j,10))*Cell_size(1,2),1+floor(j/10)*Cell_size(1,1):(floor(j/10)+1)*Cell_size(1,1),1:size(crop_im,3))=crop_im_reg;
           end
           
end
for k=1:Cell_size(1,3)
    imwrite(uint16(map_lip(:,:,k)),[pth,'\crop\Reg',char(sp_name(6,1)),'_',char(sp_name(7,1)),'.tif'],'Writemode','Append','Compression', 'none')
    %imwrite(uint16(map_filt(:,:,k)),[pth,'\crop\','Filt_',char(sp_name(6,1)),'_',char(sp_name(7,1)),'.tif'],'Writemode','Append','Compression', 'none')
end
end
