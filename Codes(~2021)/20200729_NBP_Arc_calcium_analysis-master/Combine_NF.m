
%%
% MODIFICATION HISTORY : 
%           Written by Byung Hun Lee, Deptartment of Physics and Astronomy,
%           Seoul National University, 2020/08/13

%% Load files : Image files, Cell list, Blinder 1, Blinder 2.
clear
[fnmb pthb]= uigetfile('.mat','Select the 1st binder"s NF_map','multiselect','on');  % Cropped image folder
[fnmb2 pthb2]= uigetfile('.mat','Select the 2st binder"s NF_map','multiselect','on');
[fnmim pthim]= uigetfile('.mat','Select the images','multiselect','on');

Class_folder='E:\BACKUP\대학원\연구실\MY_Projects\In_vivo_imaging\Cell_detection\TXN_detection_deep_learning\[Final]Classification_individual_cell_images';
Cell_size=[71 71 71];
Result_pth='\\NEUROBIOPHYSICS\Virtual Reality\VR_image_data';
ticks={'Day1Af','Day1Bf','Day2Af','Day2Bf','Day3Af','Day3Bf'};
%%
tic;load(fullfile(pthim,fnmim))
toc;
disp('Images and cell list are loaded')
for i=1:size(im_aligned,2)
    sp=split(fnmim,'.');
    sp_name=split(sp{1,1},'_');
    mkdir([Class_folder '\' sp_name{1,1} '\'])
    mkdir([Class_folder '\' sp_name{1,1} '\' ticks{1,i} '\TXN'])
    mkdir([Class_folder '\' sp_name{1,1} '\' ticks{1,i} '\Non_TXN'])
    mkdir([Class_folder '\' sp_name{1,1} '\' ticks{1,i} '\NA'])
end


    B1=importdata([pthb char(fnmb)]);
    B2=importdata([pthb2 char(fnmb2)]);
    for i=1:size(im_aligned,2)
        BL_list{1}(:,i)=B1.NF_map.list(:,B1.rewnsa(1,i)+3);
        BL_list{2}(:,i)=B2.NF_map.list(:,B2.rewnsa(1,i)+3);
    end
    prod=BL_list{1}.*BL_list{2};
    if ~isempty(find(prod==0))
         error('Check the list !!!! You might have 0')
    end
    prod(prod==2)=777;
    prod(prod==4)=2;
    prod(prod==9 | prod==6 | prod==3)=3;
    
    Data=[B1.NF_map.Cell_position prod];

%% Save consistent TXN cells
 clear im_exp
load(fullfile(pthim,fnmim),'cell_list')
for i=1:size(im_aligned,2)
    im_exp{i}=zeros(size(im_aligned{1},1)+2*floor(Cell_size(1,1)/2),size(im_aligned{1},2)+2*floor(Cell_size(1,2)/2),size(im_aligned{1},3)+2*floor(Cell_size(1,3)/2));
    im_exp{i}(floor(Cell_size(1,2)/2)+1:floor(Cell_size(1,2)/2)+size(im_aligned{1},1),floor(Cell_size(1,1)/2)+1:floor(Cell_size(1,1)/2)+size(im_aligned{1},2),...
           floor(Cell_size(1,3)/2)+1:floor(Cell_size(1,3)/2)+size(im_aligned{1},3))=im_aligned{i};
end

 for k=1:size(Cell_size,2)
        cell_list(:,k)=cell_list(:,k)+floor(Cell_size(1,k)/2);
 end
    OK_mask=Data(:,4:end)~=777;
for i=1:size(im_exp,2)
for j=find(OK_mask(:,i)==1)'
    clear crop_im max_crop_im mont filt_mont;
    
    crop_im=im_exp{i}(round(cell_list(j,2))-floor(Cell_size(1,2)/2):round(cell_list(j,2))+floor(Cell_size(1,2)/2),...
               round(cell_list(j,1))-floor(Cell_size(1,1)/2):round(cell_list(j,1))+floor(Cell_size(1,1)/2),...
               round(cell_list(j,3))-floor(Cell_size(1,3)/2):round(cell_list(j,3))+floor(Cell_size(1,3)/2));

max_crop_im=max(crop_im,[],3);

sw=Data(j,i+3);
if sw==1
    for jj=1:size(crop_im,3)
        if jj==1
   imwrite(uint16(crop_im(:,:,jj)),[Class_folder '\' sp_name{1,1} '\' ticks{1,i},'\Non_TXN\',num2str(j),'.tif'],'Compression','none');
        else
   imwrite(uint16(crop_im(:,:,jj)),[Class_folder '\' sp_name{1,1} '\' ticks{1,i},'\Non_TXN\',num2str(j),'.tif'],'Writemode','Append','Compression','none');
        end
    end
end

if sw==2
   for jj=1:size(crop_im,3)
       if jj==1
   imwrite(uint16(crop_im(:,:,jj)),[Class_folder '\' sp_name{1,1} '\' ticks{1,i},'\TXN\',num2str(j),'.tif'],'Compression','none');
       else
   imwrite(uint16(crop_im(:,:,jj)),[Class_folder '\' sp_name{1,1} '\' ticks{1,i},'\TXN\',num2str(j),'.tif'],'Writemode','Append','Compression','none');
       end
   end
end

if sw==3
   for jj=1:size(crop_im,3)
       if jj==1
   imwrite(uint16(crop_im(:,:,jj)),[Class_folder '\' sp_name{1,1} '\' ticks{1,i},'\NA\',num2str(j),'.tif'],'Compression','none');        
       else
   imwrite(uint16(crop_im(:,:,jj)),[Class_folder '\' sp_name{1,1} '\' ticks{1,i},'\NA\',num2str(j),'.tif'],'Writemode','Append','Compression','none');
       end
   end
end
if sw==4
    close(figure(1))
    close(figure(2))
    break
end
end
end
    %% Classify again
 load(fullfile(pthim,fnmim),'cell_list')
 OK_mask=Data(:,4:end)~=777;
 recheck_mask=1-OK_mask;
    n=10;
 xx={'TXN','No TXN'};
     handles.fig = figure(1);
set(gcf, 'Position',  [10, 10, 1800 950])

 for k=1:size(Cell_size,2)
        cell_list(:,k)=cell_list(:,k)+floor(Cell_size(1,k)/2);
 end
    for i=2:size(im_aligned,2)%%%%
        i
    list_recheck=find(recheck_mask(:,i)==1);
    j=1;
while j<=size(list_recheck,1)
    
    clear crop_im max_crop_im mont filt_mont;
clf('reset')
    crop_im=im_exp{i}(round(cell_list(list_recheck(j,1),2))-floor(Cell_size(1,2)/2):round(cell_list(list_recheck(j,1),2))+floor(Cell_size(1,2)/2),...
               round(cell_list(list_recheck(j,1),1))-floor(Cell_size(1,1)/2):round(cell_list(list_recheck(j,1),1))+floor(Cell_size(1,1)/2),...
               round(cell_list(list_recheck(j,1),3))-floor(Cell_size(1,3)/2):round(cell_list(list_recheck(j,1),3))+floor(Cell_size(1,3)/2));
    max_crop_im=max(crop_im,[],3);

for jj=1:size(crop_im,3)
Bgimfilt=imgaussfilt(crop_im(:,:,jj), 6);
subim=abs(crop_im(:,:,jj)-Bgimfilt);
Filt_im(:,:,jj)=imgaussfilt(subim,1.5);
end %make filtered image

for jj=1:size(crop_im,3)  % make montage image
              if mod(jj,n)==0
mont(1+(floor(jj/n)-1)*Cell_size(1,1):(floor(jj/n))*Cell_size(1,1),1+(n-1)*Cell_size(1,2):(n)*Cell_size(1,2))=crop_im(:,:,jj);
filt_mont(1+(floor(jj/n)-1)*Cell_size(1,1):(floor(jj/n))*Cell_size(1,1),1+(n-1)*Cell_size(1,2):(n)*Cell_size(1,2))=Filt_im(:,:,jj);
           else
mont(1+floor(jj/n)*Cell_size(1,1):(floor(jj/n)+1)*Cell_size(1,1),1+(mod(jj,n)-1)*Cell_size(1,2):(mod(jj,n))*Cell_size(1,2))=crop_im(:,:,jj);
filt_mont(1+floor(jj/n)*Cell_size(1,1):(floor(jj/n)+1)*Cell_size(1,1),1+(mod(jj,n)-1)*Cell_size(1,2):(mod(jj,n))*Cell_size(1,2))=Filt_im(:,:,jj);
           end
end
clf('reset')

handles.axes2 = axes('Units','pixels','Position',[30 30 1550 900]);
mont(mont<0)=0;
imagesc(mont)
%title(['1^s^t blinder says   ' xx{1,BL_list{1}(list_recheck(j,1),i)} '  While, 2^n^d blinder says' xx{1,BL_list{2}(list_recheck(j,1),i)}])
colormap('gray')
axis equal tight off
hold on
handles.axes3 = axes('Units','pixels','Position',[1500 50 300 300]);
imagesc(max_crop_im)
axis equal tight off
colormap('gray')
handles.axes4 = axes('Units','pixels','Position',[1500 450 300 300]);
imagesc(max(Filt_im,[],3))
axis equal tight off
colormap('gray')
hold all
sw=input(['( ' num2str(list_recheck(j,1)) '/' num2str(size(cell_list,1)) ' )th cell. Is the cell has TXN site?\n']);
if sw==1 % non txn
    for jj=1:size(crop_im,3)
        if jj==1
   imwrite(uint16(crop_im(:,:,jj)),[Class_folder '\' sp_name{1,1} '\' ticks{1,i},'\Non_TXN\',num2str(list_recheck(j,1)),'.tif'],'Compression','none');
        else
   imwrite(uint16(crop_im(:,:,jj)),[Class_folder '\' sp_name{1,1} '\' ticks{1,i},'\Non_TXN\',num2str(list_recheck(j,1)),'.tif'],'Writemode','Append','Compression','none');
        end
    end
    j=j+1;
    Data(list_recheck(j-1,1),i+3)=sw;
end

if sw==2 % TXN
   for jj=1:size(crop_im,3)
       if jj==1
           imwrite(uint16(crop_im(:,:,jj)),[Class_folder '\' sp_name{1,1} '\' ticks{1,i},'\TXN\',num2str(list_recheck(j,1)),'.tif'],'Compression','none');
       else
           imwrite(uint16(crop_im(:,:,jj)),[Class_folder '\' sp_name{1,1} '\' ticks{1,i},'\TXN\',num2str(list_recheck(j,1)),'.tif'],'Writemode','Append','Compression','none');
       end
   end
   j=j+1;
   Data(list_recheck(j-1,1),i+3)=sw;
end

if sw==3 %NA
   for jj=1:size(crop_im,3)
       if jj==1
   imwrite(uint16(crop_im(:,:,jj)),[Class_folder '\' sp_name{1,1} '\' ticks{1,i},'\NA\',num2str(list_recheck(j,1)),'.tif'],'Compression','none');        
       else
   imwrite(uint16(crop_im(:,:,jj)),[Class_folder '\' sp_name{1,1} '\' ticks{1,i},'\NA\',num2str(list_recheck(j,1)),'.tif'],'Writemode','Append','Compression','none');
       end
   end
   j=j+1;
   Data(list_recheck(j-1,1),i+3)=sw;
end
if sw==4 %Break
    close(figure(1))
    close(figure(2))
    break
end

if sw==5 %get back
    j=j-1;
     if Data(list_recheck(j,1),i+3)==1
        delete([Class_folder '\' sp_name{1,1} '\' ticks{1,i},'\Non_TXN\',num2str(list_recheck(j,1)),'.tif'])
    end
    
    if Data(list_recheck(j,1),i+3)==2
        delete([Class_folder '\' sp_name{1,1} '\' ticks{1,i},'\TXN\',num2str(list_recheck(j,1)),'.tif'])
    end
    
    if Data(list_recheck(j,1),i+3)==3
        delete([Class_folder '\' sp_name{1,1} '\' ticks{1,i},'\NA\',num2str(list_recheck(j,1)),'.tif'])
    end
end
end
if sw==4
    break
end
    end
save([Result_pth '\' sp_name{1,1} '_Arc-TXN.mat'],'Data')

%%
