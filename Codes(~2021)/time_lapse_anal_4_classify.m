%% Classify the cell in the first and second frame
clear
[Filename PathName] = uigetfile('.tif','Select the time-lapse image files to analysis');  % Cropped image folder
Class_folder='E:\BACKUP\대학원\연구실\MY_Projects\In_vivo_imaging\Cell_detection\TXN_detection_deep_learning\Classification\OE9\';
Cell_list_file='H:\Image_data\OlympusTPM\20181008-20181012\OE9\Mis\Cell_list_subreg_20181123_OE9_ST_.tif.mat';
load(Cell_list_file)
Cell_size=[50 50 ];

%% 

    sp_name=split(Filename,'_');
    mkdir([Class_folder sp_name{4,1}])
    mkdir([Class_folder sp_name{4,1} '\TXN'])
    mkdir([Class_folder sp_name{4,1} '\Non_TXN'])
    mkdir([Class_folder sp_name{4,1} '\NA'])

%%
dir_ts=dir([Class_folder sp_name{4,1} '\TXN']);
dir_no=dir([Class_folder sp_name{4,1} '\Non_TXN']);
dir_NA=dir([Class_folder sp_name{4,1} '\NA']);

dir_pre_cells=[dir_no(3:size(dir_no,1)); dir_ts(3:size(dir_ts,1)) ; dir_NA(3:size(dir_NA,1))];
fnd_result_file=dir(Class_folder);
if ~isempty(cell2mat(strfind({fnd_result_file.name},['NF_map_',sp_name{5,1},'.mat'])))
load([Class_folder,'NF_map_',sp_name{5,1},'.mat'])
else
    NF_map.list=[];
end
if ~size(dir_pre_cells,1)==0  % Find the date in the categorized folder
    i=1;
    while ~isempty(cell2mat(strfind({dir_pre_cells.name},[num2str(i) ,'.tif']))) && i<=size(Cells.list,1)
        i=i+1;
        start_here=i;
    end
else
    start_here=1;
end
%%
% engram_class.filename=img_name;
% if strfind(img_info{4,1},'HC')
% engram_class.train=0;
% else if strfind(img_info{4,1},'CFC') engram_class.train=1;
% else if strfind(img_info{4,1},'Ret') engram_class.train=2;
% end
% end
% end
% clear tar_im
% % tar_im(:,:,1)=imread(target_RGB_Image,1);
% % tar_im(:,:,2)=imread(target_RGB_Image,2);
% % tar_im(:,:,3)=imread(target_RGB_Image,3);

%%
    clear im_exp
    load(Cell_list_file)
    NF_map.Cell_position=Cells.list;
    for i=1:Cells.time
    for k=1:Cells.zstack
    im{i}(:,:,k)=imread(fullfile(PathName,Filename),k+(i-1)*Cells.zstack);
    end
    im_exp{i}=zeros(size(im{i},1)+2*floor(Cell_size(1,1)/2),size(im{i},2)+2*floor(Cell_size(1,2)/2),size(im{i},3)+2*floor(Cell_size(1,3)/2));
    im_exp{i}(floor(Cell_size(1,2)/2)+1:floor(Cell_size(1,2)/2)+size(im{i},1),floor(Cell_size(1,1)/2)+1:floor(Cell_size(1,1)/2)+size(im{i},2),floor(Cell_size(1,3)/2)+1:floor(Cell_size(1,3)/2)+size(im{i},3))=im{i};
    end
    %%
     load(Cell_list_file)
    NF_map.Cell_position=Cells.list;
    for k=1:size(Cell_size,2)
        Cells.list(:,k)=Cells.list(:,k)+floor(Cell_size(1,k)/2);
    end
for j=start_here:size(Cells.list,1)    
    clear crop_im max_crop_im;
for i=1:Cells.time
    crop_im{i}=im_exp{i}(round(Cells.list(j,2))-floor(Cell_size(1,2)/2):round(Cells.list(j,2))+floor(Cell_size(1,2)/2),...
               round(Cells.list(j,1))-floor(Cell_size(1,1)/2):round(Cells.list(j,1))+floor(Cell_size(1,1)/2),...
               round(Cells.list(j,3))-floor(Cell_size(1,3)/2):round(Cells.list(j,3))+floor(Cell_size(1,3)/2));

           max_crop_im(:,:,i)=max(crop_im{i},[],3);
           if i>1
%max_crop_im(:,:,i)=max_crop_im(:,:,i)*median(reshape(max_crop_im(:,:,1),size(max_crop_im,1)*size(max_crop_im,2),1))/median(reshape(max_crop_im(:,:,i),size(max_crop_im,1)*size(max_crop_im,2),1));
max_crop_im(:,:,i)=imhistmatch(uint16(max_crop_im(:,:,i))*10,uint16(max_crop_im(:,:,1)*10))/10;
           end
           
Bgimfilt=imgaussfilt(max_crop_im(:,:,i), 4);
subim=max_crop_im(:,:,i)-Bgimfilt;
Filt_im(:,:,i)=imgaussfilt(subim,0.7);
end           
figure(1)
montage(max_crop_im,'DisplayRange', [0 2500])
colormap('gray')
axis equal
hold on
title('Subtracted Image')

figure(2)
montage(Filt_im,'DisplayRange', [0 1000])
colormap('gray')
axis equal
title('Filtered Image')
figure(3)
subplot(1,1,1)
imagesc(Cells.max_image)
colormap('gray')
axis equal
hold all
plot(NF_map.Cell_position(j,1),NF_map.Cell_position(j,2),'ro','markersize',50)

sw=menu([num2str(j),'th cell, Is the cell has TXN site?'],'Non TS','TS','N/A Cell','Break');

if sw==1
  for t=1:Cells.time
   if t==1
          imwrite(uint16(max_crop_im(:,:,t)),[Class_folder,char(sp_name(4,1)),'\Non_TXN\',num2str(j),'.tif'],'Compression','none');
      else
   imwrite(uint16(max_crop_im(:,:,t)),[Class_folder,char(sp_name(4,1)),'\Non_TXN\',num2str(j),'.tif'],'Writemode','Append','Compression','none');
      end
  end
end

if sw==2
     for t=1:Cells.time
   if t==1
          imwrite(uint16(max_crop_im(:,:,t)),[Class_folder,char(sp_name(4,1)),'\TXN\',num2str(j),'.tif'],'Compression','none');
      else
   imwrite(uint16(max_crop_im(:,:,t)),[Class_folder,char(sp_name(4,1)),'\TXN\',num2str(j),'.tif'],'Writemode','Append','Compression','none');
      end
  end
end

if sw==3
  for t=1:Cells.time
      if t==1
          imwrite(uint16(max_crop_im(:,:,t)),[Class_folder,char(sp_name(4,1)),'\NA\',num2str(j),'.tif'],'Compression','none');
      else
   imwrite(uint16(max_crop_im(:,:,t)),[Class_folder,char(sp_name(4,1)),'\NA\',num2str(j),'.tif'],'Writemode','Append','Compression','none');
      end
  end
end
if sw==4
    close(figure(1))
    close(figure(2))
    break
end
NF_map.list(j,:)=[NF_map.Cell_position(j,1) NF_map.Cell_position(j,2) Cells.list(j,3)-floor(Cell_size(1,3)/2) sw];
end
NF_map.behavior=char(sp_name{4,1});
NF_map.zstack=Cells.zstack;
NF_map.time=Cells.time;
NF_map.filename=Cells.filename;

save([Class_folder,'NF_map_',sp_name{4,1},'.mat'],'NF_map');
