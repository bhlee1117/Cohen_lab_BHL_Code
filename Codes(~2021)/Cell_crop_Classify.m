%  Manually classfication code 
%
%     Website : http://neurobiophysics.snu.ac.kr/
%
% This code shows cell image and let user classify the TS cell and no TS
% cell.
% The analysis is proceed as following
% 1. Load the cropped cell image.
% 2. Show cell image to the user.
% 3. User classify the cell.
% 4. Save the subtracted image into different folder.

% INPUTS 

% The folder with images (3 stack, Reg, Green, Subtracted).
% 

% OUTPUTS


% MODIFICATION HISTORY : 
%           Written by Byung Hun Lee, Deptartment of Physics and Astronomy, Seoul National University, 6/11/2018.

%% Parameter setting
clear
[Filename PathName] = uigetfile('.tif','Pick cropped image files to analysis');  % Cropped image folder
Class_folder='E:\BACKUP\대학원\연구실\MY_Projects\In_vivo_imaging\Cell_detection\TXN_detection_deep_learning\Classification\OE9\';
Cell_list_file='H:\Image_data\OlympusTPM\20181008-20181012\OE9\Transformed\Cell_list_T_20181008_OE9_HC_.tif.mat';
load(Cell_list_file)
Cell_size=[60 60 49];
% target_RGB_Image='H:\Image_data\OlympusTPM\20180607\20180607_Arc_oe1_CFC_crop\t01_z01_cell_num0019.tif';
% target_Image='H:\Image_data\OlympusTPM\20180607\20180607_Arc_oe1_CFC_crop_sub\t01_z01_cell_num0019.tif';
%%
    sp_name=split(Filename,'_');
    mkdir([Class_folder sp_name{5,1}])
    mkdir([Class_folder sp_name{5,1} '\TXN'])
    mkdir([Class_folder sp_name{5,1} '\Non_TXN'])
    mkdir([Class_folder sp_name{5,1} '\NA'])

%%
dir_ts=dir([Class_folder sp_name{5,1} '\TXN']);
dir_no=dir([Class_folder sp_name{5,1} '\Non_TXN']);
dir_NA=dir([Class_folder sp_name{5,1} '\NA']);

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
    for j=1:Cells.zstack
    im(:,:,j)=imread(fullfile(PathName,Filename),j);
    end
    im_exp=zeros(size(im,1)+2*floor(Cell_size(1,1)/2),size(im,2)+2*floor(Cell_size(1,2)/2),size(im,3)+2*floor(Cell_size(1,3)/2));
    im_exp(floor(Cell_size(1,2)/2)+1:floor(Cell_size(1,2)/2)+size(im,1),floor(Cell_size(1,1)/2)+1:floor(Cell_size(1,1)/2)+size(im,2),floor(Cell_size(1,3)/2)+1:floor(Cell_size(1,3)/2)+size(im,3))=im;
    for k=1:size(Cell_size,2)
        Cells.list(:,k)=Cells.list(:,k)+floor(Cell_size(1,k)/2);
    end
for j=start_here:size(Cells.list,1)    
    clear crop_im max_crop_im;

    crop_im=im_exp(round(Cells.list(j,2))-floor(Cell_size(1,2)/2):round(Cells.list(j,2))+floor(Cell_size(1,2)/2),...
               round(Cells.list(j,1))-floor(Cell_size(1,1)/2):round(Cells.list(j,1))+floor(Cell_size(1,1)/2),...
               round(Cells.list(j,3))-floor(Cell_size(1,3)/2):round(Cells.list(j,3))+floor(Cell_size(1,3)/2));

           max_crop_im=max(crop_im,[],3);
Bgimfilt=imgaussfilt(max_crop_im, 6);
subim=max_crop_im-Bgimfilt;
Filt_im=imgaussfilt(subim,1.5);
           
figure(1)
subplot(1,2,1)
imagesc(max_crop_im)
colormap('gray')
axis equal
hold on
title('Subtracted Image')

subplot(1,2,2)
imagesc(Filt_im)
colormap('gray')
axis equal
title('Filtered Image')
figure(2)
subplot(1,1,1)
imagesc(im(:,:,round(Cells.list(j,3))-+floor(Cell_size(1,k)/2)))
colormap('gray')
axis equal
hold all
plot(NF_map.Cell_position(j,1),NF_map.Cell_position(j,2),'ro','markersize',50)
sw=menu('Is the cell has TXN site?','Non TS','TS','N/A Cell','Break');
if sw==1
   imwrite(uint16(max_crop_im),[Class_folder,char(sp_name(5,1)),'\Non_TXN\',num2str(j),'.tif'],'Compression','none');
end

if sw==2
   imwrite(uint16(max_crop_im),[Class_folder,char(sp_name(5,1)),'\TXN\',num2str(j),'.tif'],'Compression','none');
end

if sw==3
imwrite(uint16(max_crop_im),[Class_folder,char(sp_name(5,1)),'\NA\',num2str(j),'.tif'],'Compression','none');
end
if sw==4
    close(figure(1))
    break
end
NF_map.list(j,:)=[NF_map.Cell_position(j,1) NF_map.Cell_position(j,2) Cells.list(j,3)-floor(Cell_size(1,3)/2) sw];
end
NF_map.behavior=char(sp_name{5,1});
NF_map.zstack=Cells.zstack;
NF_map.filename=Cells.filename;

save([Class_folder,'NF_map_',sp_name{5,1},'.mat'],'NF_map');
