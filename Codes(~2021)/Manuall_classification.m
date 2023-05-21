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
[PathName1] = uigetdir('Pick cropped image folder to analysis');  % Cropped image folder
[PathName2] = uigetdir('Pick cropped_sub image folder to analysis');
file_name=dir(PathName1);
Class_folder='E:\BACKUP\대학원\연구실\MY_Projects\In_vivo_imaging\Cell_detection\TXN_detection_deep_learning\Subtraction_0629\';
Class_RGBfolder='E:\BACKUP\대학원\연구실\MY_Projects\In_vivo_imaging\Cell_detection\TXN_detection_deep_learning\RGB_0629\';

Result_folder='E:\BACKUP\대학원\연구실\MY_Projects\In_vivo_imaging\Cell_detection\TXN_detection_deep_learning\Result\';
Cell_cent_folder='H:\Image_data\OlympusTPM\20180626_28_In_vivo\';

% target_RGB_Image='H:\Image_data\OlympusTPM\20180607\20180607_Arc_oe1_CFC_crop\t01_z01_cell_num0019.tif';
% target_Image='H:\Image_data\OlympusTPM\20180607\20180607_Arc_oe1_CFC_crop_sub\t01_z01_cell_num0019.tif';
%%
sp_pt_name=split(PathName1,'\');
dir_no=dir([Class_folder,'Non_TXN']);
dir_ts=dir([Class_folder,'TXN']);

dir_pre_cells=[dir_no(3:size(dir_no,1)); dir_ts(3:size(dir_ts,1))];%; dir_NA(3:size(dir_NA,1))];
sp_pth_name=split(PathName1,'\');
folder_name=sp_pth_name(end,1);
img_info=split(folder_name,'_');
img_name=[];

for i=1:size(img_info,1)-1
img_name=[img_name img_info{i,1} '_'];
end
load([Cell_cent_folder,img_name])
fnd_result_file=dir(Result_folder);
if ~isempty(cell2mat(strfind({fnd_result_file.name},['Class_',img_name])))
load([Result_folder,'Class_',img_name])
else
    
end
if ~isempty(cell2mat(strfind({dir_pre_cells.name},img_info{1,1})))  % Find the date in the categorized folder
    i=3;
    while ~isempty(cell2mat(strfind({dir_pre_cells.name},[img_name,'@',file_name(i).name]))) && i<=size(file_name,1)
        i=i+1;
        start_here=i;
    end
else
    start_here=3;
end
%%
engram_class.filename=img_name;
if strfind(img_info{4,1},'HC')
engram_class.train=0;
else if strfind(img_info{4,1},'CFC') engram_class.train=1;
else if strfind(img_info{4,1},'Ret') engram_class.train=2;
end
end
end
clear tar_im
% tar_im(:,:,1)=imread(target_RGB_Image,1);
% tar_im(:,:,2)=imread(target_RGB_Image,2);
% tar_im(:,:,3)=imread(target_RGB_Image,3);

%%
for i=start_here:length(file_name)-2
    fn_sp=split(file_name(i).name,'.');
    clear I im
I(:,:,1)=double(imread(fullfile(PathName1,file_name(i).name),1));
I(:,:,2)=double(imread(fullfile(PathName1,file_name(i).name),2));
I(:,:,3)=zeros;

im=uint8(I);
% im(:,:,1)=imhistmatch(uint8(I(:,:,1)),tar_im(:,:,1));
% im(:,:,2)=imhistmatch(uint8(I(:,:,2)),tar_im(:,:,2));
% im(:,:,3)=zeros;

% tar_im_sub=imread(target_Image,1);

sub=imread(fullfile(PathName2,file_name(i).name));
%sub=imhistmatch(sub,tar_im_sub);
% Information of cropped image
crop_info=split(file_name(i).name,'_');
tzc=[str2num(crop_info{1,1}(2:3)),str2num(crop_info{2,1}(2:3)),str2num(crop_info{4,1}(4:7))];
figure(1)
subplot(1,2,1)
imagesc(im)
axis equal
hold on
subplot(1,2,2)
imagesc(sub)
colormap('gray')
axis equal
sw=menu('Is the cell has TXN site?','Non TS','TS','N/A Cell','Break');
if sw==1
    
   imwrite(im(:,:,1),[Class_RGBfolder,'Non_TXN\',img_name,'@',file_name(i).name],'Compression','none');
   imwrite(im(:,:,2),[Class_RGBfolder,'Non_TXN\',img_name,'@',file_name(i).name],'Writemode','Append','Compression','none');
   imwrite(im(:,:,3),[Class_RGBfolder,'Non_TXN\',img_name,'@',file_name(i).name],'Writemode','Append','Compression','none');
   
   imwrite(sub,[Class_folder,'Non_TXN\',img_name,'@',file_name(i).name],'Compression','none');
end

if sw==2
   imwrite(im(:,:,1),[Class_RGBfolder,'TXN\',img_name,'@',file_name(i).name],'Compression','none');
   imwrite(im(:,:,2),[Class_RGBfolder,'TXN\',img_name,'@',file_name(i).name],'Writemode','Append','Compression','none');
   imwrite(im(:,:,3),[Class_RGBfolder,'TXN\',img_name,'@',file_name(i).name],'Writemode','Append','Compression','none');
   
   imwrite(sub,[Class_folder,'TXN\',img_name,'@',file_name(i).name],'Compression','none');
    %imwrite(uint16(sub),[Class_folder,'TXN\',img_name,'@',file_name(i).name],'Compression','none');
end

% if sw==3
%    imwrite(uint16(im(:,:,1)),[Class_folder,'NA\',img_name,'@',file_name(i).name],'Compression','none');
%    imwrite(uint16(im(:,:,2)),[Class_folder,'NA\',img_name,'@',file_name(i).name],'Writemode','Append','Compression','none');
%    imwrite(uint16(sub),[Class_folder,'NA\',img_name,'@',file_name(i).name],'Writemode','Append','Compression','none');
%     %imwrite(uint16(sub),[Class_folder,'NA\',img_name,'@',file_name(i).name],'Compression','none');
% end
if sw==4
    close(figure(1))
    break
end
engram_class.result(i-2,1:6)=[tzc(1,1) Cell_cent{tzc(1,1),tzc(1,2)}(tzc(1,3),1) Cell_cent{tzc(1,1),tzc(1,2)}(tzc(1,3),2) tzc(1,2) sw tzc(1,3)];

end
save([Result_folder,'Class_',img_name],'engram_class');
%% 2. Cell segmentation of each z-stack & Cell centroid detection % Crop the Cell image and save

for i=2:length(FileName1)
    fo_na=split(char(FileName1{i}),'.');
    mkdir([PathName1,char(fo_na{1,1}),'_crop'])
    for j=1:stack{i}.zstack
        for t=1:stack{i}.time
im=imread([PathName1,char(FileName1{i})],stack{i}.time*stack{i}.zstack+stack{i}.time*(j-1)+t);
im_R=imread([PathName1,char(FileName1{i})],stack{i}.time*(j-1)+t);
Cell_cent{t,j}=Cell_segment(im,250,3,0.7);
% figure(1)
% imagesc(im(:,:,j))
% colormap('gray')
% hold all
% for k=1:size(Cell_cent{j},1)
% plot(Cell_cent{j}(k,1),Cell_cent{j}(k,2),'ro')
%      text(Cell_cent{j}(k,1),Cell_cent{j}(k,2)+10,num2str(k),'color','r',...
%  'fontweight','bold')
% end
% close(figure(1))
for k=1:size(Cell_cent{j},1)
    if round(Cell_cent{j}(k,2))>28 && round(Cell_cent{j}(k,2))+27<=size(im,1) && round(Cell_cent{j}(k,1))>26 && round(Cell_cent{j}(k,1))+25<=size(im,2) 
    im_crop(:,:,2)=im(round(Cell_cent{j}(k,2))-28:round(Cell_cent{j}(k,2))+27,round(Cell_cent{j}(k,1))-26:round(Cell_cent{j}(k,1))+25);
    im_crop(:,:,1)=im_R(round(Cell_cent{j}(k,2))-28:round(Cell_cent{j}(k,2))+27,round(Cell_cent{j}(k,1))-26:round(Cell_cent{j}(k,1))+25);
    
    subtraction=subtraction_image_ftn(im_crop); 
    imwrite(im_crop(:,:,1),[PathName1,char(fo_na{1,1}),'_crop\','t',num2str(t,'%02d'),'_z',num2str(j,'%02d'),'_cell_num',num2str(k,'%02d'),'.tif'],'Compression','none');
    imwrite(im_crop(:,:,2),[PathName1,char(fo_na{1,1}),'_crop\','t',num2str(t,'%02d'),'_z',num2str(j,'%02d'),'_cell_num',num2str(k,'%02d'),'.tif'],'WriteMode','append','Compression','none');
    imwrite(subtraction,[PathName1,char(fo_na{1,1}),'_crop\','t',num2str(t,'%02d'),'_z',num2str(j,'%02d'),'_cell_num',num2str(k,'%02d'),'.tif'],'WriteMode','append','Compression','none');
    end
end
        end
    end
    save([PathName1,char(fo_na(1,1)),'.mat'],'Cell_cent')
end

%% 6. Classify the cells manually : TS cell, No TS cell

%% 7. Train the Machine until the accuracy exceed 90%.

