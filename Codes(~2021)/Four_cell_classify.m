%  Image enhancement code
%
%     Website : http://neurobiophysics.snu.ac.kr/
%
% This code enhances the signal with specific PSF size.

% INPUTS 

% Image stack
% 

% OUTPUTS
% Filtered image stack

% MODIFICATION HISTORY : 
%     2019.01.23.
% Modification on filtered image (3D Gaussian -> 2D).
% Add close figure at the last part.
% Show the total cell number in the menu function.
%  2019.07.01 
%  Modification on zero column error.
%  2019.08.07
%  Modification on showing pop up window. 
%           Written by Byung Hun Lee, Deptartment of Physics and Astronomy, Seoul National University, 6/11/2018.
%% Classify the cell in the first and second frame
clear
[Filename PathName] = uigetfile('.tif','Select the subtracted registrated cropped image files to analysis');  % Cropped image folder
Class_folder='E:\BACKUP\대학원\연구실\MY_Projects\In_vivo_imaging\Cell_detection\TXN_detection_deep_learning\Classification_Blind\OE41\';
[fnm_cell pth_cell]=uigetfile('.mat','Select the cell coordinate');
Cell_list_file=fullfile(pth_cell,fnm_cell);
load(Cell_list_file)
Cell_size=[61 61 61];

%% 
    sp=split(Filename,'.');
    sp_name=split(sp{1,1},'_');
    mkdir([Class_folder sp_name{2,1}])
    mkdir([Class_folder sp_name{2,1} '\TXN'])
    mkdir([Class_folder sp_name{2,1} '\Non_TXN'])
    mkdir([Class_folder sp_name{2,1} '\NA'])

%%
dir_ts=dir([Class_folder sp_name{2,1} '\TXN']);
dir_no=dir([Class_folder sp_name{2,1} '\Non_TXN']);
dir_NA=dir([Class_folder sp_name{2,1} '\NA']);

dir_pre_cells=[dir_no(3:size(dir_no,1)); dir_ts(3:size(dir_ts,1)) ; dir_NA(3:size(dir_NA,1))];
fnd_result_file=dir(Class_folder);
if ~isempty(cell2mat(strfind({fnd_result_file.name},['NF_map_',sp_name{2,1},'.mat'])))
load([Class_folder,'NF_map_',sp_name{2,1},'.mat'])
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
    clear im_exp
    load(Cell_list_file)
    NF_map.Cell_position=Cells.list;

    for k=1:Cells.zstack
    im(:,:,k)=imread(fullfile(PathName,Filename),k);
    end
    im_exp=zeros(size(im,1)+2*floor(Cell_size(1,1)/2),size(im,2)+2*floor(Cell_size(1,2)/2),size(im,3)+2*floor(Cell_size(1,3)/2));
    im_exp(floor(Cell_size(1,2)/2)+1:floor(Cell_size(1,2)/2)+size(im,1),floor(Cell_size(1,1)/2)+1:floor(Cell_size(1,1)/2)+size(im,2),...
           floor(Cell_size(1,3)/2)+1:floor(Cell_size(1,3)/2)+size(im,3))=im;
  
    %%
     load(Cell_list_file)
    NF_map.Cell_position=Cells.list;
    n=12;
    for k=1:size(Cell_size,2)
        Cells.list(:,k)=Cells.list(:,k)+floor(Cell_size(1,k)/2);
    end
for j=start_here:size(Cells.list,1)    
    clear crop_im max_crop_im mont filt_mont;

    crop_im=im_exp(round(Cells.list(j,2))-floor(Cell_size(1,2)/2):round(Cells.list(j,2))+floor(Cell_size(1,2)/2),...
               round(Cells.list(j,1))-floor(Cell_size(1,1)/2):round(Cells.list(j,1))+floor(Cell_size(1,1)/2),...
               round(Cells.list(j,3))-floor(Cell_size(1,3)/2):round(Cells.list(j,3))+floor(Cell_size(1,3)/2));

max_crop_im=max(crop_im,[],3);

for jj=1:size(crop_im,3)
Bgimfilt=imgaussfilt(crop_im(:,:,jj), 6);
subim=abs(crop_im(:,:,jj)-Bgimfilt);
Filt_im(:,:,jj)=imgaussfilt(subim,1.5);
end

for jj=1:size(crop_im,3)
              if mod(jj,n)==0
mont(1+(floor(jj/n)-1)*Cell_size(1,1):(floor(jj/n))*Cell_size(1,1),1+(n-1)*Cell_size(1,2):(n)*Cell_size(1,2))=crop_im(:,:,jj);
filt_mont(1+(floor(jj/n)-1)*Cell_size(1,1):(floor(jj/n))*Cell_size(1,1),1+(n-1)*Cell_size(1,2):(n)*Cell_size(1,2))=Filt_im(:,:,jj);
           else
mont(1+floor(jj/n)*Cell_size(1,1):(floor(jj/n)+1)*Cell_size(1,1),1+(mod(jj,n)-1)*Cell_size(1,2):(mod(jj,n))*Cell_size(1,2))=crop_im(:,:,jj);
filt_mont(1+floor(jj/n)*Cell_size(1,1):(floor(jj/n)+1)*Cell_size(1,1),1+(mod(jj,n)-1)*Cell_size(1,2):(mod(jj,n))*Cell_size(1,2))=Filt_im(:,:,jj);
           end
end


handles.fig = figure(1);
set(gcf, 'Position',  [50, 30, 1500, 950])
handles.axes1 = axes('Units','pixels','Position',[50 520 1100 450]);   
imshow(mont,[min(min(mont)) 1.1*max(max(mont))])
title(['( ' num2str(j) '/' num2str(size(Cells.list,1)) ' )th cell.'])
% imshow(mont,[0 max(max(mont))])
colormap('gray')
axis equal tight off
hold on
handles.axes2 = axes('Units','pixels','Position',[1150 520 400 450]);
imagesc(max_crop_im)
axis equal tight off
colormap('gray')

handles.axes3 = axes('Units','pixels','Position',[50 50 1100 450]);
imshow(filt_mont,[min(min(filt_mont)) 1.3*max(max(filt_mont))])
colormap('gray')
axis equal tight off
title('Filtered Image')
handles.axes4 = axes('Units','pixels','Position',[1150 50 400 450]);
imagesc(max(Filt_im,[],3))
axis equal tight off
colormap('gray')

% figure(2)
% subplot(1,1,1)
% imagesc(im(:,:,round(Cells.list(j,3))-+floor(Cell_size(1,k)/2)))
% colormap('gray')
% axis equal
hold all
%plot(NF_map.Cell_position(j,1),NF_map.Cell_position(j,2),'ro','markersize',50)
%sw=menu(['( ' num2str(j) '/' num2str(size(Cells.list,1)) ' )th cell. Is the cell has TXN site?'],'Non TS','TS','N/A Cell','Break');
sw=input(['( ' num2str(j) '/' num2str(size(Cells.list,1)) ' )th cell. Is the cell has TXN site?\n']);
if sw==1
    for jj=1:size(crop_im,3)
        if jj==1
            imwrite(uint16(crop_im(:,:,jj)),[Class_folder,char(sp_name(2,1)),'\Non_TXN\',num2str(j),'.tif'],'Compression','none');
        else
   imwrite(uint16(crop_im(:,:,jj)),[Class_folder,char(sp_name(2,1)),'\Non_TXN\',num2str(j),'.tif'],'Writemode','Append','Compression','none');
        end
    end
end

if sw==2
   for jj=1:size(crop_im,3)
       if jj==1
           imwrite(uint16(crop_im(:,:,jj)),[Class_folder,char(sp_name(2,1)),'\TXN\',num2str(j),'.tif'],'Compression','none');
       else
   imwrite(uint16(crop_im(:,:,jj)),[Class_folder,char(sp_name(2,1)),'\TXN\',num2str(j),'.tif'],'Writemode','Append','Compression','none');
       end
   end
end

if sw==3
   for jj=1:size(crop_im,3)
       if jj==1
   imwrite(uint16(crop_im(:,:,jj)),[Class_folder,char(sp_name(2,1)),'\NA\',num2str(j),'.tif'],'Compression','none');        
       else
   imwrite(uint16(crop_im(:,:,jj)),[Class_folder,char(sp_name(2,1)),'\NA\',num2str(j),'.tif'],'Writemode','Append','Compression','none');
       end
   end
end
if sw==4
    close(figure(1))
    close(figure(2))
    break
end
NF_map.list(j,:)=[NF_map.Cell_position(j,1) NF_map.Cell_position(j,2) Cells.list(j,3)-floor(Cell_size(1,3)/2) sw];
if mod(j,50)==0
close(figure(1))
    close(figure(2))
end
end
  close(figure(1))
    close(figure(2))

if j==size(Cells.list,1)
pth=[Class_folder char(sp_name(2,1))];
imds = imageDatastore(pth,'IncludeSubfolders',true,'LabelSource',...
    'foldernames');
sw={'Non_TXN','TXN','NA'};

tmp=split(imds.Files,'.');
for i=1:size(tmp,1)
spp_name{i,1}=split(tmp{i,1},'\');
end
for i=1:size(spp_name{1,1},1)
    if ~isempty(findstr(char(spp_name{1,1}(i,1)),'NA')) || ~isempty(findstr(char(spp_name{1,1}(i,1)),'Non_TXN')) || ~isempty(findstr(char(spp_name{1,1}(i,1)),'TXN'))
        ref=i;
    end
end
clear tmp

NF=zeros(size(imds.Files,1),1);
for i=1:size(imds.Files,1)
    for j=1:size(sw,2)
       if ~isempty(strfind(sw{1,j},char(imds.Labels(i,1)))) && ~isempty(strfind(char(imds.Labels(i,1)),sw{1,j}))
           g=str2num(cell2mat(spp_name{i,1}(ref+1,1)));
           if NF(g,1)==0
    NF(g,1)=j;
           else
               error(['There are two' num2str(g) 'cell in different folder !!!!!!'])
           end
           end
    end
    
end
NF_map.list(:,4)=NF;
NF_map.list(:,1:3)=NF_map.Cell_position;
end

NF_map.behavior=char(sp_name{2,1});
NF_map.zstack=Cells.zstack;
NF_map.filename=Cells.filename;

save([Class_folder,'NF_map_',char(sp_name{2,1}),'.mat'],'NF_map');
