%  Classification code for RNA scope image.
%
%     Website : http://neurobiophysics.snu.ac.kr/
%
% INPUTS 

% RNA scope image

% OUTPUTS
% 1. Cropped cell image
% 2. TXN list
%
%           Written by Byung Hun Lee, Deptartment of Physics and Astronomy, Seoul National University, 03/04/2018.
%% Classify the cell in the first and second frame
clear
[fnm pth] = uigetfile('.tiff','Select the subtracted registrated cropped image files to analysis','Multiselect','on'); % Cropped image folder

save_dir='H:\Image_data\OlympusTPM\RNA_scope\Cropped\';

Class_folder='H:\Image_data\OlympusTPM\RNA_scope\Classification\';
Cell_size=[51 51];
expander=201;
%% 

clear start_here
    sp=split(fnm,'.');
    list_file=[save_dir 'list' sp{1,1} '.mat'];
    load(list_file)
    
    mkdir([Class_folder sp{1,1}])
    mkdir([Class_folder sp{1,1} '\TXN'])
    mkdir([Class_folder sp{1,1} '\Non_TXN'])
    mkdir([Class_folder sp{1,1} '\NA'])

dir_ts=dir([Class_folder sp{1,1} '\TXN']);
dir_no=dir([Class_folder sp{1,1} '\Non_TXN']);
dir_NA=dir([Class_folder sp{1,1} '\NA']);

dir_pre_cells=[dir_no(3:size(dir_no,1)); dir_ts(3:size(dir_ts,1)) ; dir_NA(3:size(dir_NA,1))];
fnd_result_file=dir(Class_folder);
if ~isempty(cell2mat(strfind({fnd_result_file.name},['NF_map_',sp{1,1},'.mat'])))
load([Class_folder,'NF_map_',sp{1,1},'.mat'])
else
    NF_map.list=[];
end
if ~size(dir_pre_cells,1)==0  % Find the date in the categorized folder
    g=1;
    while ~isempty(cell2mat(strfind({dir_pre_cells.name},[num2str(g) ,'.tif']))) && g<=size(Cell_list,1)
        g=g+1;
        start_here=g;
    end
else
    start_here=1;
end

%%
    clear im im_exp
    NF_map.Cell_position=Cell_list;

    for k=1
    im{k}(:,:,1)=imread(fullfile(pth,fnm),2);
    im{k}(:,:,2)=imread(fullfile(pth,fnm),4);
    im{k}(:,:,3)=imread(fullfile(pth,fnm),5);
    im_exp{k}=zeros(size(im{k},1)+2*floor(Cell_size(1,1)/2),size(im{k},2)+2*floor(Cell_size(1,2)/2),3);
    im_exp{k}(floor(Cell_size(1,2)/2)+1:floor(Cell_size(1,2)/2)+size(im{k},1),floor(Cell_size(1,1)/2)+1:floor(Cell_size(1,1)/2)+size(im{k},2),1:3)=im{k};
    
    im_tre_exp{k}=zeros(size(im{k},1)+2*floor(expander/2),size(im{k},2)+2*floor(expander/2),3);
    im_tre_exp{k}(floor(expander/2)+1:floor(expander/2)+size(im{k},1),floor(expander/2)+1:floor(expander/2)+size(im{k},2),1:3)=im{k};
    
    end
  
  
    %%

 sw=1;
    for i=1
        load(list_file)
        if sw==4
            break
        end
          sp=split(fnm,'.');
          clear Cells.list
           for k=1:size(Cell_size,2)
        Cells.list(:,k)=Cell_list(:,k)+floor(Cell_size(1,k)/2);
        Cells.list_exp(:,k)=Cell_list(:,k)+floor(expander/2);
           end
    
for j=start_here:size(Cell_list,1)    
    clear crop_im max_crop_im mont filt_mont RGB_im;

    crop_im=im_exp{i}(round(Cells.list(j,2))-floor(Cell_size(1,2)/2):round(Cells.list(j,2))+floor(Cell_size(1,2)/2),...
               round(Cells.list(j,1))-floor(Cell_size(1,1)/2):round(Cells.list(j,1))+floor(Cell_size(1,1)/2),1:3);
           
    crop_exp=im_tre_exp{i}(round(Cells.list_exp(j,2))-floor(expander/2):round(Cells.list_exp(j,2))+floor(expander/2),...
               round(Cells.list_exp(j,1))-floor(expander/2):round(Cells.list_exp(j,1))+floor(expander/2),1:3);
    RGB_im(:,:,3)=crop_im(:,:,1)/max(max(im_exp{1}(:,:,1)));
    RGB_im(:,:,1)=crop_im(:,:,2)/max(max(im_exp{1}(:,:,2)))*1.1;
    RGB_im(:,:,2)=crop_im(:,:,3)/max(max(im_exp{1}(:,:,3)))*1.1;
    
    RGB_im_exp(:,:,3)=crop_exp(:,:,1)/max(max(im_exp{1}(:,:,1)))*4.5;
    RGB_im_exp(:,:,1)=crop_exp(:,:,2)/max(max(im_exp{1}(:,:,2)))*1.1;
    RGB_im_exp(:,:,2)=crop_exp(:,:,3)/max(max(im_exp{1}(:,:,3)))*1.1;
    
RGB_im=flipud(RGB_im);
    figure(1)

subplot(2,2,1)

imagesc(crop_im(:,:,1))
colormap('gray')
axis equal
subplot(2,2,2)
imagesc(crop_im(:,:,2))
colormap('gray')
axis equal
subplot(2,2,3)
imagesc(crop_im(:,:,3))
colormap('gray')
axis equal
subplot(2,2,4)
hold all
imagesc(RGB_im)
axis equal
nTXN=0;

figure(2)
imagesc(RGB_im_exp)
hold all
axis equal
plot(expander/2,expander/2,'go','markersize',70)
sw=input(['( ' num2str(j) '/' num2str(size(Cells.list,1)) ' )th cell. Is the cell has TXN site?\n']);
if sw==1
    for jj=1:size(crop_im,3)
        if jj==1
            imwrite(uint16(crop_im(:,:,jj)),[Class_folder,char(sp(1,1)),'\Non_TXN\',num2str(j),'.tif'],'Compression','none');
        else
   imwrite(uint16(crop_im(:,:,jj)),[Class_folder,char(sp(1,1)),'\Non_TXN\',num2str(j),'.tif'],'Writemode','Append','Compression','none');
        end
    end
end

if sw==2
   for jj=1:size(crop_im,3)
       if jj==1
           imwrite(uint16(crop_im(:,:,jj)),[Class_folder,char(sp(1,1)),'\TXN\',num2str(j),'.tif'],'Compression','none');
       else
   imwrite(uint16(crop_im(:,:,jj)),[Class_folder,char(sp(1,1)),'\TXN\',num2str(j),'.tif'],'Writemode','Append','Compression','none');
       end
   end
   nTXN=input(['How many TXN site?\n']);
end

if sw==3
   for jj=1:size(crop_im,3)
       if jj==1
   imwrite(uint16(crop_im(:,:,jj)),[Class_folder,char(sp(1,1)),'\NA\',num2str(j),'.tif'],'Compression','none');        
       else
   imwrite(uint16(crop_im(:,:,jj)),[Class_folder,char(sp(1,1)),'\NA\',num2str(j),'.tif'],'Writemode','Append','Compression','none');
       end
   end
end
if sw==4
    close(figure(1))
    close(figure(2))
    break
end
NF_map.list(j,:)=[Cell_list(j,1) Cell_list(j,2) sw nTXN];
end
  close(figure(1))
    close(figure(2))
save([Class_folder,'NF_map_',sp{1,1},'.mat'],'NF_map');

    end


