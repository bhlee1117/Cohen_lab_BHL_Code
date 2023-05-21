%  Analysis protocol for In vivo Transcription data
%
%     Website : http://neurobiophysics.snu.ac.kr/
%
% This code is made on the purpose of detecting the engram cells which have
% Arc transcription site. 
% The analysis is proceed as following
% 1. Load the Raw data (2 channel, zstack, timelapse). The data should be
%    sequenced as ( C1Z1T1 -> C1Z1T2 -> C1Z1T3 -> C1Z2T1 -> C1Z2T2 ...)
% 2. Subtract the auto-fluorescent particles. G-R*max(G)/max(R)
% 3. Cell segmentation of each z-stack
% 4. Cell centroid detection
% 5. Crop each cells, 52 X 56 pixel
% 6. Classify the cells manually : TS cell, No TS cell
% 7. Train the Machine until the accuracy exceed 90%.

% INPUTS 

% Two channel, zstack and timelapse image from Two-photon microscope.
% 

% OUTPUTS
% Cell data containing file name, z-stack, x position, y position and
% existence of TS.

% MODIFICATION HISTORY : 
%           Written by Byung Hun Lee, Deptartment of Physics and Astronomy, Seoul National University, 6/11/2018.

%% Parameter setting
clear
z=81;
c=2;
% t= various;
cropxy=[41 41];
class_folder='E:\BACKUP\대학원\연구실\MY_Projects\In vivo imaging\Cell_detection\TXN_detection_deep_learning\Subtraction\'; % manually classified foler


%% 1. Load the Raw data (2 channel, zstack, timelapse). The data should be
%    sequenced as ( C1Z1T1 -> C1Z1T2 -> C1Z1T3 -> C1Z2T1 -> C1Z2T2 ...)

[FileName1, PathName1] = uigetfile('*.tif', 'Pick image files to analysis','Multiselect','on');
if ischar(FileName1)    tmp=FileName1; clear FileName1
    FileName1{1}=tmp; end

for i=1:length(FileName1)
imginfo=imfinfo([PathName1 char(FileName1{i})]);
frmmax=numel(imginfo);
stack{i}.Name=[PathName1 char(FileName1{i})];
stack{i}.zstack=z;
stack{i}.channel=c;
stack{i}.time=frmmax/(z*c);
end

%% 2. Cell segmentation of each z-stack & Cell centroid detection % Crop the Cell image and save

for i=1:length(FileName1)
    fo_na=split(char(FileName1{i}),'.');
    mkdir([PathName1,char(fo_na{1,1}),'crop'])
    mkdir([PathName1,char(fo_na{1,1}),'crop_sub'])
    for j=1:3%stack{i}.zstack
        for t=1%:stack{i}.time
            clear mod_Cell_cent
im=imread([PathName1,char(FileName1{i})],stack{i}.time*stack{i}.zstack+stack{i}.time*(j-1)+t);
im_R=imread([PathName1,char(FileName1{i})],stack{i}.time*(j-1)+t);
%Cell_cent{t,j}=Cell_segment2(im,200,3,0.6);  %Cell_cent=Cell_segment2(I,Cell_size,kernel,vac_ratio)
[centers, radii]=Cell_segment_circle(im);
Cell_cent{t,j}=centers;
mod_Cell_cent=mod_cell_cent(Cell_cent{t,j},10);
Cell_cent{t,j}=mod_Cell_cent;
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
for k=1:size(mod_Cell_cent,1)
    clear im_crop
    if round(mod_Cell_cent(k,2))>floor(cropxy(1,2)/2) && round(mod_Cell_cent(k,2))+floor(cropxy(1,2)/2)<=size(im,1) && round(mod_Cell_cent(k,1))>floor(cropxy(1,1)/2) && round(mod_Cell_cent(k,1))+floor(cropxy(1,1)/2)<=size(im,2) 
    im_crop(:,:,2)=im(round(mod_Cell_cent(k,2))-floor(cropxy(1,2)/2):round(mod_Cell_cent(k,2))+floor(cropxy(1,2)/2),round(mod_Cell_cent(k,1))-floor(cropxy(1,1)/2):round(mod_Cell_cent(k,1))+floor(cropxy(1,1)/2));
    im_crop(:,:,1)=im_R(round(mod_Cell_cent(k,2))-floor(cropxy(1,2)/2):round(mod_Cell_cent(k,2))+floor(cropxy(1,2)/2),round(mod_Cell_cent(k,1))-floor(cropxy(1,1)/2):round(mod_Cell_cent(k,1))+floor(cropxy(1,1)/2));
    im_crop(:,:,3)=zeros;
    im_crop=uint8(double(im_crop)/4095*255);
    
    subtraction=subtraction_image_ftn(im_crop); 
    
    
    imwrite(im_crop(:,:,1),[PathName1,char(fo_na{1,1}),'crop\','t',num2str(t,'%02d'),'_z',num2str(j,'%02d'),'_cell_num',num2str(k,'%04d'),'.tif'],'Compression','none');
    imwrite(im_crop(:,:,2),[PathName1,char(fo_na{1,1}),'crop\','t',num2str(t,'%02d'),'_z',num2str(j,'%02d'),'_cell_num',num2str(k,'%04d'),'.tif'],'WriteMode','append','Compression','none');
    imwrite(im_crop(:,:,3),[PathName1,char(fo_na{1,1}),'crop\','t',num2str(t,'%02d'),'_z',num2str(j,'%02d'),'_cell_num',num2str(k,'%04d'),'.tif'],'WriteMode','append','Compression','none');
    
    imwrite(subtraction,[PathName1,char(fo_na{1,1}),'crop_sub\','t',num2str(t,'%02d'),'_z',num2str(j,'%02d'),'_cell_num',num2str(k,'%04d'),'.tif'],'Compression','none');
    end
end
        end
    end
    save([PathName1,char(fo_na(1,1)),'.mat'],'Cell_cent')
end

%% 6. Classify the cells manually : TS cell, No TS cell

%% 7. Train the Machine until the accuracy exceed 90%.

