%%%  TXN analysis protocol : Subtraction & Max projection & Cell Segmentation 
%
%     Website : http://neurobiophysics.snu.ac.kr/
%
% This code shows cell image and let user classify the TS cell and no TS
% cell.
% The analysis is proceed as following
% 1. Load the registrated images
% 2. Do subtraction in each z stack
% 3. Max projection
% 4. Cell Segment with Circle finder.

% INPUTS 

% The folder with images (3 stack, Reg, Green, Subtracted).
% 

% OUTPUTS


% MODIFICATION HISTORY : 
%           Written by Byung Hun Lee, Deptartment of Physics and Astronomy, Seoul National University, 6/11/2018.

%% Parameter setting
clear
[fnm,pth]=uigetfile('*.tif','Select the registerated HC, Train and retrieval images','Multiselect','on');
channel=2;
timeline={'HC','Training','Retrieval'};
% target_RGB_Image='H:\Image_data\OlympusTPM\20180607\20180607_Arc_oe1_CFC_crop\t01_z01_cell_num0019.tif';
% target_Image='H:\Image_data\OlympusTPM\20180607\20180607_Arc_oe1_CFC_crop_sub\t01_z01_cell_num0019.tif';
%% Subtraction & Max Projection

for i=1:length(fnm)
    
    
imginf=imfinfo([pth char(fnm{1,i})]);
numstack=size(imginf,1);    
zstack=numstack/channel;
for j=1:zstack
    im(:,:,1)=imread([pth char(fnm{1,i})],j);
    imG=imread([pth char(fnm{1,i})],j+zstack);
    
    pseudo=imG>2500;
    imG(pseudo)=median(reshape(imG,size(imG,1)*size(imG,2),1));
    im(:,:,2)=imG;
    
subtraction(:,:,j)=subtraction_image_ftn(im);
end
max_sub(:,:,i)=max(subtraction,[],3);
if i==1
    imwrite(max_sub(:,:,i),[pth 'max_sub.tif'],'Compression','none');
else
imwrite(max_sub(:,:,i),[pth 'max_sub.tif'],'Writemode','Append','Compression','none');
end
end
%% Cell Segmentation in each image and match
for i=1:3
    %[centers, radii]=Cell_segment_circle(max_sub(:,:,i));
    cir
    %centers=Cell_segment2(max_sub(:,:,i),200,3,0.5);
mod_Cell_cent=mod_cell_cent(centers,15);
    Cell_cent{1,i}=mod_Cell_cent;
subplot(2,2,i)
imagesc(max_sub(:,:,i))
colormap('gray')
hold all
plot(Cell_cent{1,i}(:,1),Cell_cent{1,i}(:,2),'ro','linewidth',1)
title([char(timeline{1,i}),'(n=',num2str(size(Cell_cent{1,i},1)),')'])
end