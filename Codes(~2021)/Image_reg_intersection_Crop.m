%%%  TXN Analysis Protocol : 1. Image Registration and make the
%%%  intersection map of three day and crop the images.

%
%     Website : http://neurobiophysics.snu.ac.kr/
%
% The analysis is proceed as following
% 1. Load the image of three day
% 2. Image registeration. Reference : HC
% 3. Crop the image and save.

% INPUTS 

% The folder with images (3 stack, Reg, Green, Subtracted).
% 

% OUTPUTS


% MODIFICATION HISTORY : 
%           Written by Byung Hun Lee, Deptartment of Physics and Astronomy, Seoul National University, 6/11/2018.

%% Parameter setting
clear
[fnm PathName1] = uigetfile('*.tif','Pick cropped image folder to analysis','Multiselect','on');  % Cropped image folder
% imds = imageDatastore(PathName1,'IncludeSubfolders',true,'LabelSource',...
%     'foldernames');       %^%      %#ok
zstack_inter=[26 56; 26 56;];
zstack=81;
trans_stack=25;
save_path='H:\Image_data\OlympusTPM\20181008-20181012\OE9\Transformed\';
mkdir(save_path);
% target_RGB_Image='H:\Image_data\OlympusTPM\20180607\20180607_Arc_oe1_CFC_crop\t01_z01_cell_num0019.tif';
% target_Image='H:\Image_data\OlympusTPM\20180607\20180607_Arc_oe1_CFC_crop_sub\t01_z01_cell_num0019.tif';
%%
ref=1;
 sp_name=split([PathName1 fnm{1,ref}],'\');
for i=zstack_inter(ref,1):zstack_inter(ref,2)
    im=imread([PathName1 fnm{1,ref}],i);
    if i==zstack_inter(ref,1)
    imwrite(im,[save_path 'T_' sp_name{6,1}],'tif','Compression', 'none');
    else
    imwrite(im,[save_path 'T_' sp_name{6,1}],'tif','Writemode','Append','Compression', 'none');
    end
end
for i=zstack_inter(ref,1):zstack_inter(ref,2)
    im=imread([PathName1 fnm{1,ref}],i+zstack);
    imwrite(im,[save_path 'T_' sp_name{6,1}],'tif','Writemode','Append','Compression', 'none');
end
%%
i=1;
while i<=size(fnm,2)
    if i~=ref
BaseFile = [PathName1 fnm{1,ref}];
FileToWarp = [PathName1 fnm{1,i}];
sp_name=split(FileToWarp,'\');
WarpImageBaseFile = sscanf(FileToWarp,'%c',size(FileToWarp,2)-4);
imginf=imfinfo(FileToWarp);
numstack=size(imginf,1);
timestack=numstack/(2*zstack);

ImageToWarp = imread(FileToWarp,'tif',zstack+zstack_inter(i,1)+trans_stack);
BaseImage = imread(BaseFile,'tif',zstack+zstack_inter(ref,1)+trans_stack);

tempWarp = double(ImageToWarp);
tempBase = double(BaseImage);

tempWarp = tempWarp./max(max(tempWarp))*3;
tempBase = tempBase./max(max(tempBase))*3;

% Select the control points.
h = cpselect(tempWarp, tempBase);
uiwait(msgbox('Click OK after closing the CPSELECT window.','Waiting...'))

movingPoints = cpcorr(movingPoints, fixedPoints,tempWarp,tempBase);
tform2 = maketform('projective',movingPoints, fixedPoints);
for c=1:2
for z=zstack_inter(i,1):zstack_inter(i,2)
    for t=1

    ImageToWarp = imread(FileToWarp,'tif',(z-1)*timestack+t+(c-1)*numstack/2);
[transVis xdata ydata] = imtransform(ImageToWarp, tform2);
transVisWithOffset = imtransform(ImageToWarp, tform2, 'XData', [1 (size(ImageToWarp,2)+tform2.tdata.T(3,1))],'YData', [1 (size(ImageToWarp,1)+ tform2.tdata.T(3,2))]);
transVisBuffrd = uint16(zeros(size(BaseImage)));
transVisBuffrd(1:size(transVisWithOffset,1),1:size(transVisWithOffset,2)) = transVisWithOffset;
Transformation_Matrix=tform2.tdata.T(1:2,1:2);
Translation_Vector=tform2.tdata.T(3,1:2)';
transVisBuffrd = transVisBuffrd((1:size(ImageToWarp,1)), (1:size(ImageToWarp,2)));

imwrite(transVisBuffrd,[save_path 'T_' char(sp_name{6,1})],'tif','Writemode','Append','Compression', 'none');
    end
end
end
save([WarpImageBaseFile,'T.mat'], 'tform2');
    end
    i=i+1;
end