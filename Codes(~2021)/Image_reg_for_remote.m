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
[fnmH PathName1] = uigetfile('.tif','Pick Homecage image');  % Cropped image folder
[fnmR PathNameR] = uigetfile('.tif','Pick Remote memory image');  % Cropped image folder

zstack_inter=[1 61; 21 81];
zstack=81;
trans_stack=10;
save_path='H:\Image_data\OlympusTPM\20180626_28_In_vivo\OE5\Transformed\';

% target_RGB_Image='H:\Image_data\OlympusTPM\20180607\20180607_Arc_oe1_CFC_crop\t01_z01_cell_num0019.tif';
% target_Image='H:\Image_data\OlympusTPM\20180607\20180607_Arc_oe1_CFC_crop_sub\t01_z01_cell_num0019.tif';

%%
i=2;
ref=1;

BaseFile = fullfile(PathName1,fnmH);
FileToWarp = fullfile(PathNameR,fnmR);
sp_name=split(FileToWarp,'\');
WarpImageBaseFile = sscanf(FileToWarp,'%c',size(FileToWarp,2)-4);
imginf=imfinfo(FileToWarp);
numstack=size(imginf,1);
timestack=numstack/(2*zstack);

ImageToWarp = imread(FileToWarp,'tif',zstack+zstack_inter(i,1)+trans_stack);
BaseImage = imread(BaseFile,'tif',zstack+zstack_inter(ref,1)+trans_stack);

tempWarp = double(ImageToWarp);
tempBase = double(BaseImage);

tempWarp = tempWarp./max(max(tempWarp));
tempBase = tempBase./max(max(tempBase));

% Select the control points.
h = cpselect(tempWarp, tempBase);
uiwait(msgbox('Click OK after closing the CPSELECT window.','Waiting...'))

movingPoints = cpcorr(movingPoints, fixedPoints,tempWarp,tempBase);
tform2 = maketform('projective',movingPoints, fixedPoints);
%%
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

imwrite(transVisBuffrd,[save_path 'T_' char(sp_name{5,1})],'tif','Writemode','Append','Compression', 'none');
    end
end
end
save([WarpImageBaseFile,'T.mat'], 'tform2');

