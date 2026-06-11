%% Load file path
clear
clc;
[~, ~, raw] = xlsread(['/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/Prism_V2+Glu_Data_Arrangement.xlsx'], 'Sheet1', 'C5:AA31');
load('/Volumes/BHL18TB_D2/20260203_SD_V2+iGluSNFR4/25X_transformationMatrix.mat');
fpath=raw(:,1)';
V2moviemaxTime=15000;
GlumoviemaxTime=5000;
Endframe=cell2mat(raw(:,5));
foi=[3 5 6 8 9 10 12 13];
StructureData=raw(:,6);
%% Map Glu Footprint onto morphology
% Voltage image >> Structure image >> Glu

f=13;
GluResult=importdata(fullfile(fpath{f},"Glu_Result.mat"));
VoltResult=importdata(fullfile(fpath{f},"Volt_Result.mat"));
load(fullfile(fpath{f},"output_data.mat"));
swcfiles=dir(fullfile(fpath{f}, 'Tracing*.swc'));
sz=size(GluResult.AvgGluImg);
if isnan(StructureData{f})
    StructureStack=VoltResult.ref_im;
else
    StructureStack=read_tiff(StructureData{f});
end

% Registration
figure(5); clf;
imshow2(VoltResult.ref_im,[]);
figure(6); clf;
moviefixsc(StructureStack);
Z_in=input('Putative stack: \n');
[RegStrImg, tformStrReg] = interactiveImageRegistration(StructureStack(:,:,Z_in), VoltResult.ref_im);
RegStrImg = imwarp(StructureStack(:,:,Z_in), tformStrReg, 'OutputView', imref2d(size(VoltResult.ref_im)));
VoltResult.tform_Str2Volt=tformStrReg;

alignedVolt = transformCamera_O2B(Device_Data, tformReg, VoltResult.ref_im, GluResult.AvgGluImg);
alignedStructure = transformCamera_O2B(Device_Data, tformReg, RegStrImg, GluResult.AvgGluImg);

Vmeta = Device_Data{1,3};  Gmeta = Device_Data{1,4};
Rglu = refFromROI(size(GluResult.AvgGluImg), double(Gmeta.ROI));
xV = double(Vmeta.ROI(1));  yV = double(Vmeta.ROI(3));

% Voltron cropped-pixel -> Voltron world
tformVoltPix2World = affine2d([1 0 0; 0 1 0; xV yV 1]);
tformVolt2Glu = invert(tformReg);

% Compose Structure -> Glu
T = tformStrReg.T * tformVoltPix2World.T * tformVolt2Glu.T;
tformStr2Glu = projective2d(T);
RegStrImg_Glu = imwarp(StructureStack(:,:,Z_in), imref2d(size(StructureStack(:,:,Z_in))), tformStr2Glu, 'OutputView', Rglu);
VoltResult.tform_Str2Glu=tformStr2Glu;
VoltResult.Rglu=Rglu;
VoltResult.Z_step=Z_in;

figure(7); clf; tiledlayout(3,1); ax1=[];
ax1=[ax1 nexttile([1 1])];
imshow2(RegStrImg_Glu,[])
ax1=[ax1 nexttile([1 1])];
imshow2(alignedVolt,[])
ax1=[ax1 nexttile([1 1])];
imshow2(GluResult.AvgGluImg,[])
linkaxes([ax1],'xy');

save(fullfile(fpath{f}, 'Volt_Result.mat'), '-struct', 'VoltResult', '-v7.3');
disp('Volt result with Structure registration is saved.')

