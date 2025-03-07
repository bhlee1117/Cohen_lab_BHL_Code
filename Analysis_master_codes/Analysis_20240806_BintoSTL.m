clear
clc;
cd '/Volumes/BHL18TB_D2/Arranged_Data/Prism_OptopatchResult';
[~, ~, raw] = xlsread(['/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/' ...
    'Prism_OptopatchData_Arrangement.xlsx'], 'Sheet1', 'B5:K154');

save_to='/Volumes/BHL18TB_D2/Arranged_Data/Prism_OptopatchResult';
fpath=raw(:,1);
Mouse=cell2mat(raw(:,2));
NeuronInd=cell2mat(raw(:,5));
CamType=raw(:,3);
StructureData=raw(:,10);

%%
i=72;
StructureStack=mat2gray(double(tiffreadVolume([StructureData{i} '_Structure.tiff'])));

[x, y, z] = meshgrid(1:xSize, 1:ySize, 1:zSize);
[newX, newY, newZ] = meshgrid(1:2:ySize, 1:2:xSize, 1:4:zSize);
StructureStack_rsz = interp3(double(StructureStack), newX, newY, newZ);
StructureStack_rsz=StructureStack_rsz>0;

iso = isosurface(StructureStack_rsz, 0.5);
faces = iso.faces;
vertices = iso.vertices;
tri = triangulation(faces, vertices);
stlwrite(tri,stlFilename);
