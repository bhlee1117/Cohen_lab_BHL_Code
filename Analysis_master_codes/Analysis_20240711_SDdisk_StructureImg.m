clear; clc;
cd '/Volumes/BHL18TB_D2/20240730_Structure/20240731_structural_SD'
folder='/Volumes/BHL18TB_D2/20240730_Structure/20240731_structural_SD';
 binFilePaths = findBinFiles(folder);
 [folderPart, fileName, fileExt] = fileparts(binFilePaths);
 %%
 for f=1:length(binFilePaths)
     load(fullfile(folderPart{f},"output_data.mat"));
     sz=double(Device_Data{1, 3}.ROI([2 4]));
     imStack=readBinMov([binFilePaths{f}],sz(2),sz(1));
     writeTiff(imStack,[folderPart{f} '.tiff'])
 end