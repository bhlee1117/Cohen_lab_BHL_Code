clear
clc;
cd '/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/FromBackup/PP72_PlaceCellResults';
[~, ~, raw] = xlsread(['/Volumes/cohen_lab/Lab/Labmembers' ...
    '/Byung Hun Lee/Data/PrismPCdata_Arrangement.xlsx'], 'Sheet1', 'C5:S31');

% [~, ~, NeuronsToUse]=xlsread(['/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/' ...
%     'PlaceCellData_Arrangement.xlsx'], 'Sheet1', 'L8:M46');
%
% NeuronsToUse=cellfun(@(x) (str2num(num2str(x))),NeuronsToUse,'UniformOutput',false);
ref_ROI=cellfun(@(x) (str2num(num2str(x))),raw(:,10),'UniformOutput',false);
fpath=raw(:,1)';
StructureData=raw(:,9);
EndFrame=cell2mat(raw(:,13));
ifmotionReject=cell2mat(raw(:,14));
ifregressROI=cell2mat(raw(:,17));
save_figto='/Volumes/BHL_WD18TB/PP72_PlaceCellResults';
place_bin=150; time_segment=15000; overlap=200;
alignedMovFN = {'STA_Mat_SS','STA_Mat_CS','STA_Mat_dSP'};
set(0,'DefaultFigureWindowStyle','docked')

%%
 for i = 17%:length(fpath)
        sourceFile = fullfile(fpath{i},'PC_Result.mat');

        % Check if the source file exists
        if exist(sourceFile)
            % Replace the old segment with the new segment in the path
            targetFile = strrep(sourceFile, 'BHL18TB_D2', 'cohen_lab/Lab/Labmembers/Byung Hun Lee/Data');
            
            % Ensure the directory for the target file exists; if not, create it
            targetDir = fileparts(targetFile);
            if ~exist(targetDir, 'dir')
                mkdir(targetDir);
                fprintf('Created target directory: %s\n', targetDir);
            end
            
            % Copy the file to the updated target path
            copyfile(sourceFile, targetFile);
            fprintf('Copied: %s to %s\n', sourceFile, targetFile);
        else
            fprintf('File not found: %s\n', sourceFile);
        end
 end