function pathOut = convertDrivePath(pathIn, targetDrive)
% convertDrivePath  Convert a data path to point to a different drive
%
%   pathOut = convertDrivePath(pathIn, targetDrive)
%
%   targetDrive options:
%     'cohen_lab'  -> /Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data
%     'BHL18TB_D2' -> /Volumes/BHL18TB_D2
%     'BHL18TB_D1' -> /Volumes/BHL18TB_D1        % add more as needed
%
%   Example:
%     path2 = convertDrivePath(path1, 'BHL18TB_D2')
%     path1 = convertDrivePath(path2, 'cohen_lab')

% --- Define all known drive roots here ---
drives = struct();
drives.cohen_lab  = '/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data';
drives.BHL18TB_D2 = '/Volumes/BHL18TB_D2';
drives.BHL18TB_D1 = '/Volumes/BHL18TB_D1';

% --- Validate target drive ---
driveNames = fieldnames(drives);
if ~ismember(targetDrive, driveNames)
    error('Unknown drive "%s". Available options: %s', ...
          targetDrive, strjoin(driveNames, ', '));
end

% --- Detect which drive the input path is currently on ---
sourceRoot = '';
for i = 1:length(driveNames)
    root = drives.(driveNames{i});
    if startsWith(pathIn, root)
        sourceRoot = root;
        break;
    end
end

if isempty(sourceRoot)
    error('Could not detect source drive from path:\n  %s', pathIn);
end

% --- Replace source root with target root ---
targetRoot = drives.(targetDrive);
pathOut = strrep(pathIn, sourceRoot, targetRoot);

end