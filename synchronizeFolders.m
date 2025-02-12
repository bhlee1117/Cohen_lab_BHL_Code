function synchronizeFolders(sourceFolder, targetFolder)
    % Synchronize contents of sourceFolder with targetFolder
    
    % Validate folders
    if ~isfolder(sourceFolder)
        error('Source folder does not exist: %s', sourceFolder);
    end
    if ~isfolder(targetFolder)
        mkdir(targetFolder); % Create target folder if it doesn't exist
    end

    % Get list of files in source and target folders
    sourceFiles = dir(fullfile(sourceFolder, '**', '*')); % Include subfolders
    targetFiles = dir(fullfile(targetFolder, '**', '*'));

    % Convert file lists to containers for quick lookup
    sourceFileMap = createFileMap(sourceFiles, sourceFolder);
    targetFileMap = createFileMap(targetFiles, targetFolder);

    % Synchronize files
    for sourceFile = sourceFileMap.keys
        sourceFilePath = sourceFileMap(sourceFile{1});
        targetFilePath = fullfile(targetFolder, sourceFile{1});

        % Check if the file exists in the target folder and is identical
        if ~isKey(targetFileMap, sourceFile{1}) || ...
           ~filesAreEqual(sourceFilePath, targetFilePath)
            % Copy file if it doesn't exist or is different
            targetDir = fileparts(targetFilePath);
            if ~isfolder(targetDir)
                mkdir(targetDir); % Create directory if needed
            end
            copyfile(sourceFilePath, targetFilePath);
            fprintf('Copied: %s -> %s\n', sourceFilePath, targetFilePath);
        end
    end

    % Remove extra files in the target folder
    for targetFile = targetFileMap.keys
        if ~isKey(sourceFileMap, targetFile{1})
            delete(targetFileMap(targetFile{1}));
            fprintf('Deleted: %s\n', targetFileMap(targetFile{1}));
        end
    end
end

function fileMap = createFileMap(files, rootFolder)
    % Helper function to create a map of relative paths to full paths
    fileMap = containers.Map;
    for i = 1:numel(files)
        if ~files(i).isdir
            % Get the relative path to the root folder
            relativePath = strrep(fullfile(files(i).folder, files(i).name), ...
                                  [rootFolder, filesep], '');
            fileMap(relativePath) = fullfile(files(i).folder, files(i).name);
        end
    end
end

function isEqual = filesAreEqual(file1, file2)
    % Helper function to compare two files based on size and modification time
    if ~exist(file1, 'file') || ~exist(file2, 'file')
        isEqual = false;
        return;
    end

    % Get file information
    fileInfo1 = dir(file1);
    fileInfo2 = dir(file2);

    % Compare size and modification time
    isEqual = (fileInfo1.bytes == fileInfo2.bytes) && ...
              abs(fileInfo1.datenum - fileInfo2.datenum) < 1e-6; % Allow minor time differences
end