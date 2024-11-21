function backupServer(fpath,ChangeStr,targetStr,Name)
        sourceFile = fullfile(fpath,Name);

        % Check if the source file exists
        if exist(sourceFile)
            % Replace the old segment with the new segment in the path
            targetFile = strrep(sourceFile, ChangeStr, targetStr);
            
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