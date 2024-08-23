function binFilePaths = findBinFiles(folder)
    % This function finds all .bin files in the specified folder and its subfolders.
    % Input:
    %   folder - The root folder where the search will begin.
    % Output:
    %   binFilePaths - A cell array containing the full paths of all .bin files found.
    
    % Initialize cell array to hold file paths
    binFilePaths = {};
    
    % Get list of all files and folders in the folder
    filesAndFolders = dir(folder);
    
    % Loop through each item in the directory
    for i = 1:length(filesAndFolders)
        % Get the name of the item
        itemName = filesAndFolders(i).name;
        
        % Skip '.' and '..' directories
        if strcmp(itemName, '.') || strcmp(itemName, '..')
            continue;
        end
        
        % Get the full path of the item
        itemFullPath = fullfile(folder, itemName);
        
        % Check if the item is a folder
        if filesAndFolders(i).isdir
            % Recursively search in the subfolder
            subfolderBinFiles = findBinFiles(itemFullPath);
            % Append the results to the main list
            binFilePaths = [binFilePaths; subfolderBinFiles];
        else
            % Check if the file has a .bin extension
            [~, ~, ext] = fileparts(itemFullPath);
            if strcmp(ext, '.bin')
                % Add the file path to the list
                binFilePaths{end+1} = itemFullPath;
            end
        end
    end
    
    % Display the results
    if isempty(binFilePaths)
        disp('No .bin files found.');
    else
        disp('Found .bin files:');
        for i = 1:length(binFilePaths)
            disp(binFilePaths{i});
        end
    end
end