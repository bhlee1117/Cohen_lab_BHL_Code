% Get workspace variables
S = whos;

% Sort by memory (bytes, descending)
[~, idx] = sort([S.bytes], 'descend');
S = S(idx);

% Convert to table for nicer display
T = struct2table(S);
T_MB = T.bytes / 1e6;

% Display variable list
disp('--- Variables sorted by memory (MB) ---')
disp(table(T.name, T.size, T_MB, ...
    'VariableNames', {'Name','Size','Memory_MB'}))

% Total workspace memory
totalBytes = sum([S.bytes]);
totalGB = totalBytes / 1e9;

fprintf('\nTotal workspace memory: %.3f GB\n', totalGB);
