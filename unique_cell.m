function result = unique_cell(C)
% UNIQUE_CELL  Find unique and repeated entries in a cell array of strings.
%
%   result = unique_cell(C)
%
%   INPUT
%     C            : cell array of strings or NaN values
%
%   OUTPUT  (struct)
%     result.labels          : N x 1 cell, NaN values converted to 'NaN'
%     result.uniqueLabels    : unique strings (stable/first-appearance order)
%     result.groupIdx        : N x 1 index — same number = same label
%     result.uniqueOnly      : indices whose label appears exactly once
%     result.repeatedAll     : indices that share a label with at least one other
%     result.repeatedGroups  : cell array, each entry = indices sharing a label
%
%   EXAMPLE
%     result = unique_cell(StructureData);
%     disp(result.uniqueOnly)
%     disp(result.repeatedGroups)

N = numel(C);

% Normalise: convert NaN entries to string 'NaN'
labels = cell(N, 1);
for f = 1:N
    if ischar(C{f}) || isstring(C{f})
        labels{f} = char(C{f});
    else
        labels{f} = 'NaN';
    end
end

% Unique labels and group membership
[uniqueLabels, ~, groupIdx] = unique(labels, 'stable');

% Count occurrences per label
counts = accumarray(groupIdx, 1);

% Indices with a one-off label vs shared label
uniqueOnly  = find(ismember(groupIdx, find(counts == 1)));
repeatedAll = find(ismember(groupIdx, find(counts  > 1)));

% Group repeated indices together
repeatedGroups = {};
for u = find(counts > 1)'
    repeatedGroups{end+1} = find(groupIdx == u);  %#ok<AGROW>
end

% Print summary
fprintf('=== unique_cell: %d entries, %d unique labels ===\n', N, numel(uniqueLabels));
for u = 1:numel(uniqueLabels)
    memberIdx = find(groupIdx == u);
    if numel(memberIdx) == 1
        fprintf('  [UNIQUE]   "%s"  → index: %s\n', uniqueLabels{u}, mat2str(memberIdx'));
    else
        fprintf('  [REPEATED] "%s"  → indices: %s\n', uniqueLabels{u}, mat2str(memberIdx'));
    end
end
fprintf('  Unique-only indices:   %s\n', mat2str(uniqueOnly'));
fprintf('  Repeated indices:      %s\n', mat2str(repeatedAll'));

% Pack output
result.labels         = labels;
result.uniqueLabels   = uniqueLabels;
result.groupIdx       = groupIdx;
result.uniqueOnly     = uniqueOnly;
result.repeatedAll    = repeatedAll;
result.repeatedGroups = repeatedGroups;

end
