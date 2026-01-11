function fpath_mac = fpath_window2mac(fpath)
%FPATH_WINDOW2MAC Convert Windows path to Mac path.
%   If the path contains a folder named 'Lab', it maps:
%     X:\Lab\...  ->  /Volumes/Lab/...
%   Otherwise it just converts \ to /.

sp = split(fpath, '\');
i = 1;

try
    while isempty(strfind(sp{i}, 'Lab'))
        i = i + 1;
    end
    ref = i;

    fpath_mac = '/Volumes/Lab';
    for j = ref+1:length(sp)
        if ~isempty(sp{j})
            fpath_mac = [fpath_mac '/' sp{j}];
        end
    end

catch
    % Fallback: just replace separators
    fpath_mac = '';
    for j = 1:length(sp)
        if ~isempty(sp{j})
            fpath_mac = [fpath_mac '/' sp{j}];
        end
    end
end
end
