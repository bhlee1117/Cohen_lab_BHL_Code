function s = flow_val2str(v)
% s = flow_val2str(v)
%
% A parameter value as a MATLAB literal the user can edit in the GUI table, and that
% flow_str2val reads back exactly. Round-trip is the contract: flow_str2val(flow_val2str(v))
% must equal v for every type flow_defaults uses.
%
% 2026.08.28, Byung Hun Lee

if ischar(v)
    s = ['''' strrep(v, '''', '''''') ''''];
elseif isstring(v)
    s = ['"' char(v) '"'];
elseif islogical(v) || isnumeric(v)
    if isempty(v), s = '[]'; else, s = mat2str(v, 6); end
elseif iscell(v)
    s = ['{' strjoin(cellfun(@flow_val2str, v, 'uni', 0), ', ') '}'];
else
    s = ['<' class(v) '>'];
end
end
