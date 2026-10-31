function [v, ok, msg] = flow_str2val(s)
% [v, ok, msg] = flow_str2val(s)
%
% Read a parameter literal typed into the GUI table.
%
% Only LITERALS are accepted -- numbers, strings, logicals, and bracketed/braced
% combinations of them. Anything that looks like a function call, an assignment, a
% statement separator or a workspace reference is rejected before evaluation. The input
% comes from the user's own keyboard in their own session, so this is not a security
% boundary; it is there so that a typo in a parameter box reports "could not parse" in
% the log instead of executing something surprising or throwing out of a callback.
%
% ok = false and msg explains why; v is [] in that case.
%
% 2026.08.28, Byung Hun Lee

v = [];  ok = false;  msg = '';
s = strtrim(char(s));
if isempty(s), msg = 'empty'; return, end
if numel(s) > 500,  msg = 'too long for a parameter'; return, end

% reject anything that is not a bare literal
if ~isempty(regexp(s, '[;=@]', 'once'))
    msg = 'no semicolons, assignments or function handles'; return
end
if ~isempty(regexp(s, '[A-Za-z_]\w*\s*\(', 'once'))
    msg = 'function calls are not allowed'; return
end
% the only bare words permitted are the logical constants
words = regexp(s, '[A-Za-z_]\w*', 'match');
inStr = regexp(s, '''[^'']*''|"[^"]*"', 'match');
for k = 1:numel(words)
    w = words{k};
    if any(strcmp(w, {'true','false','Inf','inf','NaN','nan','pi','eps'})), continue, end
    if any(cellfun(@(q) contains(q, w), inStr)), continue, end   % inside a quoted string
    msg = sprintf('"%s" is not a literal', w); return
end

try
    v = eval(s);
    ok = true;
catch ME
    msg = ME.message;
end
end
