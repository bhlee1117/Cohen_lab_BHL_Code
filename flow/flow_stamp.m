function h = flow_stamp(varargin)
% h = flow_stamp(cfgSubStruct, upstreamStamp, ...)
%
% A short hash of everything a stage depends on: its own parameters, plus the stamps of
% the stages it consumes. Two runs share a stamp exactly when nothing that could change
% the result has changed.
%
% This is what lets flow_pipeline skip work safely. The script version had no way to know
% whether a cached variable was still valid, so the only trustworthy move was to re-run
% everything -- and the tempting move was to re-run one cell and hope, which is how a
% stale variable silently produces a wrong answer. Chaining the upstream stamp in means a
% change to the band-pass invalidates the tracking, the segments, and everything after,
% without anyone having to remember that.
%
% NOT a checksum of the DATA. If you point the pipeline at a different movie with the same
% settings the stamp is unchanged -- flow_pipeline stores a data tag separately.
%
% INPUT   any mix of structs, char, numeric, logical, cell
% OUTPUT  h   16-character hex string
%
% 2026.08.28, Byung Hun Lee

s = '';
for k = 1:numel(varargin)
    s = [s '|' ser(varargin{k})];  %#ok<AGROW>
end
% Two independent polynomial hashes, 8 hex digits each. NOT uint64 FNV: MATLAB's integer
% arithmetic SATURATES instead of wrapping, so the usual multiply-and-overflow hash pins
% at intmax and every input collapses to the same value. These stay inside double's exact
% integer range (2^31 * 131 < 2^53), so the modular arithmetic is exact.
h = lower([dec2hex(polyHash(s, 2166136261, 131), 8), ...
           dec2hex(polyHash(s,      5381, 33), 8)]);
end

function s = ser(v)
% deterministic text form of a value; struct fields are sorted so field ORDER never
% changes the stamp
if isstruct(v)
    if numel(v) ~= 1
        s = ''; for k = 1:numel(v), s = [s ';' ser(v(k))]; end %#ok<AGROW>
        return
    end
    f = sort(fieldnames(v));
    s = '';
    for k = 1:numel(f), s = [s f{k} '=' ser(v.(f{k})) ';']; end %#ok<AGROW>
elseif iscell(v)
    s = '{'; for k = 1:numel(v), s = [s ser(v{k}) ',']; end, s = [s '}']; %#ok<AGROW>
elseif ischar(v) || isstring(v)
    s = char(v);
elseif islogical(v)
    s = sprintf('%d', v(:));
elseif isnumeric(v)
    if isempty(v), s = '[]';
    else,          s = sprintf('%.10g,', double(v(:)));
    end
else
    s = class(v);              % function handles etc: identity by type only
end
end

function u = polyHash(s, seed, mult)
% Horner polynomial hash mod (2^31 - 1), in double so nothing saturates. Self-contained,
% so the pipeline does not depend on a Java or toolbox hash being available.
M = 2147483647;                      % 2^31 - 1, prime
u = mod(seed, M);
b = double(uint8(s));
for k = 1:numel(b)
    u = mod(u*mult + b(k), M);       % u*mult < 2^31 * 131 < 2^53: exact in double
end
end
