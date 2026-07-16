function G = loadGlu(fn, varargin)
% LOADGLU  Load Glu_Result.mat but skip the large raw-trace matrices that the
% STA/place-tuning pipeline never uses (F0_glu, T_glu ~ 1.9 GB each), so loading
% is ~3x faster. Drop-in replacement for importdata(fn).
%   G = loadGlu(fn)        % loads everything except F0_glu, T_glu
%   G = loadGlu(fn,'all')  % loads everything
skip = {'F0_glu','T_glu'};
if any(strcmpi(varargin,'all')), skip = {}; end
vars = who('-file', fn);              % top-level variables (fields, saved with -struct)
vars = setdiff(vars, skip);
G = load(fn, vars{:});
end
