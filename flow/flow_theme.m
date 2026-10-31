function flow_theme(fig)
% flow_theme(fig)
% flow_theme()          % applies to every open figure
%
% Force the light theme. MATLAB R2025a+ figures follow the desktop theme by default,
% which turns these diagnostic figures dark and makes the blue-white-red correlation maps
% and the white-on-vessel contours hard to read -- and unlike the axes colours, the
% theme is not something the plotting code sets.
%
% A per-figure setting is the only way: MATLAB explicitly refuses a groot default for
% Theme ("SET or GET of default values is not permitted on the Theme property"). Silently
% does nothing on releases without the property.
%
% 2026.08.28, Byung Hun Lee

if nargin < 1 || isempty(fig), fig = findobj(groot, 'Type', 'figure'); end
for k = 1:numel(fig)
    f = fig(k);
    if ~isvalid(f), continue, end
    try
        if isprop(f, 'Theme'), set(f, 'Theme', 'light'); end
    catch
        % older release, or a figure type that has no theme: nothing to do
    end
end
end
