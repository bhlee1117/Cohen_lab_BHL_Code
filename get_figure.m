function fh = get_figure(name, varargin)
% fh = get_figure(name, ...)
%
% Get the figure called NAME, or create it if it does not exist yet, bring it to
% front and clear it.  Re-running a cell then UPDATES the existing window instead
% of piling up a new one every time -- which matters here because the script sets
% DefaultFigureWindowStyle to 'docked'.
%
%   figure('Name','X','NumberTitle','off'); clf     % new window on every run
%   get_figure('X')                                 % same window, redrawn
%
% Any extra Name/Value pairs are passed to figure() when the window is created
% (they are ignored when an existing window is reused), e.g.
%   get_figure('R3-seedmovie', 'Color','w', 'WindowStyle','normal')
%
% See also: figure, clf
%
% 2026.08.26, Byung Hun Lee

fh = findobj(groot, 'Type', 'figure', 'Name', name);
if isempty(fh)
    fh = figure('Name', name, 'NumberTitle', 'off', varargin{:});
else
    fh = fh(1);                 % oldest match if the name was somehow duplicated
    figure(fh);
end
clf(fh);
end
