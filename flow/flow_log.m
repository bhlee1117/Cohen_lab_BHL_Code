function lines = flow_log(ax, lines, action, txt, isErr)
% lines = flow_log(ax, [],    'init')               prepare an axes to be used as a log
% lines = flow_log(ax, lines, 'add', txt, isErr)    append and redraw
%
% A colour-capable log panel. A multiline uicontrol edit box cannot colour individual
% lines, and a classic-figure uitable cannot take addStyle, so the log is drawn as text
% objects in an axes -- which gives per-line colour for free.
%
% Errors are RED and prefixed, ordinary messages black, older ones fade to grey. Without
% that, a failure scrolls past looking exactly like a progress message.
%
% lines is an [n x 1] struct array with fields .txt and .err; the caller owns it and
% passes it back in, so the log survives guidata round-trips.
%
% 2026.08.28, Byung Hun Lee

MAXKEEP = 300;      % lines retained
MAXSHOW = 16;       % lines drawn (newest first)

if strcmp(action, 'init')
    cla(ax, 'reset');
    set(ax, 'XLim',[0 1], 'YLim',[0 1], 'XTick',[], 'YTick',[], ...
            'Color',[1 1 1], 'Box','on');
    lines = struct('txt', {}, 'err', {});
    return
end

new = cellstr(splitlines(char(txt)));
for k = 1:numel(new)
    lines(end+1).txt = new{k};   %#ok<AGROW>
    lines(end).err   = isErr;
end
if numel(lines) > MAXKEEP, lines = lines(end-MAXKEEP+1:end); end

cla(ax);
set(ax, 'XLim',[0 1], 'YLim',[0 1], 'XTick',[], 'YTick',[], 'Color',[1 1 1], 'Box','on');
n = min(MAXSHOW, numel(lines));
for k = 1:n
    L = lines(end-k+1);                      % newest at the top
    if L.err
        col = [0.80 0 0];  w = 'bold';
    elseif k <= 3
        col = [0 0 0];     w = 'normal';     % the most recent few stay black
    else
        col = [0.45 0.45 0.45];  w = 'normal';
    end
    text(ax, 0.012, 1 - (k-0.5)/MAXSHOW, L.txt, 'Units','normalized', ...
        'FontName','Courier', 'FontSize', 8.5, 'Color', col, 'FontWeight', w, ...
        'VerticalAlignment','middle', 'Interpreter','none', 'Clipping','on');
end
end
