function interactive_scatter_multiselectViewer(Mat1, A, B)
% interactive_scatter_multiselectViewer - Visualize average image and trace for brushed points.
%   Mat1: N x T x R (data)
%   A, B: 1 x R vectors for scatter plot

assert(isvector(A) && isvector(B) && numel(A) == numel(B), 'A and B must be vectors of equal length');
assert(ndims(Mat1) == 3 && size(Mat1,3) == numel(A), 'Mat1 must be N x T x R matching A and B');

[N, T, R] = size(Mat1);

fig = figure('Name', 'Multi-selection Scatter Viewer', 'NumberTitle', 'off');

% -- SCATTER PLOT --
ax1 = subplot(2,2,[1 3]);
scatterPlot = scatter(A, B, 40, 'k', 'filled');
xlabel('A'); ylabel('B'); title('Brush points and press Confirm');
brush on;
hold on;
highlightPlot = scatter(nan, nan, 60, 'ro', 'LineWidth', 1.5);
hold off;
grid on;

% -- AVERAGE IMAGE PANEL --
ax2 = subplot(2,2,2);
avgImg = imagesc(zeros(N, T));
colormap(ax2, parula);
colorbar;
axis tight off;
title('Average Matrix Slice');

% -- AVERAGE TRACE PANEL --
ax3 = subplot(2,2,4);
avgTrace = plot(nan, nan, 'k-', 'LineWidth', 1.5);
grid on;
xlabel('Time'); ylabel('Signal');
title('Average Trace');
axis tight;

% -- BUTTON CONTROL --
uicontrol('Style', 'pushbutton', 'String', 'Confirm Selection', ...
    'Position', [20 20 140 30], 'Callback', @confirmSelection);

% -- Selection logic --
    function confirmSelection(~, ~)
        bd = get(scatterPlot, 'BrushData');
        if isempty(bd) || ~any(bd)
            warning('No points selected. Use brush tool to select.');
            return;
        end
        idx = find(logical(bd(:)));
        set(highlightPlot, 'XData', A(idx), 'YData', B(idx));
        updatePanels(idx);
    end

    function updatePanels(idx)
        if isempty(idx), return; end

        mat = Mat1(:,:,idx);                  % N x T x R'
        avgMat = mean(mat, 3,'omitnan');                % N x T
        %avgTr = mean(avgMat, 1);              % 1 x T

        % Update image (top-right)
        set(avgImg, 'CData', avgMat);
        ax2.CLim = prctile(avgMat(:), [5 95]);

        % Update trace plot (bottom-right)
        cla(ax3); hold(ax3, 'on');
        traces = squeeze(mean(mat,3,'omitnan'));        % N x T'
        l=plot(ax3, traces');  % All traces in light gray
        arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(turbo(N),2))
        %plot(ax3, 1:T, avgTr, 'k-', 'LineWidth', 1.5);  % Average in black
        hold(ax3, 'off');
        grid(ax3, 'on');
        xlabel(ax3, 'Time'); ylabel(ax3, 'Signal');
        title(ax3, 'Average of red points');
        xlim(ax3, [1 T]);
        ylim(ax3, [min(traces(:)), max(traces(:))]);
    end
end
