function interactive_scatter_matrix_viewer(Mat1, A, B)
% INTERACTIVE_SCATTER_MATRIX_VIEWER Visualizes data by linking scatter plot to matrix slices.
%
%   interactive_scatter_matrix_viewer(Mat1, A, B) creates an interactive
%   figure where A and B define the scatter plot, and Mat1 can be either
%   a 3D numeric array (N x T x R) or a cell array where each cell is N_i x T_i.

    % Validate A/B
    assert(isvector(A) && isvector(B) && numel(A) == numel(B), ...
        'A and B must be vectors of equal length');

    useCell = iscell(Mat1);
    if useCell
        assert(numel(Mat1) == numel(A), 'Length of cell array Mat1 must match A/B');
    else
        assert(ndims(Mat1) == 3, 'Mat1 must be N x T x R if numeric');
        assert(size(Mat1, 3) == numel(A), 'Third dimension of Mat1 must match A/B');
    end

    % Precompute color axis
    if useCell
        allvals=cell2mat(cellfun(@(x) x(:),Mat1,'UniformOutput',false)');
        %allvals = vertcat(Mat1{:});
        allvals = allvals(:);
    else
        allvals = Mat1(:);
    end
    cax = [prctile(allvals, 5), prctile(allvals, 99)];

    % Create figure
    fig = figure('Name', 'Interactive Matrix Viewer', 'NumberTitle', 'off', ...
                 'WindowScrollWheelFcn', @scroll_callback);

    % Left panel: Scatter plot
    ax1 = subplot(2,2,[1 3]);
    scatterPlot = scatter(A, B, 'filled');
    hold(ax1, 'on'); grid on;
    highlightHandle = plot(ax1, NaN, NaN, 'ro', 'LineWidth', 2, 'MarkerSize', 10);
    hold(ax1, 'off');
    title('Click a point');
    xlabel('A'); ylabel('B');
    axis tight;

    % Top-right panel: Image
    ax2 = subplot(2,2,2);
    imgHandle = imagesc(nan, cax);
    colormap(ax2, parula);
    axis tight;
    title('Matrix Slice (N x T)');

    % Bottom-right panel: Traces
    ax3 = subplot(2,2,4);
    traceHandles = gobjects(0,1);
    hold(ax3, 'on');
    hold(ax3, 'off');
    title('All Traces');
    xlabel('Time'); ylabel('Signal');
    grid on;

    % Link axes
    linkaxes([ax2, ax3], 'x');

    % State
    selectedIdx = NaN;

    % Click callback
    scatterPlot.ButtonDownFcn = @click_callback;

    function click_callback(~, event)
        pt = event.IntersectionPoint(1:2);
        dist = vecnorm([A(:) B(:)]' - pt(:), 2, 1);
        [~, selectedIdx] = min(dist);
        update_panels();
    end

    function scroll_callback(~, ~)
        % Reserved for future use
    end

    function update_panels()
        if isnan(selectedIdx), return; end

        % Update red circle on selected point
        set(highlightHandle, 'XData', A(selectedIdx), 'YData', B(selectedIdx));

        % Update matrix slice
        if useCell
            matSlice = Mat1{selectedIdx};
        else
            matSlice = Mat1(:,:,selectedIdx);
        end

        [N, T] = size(matSlice);
        set(imgHandle, 'CData', matSlice);
        set(ax2, 'XLim', [0.5, T+0.5], 'YLim', [0.5, N+0.5]);

        % Clear and replot traces
        cla(ax3);
        hold(ax3, 'on');
        cmap = turbo(N);
        traceHandles = gobjects(N,1);
        for i = 1:N
            traceHandles(i) = plot(ax3, 1:T, matSlice(i,:), 'Color', cmap(i,:));
        end
        hold(ax3, 'off');
        xlim(ax3, [1 T]);
        ylim(ax3, [min(matSlice(:)), max(matSlice(:))]);

        title(ax2, sprintf('Slice #%d (N x T)', selectedIdx));
        title(ax3, sprintf('All Traces (Slice #%d)', selectedIdx));
    end
end
