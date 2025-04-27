function interactive_multi_select_viewer(Mat1, Mat2, A, B)
    % Validate input dimensions
    assert(size(A,1) == 2, 'A must be 2 x S1');
    assert(size(B,1) == 2, 'B must be 2 x S2');
    assert(size(Mat1,1) == size(Mat2,1), 'N must be the same for both matrices');
    
    N = size(Mat1,1);
    S1 = size(Mat1,2);
    S2 = size(Mat2,2);
    T = size(Mat1,3);
    T2 = size(Mat2,3);

    % Create figure
    fig = figure('Name', 'Multi-Select Scatter Plot', 'NumberTitle', 'off');

    % Left Panel: Scatter plot
    ax1 = subplot(1,2,1);
    hold on;
    scatterA = scatter(A(1,:), A(2,:), 'bo', 'DisplayName', 'A points');
    scatterB = scatter(B(1,:), B(2,:), 'r.', 'DisplayName', 'B points');
    legend;
    title('Click & Drag to Select Points');
    xlabel('X'); ylabel('Y');

    % Right Panel: Image display
    ax2 = subplot(1,2,2);
    img = imagesc(zeros(N, T)); % Placeholder image
    colormap(ax2, gray);
    axis off;
    title('Averaged Matrix');

    % Store selection state
    selectedIndicesA = [];
    selectedIndicesB = [];
    selectedSetA = false;
    selectedSetB = false;

    % Set callbacks for selection
    rect = drawrectangle(ax1, 'Visible', 'off', 'Color', 'g');
    fig.WindowButtonDownFcn = @start_select;
    fig.WindowButtonUpFcn = @end_select;

    function start_select(~, ~)
        rect.Position = [0, 0, 0, 0];
        rect.Visible = 'on';
        rect.Position = [ax1.CurrentPoint(1,1:2), 0, 0];
        fig.WindowButtonMotionFcn = @update_rect;
    end

    function update_rect(~, ~)
        pos = rect.Position;
        cp = ax1.CurrentPoint(1,1:2);
        rect.Position = [pos(1), pos(2), cp(1)-pos(1), cp(2)-pos(2)];
    end

    function end_select(~, ~)
        fig.WindowButtonMotionFcn = [];
        rect.Visible = 'off';
        pos = rect.Position;

        % Find selected points within the rectangle
        selectedIndicesA = find(A(1,:) >= pos(1) & A(1,:) <= pos(1) + pos(3) & ...
                                A(2,:) >= pos(2) & A(2,:) <= pos(2) + pos(4));
        selectedIndicesB = find(B(1,:) >= pos(1) & B(1,:) <= pos(1) + pos(3) & ...
                                B(2,:) >= pos(2) & B(2,:) <= pos(2) + pos(4));
        
        selectedSetA = ~isempty(selectedIndicesA);
        selectedSetB = ~isempty(selectedIndicesB);

        update_image();
    end

    function update_image()
        if selectedSetA && selectedSetB
            % If points from both sets are selected, compute mean of both
            avgMat = (mean(Mat1(:, selectedIndicesA, :), 2,'omitnan') + ...
                      mean(Mat2(:, selectedIndicesB, :), 2,'omitnan')) / 2;
        elseif selectedSetA
            avgMat = mean(Mat1(:, selectedIndicesA, :), 2,'omitnan');
        elseif selectedSetB
            avgMat = mean(Mat2(:, selectedIndicesB, :), 2,'omitnan');
        else
            avgMat = zeros(N, T);
        end

        % Reshape to N x T and update image
        img.CData = squeeze(avgMat);
        title(ax2, sprintf('Average of %d (A) and %d (B) points', length(selectedIndicesA), length(selectedIndicesB)));
    end
end
