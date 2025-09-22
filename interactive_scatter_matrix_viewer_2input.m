function interactive_scatter_matrix_viewer_2input(Mat1, Mat2, A, B)
%
% This function displays two sets of scatter points (A and B) in the left panel.
% Clicking a point interactively selects the closest point from either set and updates
% the right panel to display a corresponding matrix slice from Mat1 or Mat2.
% Users can navigate through the third dimension of the selected matrix slice using
% the mouse scroll wheel.
%
% USAGE:
%   interactive_scatter_matrix_viewer_2input(Mat1, Mat2, A, B)
%
% INPUTS:
%   Mat1 : N x S1 x T numeric matrix, where N is spatial dimension, S1 is the number
%          of points in A, and T is the number of time slices.
%   Mat2 : N x S2 x T2 numeric matrix, analogous to Mat1 but corresponding to points in B.
%   A    : 2 x S1 numeric array containing [X; Y] coordinates of scatter points for set A.
%   B    : 2 x S2 numeric array containing [X; Y] coordinates of scatter points for set B.
%
% INTERACTIONS:
%   - Left-click a scatter point to display the corresponding matrix slice.
%   - Scroll mouse wheel to navigate through different time slices of the selected point.
%
% EXAMPLE:
%   interactive_scatter_matrix_viewer_2input(rand(100,10,20), rand(100,15,30), rand(2,10), rand(2,15));    

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
    fig = figure('Name', 'Interactive Scatter Plot', 'NumberTitle', 'off', ...
        'WindowScrollWheelFcn', @scroll_callback);

    % Left Panel: Scatter plot
    ax1 = subplot(1,2,1);
    hold on;
    scatterA = scatter(A(1,:), A(2,:), 'bo', 'DisplayName', 'A points');
    scatterB = scatter(B(1,:), B(2,:), 'r.', 'DisplayName', 'B points');
    legend;
    title('Select a point');
    xlabel('X'); ylabel('Y');

    % Right Panel: Image display
    ax2 = subplot(1,2,2);
    img = imagesc(zeros(N, N)); % Placeholder image
    colormap(ax2, gray);
    axis tight off;
    title('Selected Matrix Slice');

    % State tracking
    selectedIdx = NaN;
    selectedSet = NaN; % 1 for Mat1, 2 for Mat2
    currentTime = 1;

    % Set callbacks
    scatterA.ButtonDownFcn = @point_callback;
    scatterB.ButtonDownFcn = @point_callback;

    function point_callback(~, event)
        clickPoint = event.IntersectionPoint(1:2)';
        
        % Find nearest point in A
        [~, idxA] = min(vecnorm(A - clickPoint, 2, 1));
        distA = norm(A(:, idxA) - clickPoint);
        
        % Find nearest point in B
        [~, idxB] = min(vecnorm(B - clickPoint, 2, 1));
        distB = norm(B(:, idxB) - clickPoint);
        
        % Choose the closest point
        if distA < distB
            selectedIdx = idxA;
            selectedSet = 1;
            currentTime = min(currentTime, T);
        else
            selectedIdx = idxB;
            selectedSet = 2;
            currentTime = min(currentTime, T2);
        end
        
        % Update plot
        update_image();
    end

    function scroll_callback(~, event)
        if isnan(selectedIdx), return; end
        % Scroll up/down
        if event.VerticalScrollCount > 0
            currentTime = max(1, currentTime - 1);
        else
            if selectedSet == 1
                currentTime = min(T, currentTime + 1);
            else
                currentTime = min(T2, currentTime + 1);
            end
        end
        update_image();
    end

    function update_image()
        if selectedSet == 1
            img.CData = squeeze(Mat1(:, selectedIdx, :));
        else
            img.CData = squeeze(Mat2(:, selectedIdx, :));
        end
        title(ax2, sprintf('Selected: %s(%d) at T=%d', ...
            char(64+selectedSet), selectedIdx, currentTime));
    end
end
