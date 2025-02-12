function modifiedMarkerMatrix = interactiveFrameCheck(dataMatrix, markerMatrix, windowSize)
    % dataMatrix: N x T matrix (data)
    % markerMatrix: k x T binary matrix (markers for k features)
    % windowSize: Half-width of the window around each marked frame to display
    % lineAlpha: Alpha (transparency) value for xlines (0 to 1)

    lineAlpha=0.4;
    % Validate inputs
    [N, T] = size(dataMatrix);
    k = size(markerMatrix, 1);
    assert(size(markerMatrix, 2) == T, ...
        'Marker matrix must have dimensions k x T (k features, T time points).');

    % Generate unique colors for k features
    colors = lines(k); % Generates k unique colors

    % Get marked frames for navigation
    allFrames = getAllMarkedFrames();

    % Initialize current index and color scale
    currentIndex = 1;
    colorLimits = [prctile(dataMatrix(:),2), prctile(dataMatrix(:),99.5)]; % Set consistent color scale

    % Create a figure for visualization
    hFig = figure('Name', 'Interactive Frame Checker (k Features)', ...
                  'NumberTitle', 'off', ...
                  'KeyPressFcn', @keyPressHandler, ...
                  'WindowButtonDownFcn', @mouseClickHandler);

    % Display the first frame
    displayFrame();

    % Wait until the figure is closed
    uiwait(hFig);

    % Output the modified markerMatrix
    modifiedMarkerMatrix = markerMatrix;

    function allFrames = getAllMarkedFrames()
        % Combine and sort all marked frames
        markedFrames = [];
        for featureIdx = 1:k
            markedFrames = [markedFrames, find(markerMatrix(featureIdx, :) == 1)];
        end
        allFrames = unique(markedFrames); % Sort and remove duplicates
    end

    function displayFrame()
        % Ensure the current index is valid
        if currentIndex > 0 && currentIndex <= numel(allFrames)
            frame = allFrames(currentIndex);
            startIdx = max(1, frame - windowSize);
            endIdx = min(T, frame + windowSize);

            % Extract and display the data segment
            dataSegment = dataMatrix(:, startIdx:endIdx);

            clf;
            imagesc(dataSegment, colorLimits);
            colormap('jet');
            colorbar;
            xlabel('Time');
            ylabel('Channels');
            title(sprintf('Frame %d (Features: %s)', frame, getFeatureTypes(frame)));
            hold on;

            % Highlight the current marked frame in red
            plotXLine(frame - startIdx + 1, [1, 0, 0], lineAlpha);

            % Highlight other marked frames within the window
            highlightOtherMarkedFrames(startIdx, endIdx, frame, startIdx);
        else
            disp('No more frames to display.');
        end
    end

    function featureTypes = getFeatureTypes(frame)
        % Determine the feature types for the given frame
        featureTypes = find(markerMatrix(:, frame) == 1);
        featureTypes = sprintf('%d ', featureTypes);
    end

    function highlightOtherMarkedFrames(startIdx, endIdx, currentFrame, displayStart)
        % Highlight other marked frames in different colors
        for featureIdx = 1:k
            featureFrames = find(markerMatrix(featureIdx, startIdx:endIdx) == 1) + startIdx - 1;

            % Exclude the current frame
            featureFrames(featureFrames == currentFrame) = [];

            % Plot frames in their corresponding feature color
            for f = featureFrames
                plotXLine(f - displayStart + 1, colors(featureIdx, :), lineAlpha);
            end
        end
    end

    function plotXLine(x, color, alpha)
        % Helper function to plot xline with transparency
        line([x, x], ylim, 'Color', [color, alpha], 'LineWidth', 2);
    end

    function keyPressHandler(~, event)
        % Handle keyboard input
        frame = allFrames(currentIndex);
        switch event.Key
            case 'rightarrow'
                currentIndex = min(currentIndex + 1, numel(allFrames));
                displayFrame();
            case 'leftarrow'
                currentIndex = max(currentIndex - 1, 1);
                displayFrame();
            case 'u' % Unmark the current frame
                if ~isempty(allFrames)
                    markerMatrix(:, frame) = 0; % Unmark all features for this frame
                    disp(['Unmarked frame: ', num2str(frame)]);
                    allFrames = getAllMarkedFrames();
                    currentIndex = min(currentIndex, numel(allFrames));
                end
            case 'q' % Quit
                disp('Exiting interactive frame checker...');
                close(hFig);
            otherwise
                % Handle number keys to mark as specific features
                numKey = str2double(event.Key);
                if ~isnan(numKey) && numKey >= 1 && numKey <= k
                    markerMatrix(:, frame) = 0; % Clear all other markings for this frame
                    markerMatrix(numKey, frame) = 1; % Mark this frame as the selected feature
                    disp(['Marked frame ', num2str(frame), ' as Feature ', num2str(numKey)]);
                    allFrames = getAllMarkedFrames();
                end
        end
    end

    function mouseClickHandler(~, ~)
        % Handle mouse clicks
        clickPoint = get(gca, 'CurrentPoint');
        clickedTimeIndex = round(clickPoint(1, 1)); % Time index clicked
        frame = allFrames(currentIndex);
        startIdx = max(1, frame - windowSize);

        % Map clicked point to the original data matrix
        actualIndex = startIdx + clickedTimeIndex - 1;

        if actualIndex >= 1 && actualIndex <= T
            % Display dialog for marking
            featureList = arrayfun(@(idx) sprintf('Feature %d', idx), 1:k, 'UniformOutput', false);
            choice = listdlg('PromptString', 'Select feature to mark:', ...
                             'SelectionMode', 'single', ...
                             'ListString', featureList, ...
                             'CancelString', 'Cancel');
            if ~isempty(choice)
                markerMatrix(choice, actualIndex) = 1; % Mark selected feature
                markerMatrix(setdiff(1:k, choice), actualIndex) = 0; % Ensure exclusivity
                disp(['Marked frame ', num2str(actualIndex), ' as Feature ', num2str(choice)]);
                allFrames = getAllMarkedFrames();
                currentIndex = find(allFrames == actualIndex, 1); % Jump to the new frame
            end
        else
            disp('Click outside the valid range.');
        end
    end
end
