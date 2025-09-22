function modifiedMarkerMatrix = interactiveFrameCheck(dataMatrix, markerMatrix, windowSize, plotMode)
% dataMatrix: N x T matrix (data)
% markerMatrix: k x T binary matrix (markers for k features)
% windowSize: Half-width of the window around each marked frame to display
% plotMode: 'image' or 'trace'

if nargin < 4
    plotMode = 'image'; % Default mode
end
lineAlpha = 0.7;

[N, T] = size(dataMatrix);
k = size(markerMatrix, 1);
assert(size(markerMatrix, 2) == T, 'Marker matrix must be k x T');

disp('press up arrow to unmark, q to quit, number to mark');
colors = lines(k);
allFrames = getAllMarkedFrames();

currentIndex = 1;
colorLimits = [prctile(dataMatrix(:),2), prctile(dataMatrix(:),99.9)];

hFig = figure('Name', 'Interactive Frame Checker', ...
    'NumberTitle', 'off', ...
    'KeyPressFcn', @keyPressHandler, ...
    'WindowButtonDownFcn', @mouseClickHandler);
displayFrame();
uiwait(hFig);
modifiedMarkerMatrix = markerMatrix;

    function allFrames = getAllMarkedFrames()
        markedFrames = [];
        for featureIdx = 1:k
            markedFrames = [markedFrames, find(markerMatrix(featureIdx, :) == 1)];
        end
        allFrames = unique(markedFrames);
    end

    function displayFrame()
        if currentIndex < 1 || currentIndex > numel(allFrames)
            disp('No more frames to display.'); return;
        end

        frame = allFrames(currentIndex);
        startIdx = max(1, frame - windowSize);
        endIdx = min(T, frame + windowSize);
        dataSegment = dataMatrix(:, startIdx:endIdx);
        clf;

        switch lower(plotMode)
            case 'image'
                imagesc(dataSegment, colorLimits);
                colormap('turbo'); colorbar;
                xlabel('Time'); ylabel('Channels');
                title(sprintf('Frame %d (Features: %s)', frame, getFeatureTypes(frame)));
                hold on;
                plotXLine(frame - startIdx + 1, [1, 0, 0], lineAlpha);
                highlightOtherMarkedFrames(startIdx, endIdx, frame, startIdx);

            case 'trace'
                offset = max(abs(dataSegment(:))) * 1.5;
                yOffsets = offset * (N-1:-1:0)';
                stackedTraces = dataSegment + yOffsets;
                hold on;

                x = 1:size(dataSegment,2);
                plot(x, stackedTraces', 'LineWidth', 1.2); grid on;
                axChildren = get(gca, 'Children');
                for ch = 1:N
                    set(axChildren(N - ch + 1), 'Color', colors(mod(ch-1, k)+1, :));
                end

                xlabel('Time');
                xlim([0, 2*windowSize+1]);
                yticks(yOffsets);
                yticklabels(arrayfun(@(ch) sprintf('Ch %d', ch), 1:N, 'UniformOutput', false));
                title(sprintf('Traces at Frame %d (Features: %s)', frame, getFeatureTypes(frame)));

                t = frame - startIdx + 1;
                for f = startIdx:endIdx
                    tIdx = f - startIdx + 1;
                    for featureIdx = 1:k
                        if markerMatrix(featureIdx, f)>0
                            ch = featureIdx;
                            if f > T || tIdx < 1 || tIdx > size(dataSegment,2)
                                continue;
                            end                           
                            yBase = dataMatrix(1, f);
                            isCurrent = (f == frame);
                            sz = 100 * isCurrent + 60 * ~isCurrent;
                            filled = isCurrent;
                            if filled
                                scatter(tIdx, yBase + offset * 0.2, sz, colors(ch,:), 'filled', 'v');
                            else
                                scatter(tIdx, yBase + offset * 0.2, sz, colors(ch,:), 'v');
                            end
                        end
                    end
                end
                ylim([min(yOffsets) - offset/2, max(yOffsets) + offset]);
        end
    end

    function featureTypes = getFeatureTypes(frame)
        featureTypes = find(markerMatrix(:, frame) == 1);
        featureTypes = sprintf('%d ', featureTypes);
    end

    function highlightOtherMarkedFrames(startIdx, endIdx, currentFrame, displayStart)
        for featureIdx = 1:k
            frames = find(markerMatrix(featureIdx, startIdx:endIdx)) + startIdx - 1;
            frames(frames == currentFrame) = [];
            for f = frames
                plotXLine(f - displayStart + 1, colors(featureIdx,:), lineAlpha);
            end
        end
    end

    function plotXLine(x, color, alpha)
        yl = ylim;
        line([x x], yl, 'Color', [color, alpha], 'LineWidth', 2);
    end

    function keyPressHandler(~, event)
        frame = allFrames(currentIndex);
        switch event.Key
            case 'rightarrow'
                currentIndex = min(currentIndex + 1, numel(allFrames));
                displayFrame();
            case 'leftarrow'
                currentIndex = max(currentIndex - 1, 1);
                displayFrame();
            case 'uparrow'
                markerMatrix(:, frame) = 0;
                disp(['Unmarked frame ', num2str(frame)]);
                allFrames = getAllMarkedFrames();
                currentIndex = min(currentIndex, numel(allFrames));
                displayFrame();
            case 'q'
                disp('Exiting interactive viewer...');
                close(hFig);
            otherwise
                numKey = str2double(event.Key);
                if ~isnan(numKey) && numKey >= 1 && numKey <= k
                    markerMatrix(:, frame) = 0;
                    markerMatrix(numKey, frame) = 1;
                    disp(['Marked frame ', num2str(frame), ' as Feature ', num2str(numKey)]);
                    allFrames = getAllMarkedFrames();
                end
        end
    end

    function mouseClickHandler(~, ~)
        clickPoint = get(gca, 'CurrentPoint');
        clickedTimeIndex = round(clickPoint(1, 1));
        frame = allFrames(currentIndex);
        startIdx = max(1, frame - windowSize);
        actualIndex = startIdx + clickedTimeIndex - 1;
        if actualIndex >= 1 && actualIndex <= T
            featureList = arrayfun(@(i) sprintf('Feature %d', i), 1:k, 'UniformOutput', false);
            choice = listdlg('PromptString', 'Select feature to mark:', ...
                'SelectionMode', 'single', ...
                'ListString', featureList, 'CancelString', 'Cancel');
            if ~isempty(choice)
                markerMatrix(choice, actualIndex) = 1;
                markerMatrix(setdiff(1:k, choice), actualIndex) = 0;
                disp(['Marked frame ', num2str(actualIndex), ' as Feature ', num2str(choice)]);
                allFrames = getAllMarkedFrames();
                currentIndex = find(allFrames == actualIndex, 1);
                displayFrame();
            end
        else
            disp('Click outside valid range.');
        end
    end
end

% function modifiedMarkerMatrix = interactiveFrameCheck(dataMatrix, markerMatrix, windowSize, plotMode)
% % dataMatrix: N x T matrix (data)
% % markerMatrix: k x T binary matrix (markers for k features)
% % windowSize: Half-width of the window around each marked frame to display
% % plotMode: 'image' or 'trace'
% 
% if nargin < 4
%     plotMode = 'image'; % Default mode
% end
% lineAlpha = 0.7;
% 
% [N, T] = size(dataMatrix);
% k = size(markerMatrix, 1);
% assert(size(markerMatrix, 2) == T, 'Marker matrix must be k x T');
% 
% disp('press up arrow to unmark, q to quit, number to mark');
% colors = lines(k);
% allFrames = getAllMarkedFrames();
% 
% currentIndex = 1;
% colorLimits = [prctile(dataMatrix(:),2), prctile(dataMatrix(:),99.9)];
% 
% hFig = figure('Name', 'Interactive Frame Checker', ...
%     'NumberTitle', 'off', ...
%     'KeyPressFcn', @keyPressHandler, ...
%     'WindowButtonDownFcn', @mouseClickHandler);
% displayFrame();
% uiwait(hFig);
% modifiedMarkerMatrix = markerMatrix;
% 
%     function allFrames = getAllMarkedFrames()
%         markedFrames = [];
%         for featureIdx = 1:k
%             markedFrames = [markedFrames, find(markerMatrix(featureIdx, :) == 1)];
%         end
%         allFrames = unique(markedFrames);
%     end
% 
%     function displayFrame()
%         if currentIndex < 1 || currentIndex > numel(allFrames)
%             disp('No more frames to display.'); return;
%         end
% 
%         frame = allFrames(currentIndex);
%         startIdx = max(1, frame - windowSize);
%         endIdx = min(T, frame + windowSize);
%         dataSegment = dataMatrix(:, startIdx:endIdx);
%         clf;
% 
%         switch lower(plotMode)
%             case 'image'
%                 imagesc(dataSegment, colorLimits);
%                 colormap('turbo'); colorbar;
%                 xlabel('Time'); ylabel('Channels');
%                 title(sprintf('Frame %d (Features: %s)', frame, getFeatureTypes(frame)));
%                 hold on;
%                 plotXLine(frame - startIdx + 1, [1, 0, 0], lineAlpha);
%                 highlightOtherMarkedFrames(startIdx, endIdx, frame, startIdx);
% 
%             case 'trace'
%                 offset = max(abs(dataSegment(:))) * 1.5;
%                 yOffsets = offset * (N-1:-1:0)';
%                 stackedTraces = dataSegment + yOffsets;
%                 hold on;
% 
%                 % Plot all traces at once using line()
%                 x = 1:size(dataSegment,2);
%                 plot(x, stackedTraces', 'LineWidth', 1.2);
% 
%                 % Set colors for each trace
%                 axChildren = get(gca, 'Children');
%                 for ch = 1:N
%                     set(axChildren(N - ch + 1), 'Color', colors(mod(ch-1, k)+1, :));
%                 end
% 
%                 % Axes and labels
%                 xlabel('Time');
%                 xlim([0, 2*windowSize+1]);
%                 yticks(yOffsets);
%                 yticklabels(arrayfun(@(ch) sprintf('Ch %d', ch), 1:N, 'UniformOutput', false));
%                 title(sprintf('Traces at Frame %d (Features: %s)', frame, getFeatureTypes(frame)));
% 
%                 % Mark current feature
%                 t = frame - startIdx + 1;
%                 activeFeatures = find(markerMatrix(:, frame));
%                 for idx = 1:length(activeFeatures)
%                     ch = activeFeatures(idx);
%                     yBase = dataSegment(ch, t) + offset * (N - ch);
%                     scatter(t, yBase + 0.2 * offset, 100, colors(ch,:), 'filled', 'v');
%                 end
% 
%                 ylim([min(yOffsets) - 2*offset, max(yOffsets) + 2*offset]);
% 
%         end
%     end
% 
%     function featureTypes = getFeatureTypes(frame)
%         featureTypes = find(markerMatrix(:, frame) == 1);
%         featureTypes = sprintf('%d ', featureTypes);
%     end
% 
%     function highlightOtherMarkedFrames(startIdx, endIdx, currentFrame, displayStart)
%         for featureIdx = 1:k
%             frames = find(markerMatrix(featureIdx, startIdx:endIdx)) + startIdx - 1;
%             frames(frames == currentFrame) = [];
%             for f = frames
%                 plotXLine(f - displayStart + 1, colors(featureIdx,:), lineAlpha);
%             end
%         end
%     end
% 
%     function plotXLine(x, color, alpha)
%         yl = ylim;
%         line([x x], yl, 'Color', [color, alpha], 'LineWidth', 2);
%     end
% 
%     function keyPressHandler(~, event)
%         frame = allFrames(currentIndex);
%         switch event.Key
%             case 'rightarrow'
%                 currentIndex = min(currentIndex + 1, numel(allFrames));
%                 displayFrame();
%             case 'leftarrow'
%                 currentIndex = max(currentIndex - 1, 1);
%                 displayFrame();
%             case 'uparrow'
%                 markerMatrix(:, frame) = 0;
%                 disp(['Unmarked frame ', num2str(frame)]);
%                 allFrames = getAllMarkedFrames();
%                 currentIndex = min(currentIndex, numel(allFrames));
%                 displayFrame();
%             case 'q'
%                 disp('Exiting interactive viewer...');
%                 close(hFig);
%             otherwise
%                 numKey = str2double(event.Key);
%                 if ~isnan(numKey) && numKey >= 1 && numKey <= k
%                     markerMatrix(:, frame) = 0;
%                     markerMatrix(numKey, frame) = 1;
%                     disp(['Marked frame ', num2str(frame), ' as Feature ', num2str(numKey)]);
%                     allFrames = getAllMarkedFrames();
%                 end
%         end
%     end
% 
%     function mouseClickHandler(~, ~)
%         clickPoint = get(gca, 'CurrentPoint');
%         clickedTimeIndex = round(clickPoint(1, 1));
%         frame = allFrames(currentIndex);
%         startIdx = max(1, frame - windowSize);
%         actualIndex = startIdx + clickedTimeIndex - 1;
%         if actualIndex >= 1 && actualIndex <= T
%             featureList = arrayfun(@(i) sprintf('Feature %d', i), 1:k, 'UniformOutput', false);
%             choice = listdlg('PromptString', 'Select feature to mark:', ...
%                 'SelectionMode', 'single', ...
%                 'ListString', featureList, 'CancelString', 'Cancel');
%             if ~isempty(choice)
%                 markerMatrix(choice, actualIndex) = 1;
%                 markerMatrix(setdiff(1:k, choice), actualIndex) = 0;
%                 disp(['Marked frame ', num2str(actualIndex), ' as Feature ', num2str(choice)]);
%                 allFrames = getAllMarkedFrames();
%                 currentIndex = find(allFrames == actualIndex, 1);
%                 displayFrame();
%             end
%         else
%             disp('Click outside valid range.');
%         end
%     end
% end