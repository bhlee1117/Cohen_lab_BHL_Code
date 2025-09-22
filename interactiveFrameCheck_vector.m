function modifiedMarkerMatrix = interactiveFrameCheck_vector(dataVector, markerMatrix, windowSize)
% interactiveFrameCheck_vector - Interactive viewer for 1xT vector data with N-class markers.
%
% Inputs:
%   dataVector   : 1 x T vector of data
%   markerMatrix : N x T binary matrix (N classes)
%   windowSize   : half-width of time window to display around each mark

assert(isvector(dataVector), 'dataVector must be a vector');
dataVector = dataVector(:)';
[N, T] = size(markerMatrix);
assert(length(dataVector) == T, 'markerMatrix must match dataVector length');

colors = lines(N);
allFrames = getAllMarkedFrames();
currentIndex = 1;

hFig = figure('Name', 'Interactive Frame Check (Vector)', ...
              'KeyPressFcn', @keyPressHandler);
plotFrame();
uiwait(hFig);
modifiedMarkerMatrix = markerMatrix;

    function allFrames = getAllMarkedFrames()
        frames = [];
        for i = 1:N
            frames = [frames, find(markerMatrix(i, :) == 1)];
        end
        allFrames = unique(frames);
    end

    function plotFrame()
        if isempty(allFrames), return; end
        frame = allFrames(currentIndex);
        startIdx = max(1, frame - windowSize);
        endIdx = min(T, frame + windowSize);

        clf;
        plot(startIdx:endIdx, dataVector(startIdx:endIdx), 'w');
        hold on;

        for cls = 1:N
            classMarkers = find(markerMatrix(cls, startIdx:endIdx)) + startIdx - 1;
            for i = 1:length(classMarkers)
                f = classMarkers(i);
                if f == frame
                    scatter(f, dataVector(f), 100, colors(cls,:), 'filled', 'v');
                else
                    scatter(f, dataVector(f), 60, colors(cls,:), 'v');
                end
            end
        end

        title(sprintf('Frame %d', frame));
        xlabel('Time'); ylabel('Signal');
        xlim([startIdx endIdx]);
    end

    function keyPressHandler(~, event)
        frame = allFrames(currentIndex);
        switch event.Key
            case 'rightarrow'
                currentIndex = min(currentIndex + 1, numel(allFrames));
                plotFrame();
            case 'leftarrow'
                currentIndex = max(currentIndex - 1, 1);
                plotFrame();
            case 'uparrow'
                markerMatrix(:, frame) = 0;
                disp(['Unmarked frame ', num2str(frame)]);
                allFrames = getAllMarkedFrames();
                currentIndex = min(currentIndex, numel(allFrames));
                if isempty(allFrames), close(hFig); return; end
                plotFrame();
            case 'q'
                disp('Exiting viewer.');
                close(hFig);
            otherwise
                numKey = str2double(event.Key);
                if ~isnan(numKey) && numKey >= 1 && numKey <= N
                    markerMatrix(:, frame) = 0;
                    markerMatrix(numKey, frame) = 1;
                    disp(['Marked frame ', num2str(frame), ' as Class ', num2str(numKey)]);
                    allFrames = getAllMarkedFrames();
                    plotFrame();
                end
        end
    end
end
