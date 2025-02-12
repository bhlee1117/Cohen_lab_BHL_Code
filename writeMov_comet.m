function writeMov_comet(vector1, vector2, windowLength, filename , labels, skipFrm)
    % Validate input
    if length(vector1) ~= length(vector2)
        error('The input vectors must have the same length.');
    end

    if nargin<5
labels={'Vector 1','Vector 2'};
    end

        if nargin<6
skipFrm=1;
    end

    numPoints = length(vector1);
    colormap('turbo'); % Set colormap
    cmap = colormap; % Get the colormap matrix
    numColors = size(cmap, 1);

    % Prepare video writer
    videoWriter = VideoWriter(filename, 'MPEG-4'); % Save as MP4 file
    videoWriter.FrameRate = 25; % Adjust frame rate if needed
    open(videoWriter);

    % Prepare figure
    fig=figure;
    set(gcf, 'Position', [100, 100, 1000, 800]); % Set figure window size [left, bottom, width, height]
    hold on;
    grid on;
    %axis equal tight;
    xlim([min(vector1), max(vector1)]);
    ylim([min(vector2), max(vector2)]);
    xlabel(labels{1});
    ylabel(labels{2});
    title('Moving Window Animation with Chronological Color');
    colorbar;
    colormap(turbo);
    caxis([0, windowLength]); % Set colorbar range
    cb = colorbar;
    cb.Ticks = linspace(0, windowLength, 5); % Set ticks from 0 to window size
    ylabel(cb, 'Time (ms)');

    set(gca, 'FontSize', 20); % Set font size for axes
    box on;

    % Animation
    for i = 1:skipFrm:(numPoints - windowLength + 1)
        % Define the window
        windowStart = i;
        windowEnd = i + windowLength - 1;

        % Map the window indices to colormap
        colorIndices = linspace(1, numColors, windowLength);
        colors = cmap(round(colorIndices), :);

        l=plot([vector1(windowStart:windowEnd-1); vector1(windowStart+1:windowEnd)],[vector2(windowStart:windowEnd-1); vector2(windowStart+1:windowEnd)],'LineWidth',2);
        arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(turbo(windowLength-1),2));
        
        % Capture the frame
        frame = getframe(gcf);
        writeVideo(videoWriter, frame);
        title(['Leading point is frame # ' num2str(i)])

        % Clear the figure for the next frame (comment this if you want to keep all frames)
        cla;
    end

    % Finalize
    hold off;
    close(videoWriter);
    close(fig);

    disp(['Animation saved as ', filename]);
end