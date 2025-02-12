function [binnedZ, centerX, centerY] = show3Dbinning(cellArray, var2, var3, plotType)
    % Function to bin, concatenate, and display 3D data from a cell array
    % Input:
    %   cellArray - N x M cell array containing data
    %   var2      - Number of bins or custom edges for X-axis
    %   var3      - Number of bins or custom edges for Y-axis
    %   plotType  - Type of plot: 'errorbar', 'surf', or 'image'

    % Check inputs
    if nargin < 4
        error('You must specify the cell array, binning parameters, and plot type.');
    end

     % Initialize arrays to store concatenated data
    allX = [];
    allY = [];
    allZ = [];

    % Process each cell in the array
    for i = 1:numel(cellArray)
        % Skip if cell is empty
        if isempty(cellArray{i})
            continue;
        end

        % Extract data from the cell
        data = cellArray{i};

        % First row and column represent the axes
        timeAxis = data(1, 2:end);
        xAxis = data(2:end, 1);
        zData = data(2:end, 2:end);

        % Flatten data for concatenation
        [X, Y] = meshgrid(timeAxis, xAxis);
        allX = [allX; X(:)];
        allY = [allY; Y(:)];
        allZ = [allZ; zData(:)];
    end

    
    % Determine if var2 and var3 are bin counts or edges
    if isscalar(var2)
        numBinsX = var2;
        edgesX = linspace(min(allX) - abs(min(allX)) * 0.1, max(allX), numBinsX + 1);
    else
        edgesX = var2;
    end

    if isscalar(var3)
        numBinsY = var3;
        edgesY = linspace(min(allY) - abs(min(allY)) * 0.1, max(allY), numBinsY + 1);
    else
        edgesY = var3;
    end

   

    % Compute bin centers
    centerX = (edgesX(1:end-1) + edgesX(2:end)) / 2;
    centerY = (edgesY(1:end-1) + edgesY(2:end)) / 2;

    % Use histcounts2 to bin data
    [~, ~, ~, binX, binY] = histcounts2(allX, allY, edgesX, edgesY);

    % Compute mean and std for Z values in each bin
    binnedZ = zeros(numel(edgesX) - 1, numel(edgesY) - 1);
    binnedStd = zeros(size(binnedZ));
    for i = 1:numel(edgesX) - 1
        for j = 1:numel(edgesY) - 1
            % Find points in the current bin
            binMask = (binX == i & binY == j);
            if any(binMask)
                binnedZ(i, j) = mean(allZ(binMask), 'omitnan');
                binnedStd(i, j) = std(allZ(binMask), 'omitnan');
            else
                binnedZ(i, j) = NaN;
                binnedStd(i, j) = NaN;
            end
        end
    end

    % Plot data
    switch plotType
        case 'surf'
            % Surface plot
            h = surf(centerX, centerY, binnedZ');
            h.FaceAlpha = 0.7; hold all
            %scatter3(allX,allY,allZ,20,[0.7 0.7 0.7],'filled','o')
            xlabel('X Axis');
            ylabel('Y Axis');
            zlabel('Z Values');
            colorbar;
        case 'errorbar'
            % 3D Error bar plot with colored points
            [meshX, meshY] = meshgrid(centerX, centerY);
            errorbar3(meshX(:), meshY(:), binnedZ(:), binnedStd(:));
            xlabel('X Axis');
            ylabel('Y Axis');
            zlabel('Z Values');
            colorbar;
        case 'image'
            % 2D image plot
            pcolor(centerX, centerY, binnedZ');
            shading interp;
            colorbar;
            xlabel('X Axis');
            ylabel('Y Axis');
    end
end

function errorbar3(X, Y, Z, E)
    % Helper function for 3D error bars with colored points
    hold on;
    scatter3(X, Y, Z, 50, Z, 'filled'); % Colored points
    colormap('jet'); % Set colormap for Z values
    for i = 1:numel(X)
        if ~isnan(Z(i)) && ~isnan(E(i))
            % Error bar lines
            plot3([X(i), X(i)], [Y(i), Y(i)], [Z(i)-E(i), Z(i)+E(i)], 'k');
        end
    end
    grid on;
    hold off;
end