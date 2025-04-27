function [length, width, rotatedBox] = measureRectangle(mask, visualize)
    % Ensure the mask is binary
    if size(mask, 3) > 1
        mask = rgb2gray(mask);
    end
    mask = imbinarize(mask);

    % Find object boundaries
    [B, ~] = bwboundaries(mask, 'noholes');

    if isempty(B)
        error('No objects detected in the mask.');
    end

    % Get largest object boundary
    largestBoundary = B{1};

    % Compute minimum area bounding box using PCA
    covMatrix = cov(largestBoundary);
    [~, S, V] = svd(covMatrix); % Singular value decomposition

    % Extract principal axes
    majorAxis = V(:,1);
    minorAxis = V(:,2);

    % Project boundary points onto principal axes
    projMajor = largestBoundary * majorAxis;
    projMinor = largestBoundary * minorAxis;

    % Compute length and width (range along axes)
    length = max(projMajor) - min(projMajor);
    width = max(projMinor) - min(projMinor);

    % Ensure length is the longer side
    if width > length
        temp = length;
        length = width;
        width = temp;
    end

    % Display results
    fprintf('Measured Length: %.2f pixels\n', length);
    fprintf('Measured Width: %.2f pixels\n', width);

    % Draw rotated bounding box
        cornerPoints = [min(projMajor), min(projMinor); 
                        max(projMajor), min(projMinor);
                        max(projMajor), max(projMinor);
                        min(projMajor), max(projMinor)];

        % Rotate back to image coordinates
        rotatedBox = cornerPoints * [majorAxis'; minorAxis'];

        % Connect points to form the rectangle
        rotatedBox = [rotatedBox; rotatedBox(1,:)]; % Close the rectangle
        
    % Visualization
    if visualize
        figure;
        imshow(mask);
        hold on;
        plot(largestBoundary(:,2), largestBoundary(:,1), 'g', 'LineWidth', 2); % Object boundary

        
        plot(rotatedBox(:,2), rotatedBox(:,1), 'r', 'LineWidth', 2);

        title('Rotated Bounding Box');
        hold off;
    end
end
