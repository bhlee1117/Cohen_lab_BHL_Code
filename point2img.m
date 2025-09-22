function binaryImage = point2img(points, radius, imageSize)
    % point2img: Generate a binary image with pixels near (x, y) points set to 1
    % Inputs:
    %   points: Nx2 matrix of [x, y] coordinates
    %   radius: Scalar, radius (in pixels) around points to set to 1
    %   imageSize: Scalar, size of square image (e.g., 50 for 50x50)
    % Output:
    %   binaryImage: Logical matrix of size [imageSize, imageSize]

    % Pre-allocate binary image as logical matrix
    binaryImage = false(imageSize(1), imageSize(2));

    % Round coordinates and check bounds
    x = round(points(:, 2));
    y = round(points(:, 1));
    valid = x >= 1 & x <= imageSize(1) & y >= 1 & y <= imageSize(2);
    x = x(valid);
    y = y(valid);

    % If no valid points, return empty image
    if isempty(x)
        return;
    end

    % Create grid for distance calculation
    [xGrid, yGrid] = ndgrid(1:imageSize(1), 1:imageSize(2));

    % Vectorized assignment of pixels within radius
    for i = 1:length(x)
        if i==1
        distance = sqrt((xGrid - x(i)).^2 + (yGrid - y(i)).^2);
        binaryImage = binaryImage | (distance <= radius*2 & distance > 0);    
        else
        distance = sqrt((xGrid - x(i)).^2 + (yGrid - y(i)).^2);
        binaryImage = binaryImage | (distance <= radius & distance > 0);
        end
    end
end