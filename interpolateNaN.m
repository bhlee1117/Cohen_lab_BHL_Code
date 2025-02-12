function interpolatedMatrix = interpolateNaN(matrix)
    % interpolateNaN - Interpolates NaN values along the T direction (columns)
    %
    % Syntax: interpolatedMatrix = interpolateNaN(matrix)
    %
    % Inputs:
    %   matrix - N x T matrix containing NaN values
    %
    % Output:
    %   interpolatedMatrix - N x T matrix with NaN values interpolated along columns
    
    % Get the size of the input matrix
    [N, T] = size(matrix);
    
    % Initialize the output matrix
    interpolatedMatrix = matrix;

    % Loop through each row
    for i = 1:N
        row = matrix(i, :); % Extract the current row
        nanIdx = isnan(row); % Identify NaN indices
        if any(nanIdx)
            validIdx = ~nanIdx; % Identify valid indices
            % Interpolate NaN values using linear interpolation
            interpolatedMatrix(i, nanIdx) = interp1(find(validIdx), row(validIdx), find(nanIdx), 'linear', 'extrap');
        end
    end
end