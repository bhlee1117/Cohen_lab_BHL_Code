function smoothed_data = imaveragefilt(data, window_size)
    % MOVING_AVERAGE_2D Applies a 2D moving average filter to a matrix.
    %
    %   smoothed_data = moving_average_2d(data, window_size)
    %
    %   data        - Input 2D matrix (image or data grid)
    %   window_size - Scalar or [rows, cols] specifying the window size
    %
    %   smoothed_data - Output matrix with moving average applied

    kernel = ones(window_size) / window_size^2;
smoothed_data = imfilter(data, kernel, 'replicate');
end
