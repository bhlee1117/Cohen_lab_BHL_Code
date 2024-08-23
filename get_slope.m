function slopes=get_slope(data,window_size)
half_window = floor(window_size / 2);

% Initialize the slope vector
slopes = nan(size(data));

% Calculate the slope for each window position
for i = 1:length(data)
    % Determine the start and end indices of the window
    start_idx = max(1, i - half_window);
    end_idx = min(length(data), i + half_window);
    
    % Extract the data points within the window
    window_data = data(start_idx:end_idx);
    
    % Generate the x-values for the window
    x = start_idx:end_idx;
    
    % Fit a linear polynomial to the data points
    if length(window_data) >= 2 % Ensure there are at least 2 points to fit
        p = polyfit(x, window_data, 1);
        slopes(i) = p(1); % The slope is the first coefficient
    end
end
