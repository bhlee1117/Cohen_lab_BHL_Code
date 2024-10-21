function [transients_final] = detect_transientSimple(traces_hi, thres)

thres_on = thres(1);
thres_off = thres(2);
[num_neurons, num_points] = size(traces_hi);

transients_final = false(num_neurons, num_points); % Initialize as logical for better performance

% Precompute the threshold-based masks
above_thres_on = traces_hi > thres_on;
above_thres_off = traces_hi > thres_off;

% Moving average with a window size of 3
smoothed_traces = movmean(traces_hi, 3, 2);
above_smoothed_thres_on = smoothed_traces > thres_on;

for i = 1:num_neurons
    % Find initial transients and background transients
    [transients, n] = bwlabel(above_thres_on(i, :));
    [transients_back, ~] = bwlabel(above_thres_off(i, :));
    
    % Update the final transients based on the smoothed trace
    transients_final(i, :) = above_smoothed_thres_on(i, :);
    
    for t = 1:n
        tr_ind = find(transients == t);
        if isempty(tr_ind)
            continue;
        end
        t2 = find(transients_back == transients_back(tr_ind(end)));
        
        % Update transients_final within the identified range
        transients_final(i, t2(1):t2(end)) = true;
    end
end

% Convert logical array to double
transients_final = double(transients_final);

end