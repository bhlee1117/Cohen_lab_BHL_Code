function shrunkBinary = shrinkToMax(binaryVector, A)
    % binaryVector: Input binary vector (1D array of 0s and 1s)
    % A: Vector of the same size, containing values to evaluate
    % shrunkBinary: Output binary vector with consecutive trues reduced to one

    % Validate inputs
    assert(isvector(binaryVector) && isvector(A), 'Inputs must be vectors.');
    assert(length(binaryVector) == length(A), 'Vectors must have the same length.');
    
    % Initialize the output binary vector
    shrunkBinary = zeros(size(binaryVector));
    
    % Find indices of consecutive true regions
    starts = find(diff([0, binaryVector]) == 1); % Start of a true segment
    ends = find(diff([binaryVector, 0]) == -1);  % End of a true segment
    
    % Loop through each segment
    for i = 1:length(starts)
        segment = starts(i):ends(i); % Indices of the segment
        [~, maxIdx] = max(A(segment)); % Find index of max value within the segment
        shrunkBinary(segment(maxIdx)) = 1; % Set the max index to true
    end
    disp('Consecutive vector squeezed')
end