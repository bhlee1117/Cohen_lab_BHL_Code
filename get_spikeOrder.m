function orderVector=get_spikeOrder(A,K)
 % A,BinaryCondition: Binary trace indicating intervals
    % K,spike: Binary trace to be converted to order
    % orderVector: Output vector with order of true values in K within the 'on' intervals of A

    % Validate inputs
    assert(isvector(A) && isvector(K), 'Inputs A and K must be vectors.');
    assert(length(A) == length(K), 'Vectors A and K must be of the same length.');

    % Initialize the output vector
    orderVector = zeros(size(K));
    
    % Find the "on" intervals in A
    starts = find(diff([0, A]) == 1); % Start of an "on" interval
    ends = find(diff([A, 0]) == -1);  % End of an "on" interval

    % Loop through each "on" interval
    for i = 1:length(starts)
        % Get indices for the current interval
        intervalIndices = starts(i):ends(i);
        
        % Find true values in K within this interval
        trueIndices = find(K(intervalIndices));
        
        % Assign order numbers to the true values
        orderVector(intervalIndices(trueIndices)) = 1:length(trueIndices);
    end
end