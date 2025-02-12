function orderVector = get_BurstOrder(binaryMatrix, maxInterval)
    % binaryMatrix: Input binary matrix (1D or 2D array)
    % maxInterval: Maximum interval (in frames) between true points to be in the same train
    % orderVector: Output matrix with order numbers representing the position within each train

    % Ensure binaryMatrix is a row vector for processing
    isRowVector = isrow(binaryMatrix);
    if isRowVector
        binaryMatrix = binaryMatrix(:);
    end

    % Find indices of true values
    trueIndices = find(binaryMatrix);
    
    % Handle edge cases: no true values
    if isempty(trueIndices)
        orderVector = zeros(size(binaryMatrix));
        if isRowVector
            orderVector = orderVector';
        end
        return;
    end

    % Identify start of new trains
    diffIndices = [Inf; diff(trueIndices)];
    newTrainStart = diffIndices > maxInterval;

    % Assign a train number to each true index
    trainNumbers = cumsum(newTrainStart);

    % Create the order within each train
    orderInTrain = zeros(size(trueIndices));
    for train = unique(trainNumbers)'
        trainIndices = trainNumbers == train;
        orderInTrain(trainIndices) = 1:sum(trainIndices);
    end

    % Create the output vector
    orderVector = zeros(size(binaryMatrix));
    orderVector(trueIndices) = orderInTrain;

    % Ensure the output has the same orientation as the input
    if isRowVector
        orderVector = orderVector';
    end
end