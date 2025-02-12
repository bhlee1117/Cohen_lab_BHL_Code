function B = get_Class2index(A)
    % A is the input binary matrix (3 x T)
    % The row of output index is true until there is another true frame in
    % the other class.
    [numRows, numCols] = size(A); % Get the size of the matrix
    B = zeros(numRows, numCols);  % Initialize the output matrix
    
    currentRow = 0; % Start with no active row
    for col = 1:numCols
        % Find all rows with true values in the current column
        rowIdx = find(A(:, col));
        
        % Check if there is more than one true value
        if length(rowIdx) > 1
            error('Column %d contains more than one true value.', col);
        end
        
        % If there's exactly one true value, update the current active row
        if ~isempty(rowIdx)
            currentRow = rowIdx;
        end
        
        % Update the output matrix based on the current active row
        if currentRow > 0
            B(currentRow, col) = 1;
        end
    end
end