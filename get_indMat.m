function binary_matrix=get_indMat(vector)
% This function convert a integer vector (1xk) with N element to (kxN)
% binary matrix
% 2024/09/20 Byung Hun Lee
numbers=unique(vector);
N=length(numbers);

    % Get the length of the vector
    k = length(vector);
    
    % Initialize a k x N binary matrix of zeros
    binary_matrix = zeros(k, N);
    
    % Create a column index array from 1 to k
    row_indices = 1:k;
    
    % Use linear indexing to assign 1s at the appropriate positions
    binary_matrix(sub2ind([k, N], row_indices, vector')) = 1;

end