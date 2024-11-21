function binaryVector=ind2binary(indices,n,dat)
% Convert indicies to a binary vector which defined indicies are true.
% 2024.11.10 Byung Hun Lee, Cohen Lab



% Determine the length of the binary vector
% This should be at least as large as the maximum index
if nargin<2
n = max(indices);
dat=ones(1,length(indices));
end

% Initialize a logical vector of length n with all false (0)
binaryVector = zeros(1, n);

% Set the specified indices to true (1)
binaryVector(indices) = dat;
end