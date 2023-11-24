%% Vectorize Example
% This example shows the speed up of vectorized code.
% However, also note that this vectorized code requires more memory to
% create temporary data (use OS Task Manager to view memory use)

% Copyright 2007 The MathWorks, Inc.

clear
a = rand(4000);
b = zeros(4000);

%% Vectorized
tic;
b = a.^a;
toc;

%% De-Vectorized

tic;
for id1 = 1:size(a,1)
  for id2 = 1:size(a,2)
    b(id1, id2) = a(id1, id2)^a(id1, id2);
  end
end
toc;






%% De-Vectorized (row-wise)

tic;
for id1 = 1:size(a,1)
  b(id1, :) = a(id1, :).^a(id1, :);
end
toc;

%% De-Vectorized (column-wise)

tic;
for id2 = 1:size(a,2)
  b(:, id2) = a(:, id2).^a(:, id2);
end
toc;

