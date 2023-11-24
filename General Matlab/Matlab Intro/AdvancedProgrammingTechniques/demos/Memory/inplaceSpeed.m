function inplaceSpeed

% This function examines the speed up of using in-place operations
%
% Copyright 2007 The MathWorks, Inc.

% Create 1000-by-1000 random variable
a = rand(1000);
b = rand(1000);

iter = 1000;

% NOT in-place operation
tic
for id = 1:iter
  b = a + id;
end
toc

% in-place operation
tic
for id = 1:iter
  a = a + id;
end
toc