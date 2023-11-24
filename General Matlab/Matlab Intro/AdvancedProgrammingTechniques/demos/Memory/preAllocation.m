%% Effect of Pre-Allocation


% Copyright 2007 The MathWorks, Inc.

clear all;clc

ln = 30000;

%% No Pre-Allocation
tic;
for id = 1:ln
  a(id) = rand();
end
t1=toc;
fprintf('Elapsed time is %g seconds.\n', t1);

%% Pre-Allocation
tic;
b = nan(1, ln);
%b = a;
for id = 1:ln
  b(id) = rand();
end
t2=toc;
fprintf('Elapsed time is %g seconds.\n', t2);

fprintf('Speedup: %d times faster\n', round(t1/t2));