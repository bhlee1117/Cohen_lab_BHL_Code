function testInPlace
 
% Copyright 2007 The MathWorks, Inc.

%% Show in-place function calling

%% Make some data
n = 20*2^20;  %n = 35*2^20;
disp 'create initial var';drawnow;
x = randn(n,1);
pause(5)

%% Call regular function with same LHS
disp 'NOT in place but SAME var';drawnow;
x = myfunc(x); 
pause(5)

%% Call inplace function with same LHS
disp 'in place and SAME var';drawnow;
x = myfuncInPlace(x);
pause(5)

%% Call regular function with different LHS
disp 'NOT in place and DIFF var';drawnow;
y = myfunc(x); 
pause(5)

%% Call inplace function with same LHS
% change to new LHS and get an error
disp 'in place but DIFF var';drawnow;
y = myfuncInPlace(x);
pause(5)

disp 'finished'


