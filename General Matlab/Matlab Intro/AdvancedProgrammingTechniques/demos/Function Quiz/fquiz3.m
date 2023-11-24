
% Copyright 2007 The MathWorks, Inc.

%% What type of function is fh a handle to?
a = 2;
b = 3;
fh = @(x) a*x + b;

%% 
% Call like any other function
fh(1:5)

%%
% Pass as input to other functions ("fun-funs")
ezplot(fh)

%%
% What if we lose or change a&b?
a = 100; b = 200;
fh(1:5)

%%
% Make a similar function
fh2 = @(x) a*x + b;
hold on;h2 = ezplot(fh2);set(h2, 'Color', 'r');
shg;

%%
% How do we figure out what's there?
f = functions(fh);
f.workspace{:}

%% (inline functions, optional)
% Create a function f,
a = 2;
b = 3;
f = inline('a*x + b','x','a','b');

%% 
% Evaluate the function
f(1:5,a,b)

%%
% What if we lose or change a&b?
a = 100; b = 200;
f(1:5,a,b)
