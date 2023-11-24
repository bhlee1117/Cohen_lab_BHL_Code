function x = myfuncInPlace(x)

% Copyright 2007 The MathWorks, Inc.

% This function can operate on data in-place if called correctly.
x = sin(2*x.^2+3*x+4);
