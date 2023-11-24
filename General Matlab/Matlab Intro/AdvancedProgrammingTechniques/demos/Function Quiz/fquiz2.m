
% Copyright 2007 The MathWorks, Inc.

%% What is fh?
fh = @sin

%%
% Call like any other function
fh(pi/4)

fh(pi/4) - sin(pi/4)

%%
% Pass as input to other functions ("fun-funs")
ezplot(fh)
