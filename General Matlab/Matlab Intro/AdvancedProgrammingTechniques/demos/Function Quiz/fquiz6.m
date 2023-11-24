function fh = fquiz6(a,b)
% Function Quiz 6
% Copyright 2007 The MathWorks, Inc.

%% We can export custom functions
fh = @makeline;

    function y = makeline(x)
        y = a*x + b;
    end

end

