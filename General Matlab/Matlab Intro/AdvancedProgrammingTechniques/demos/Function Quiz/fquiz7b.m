function countfcn = fquiz7b(initvalue)

% Copyright 2007 The MathWorks, Inc.

%% We can export custom functions with state

currentCount = initvalue; % Initial value
countfcn = @() getCounter(currentCount);   % Return handle to getCounter

end

function count = getCounter(currentCount)
% This function increments the variable 'currentCount', when it
% gets called (using its function handle) .
currentCount = currentCount + 1;
count = currentCount;
end
