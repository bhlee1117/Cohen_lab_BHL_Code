function countfcn = fquiz7(initvalue)
% Function Quiz 7
% Copyright 2007 The MathWorks, Inc.

%% We can export custom functions with state

currentCount = initvalue  ; % Initial value
countfcn     = @getCounter; % Return handle to getCounter

    function count = getCounter()
        % This function increments the variable 'currentCount', when it
        % gets called (using its function handle) .
        currentCount = currentCount + 1;
        count        = currentCount;
    end
end

