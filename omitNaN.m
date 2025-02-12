function [cleanTrace, validIdx] = omitNaN(trace)
    % omitNaN - Remove NaN values and return the cleaned trace
    %
    % Inputs:
    %   trace - 1D signal containing NaN values
    %
    % Outputs:
    %   cleanTrace - Trace with NaN values removed
    %   validIdx   - Logical index of valid (non-NaN) elements

    validIdx = ~isnan(trace); % Logical index of valid values
    cleanTrace = trace(validIdx); % Keep only non-NaN values
end