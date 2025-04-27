function R2 = get_Rsquare(Y, Y_fit)
% computeR2: Calculates coefficient of determination (R^2)
%
%   R2 = computeR2(Y, Y_fit) returns the R-squared value between the actual
%   data Y and the model-predicted values Y_fit.
%
%   Inputs:
%       Y      - Actual observed values (vector)
%       Y_fit  - Fitted or predicted values (same size as Y)
%
%   Output:
%       R2     - Coefficient of determination (scalar)

    if ~isvector(Y) || ~isvector(Y_fit)
        error('Inputs Y and Y_fit must be vectors.');
    end
    
    if length(Y) ~= length(Y_fit)
        error('Y and Y_fit must be the same length.');
    end

    Y = Y(:);       % Ensure column vector
    Y_fit = Y_fit(:);
    
    SS_res = sum((Y - Y_fit).^2);
    SS_tot = sum((Y - mean(Y)).^2);

    R2 = 1 - SS_res / SS_tot;
end
