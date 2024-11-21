function [cFitted fit_curves]=expfit_wBd(x,y,xv,initialGuess,lowerBounds,upperBounds)

% 2024.10.28 Byung Hun Lee, exponential fitting
% Define the exponential model function: y = a*exp(b*x)
modelFun = @(c, x) c(1) * exp(-x/c(2))+c(3);

% Initial guess for parameters [a, b]
%initialGuess = [1, 0.5];  % Modify this as per your estimation for a and b

% Define lower and upper bounds for [a, b]
%lowerBounds = [0, -Inf];   % Lower bound for 'a' and 'b' (e.g., a >= 0)
%upperBounds = [Inf, Inf,Inf];    % Upper bound for 'a' and 'b' (e.g., b <= 1)

% Perform fitting with lsqcurvefit
options = optimoptions('lsqcurvefit', 'Display', 'off');
cFitted = lsqcurvefit(modelFun, initialGuess, x, y, lowerBounds, upperBounds, options);
fit_curves= modelFun(cFitted,xv);
end