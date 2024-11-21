function BC = bimodality_coefficient(data)

% 2024.10.28 Byung Hun Lee, Cohen Lab,
% Calculating bimodality_coefficient
    % Calculate skewness
    skewness_val = skewness(data);
    
    % Calculate kurtosis
    kurtosis_val = kurtosis(data);
    
    % Number of samples
    n = length(data);
    
    % Bimodality Coefficient formula
    %BC = (skewness_val^2 + 1) / (kurtosis_val + 3*(n-1)^2 / ((n-2)*(n-3)));
    BC = (skewness_val^2 + 1) / kurtosis_val;

end