function customColormap = gen_colormap(keyColors, numSteps)
    if nargin < 2
        numSteps = 256;
    end

    % Define the positions of key colors in the interpolation space
    keyPositions = linspace(1, numSteps, size(keyColors, 1));

    % Create the interpolated colormap using linear interpolation
    customColormap = [interp1(keyPositions, keyColors(:, 1), 1:numSteps, 'linear')', ...
                      interp1(keyPositions, keyColors(:, 2), 1:numSteps, 'linear')', ...
                      interp1(keyPositions, keyColors(:, 3), 1:numSteps, 'linear')'];

end
