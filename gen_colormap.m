function customColormap = gen_colormap(keyColors)
numSteps = 256;

% Define key colors: start with blue, transition through white, then to red and orange

% Initialize the colormap matrix
customColormap = [];

% Generate gradients between key colors
for i = 1:size(keyColors, 1)-1
    % Calculate the gradient between two key colors
    gradient = [linspace(keyColors(i, 1), keyColors(i+1, 1), numSteps/(size(keyColors, 1)-1))' ...
                linspace(keyColors(i, 2), keyColors(i+1, 2), numSteps/(size(keyColors, 1)-1))' ...
                linspace(keyColors(i, 3), keyColors(i+1, 3), numSteps/(size(keyColors, 1)-1))'];
    
    % Append the gradient to the colormap, excluding the last color to avoid repetition
    customColormap = [customColormap; gradient(1:end-1, :)];
end

% Append the last key color
customColormap = [customColormap; keyColors(end, :)];
end