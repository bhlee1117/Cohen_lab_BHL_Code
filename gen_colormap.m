function customColormap = gen_colormap(keyColors,numSteps)
  if nargin < 2
        numSteps = 256;
    end

    % Initialize the colormap matrix
    customColormap = zeros(numSteps, 3);

    % Calculate the number of steps per gradient
    numGradients = size(keyColors, 1) - 1;
    stepsPerGradient = round(numSteps / numGradients);

    % Generate the colormap
    currentStep = 1;
    for i = 1:numGradients
        % Ensure the last gradient fills up the remaining steps
        if i == numGradients
            steps = numSteps - currentStep + 1;
        else
            steps = stepsPerGradient;
        end

        % Calculate the gradient between two key colors
        gradient = [linspace(keyColors(i, 1), keyColors(i + 1, 1), steps)', ...
                    linspace(keyColors(i, 2), keyColors(i + 1, 2), steps)', ...
                    linspace(keyColors(i, 3), keyColors(i + 1, 3), steps)'];

        % Fill the corresponding section of the colormap
        customColormap(currentStep:currentStep + steps - 1, :) = gradient;

        % Update the step index
        currentStep = currentStep + steps;
    end
end