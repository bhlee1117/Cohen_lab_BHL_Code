function myColormap = Aurora
turboColormap = turbo;
n = floor(size(turboColormap, 1)/6);  % Number of colors in 'turbo' colormap
whiteToTurbo = [linspace(1, turboColormap(20, 1), n)', ...
                linspace(1, turboColormap(20, 2), n)', ...
                linspace(1, turboColormap(20, 3), n)'];

% Combine the transition with the rest of the 'turbo' colormap
myColormap = [whiteToTurbo(1:end-1, :); turboColormap(21:end,:)];