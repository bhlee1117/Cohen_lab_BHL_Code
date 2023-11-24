function colorSpec = generateColorSpecLocal(n, mapString)
% generateColorSpecLocal(n);
% J. Hou 18 Dec 2009
% J. Hou 19 Feb 2010 modified
% J. Hou 19 March 2013 modified to allow for different color schemes (;
% J. Hou 07 July 2013 modifed such that for small n (2,3,4), give predefined
% colors
% Generates a 3*n array, representing the RGB values of a color spectrum, based off the jet color scheme.
% In the future, will expand to allow for multiple color schemes, but for
% now, jet is the only one.
% Difference between generateColorSpecLocal and generateColorSpec is that
% this one does not require access to the server (well, except to call this
% function, but does not depend on being able to find the colorJet.mat file)
% To use:
% colors = generateColorSpecLocal(nLines);
% figure(1); clf
% set(gca, 'ColorOrder', colors,  'NextPlot', 'replacechildren');
% plot(data);

if nargin == 1
    mapString = 'jet';
end
figure;
colors = colormap(mapString);
close;

switch n
    case 2
        colorSpec = [0, 0, 0.9; 0.9, 0, 0];
    case 3
        colorSpec = [0, 0, 0.9; 0, 0.9, 0; 0.9, 0, 0];
    case 4
        colorSpec = [0, 0, 0.9; 0, 0.9, 0; 0.9, 0.5, 0; 0.9, 0, 0];
    otherwise
        colorSpec = colors(round(linspace(1, size(colors, 1), n)), :);
end