function colorSpec = generateColorSpec(n, directoryLetter)
% generateColorSpec(n, directoryLetter);
% J. Hou 18 Dec 2009
% Generates a 3*n array, representing the RGB values of a color spectrum, based off the jet color scheme.
% Need to input directoryLetter to know where to find jet color scheme
% In the future, will expand to allow for multiple color schemes, but for
% now, jet is the only one.

load([directoryLetter ':\Lab\Computer Code\General Matlab\colorsJet.mat']);
colorSpec = colorsJet(round(linspace(1, size(colorsJet, 1), n)), :);