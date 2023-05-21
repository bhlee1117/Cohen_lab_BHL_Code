function [centerOffset, cor] = findCenter(img, guess, errorWin, radMax, radFlag)

% [center, cor] = findCenter(img, guess, errorWin, radMax)
% Finds the center of data that is assumed to be radially symmetric about
% an unknown point. Relies on a good guess (by the user) as to where the
% center should be, within +/- errorWin # pixels
% Created by JHH 03 August 2010
%
% Input:
% img = 2D array of data
% guess = x and y coordinates of best-guess center
% errorWin = size of half window to explore around best guess center
% radMax = maximum radius of feature expected to be radially symmetric
% radFlag = 1 if using all radii from 0 to radMax for the center-finding, or
% radFlag = 0 if just radMax for the center-finding

% Output:
% centerOffset = best fit center (column, row) of radially symmetric data
% cor = correlation function maximized to find center


nTheta = 100;
guessX = (guess(1)-errorWin):(guess(1)+errorWin);
guessY = (guess(2)-errorWin):(guess(2)+errorWin);
cor = zeros(errorWin);

scale = std(double(img(:)))*3.5; % to prevent reaching inf

if radFlag % use a range of radii
    radii = linspace(0, radMax, 20);
    radii = radii(2:end); % don't actually want r = 0
else
    radii = radMax; % just use largest radii
end

for iX = 1:2*errorWin+1
    for iY = 1:2*errorWin+1
        tmpRadInfo = radInterp(img, guessX(iX), guessY(iY), radii, nTheta);
        cor(iY, iX) = nansum(prod(tmpRadInfo.data/scale, 2));
    end
end

[bestFit, coeff, subCenter] = paraboloidFit(cor);
centerOffset = guess + subCenter - [errorWin + 1, errorWin + 1];

end
