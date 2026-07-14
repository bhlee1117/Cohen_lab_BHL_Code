function mask = make_dendrite_line_mask(imSize, centerXY, angleDeg, length_um, width_um, pix_um)
% Make a rectangular dendrite-like line mask.
%
% imSize    : [X Y]
% centerXY  : [row col] center position in pixels
% angleDeg  : dendrite angle in degrees
% length_um : dendrite length in um
% width_um  : dendrite width in um
% pix_um    : pixel size in um / pixel
%
% Output:
% mask      : X x Y binary mask

X = imSize(1);
Y = imSize(2);

[rowGrid, colGrid] = ndgrid(1:X, 1:Y);

row0 = centerXY(1);
col0 = centerXY(2);

theta = deg2rad(angleDeg);

% Coordinates relative to center
dr = rowGrid - row0;
dc = colGrid - col0;

% Coordinate along dendrite axis
s =  dc * cos(theta) + dr * sin(theta);

% Coordinate perpendicular to dendrite axis
r = -dc * sin(theta) + dr * cos(theta);

halfLength_pix = (length_um / pix_um) / 2;
halfWidth_pix  = (width_um  / pix_um) / 2;

mask = abs(s) <= halfLength_pix & abs(r) <= halfWidth_pix;
mask = double(mask);

end