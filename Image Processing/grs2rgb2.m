function res = grs2rgb2(img, map, cmin, cmax, nancolor)

%%Convert grayscale images to RGB using specified colormap.
%	IMG is the grayscale image. Must be specified as a name of the image 
%	including the directory, or the matrix.
%	MAP is the M-by-3 matrix of colors.
%   
%
%	RES = GRS2RGB(IMG) produces the RGB image RES from the grayscale image IMG 
%	using the colormap HOT with 64 colors.
%
%	RES = GRS2RGB(IMG,MAP) produces the RGB image RES from the grayscale image 
%	IMG using the colormap matrix MAP. MAP must contain 3 columns for Red, 
%	Green, and Blue components.  
%
%	RES = GRS2RGB(IMG, MAP, CMIN, CMAX) uses CMIN and CMAX to determine the input values mapped to the extrema of the
%   colormap.  Defaults are min(img(:)) and max(img(:)).
%
%	RES = GRS2RGB(IMG, MAP, CMIN, CMAX, NANCOLOR) uses CMIN and CMAX to determine the input values mapped to the extrema of the
%   colormap.  Defaults are min(img(:)) and max(img(:)). Pixels with values NaN will result in pixels in RES with RGB color NANCOLOR, which is a 1x3 array.
%   Default for NANCOLOR is MAP(1, :).
%
%
%	Example 1:
%	open 'image.tif';	
%	res = grs2rgb(image);
%
%	Example 2:
%	cmap = colormap(summer);
% 	res = grs2rgb('image.tif',cmap);
%
% 	See also COLORMAP, HOT
%
%	Written by 
%	Valeriy R. Korostyshevskiy, PhD
%	Georgetown University Medical Center
%	Washington, D.C.
%	December 2006
%
% 	vrk@georgetown.edu
%  
%   Modified to include cmin and cmax Adam Cohen, July 2011
%   Modified to include nancolor, JHH Sept 2013

% Check the arguments
if nargin<1
	error('grs2rgb:missingImage','Specify the name or the matrix of the image');
end;

if ~exist('map','var') || isempty(map)
	map = hot(64);
end;

if ~exist('cmax', 'var')
    cmin = min(img(:));
    cmax = max(img(:));
end;

[l,w] = size(map);

if w~=3
	error('grs2rgb:wrongColormap','Colormap matrix must contain 3 columns');
end;

if ischar(img)
	a = imread(img);
elseif isnumeric(img)
	a = img;
else
	error('grs2rgb:wrongImageFormat','Image format: must be name or matrix');
end;

% Calculate the indices of the colormap matrix
a = double(a);
% a(a==0) = 1; % Needed to produce nonzero index of the colormap matrix
ci = 2 + ceil((l-2)*(a-cmin)/(cmax - cmin));

ci = min(ci, l);
ci = max(ci, 1);

%%
% Colors in the new image
[il,iw] = size(a);
r = zeros(il,iw); 
g = zeros(il,iw);
b = zeros(il,iw);
r(:) = map(ci,1);
g(:) = map(ci,2);
b(:) = map(ci,3);

% New image
res = zeros(il,iw,3);
res(:,:,1) = r; 
res(:,:,2) = g; 
res(:,:,3) = b;
res2 = reshape(res, [il*iw, 3]);

%% modification by JHH to account for NaNs in image, set pixels to nancolor
%% beginning of chage by JHH
if exist('nancolor')
    for iC = 1:3
        res2(isnan(a), iC) = nancolor(iC);
    end;
    res = reshape(res2, [il, iw, 3]);
end