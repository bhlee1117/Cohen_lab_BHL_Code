function h = imshow3(varargin)
img = varargin{1};
if numel(varargin) == 1
    factor = [-3,3];
else factor = varargin{2};
    varargin(2) = [];
end
disprange = mean(img(:))+std(img(:))*factor;
varargin = [varargin, 'InitialMagnification','fit','DisplayRange',disprange];
hh = imshow(varargin{:}); 
titlestr = inputname(1);
title(titlestr);

if (nargout > 0)
% Only return handle if caller requested it.
h = hh;
end