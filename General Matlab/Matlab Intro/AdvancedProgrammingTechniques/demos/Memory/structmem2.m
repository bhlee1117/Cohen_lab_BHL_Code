%% Load an image

% Copyright 2007 The MathWorks, Inc.

X = imread('onion.png');
imshow(X)
whos  % note that X is a 3-dimensional array

%% Store one plane per field
im1.red   = X(:,:,1);                % each plane is m x n
im1.green = X(:,:,2); 
im1.blue  = X(:,:,3);     

%% Store one pixel per field
[nr,nc,np] = size(X); 
im2(nr,nc).pixel = [0 0 0];
for row = 1:nr
    for col = 1:nc
        im2(row,col).pixel = reshape(X(row,col,:),1,3);
    end
end
clear nr nc np row col

%% Which one is bigger?
whos im*