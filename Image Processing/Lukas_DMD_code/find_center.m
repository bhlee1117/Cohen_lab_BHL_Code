% this function takes three inputs: an image, a label matrix that contains
% the objects labeled in the image, and the object whose center coordinates
% you want to find

function pts = find_center(im, L, num)

% make images out of the connected regions
[y, x] = find(L == num); 

obj = zeros(2048); 
for i = 1:numel(y)
    obj(y(i), x(i)) = 1; 
end

obj = obj(min(y):max(y), min(x):max(x));

% call the image matrix A and the convolution matrix B
top_left = [min(x), min(y)]; 
btm_rgt = [max(x), max(y)]; 
A = obj; 
ycoords_vect = uint16((top_left(2) - .5*(btm_rgt(2) - top_left(2))):(btm_rgt(2) + .5*(btm_rgt(2) - top_left(2)))); 
xcoords_vect = uint16((top_left(1) - .5*(btm_rgt(1) - top_left(1))):(btm_rgt(1) + .5*(btm_rgt(1) - top_left(1)))); 
B = im(ycoords_vect, xcoords_vect); 
A = double(A); 
B = double(B); 
C = conv2(B, A, 'same'); 

% find the indices of the maximum value of the convolution
[ymax, xmax] = find(C == max(C(:)));

% convert those indices to real camera coordinates
xcam = xmax + (min(xcoords_vect) - 1); 
ycam = ymax + (min(ycoords_vect) - 1); 

pts = [xcam, ycam]; 
end