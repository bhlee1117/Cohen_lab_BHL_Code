function img = imblend(img1,img2)
%imblend    blend two images into one orange/blue map. 
% 
%   2016 Vicente Parot
%   Cohen Lab - Harvard University

img = cat(3,img1,img1/2+img2/2,img2);