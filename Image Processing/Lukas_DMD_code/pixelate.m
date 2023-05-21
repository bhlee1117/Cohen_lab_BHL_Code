% pixelate the image

function roi_points=pixelate(movie_in, pform)

refimg = mean(movie_in, 3);
[ysize, xsize] = size(refimg(:,:,1));

[y_px, x_px] = size(pform);

