function [masked_movie, masked_intens, mask] = mask_region2(movie,fig)
% function [masked_movie, masked_intens, mask] = mask_region2(movie, fig)
% Lets a user select a region of a movie to set to zero.
% Returns the movie with the indicated pixels set to zero, and the
% whole-field average intensity of the modified movie.
% fig = figure number to use;
% Adam Cohen 21 Dec. 2012

[ysize, xsize, nframes] = size(movie);
avgimg = mean(movie, 3);

[xcoords,ycoords]=meshgrid(1:xsize,1:ysize);
figure(fig)
imshow(avgimg,[],'InitialMagnification', 'fit')
title('Click around the region of interest')
[xpts, ypts] = getline(gca, 'closed');
selection = inpolygon(xcoords, ycoords, xpts, ypts);
close;

mask=zeros(size(avgimg));
mask(selection)=1;
masked_movie = movie.*repmat(mask,[1,1,nframes]);
masked_intens = squeeze(mean(mean(masked_movie)));