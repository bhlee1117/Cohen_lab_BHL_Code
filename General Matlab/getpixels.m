function [pixelsout,npixels] = getpixels(imin)
% FUNCTION [PIXELSOUT,NPIXELS] = GETPIXELS(IMIN)
% This function accepts an image as an input and presents a new window with
% a cursor.  The function will return a mask the same size as the input
% image with 1s in the selection, and 0s elsewhere.  

[ysize,xsize,nframes] = size(imin);

if nframes > 4
    avgimg = mean(imin,3);
else
    avgimg = imin;
end

hh = figure;
imshow(avgimg,[],'InitialMagnification','Fit')
[xpts, ypts] = getline(gca,'closed');
close(hh);

xcoords = ones(ysize,1)*(1:xsize);
ycoords = (1:ysize)'*ones(1,xsize);
pixelsout = inpolygon(xcoords, ycoords, xpts, ypts);
npixels = length(find(pixelsout));
