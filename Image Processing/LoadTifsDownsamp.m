function out = LoadTifsDownsamp(flist, dx, dy, dn);
% function out = LoadTifsDownsamp(flist, dx, dy, dn);
% Reads Tifs from a list and returns a downsampled 3D image stack.
% flist is a list of files generated by a dir command.  The files are
% assumed to be in order and all of the same size.
% output pixels are *averaged* over dx, dy, but only *sampled* at dn
%
% Designed to read Tif stacks from the Hamamatsu camera
% AEC 5 Aug. 2013

tmp = imread(flist(1).name);
[ysize, xsize] = size(tmp);
nframes = length(flist);

out = zeros(length(1:dy:ysize), length(1:dx:xsize), length(1:dn:nframes));
[dysize, dxsize, dnframes] = size(out);

c = 1;
for j = 1:dn:nframes;
    tmp = double(imread(flist(j).name));
    out(:,:,c) = imresize(tmp, [dysize, dxsize], 'box');
    c = c + 1
end;


