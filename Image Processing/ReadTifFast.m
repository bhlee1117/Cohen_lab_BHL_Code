function out = ReadTifFast(fname)
% function out = ReadTifFast(fname)
% reads in a Tif file exported from the Andor Solis software.  Can read
% long movies much faster than the frame-by-frame method
% The specification for the Tif file structure is in:
%
% partners.adobe.com/public/developer/en/tiff/TIFF6.pdf 
%
% Image must be > 3 frames.
%
% AEC 17 Jan. 2011.

% Read in the whole file into a UINT8 vector

fid = fopen(fname);
tmp = uint8(fread(fid, 'uint8'));
fclose(fid);
L = length(tmp);

% Get the location of the directory for the first image
ifd1 = typecast(tmp(5:8), 'uint32')+1;
% Get the number of directory entries for this image
nentries = uint32(typecast(tmp(ifd1:ifd1+1), 'uint16'));

%% This section looks up the tags for the first image.  Assuming they
%% always come in the same order, it is not necessary.
tags = zeros(nentries,1);
for j = (1:nentries);
    tags(j) = typecast(tmp(ifd1+2+12*(j-1):ifd1+3+12*(j-1)), 'uint16');
end;
xsizeloc = find(tags == 256)-1;
ysizeloc = find(tags == 257)-1;
img1OffsetLocLoc = find(tags == 273)-1;


%%
% Get the locations of the relevant tags
xsize = double(typecast(tmp(ifd1 + 2 + xsizeloc*12 + 8:ifd1 + 2 + xsizeloc*12 + 9), 'uint16'));
ysize = double(typecast(tmp(ifd1 + 2 + ysizeloc*12 + 8:ifd1 + 2 + ysizeloc*12 + 9), 'uint16'));
% The location of the vector whose first element is the image 1 offset is
% here:
img1OffsetLoc = typecast(tmp(ifd1 + 2 + img1OffsetLocLoc*12 + 8:ifd1 + 2 + img1OffsetLocLoc*12 + 11), 'uint32')+1;
% Find the offset of the first image.
img1Offset = typecast(tmp(img1OffsetLoc:img1OffsetLoc+3),'uint32')+1;

% Get the location of the directory for the second image
ifd2 = typecast(tmp(ifd1 + 2 + nentries*12:ifd1 + 2 + nentries*12 + 3), 'uint32')+1;
% directories 2:end have a different number of tags
nentries2 = uint32(typecast(tmp(ifd2:ifd2 + 1), 'uint16'));

%% Locate the tags in images 2:end.  Turns out to have the same order as in
%% image 1.
% tags2 = zeros(nentries2,1);
% for j = (1:nentries2);
%     tags2(j) = typecast(tmp(ifd2+2+12*(j-1):ifd2+3+12*(j-1)), 'uint16');
% end;

%%
% The location of the vector whose first element is the image 2 offset is
% here:
img2OffsetLoc = typecast(tmp(ifd2 + 2 + 12*7 + 8:ifd2 + 2 + 12*7 + 11), 'uint32')+1;
% Find the offset of the first image.
img2Offset = typecast(tmp(img2OffsetLoc:img2OffsetLoc+3),'uint32')+1;

%% Repeat for the third image.
ifd3 = typecast(tmp(ifd2 + 2 + nentries2*12:ifd2 + 2 + nentries2*12 + 3), 'uint32')+1;
img3OffsetLoc = typecast(tmp(ifd3 + 2 + 12*7 + 8:ifd3 + 2 + 12*7 + 11), 'uint32')+1;
img3Offset = typecast(tmp(img3OffsetLoc:img3OffsetLoc+3),'uint32')+1;

%% 
% Number of bytes per image.
bytesperimage = img3Offset - img2Offset;
nframes = length(img2Offset:bytesperimage:L) + 1;

% Put all the data in an array.
dat = zeros(bytesperimage/2, nframes);
% handle the first image separately
dat(:,1) = typecast(tmp(img1Offset:img1Offset + bytesperimage-1), 'uint16');

% Pad the data with zeros to accommodate a different number of tags in the
% last directory
Lpad = max(img2Offset:bytesperimage:L) + bytesperimage - L - 1;
tmp = [tmp; zeros(Lpad,1)];

% Line up all the images side-by-side (as vectors), with the metadata at
% the bottom.
dat(:,2:end) = reshape(typecast(tmp(img2Offset:end), 'uint16'), bytesperimage/2, nframes-1);

% Discard the metadata
dat = dat(1:xsize*ysize,:);

% Reshape into a three-dimensional array
dat2 = reshape(dat, xsize, ysize, nframes);
% adopt the conventional orientation of the axes and convert data % to Double 
out = double(permute(dat2, [2 1 3]));

