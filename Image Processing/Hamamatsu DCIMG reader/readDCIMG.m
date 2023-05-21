function out = readDCIMG(filename, framelist);
% function out = readDCIMG(filename, framelist);
% filename: string specifying dcimg file name.
% framelist: a list of frame numbers to read in.  If omitted, reads all
% frames.

[framedata,nframes] = dcimgmatlab(1, filename);
framedata = framedata';

if nargin == 1;
    framelist = 1:nframes;
end;

if max(framelist) > nframes;
    sprintf('Invalid frame list')
    out = [];
    return
end;
[ysize, xsize] = size(framedata);

readframes = length(framelist);
framelist = reshape(framelist, [1 readframes]); % convert to a row vector
out = uint16(zeros(ysize, xsize, readframes));

c = 1;
for j = framelist;
    [framedata,nframes]= dcimgmatlab(j, filename);
    out(:,:,c) = framedata';
    c = c + 1;
    j;
end;

