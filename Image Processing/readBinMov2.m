
function [mov, nframe] = readBinMov(fileName, nrow, ncol)
% [mov nframe] = readBinMov(fileName, nrow, ncol)
% If only one input: mov = readBinMov(fileName, nrow, ncol)
% The output, mov, is a 3-D array or unsigned 16-bit integers
% The binary input file is assumed to be written by Labview from a 
% Hamamatsu .dcimg file.  The binary file has 16-bit integers with a little
% endian format.
% The frame size information must be independently documented in a .txt 
% file and read into this function as the number of rows (nrow) and the 
% number of columns (ncol).

% read file into tmp vector
fProp = dir(fileName);
fSize = fProp.bytes;
fid = fopen(fileName);                  % open file

%L = int32(floor(fSize/2/(nrow*ncol)));
L = floor(fSize/2/(nrow*ncol));
mov = fread(fid,L*nrow*ncol, '*uint16', 'l');       % uint16, little endian
fclose(fid);                            % close file

% reshape vector into appropriately oriented, 3D array
% L = int32(floor(length(mov)/(nrow*ncol)));
% tmp = tmp(1:L*nrow*ncol,1);
mov = reshape(mov, [ncol nrow L]);
mov = permute(mov, [2 1 3]);

if nargout > 1
    nframe = L;
end