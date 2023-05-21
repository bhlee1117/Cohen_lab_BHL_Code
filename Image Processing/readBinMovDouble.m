
function [mov, nframe] = readBinMovDouble(fileName, nrow, ncol)
% [mov nframe] = readBinMovDouble(fileName, nrow, ncol)
% If only one input: mov = readBinMovSingle(fileName, nrow, ncol)
% The output, mov, is a 3-D array or unsigned 16-bit floats
% The binary input file is assumed to be written by Labview from a 
% Hamamatsu .dcimg file.  The binary file has 16-bit floats with a little
% endian format.
% The frame size information must be independently documented in a .txt 
% file and read into this function as the number of rows (nrow) and the 
% number of columns (ncol).

% read file into tmp vector
fid = fopen(fileName);                  % open file
tmp = fread(fid, '*double', 'l');       % uint16, little endian
fclose(fid);                            % close file

% reshape vector into appropriately oriented, 3D array
L = length(tmp)/(nrow*ncol);
mov = reshape(tmp, [nrow ncol L]);

if nargout > 1
    nframe = L;
end