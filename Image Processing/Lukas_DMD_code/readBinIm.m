function im = readBinIm(fileName, nrow, ncol)
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
fid = fopen(fileName);                  % open file
tmp = fread(fid, 'ubit1', 'b');       % uint16, little endian
fclose(fid);                            % close file

% reshape vector into appropriately oriented, 2D array
im = reshape(tmp, [ncol nrow]);
im = permute(im, [2 1 3]);
end