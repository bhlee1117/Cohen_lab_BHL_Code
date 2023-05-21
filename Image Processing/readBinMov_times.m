function [mov, nframes] = readBinMov_times(fileName, nrow, ncol, framelist);
% [mov, nframe] = readBinMov_times(fileName, nrow, ncol, framelist);
% If only one output: mov = readBinMov_times(fileName, nrow, ncol, framelist);
% The output, mov, is a 3-D array or unsigned 16-bit integers
% The binary input file is assumed to be written by Labview from a 
% Hamamatsu .dcimg file.  The binary file has 16-bit integers with a little
% endian format.
% The frame size information must be independently documented in a .txt 
% file and read into this function as the number of rows (nrow) and the 
% number of columns (ncol).
% framelist is a list of the frames to read in.


% read file into tmp vector
fid = fopen(fileName);                  % open file
nframes = length(framelist);

framesize = nrow*ncol;
startlocs = framesize*(framelist-1);

tmp = zeros(framesize, nframes,'uint16');

for j = 1:nframes;
    fseek(fid, startlocs(j)*2, 'bof');  % factor of two is for two bytes/pixel
    tmp(:,j) = fread(fid, framesize, '*uint16', 'l'); % uint16, little endian
end;
fclose(fid);                            % close file

% reshape vector into appropriately oriented, 3D array
mov = reshape(tmp, [ncol nrow nframes]);
mov = permute(mov, [2 1 3]);
