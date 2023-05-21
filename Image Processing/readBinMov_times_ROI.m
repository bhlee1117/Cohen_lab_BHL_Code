function [mov, nframes] = readBinMov_times_ROI(fileName, nrow, ncol, framelist, ROI)
% [mov nframes] = [mov, nframes] = readBinMov_times_ROI(fileName, nrow, ncol, framelist, ROI)
% framelist: a list of the frames to read in.
% ROI: format is [rmin rmax cmin cmax]
% If only one output: mov = [mov, nframes] = readBinMov_times_ROI(fileName, nrow, ncol, framelist, ROI)
% The output, mov, is a 3-D array or unsigned 16-bit integers
% The binary input file is assumed to be written by Labview from a 
% Hamamatsu .dcimg file.  The binary file has 16-bit integers with a little
% endian format.
% The frame size information must be independently documented in a .txt 
% file and read into this function as the number of rows (nrow) and the 
% number of columns (ncol).
% This script reads in data frame-by-frame.  For each frame, it first reads
% columns containing the ROI (continuous block of data), and then keep only
% the ROI region.

% define image size
rlist = (ROI(1):ROI(2));        % ROI row
clist = (ROI(3):ROI(4));        % ROI column
nrow_roi = length(rlist);       % number of rows in ROI
ncol_roi = length(clist);       % number of columns in ROI
nframes = length(framelist);    % number of frames in the frame list

% define start locations at given frames, skip to the first column of ROI
startlocs = nrow*ncol*(framelist-1) + nrow*(clist(1)-1);

% read in frame-by-frame
mov = uint16(zeros(nrow_roi,ncol_roi, nframes));
tic
fid = fopen(fileName);                  % open file
for j = 1:nframes
    fseek(fid, startlocs(j)*2, 'bof');  % factor of two is for two bytes/pixel
    % read in subframe, uint16, little endian
    tmp = fread(fid,[nrow,ncol_roi],'*uint16','l');
    mov(:,:,j) = tmp(rlist,:);
end
fclose(fid);                            % close file
toc

% notes: fread works by fill in the matrix column-wise; similarly, fwrite
% also writes matrix to binary in a column-by-column format.