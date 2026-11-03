function [mov, nframes] = readBinMovSingle_times(fileName, nrow, ncol, framelist)
% [mov, nframes] = readBinMovSingle_times(fileName, nrow, ncol, framelist)
% Reads the frames framelist of a float32 movie written frame by frame as nrow x ncol arrays
% (fwrite(fid, single(mov), 'float32'), the layout read by readBinMovSingle), e.g. the
% mc2_quicklook_fr<a>-<b>.bin of Analysis_20261025_SD_Glu_ASAP_quicklook.m.
% The output, mov, is a nrow x ncol x numel(framelist) single array.
% A contiguous framelist is read in one call.

fid = fopen(fileName);
nframes = numel(framelist);
framesize = nrow*ncol;

if nframes > 1 && all(diff(framelist) == 1)
    fseek(fid, (framelist(1)-1)*framesize*4, 'bof');   % 4 bytes/pixel
    mov = fread(fid, framesize*nframes, '*single', 'l');
else
    mov = zeros(framesize, nframes, 'single');
    for j = 1:nframes
        fseek(fid, (framelist(j)-1)*framesize*4, 'bof');
        mov(:,j) = fread(fid, framesize, '*single', 'l');
    end
end
fclose(fid);

mov = reshape(mov, [nrow ncol nframes]);
