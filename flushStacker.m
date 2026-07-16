function S = flushStacker(S, fpath_f, base, offsetVal)
% FLUSHSTACKER  Write a stacker's current raw-window buffer to a numbered .bin
% part (vm/transpose format, reloadable with readBinMov_times) and record its
% metadata.  Called by accumStack when the buffer is large, and once at the end.
if isempty(S.buf), return; end
S.part = S.part + 1;
fname = sprintf('%s_ROI%02d_%s_part%02d.bin', base, S.roi, clsName(S.cls), S.part);
MovtoWrite = vm(double(S.buf)+offsetVal);
MovtoWrite.transpose.savebin(fullfile(fpath_f, fname));
S.binFiles{end+1}  = fname;
S.evPerFile(end+1) = size(S.buf,3);
S.centerFrames     = [S.centerFrames S.centers];
S.centers = [];
S.buf = [];
end
