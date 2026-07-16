function w = loadSTAevent(I, fpath_f, k)
% LOADSTAEVENT  Reload the k-th stored RAW event window for a STAinfo entry I.
%   w = loadSTAevent(STAinfo_Glu(kk), fpath{f}, 3);   % -> H x W x winLen
% The window holds raw camera counts (motion-corrected mc/mc2); the averaged STA
% in STA_*_mov(kk).mov is the PROCESSED movie. k indexes I.centerFrames, so the
% trigger frame of the returned window is I.centerFrames(k).
cum = [0 cumsum(I.evPerFile)];
p   = find(k>cum(1:end-1) & k<=cum(2:end), 1);
if isempty(p), error('event %d out of range (have %d)', k, cum(end)); end
localIdx = k - cum(p);
v = double(readBinMov_times(fullfile(fpath_f, I.binFiles{p}), I.pixels, I.winLen, localIdx)) - I.offset;
w = reshape(v, I.H, I.W, I.winLen);
end
