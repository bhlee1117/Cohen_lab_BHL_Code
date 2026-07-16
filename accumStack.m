function S = accumStack(S, AddMovAvg, AddMovRaw, centersGlobal, fpath_f, base, offsetVal, flushBytes)
% ACCUMSTACK  Add a chunk's event windows to a stacker S (see newStacker).
%   AddMovAvg / AddMovRaw : [Npx x nEvent x winLen] from get_STA (processed / raw movie).
% The running AVERAGE uses the processed windows; the RAW stack stores the raw
% camera-count windows. The raw buffer is flushed to disk once it exceeds flushBytes.
nEv = size(AddMovAvg,2);
if nEv==0, return; end
winLen = size(AddMovAvg,3);
S.sum = S.sum + reshape(sum(AddMovAvg,2),[S.H S.W winLen]);   % pixel order matches reshape(H,W)
S.cnt = S.cnt + nEv;
S.buf = cat(3, S.buf, permute(AddMovRaw,[1 3 2]));           % -> [Npx x winLen x nEv], raw counts
S.centers = [S.centers centersGlobal(:)'];
if numel(S.buf)*8 > flushBytes                              % double bytes in the buffer
    S = flushStacker(S, fpath_f, base, offsetVal);
end
end
