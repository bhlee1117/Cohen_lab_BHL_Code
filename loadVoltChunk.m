%% ============================ Local functions ============================
function [movR, mov_mc, mtr, nT] = loadVoltChunk(fpath_f, j, sz, readFrame)
% Load one motion-corrected voltage chunk and return:
%   movR   : spike-positive motion residual (-mov_res), H x W x nT
%   mov_mc : raw motion-corrected movie,               H x W x nT
%   mtr    : motion trace (nT x k), smoothed
%   nT     : number of frames
S = load([fpath_f '/mcTrace' num2str(j,'%02d') '.mat'], 'mcTrace');
mcTrace = S.mcTrace;
if isstruct(mcTrace), mcTrace = mcTrace.xymean; end
mtr = movmean(mcTrace, 3, 1);

try
    mov_mc = double(readBinMov_times([fpath_f '/mc' num2str(j,'%02d') '.bin'], sz(2), sz(1), 1:readFrame));
catch
    mov_mc = double(readBinMov([fpath_f '/mc' num2str(j,'%02d') '.bin'], sz(2), sz(1)));
end
nT  = size(mov_mc,3);
mtr = mtr(1:nT, :);

% photobleach normalization + motion-residual movie
mean_F  = squeeze(mean(mov_mc,[1 2]));
y_fit   = expfitDM_2((1:nT)', mean_F, (1:nT)', [10000 100]);
mov_mc2 = mov_mc ./ reshape(y_fit,1,1,[]);
mov_res = mov_mc2 - mean(mov_mc2,3);
mov_res = SeeResiduals(mov_res, mtr);
mov_res = SeeResiduals(mov_res, mtr.^2);
mov_res = SeeResiduals(mov_res, mtr(:,1).*mtr(:,end));
movR = -mov_res;   % negative-going voltage -> spike-positive
end
