%% ============================ Local functions ============================
function [movR, mov_mc, mtr, nT, M] = loadVoltChunk(fpath_f, j, sz, readFrame, prefix)
% Load one motion-corrected voltage chunk and return:
%   movR   : spike-positive motion residual (-mov_res), H x W x nT
%   mov_mc : raw motion-corrected movie,               H x W x nT
%   mtr    : motion trace (nT x k), smoothed
%   nT     : number of frames
%   M      : (optional) each pixel's mean of the bleach-normalised movie, H x W. movR has every
%            pixel's chunk mean removed; M - movR is the motion-regressed movie with its baseline
%            (needed to join chunks without steps)
%
% prefix (optional): movie-file name prefix before the zero-padded segment
%   number. Default 'mc' -> mc01.bin, mc02.bin, ...
%   Use 'mc_ShutterReg' for mc_ShutterReg01.bin, mc_ShutterReg02.bin, ...
%   (The mcTrace files are always named mcTraceNN.mat.)
if nargin < 5 || isempty(prefix), prefix = 'mc'; end

S = load(fullfile(fpath_f, ['mcTrace' num2str(j,'%02d') '.mat']), 'mcTrace');
mcTrace = S.mcTrace;
if isstruct(mcTrace), mcTrace = mcTrace.xymean; end
mtr = movmean(mcTrace, 3, 1);

movfile = fullfile(fpath_f, [prefix num2str(j,'%02d') '.bin']);
try
    mov_mc = double(readBinMov_times(movfile, sz(2), sz(1), 1:readFrame));
catch
    mov_mc = double(readBinMov(movfile, sz(2), sz(1)));
end
nT  = size(mov_mc,3);
mtr = mtr(1:nT, :);

% photobleach normalization + motion-residual movie
mean_F  = squeeze(mean(mov_mc,[1 2]));
y_fit   = expfitDM_2((1:nT)', mean_F, (1:nT)', [10000 100]);
mov_mc2 = mov_mc ./ reshape(y_fit,1,1,[]);
M       = mean(mov_mc2,3);
mov_res = mov_mc2 - M;
mov_res = SeeResiduals(mov_res, mtr);
mov_res = SeeResiduals(mov_res, mtr.^2);
mov_res = SeeResiduals(mov_res, mtr(:,1).*mtr(:,end));
movR = -mov_res;   % negative-going voltage -> spike-positive
end
