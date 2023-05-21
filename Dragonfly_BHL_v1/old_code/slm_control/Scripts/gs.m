function [approxSLM, stats] = gs(commonz, zpolys, ttd, gsIter, randomInit, weightedGS, invertMatrix, useGPU, showStats, targetAmps)
if ~exist('gsIter','var'), gsIter = 2; end
if ~exist('randomInit','var'), randomInit = false; end
if ~exist('weightedGS','var'), weightedGS = true; end
if ~exist('transposedFwdProjection','var'), invertMatrix = false; end
if ~exist('useGPU','var'), useGPU = true; end
if ~exist('showStats','var'), showStats = false; end
if ~exist('targetAmps','var'), targetAmps = ones(size(ttd,2),1); end
%% GS Init
targetZernikeCoeffs = commonz*ones(1,size(ttd,2));
targetZernikeCoeffs([1 2 4],:) = targetZernikeCoeffs([1 2 4],:) + ttd;
temp = reshape(zpolys.data,prod(zpolys.imsz),[]);
if useGPU
    try
        temp = gpuArray(temp);
    catch
        warning 'gs: no gpu?, using cpu instead'
    end
end
backProjM = exp(1i*temp*targetZernikeCoeffs);
backProjM = backProjM*pinv(backProjM'*backProjM);
% normBackProjM = 1; % norm(backProjM); whos
if invertMatrix
    fwdProjM = pinv(backProjM);
end
%% GS reset
% targetAmps = ones(size(ttd,2),1); % sqrt([3; 1; 2]); % 
weightAmps = ones(size(ttd,2),1);
approxAmps = targetAmps;
if randomInit
    approxAmps = approxAmps.*exp(1i*2*pi*rand(size(targetAmps)));
end
if useGPU
    try
        backProjM = gpuArray(backProjM);
        if invertMatrix
            fwdProjM = gpuArray(fwdProjM);
        end
        targetAmps = gpuArray(targetAmps);
        approxAmps = gpuArray(approxAmps);
    catch
        warning 'gs: no gpu?, using cpu instead'
    end
end
%% GS iteration
t = tic;
for it = 1:gsIter
    % reset image space amplitudes, keep phase
    if weightedGS
        % use weighted GS: update, then apply weights
        weightAmps = weightAmps.*targetAmps.*mean(abs(approxAmps))./abs(approxAmps);
        approxAmps = weightAmps.*approxAmps./abs(approxAmps);
    else
        % alternative strategy, no weights
        approxAmps = targetAmps.*approxAmps./abs(approxAmps);
    end

    % recalculate SLM phase and amplitude
    approxSLM = backProjM*approxAmps; % ./normBackProjM.^2;
    % set SLM amplitude to 1, keep phase
    approxSLM = approxSLM./abs(approxSLM);

    % recalculate image complex amplitude and phase
    if invertMatrix
        % alternative strategy, use pseudo-inverted matrix
        approxAmps = fwdProjM*approxSLM;
    else
        % use transposed matrix, skips inversion step
        approxAmps = backProjM'*approxSLM; % ./normBackProjM.^2;
    end
end
%% wrap up
t = toc(t);
if useGPU
    approxSLM = gather(approxSLM);
    approxAmps = gather(approxAmps);
end
if showStats || nargout>1
    % imshow(sum(vm(angle(approx_cplxit),zpolys.imsz)),[])
    allim = abs(approxAmps).^2;
    en = sum(allim);
    un = 1-(max(allim)-min(allim))./(max(allim)+min(allim));
    si = std(allim,0)/mean(allim);
    stats = [en un si t];
    if showStats, disp(stats), end
end
approxSLM = toimg(approxSLM,zpolys.imsz);
