function guard = flow_noPool()
% guard = flow_noPool()
%
% Stop MATLAB opening a parallel pool for the duration of the caller. Keep the returned
% object alive (just assign it to a local) and the previous setting comes back when the
% caller returns, however it returns.
%
% WHY. nanPSD loops over its rows with parfor, which is right when it is handed hundreds
% of traces and pure overhead when it is handed one -- and this pipeline hands it one, per
% segment. Measured here, the pool start-up cost was 46 s against 2.6 s for the identical
% computation without it. On this machine the licence check during pool start-up has also
% crashed MATLAB outright, so avoiding the pool is not only a speed matter.
%
% parfor with no pool runs as an ordinary serial loop, so nothing about the result changes.
%
% 2026.08.28, Byung Hun Lee

guard = [];
try
    ps = parallel.Settings;
    prev = ps.Pool.AutoCreate;
    if prev
        ps.Pool.AutoCreate = false;
        guard = onCleanup(@() restorePool(prev));
    end
catch
    % no Parallel Computing Toolbox: parfor is already serial, nothing to do
end
end

function restorePool(prev)
try
    ps = parallel.Settings;
    ps.Pool.AutoCreate = prev;
catch
end
end
