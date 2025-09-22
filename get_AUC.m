function [AUC AUC_raw kinkAmp kink_raw] = get_AUC(trace, center, Tau_front, Tau_back)
% [AUC AUC_raw kinkAmp kink_raw] = get_AUC(trace, center, Tau_front, Tau_back)
% GET_AUC computes AUC relative to baseline around a center index for 1D, 2D, or 3D traces
%
% INPUTS:
%   trace      - 1D, 2D (N x T), or 3D (N x T x R) array of traces
%   center     - Center index (scalar for 1D, or array matching trace along dim 1 and 3)
%   Tau_front  - Samples before center to search
%   Tau_back   - Samples after center to search
%
% OUTPUT:
%   AUC        - AUC values, matching size (N x 1 x R) for 2D/3D inputs
%   AUC_raw    - AUC values but not subtract baseline, matching size (N x 1 x R) for 2D/3D inputs
%   kinkAmp    - Difference between max(segment) and baseline

    sz = size(trace);
    nd = ndims(trace);

    switch nd
        case 1  % 1D vector
            [AUC, AUC_raw, kinkAmp, kink_raw] = compute_auc_single(trace, center, Tau_front, Tau_back);
        case 2  % 2D matrix (N x T)
            N = sz(1);
            AUC = nan(N, 1);
            AUC_raw= nan(N,1);
            kinkAmp = nan(N, 1);
            kink_raw = nan(N, 1);
            for n = 1:N
                [AUC(n), AUC_raw(n), kinkAmp(n), kink_raw(n)] = compute_auc_single(trace(n,:), center(n), Tau_front, Tau_back);
            end
        case 3  % 3D matrix (N x T x R)
            N = sz(1); R = sz(3);
            AUC = nan(N, R);
            AUC_raw = nan(N, R);
            kinkAmp = nan(N, R);
            kink_raw = nan(N, R);
            for r = 1:R
                for n = 1:N
                    [AUC(n, r), AUC_raw(n,r), kinkAmp(n, r), kink_raw(n, r)] = compute_auc_single(trace(n,:,r), center(n,1,r), Tau_front, Tau_back);
                end
            end
        otherwise
            error('Trace must be 1D, 2D, or 3D.');
    end
end

function [auc, auc_raw, kink, kink_raw] = compute_auc_single(trace, center, Tau_front, Tau_back)
    N = length(trace);
    start_pre = max(center - Tau_front, 1);
    end_pre = center - floor(Tau_front/2);
    start_post = center + floor(Tau_back/2);
    end_post = min(center + Tau_back, N);

    if end_pre < start_pre || end_post < start_post
        auc = NaN;
        auc_raw= NaN;
        kink = NaN;
        return;
    end

    % Use mean of lowest 3 values instead of min
    pre_window = trace(start_pre:end_pre);
    post_window = trace(start_post:end_post);

    preVals = sort(pre_window);
    postVals = sort(post_window);

    preAmp = mean(preVals(1:min(3, numel(preVals))));
    postAmp = mean(postVals(1:min(3, numel(postVals))));

    [~, preIdxRel] = min(pre_window);
    [~, postIdxRel] = min(post_window);
    preIdx = start_pre + preIdxRel - 1;
    postIdx = start_post + postIdxRel - 1;

    segment = trace(preIdx:postIdx);
    segment_raw = trace(end_pre:start_post);
    baseline = mean([preAmp, postAmp]);
    auc = sum(segment) - baseline * numel(segment);
    auc_raw = sum(segment_raw);
    kink = max(segment) - baseline;
    kink_raw = max(segment);
end
