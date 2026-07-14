% alpha_injection_test.m
% Estimate detection limit by injecting synthetic alpha-function events
% into mov_res and running the existing detection algorithm.
%
% Requires in workspace:
%   mov_res      : base movie [Y x X x T]
%   Result.swc   : dendrite skeleton [N x 2], columns = [col(x), row(y)] in pixels

% -------------------------------------------------------------------------
% Parameter grid
% -------------------------------------------------------------------------
lengths   = [5 10 20];         % dendrite segment lengths [µm]
alphas    = [0.01 0.02 0.05];  % amplitude fractions (dF/F)
taus      = [3 5 10];          % alpha time constants [frames]

% -------------------------------------------------------------------------
% Options
% -------------------------------------------------------------------------
[nY, nX, nT] = size(mov_res);
t0         = round(nT / 3);    % injection frame
tol        = 5;                % detection tolerance [frames]
s          = -1;               % sign: -1 for Voltron2 depolarization
px_per_um  = 1;                % pixels per µm
gauss_sigma = 1.0;             % Gaussian width [px] for skeleton footprint

% seg_center: [row, col] — closest skeleton point to this is used as anchor
seg_center = [round(nY/2), round(nX/2)];   % adjust to your ROI

% -------------------------------------------------------------------------
% Step 2: Baseline image F0
% -------------------------------------------------------------------------
F0 = mean(mov_res, 3);   % [Y x X]

% -------------------------------------------------------------------------
% Step 3: Order skeleton and find anchor point closest to seg_center
% -------------------------------------------------------------------------
% Result.swc assumed [N x 2] = [col(x), row(y)]
skel_col = Result.swc(:, 1);   % x in image = column
skel_row = Result.swc(:, 2);   % y in image = row

% Find skeleton point closest to seg_center
dist_to_center = sqrt((skel_row - seg_center(1)).^2 + (skel_col - seg_center(2)).^2);
[~, anchor_idx] = min(dist_to_center);

% -------------------------------------------------------------------------
% Step 4: Build skeleton-based spatial masks for each length
% -------------------------------------------------------------------------
% For each ell, walk ell/2 px along the skeleton in each direction from anchor
% Then paint a Gaussian blob at each selected skeleton point

[xx, yy] = meshgrid(1:nX, 1:nY);   % pixel coordinate grids

% Cumulative arc-length along skeleton (to select segment by physical length)
d_arc = [0; cumsum(sqrt(diff(skel_col).^2 + diff(skel_row).^2))];

% Pre-build masks for all lengths
S_all = zeros(nY, nX, nL);

for iL = 1:nL
    ell    = lengths(iL);
    ell_px = ell * px_per_um;
    half   = ell_px / 2;

    arc_anchor = d_arc(anchor_idx);
    in_seg = abs(d_arc - arc_anchor) <= half;
    seg_rows = skel_row(in_seg);
    seg_cols = skel_col(in_seg);

    % Paint Gaussian at each skeleton point within segment
    S = zeros(nY, nX);
    for k = 1:numel(seg_rows)
        dr = yy - seg_rows(k);
        dc = xx - seg_cols(k);
        S  = S + exp(-(dr.^2 + dc.^2) / (2 * gauss_sigma^2));
    end
    if max(S(:)) > 0, S = S / max(S(:)); end

    S_all(:,:,iL) = S;
end

% -------------------------------------------------------------------------
% Footprint preview — show all length masks overlaid on F0
% -------------------------------------------------------------------------
figure('Name', 'Injection footprints');
tiledlayout(1, nL, 'TileSpacing', 'compact', 'Padding', 'compact');

for iL = 1:nL
    nexttile;
    imagesc(F0); colormap(gray); axis image off; hold on;
    % Overlay footprint in red
    S_show = S_all(:,:,iL);
    rgb = cat(3, ones(nY,nX), zeros(nY,nX), zeros(nY,nX));
    h_ov = imagesc(rgb);
    set(h_ov, 'AlphaData', S_show * 0.6);
    % Draw full skeleton
    plot(skel_col, skel_row, 'c-', 'LineWidth', 1);
    % Mark anchor
    plot(skel_col(anchor_idx), skel_row(anchor_idx), 'y*', 'MarkerSize', 8);
    title(sprintf('\\ell = %d µm', lengths(iL)));
    hold off;
end
sgtitle('Injection footprints along dendrite skeleton (red = S, cyan = full skeleton)');

% -------------------------------------------------------------------------
% Pre-allocate results
% -------------------------------------------------------------------------
nL        = numel(lengths);
nA        = numel(alphas);
nTau      = numel(taus);
P_detect  = nan(nL, nA, nTau);

% -------------------------------------------------------------------------
% Main loop
% -------------------------------------------------------------------------
nCond    = nL * nA * nTau;
cond_idx = 0;

for iL = 1:nL
    S = S_all(:,:,iL);

    for iA = 1:nA
        alpha = alphas(iA);

        for iTau = 1:nTau
            tau = taus(iTau);
            cond_idx = cond_idx + 1;
            fprintf('Condition %d/%d  ell=%.1f um  alpha=%.3f  tau=%d frames\n', ...
                cond_idx, nCond, lengths(iL), alpha, tau);

            % Step 5: Alpha waveform h(t; tau) = (t/tau)*exp(1 - t/tau)
            t_rel = (0:nT-1)' - t0;
            h = zeros(nT, 1);
            mask = t_rel > 0;
            tr = t_rel(mask);
            h(mask) = (tr / tau) .* exp(1 - tr / tau);
            if max(h) > 0, h = h / max(h); end
            h = reshape(h, 1, 1, nT);   % [1 x 1 x nT] for broadcast

            % Step 6: Build injected movie
            % I_inj(x,y,t) = I_base(x,y,t) + s*alpha*F0(x,y)*S(x,y)*h(t)
            spatial = s * alpha * F0 .* S;   % [Y x X]
            inj     = spatial .* h;          % [Y x X x nT]
            mov_inj = mov_res + inj;

            % Step 7: Run detection algorithm
            % TODO: replace with your detection function
            % e.g.: t_det = your_detection_function(mov_inj);
            t_det = [];   % <-- replace this line

            % Step 9: Match detected event to t0
            if any(abs(t_det - t0) <= tol)
                P_detect(iL, iA, iTau) = 1;
            else
                P_detect(iL, iA, iTau) = 0;
            end

        end
    end
end

% -------------------------------------------------------------------------
% Step 10: 3D scatter plot
% -------------------------------------------------------------------------
[gL, gA, gTau] = ndgrid(lengths, alphas, taus);

figure;
scatter3(gL(:), gA(:), gTau(:), 60, P_detect(:), 'filled');
colormap(jet);
cb = colorbar;
cb.Label.String = 'P_{detect}';
clim([0 1]);
xlabel('Length (\mum)');
ylabel('Amplitude \alpha (dF/F)');
zlabel('\tau (frames)');
title('Detection probability map');
grid on;
