function plot_VtrigGlu(VtrigGlu, coord, AvgImg, scoreMs, figbase, voltFtprnt, synSel)
% PLOT_VTRIGGLU  Plot SS/CS spike-triggered glutamate STA + pre/post-spike tuning maps.
%   plot_VtrigGlu(VtrigGlu, coord, AvgImg, scoreMs, figbase, voltFtprnt)
%   VtrigGlu : struct array (per good neuron) with .roi and .SS/.CS substructs, each
%              holding .STA (Nsyn x winLen), .tMs, .score (Nsyn x W), .pval (Nsyn x W),
%              .ntrig  (as produced by S4b via staGlu_prewin). Column w of score/pval
%              corresponds to row w of scoreMs.
%   coord    : Nsyn x 2 synapse locations   AvgImg : HxW cell image
%   scoreMs  : W x 2 scoring windows (ms). Row 1 = pre-spike [-50 0], row 2 = post-spike
%              [0 50] (any # of rows is accepted).
%   figbase  : figures are opened at figbase + roi (one per neuron)
%   voltFtprnt (optional) : H x W x nROI voltage footprints aligned into the glutamate
%              image space. If given, footprint(:,:,roi) is drawn as a cyan contour on
%              each spatial (right) map to mark the triggering soma location.
%   synSel (optional) : Nsyn logical (or index list) restricting the display to a subset
%              of synapses, e.g. place-tuned synapses (STA_Result.SpatialInfo_p<0.05).
%              Default = all synapses.
%
% Per neuron, for each class (SS, CS):
%   LEFT (one tall panel) : STA-glutamate heatmap of ALL synapses, ordered as
%       [ synapses tuned in window 1 (sorted by that window's AUC, desc) ;
%         synapses tuned in window 2 but not window 1 (sorted by window-2 AUC) ; ... ;
%         the remaining (untuned) synapses, sorted by mean amplitude, desc ].
%       Green dashed lines separate the blocks.
%   RIGHT (one map per window) : spatial map; marker SIZE grows with -log10(pval)
%       (lower p -> bigger dot), color = z-scored mean amplitude (score) in that window.
    if nargin<6, voltFtprnt = []; end
    if nargin<7 || isempty(synSel), synSel = true(size(coord,1),1); end
    selIdx   = find(synSel);                    % synapses to display (e.g. place-tuned)
    coord    = coord(selIdx,:);                 % restrict coordinates to the selected set
    imgRGB   = repmat(mat2gray(double(AvgImg)),1,1,3);
    classes  = {'SS','CS'};
    W        = size(scoreMs,1);                 % # score windows (pre, post, ...)
    baseEdge = min(scoreMs(:,1));               % baseline = STA before the earliest window
    for gi = 1:numel(VtrigGlu)
        roi = VtrigGlu(gi).roi;
        figure(figbase+roi); clf;
        t = tiledlayout(numel(classes)*W, 2, 'TileSpacing','compact','Padding','compact');
        for ci = 1:numel(classes)
            C   = VtrigGlu(gi).(classes{ci});
            STA = C.STA(selIdx,:); tMs = C.tMs; Nsyn = size(STA,1);
            pvalAll = C.pval(selIdx,:); scoreAll = C.score(selIdx,:);   % restricted to selected synapses
            base    = mean(STA(:, tMs < baseEdge), 2, 'omitnan');  % pre-window baseline
            STc     = STA - base;                                  % baseline-subtracted STA
            meanAmp = mean(STc, 2, 'omitnan');                     % overall mean amplitude / synapse

            % ---- combined ordering: window-1 tuned, then window-2 tuned, ... then rest ----
            placed = false(Nsyn,1);  ord = [];  blockEdge = zeros(1,W);
            for w = 1:W
                inW   = tMs>=scoreMs(w,1) & tMs<=scoreMs(w,2);
                aucW  = trapz(tMs(inW), STc(:,inW), 2);            % Nsyn x 1
                idxW  = find(pvalAll(:,w) < 0.05 & ~placed);       % tuned & not already placed
                [~,o] = sort(aucW(idxW),'descend','MissingPlacement','last');
                idxW  = idxW(o);
                ord   = [ord; idxW];  placed(idxW) = true;  blockEdge(w) = numel(ord); %#ok<AGROW>
            end
            restIdx = find(~placed);
            [~,o]   = sort(meanAmp(restIdx),'descend','MissingPlacement','last');
            ord     = [ord; restIdx(o)];

            % ---- LEFT: one tall STA heatmap (all synapses, combined order) ----
            r0 = (ci-1)*W + 1;  leftTile = (r0-1)*2 + 1;           % top-left tile of this class
            nexttile(t, leftTile, [W 1]);
            shown = STc(ord,:);
            hi = prctile(abs(shown(:)),98);                        % adaptive color scale
            if ~isfinite(hi) || hi<=0, hi = max(abs(shown(:))); end
            if isempty(hi) || ~isfinite(hi) || hi<=0, hi = 1; end
            imagesc(tMs, 1:Nsyn, shown, [-hi hi]); hold on; axis tight;
            for w = 1:W, xline(scoreMs(w,1),'w--'); xline(scoreMs(w,2),'w--'); end
            xline(0,'w-');
            for w = 1:W-1, yline(blockEdge(w)+0.5,'g-','LineWidth',1); end   % block separators
            yline(blockEdge(W)+0.5,'g-','LineWidth',1);                      % tuned | rest
            set(gca,'YDir','reverse');
            xlabel('Time from spike (ms)'); ylabel('synapse (tuned first, then by amplitude)');
            colorbar;
            title(sprintf('ROI %d  %s STA glu (n=%d)  tuned: %s', ...
                roi, classes{ci}, C.ntrig, mat2str(diff([0 blockEdge]))));

            % ---- RIGHT: one spatial map per window ----
            for w = 1:W
                win   = scoreMs(w,:);  winTag = sprintf('%g..%g ms', win(1), win(2));
                score = scoreAll(:,w);  pval = pvalAll(:,w);  sig = pval < 0.05;
                rTile = ((ci-1)*W + w - 1)*2 + 2;
                nexttile(t, rTile); image(imgRGB); hold on; axis image off;
                if ~isempty(voltFtprnt) && roi<=size(voltFtprnt,3)   % soma footprint contour
                    fp = voltFtprnt(:,:,roi); mfp = max(fp(:));
                    if isfinite(mfp) && mfp>0
                        contour(fp/mfp, [0.2 0.5], 'c', 'LineWidth', 1);
                    end
                end
                scatter(coord(:,1),coord(:,2),2,[.5 .5 .5],'filled');
                if any(sig)
                    sz  = 2 + 4*(-log10(pval(sig)));             % lower p -> bigger dot
                    zsc = (score - mean(score,'omitnan')) ./ std(score,'omitnan');   % z across synapses
                    scatter(coord(sig,1),coord(sig,2), sz, zsc(sig),'filled');
                    colormap(gca,'autumn'); mx = prctile(zsc(sig),90);
                    if isfinite(mx) && mx>0, caxis([0 mx]); end
                    cb = colorbar; cb.Label.String = 'mean amplitude (z-score)';
                end
                title(sprintf('%s tuned [%s]: %d/%d (size = -log10 p)', ...
                    classes{ci}, winTag, sum(sig), Nsyn));
            end
        end
        sgtitle(sprintf('ROI %d : spike-triggered glutamate tuning (pre & post)', roi));
    end
    drawnow;
end
