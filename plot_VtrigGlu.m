function plot_VtrigGlu(VtrigGlu, coord, AvgImg, scoreMs, figbase)
% PLOT_VTRIGGLU  Plot SS/CS spike-triggered glutamate STA + pre-spike tuning maps.
%   plot_VtrigGlu(VtrigGlu, coord, AvgImg, scoreMs, figbase)
%   VtrigGlu : struct array (per good neuron) with .roi and .SS/.CS substructs, each
%              holding .STA (Nsyn x winLen), .tMs, .score (Nsyn x1), .pval, .ntrig
%              (as produced by S4b via staGlu_prewin).
%   coord    : Nsyn x 2 synapse locations   AvgImg : HxW cell image
%   scoreMs  : [a b] pre-spike window (ms) used for the score (for labels/lines)
%   figbase  : figures are opened at figbase + roi (one per neuron)
% Per neuron (2x2): SS STA heatmap | SS pre-spike-tuned map ; CS STA heatmap | CS map.
    imgRGB  = repmat(mat2gray(double(AvgImg)),1,1,3);
    classes = {'SS','CS'};
    for gi = 1:numel(VtrigGlu)
        roi = VtrigGlu(gi).roi;
        figure(figbase+roi); clf; tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
        for ci = 1:2
            C = VtrigGlu(gi).(classes{ci});
            STA = C.STA; tMs = C.tMs; score = C.score; sig = C.pval < 0.05; Nsyn = size(STA,1);

            base = mean(STA(:, tMs < scoreMs(1)), 2, 'omitnan');   % early-pre baseline
            STc  = STA - base;
            [~,ord] = sort(score,'descend','MissingPlacement','last');
            nexttile; imagesc(tMs, 1:Nsyn, STc(ord,:), [-0.5 1]); hold on;
            xline(scoreMs(1),'w--'); xline(0,'w-');
            xlabel('Time from spike (ms)'); ylabel('synapse (sorted by pre-score)'); colorbar;
            title(sprintf('ROI %d  %s STA glu (n=%d)', roi, classes{ci}, C.ntrig));

            nexttile; image(imgRGB); hold on; axis image off;
            scatter(coord(:,1),coord(:,2),6,[.5 .5 .5],'filled');
            if any(sig), scatter(coord(sig,1),coord(sig,2),22,score(sig),'filled'); end
            colormap(gca,'turbo'); mx = max(score(sig)); if any(sig)&&mx>0, caxis([0 mx]); end; colorbar;
            title(sprintf('%s pre-spike tuned (%g..%g ms): %d/%d', classes{ci}, scoreMs(1), scoreMs(2), sum(sig), Nsyn));
        end
        sgtitle(sprintf('ROI %d : pre-spike glutamate tuning', roi));
    end
    drawnow;
end
