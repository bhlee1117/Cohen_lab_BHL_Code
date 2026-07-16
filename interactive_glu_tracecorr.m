function interactive_glu_tracecorr(AvgImg, coord, tuning, traces, maxLag, dotSize)
% INTERACTIVE_GLU_TRACECORR  Peak-position map + max-lag dF/F trace-correlation map.
%
%   interactive_glu_tracecorr(AvgImg, coord, tuning, traces, maxLag, dotSize)
%   AvgImg  : HxW average glutamate image (grayscale background)
%   coord   : N x 2   synapse locations [x y]  (pass ONLY the synapses to show)
%   tuning  : N x place_bin   place tuning curve per synapse (for the TOP peak-position color)
%   traces  : N x T   raw dF/F traces (for the BOTTOM correlation)
%   maxLag  : max frame lag for the cross-correlation (default 5)
%   dotSize : optional scatter marker size (default 14)
%
% TOP    : synapses colored by PEAK position as an ANGLE (circular hsv), from `tuning`.
% BOTTOM : click a synapse -> synapses recolored by the MAX cross-correlation of their dF/F
%          trace with the clicked synapse's trace over lags -maxLag..+maxLag.

    [H,W] = size(AvgImg);
    if nargin<5 || isempty(maxLag),  maxLag  = 5;  end
    if nargin<6 || isempty(dotSize), dotSize = 14; end
    N  = size(coord,1);  pb = size(tuning,2);
    if N==0, warning('No synapses to show.'); return; end

    % ---- top color: peak position -> angle (circular) ----
    T  = ringmovMean(tuning, 3);
    [~, pk] = max(T, [], 2);
    theta   = (pk-1)/pb * 2*pi;

    % ---- bottom metric: max-lag dF/F trace cross-correlation matrix ----
    fprintf('interactive_glu_tracecorr: computing %d x %d max-lag(%d) trace corr...\n', N, N, maxLag);
    Z  = zscore(double(traces), 0, 2);          % zero-mean/unit-var per trace (over time)
    Tn = size(Z,2);
    C  = (Z * Z') / Tn;                          % lag 0 (symmetric)
    for l = 1:maxLag
        R = (Z(:,1:Tn-l) * Z(:,1+l:Tn)') / (Tn-l);   % r_ij for j lagged by +l
        C = max(C, max(R, R'));                       % cover +/- l, keep max
    end
    C(isnan(C)) = 0;

    sel = 1;
    fig = figure('Name','Glu trace correlation (max-lag)','NumberTitle','off','Position',[120 60 780 920]);
    tl  = tiledlayout(fig,2,1,'TileSpacing','compact','Padding','compact');
    axT = nexttile(tl,1);  axB = nexttile(tl,2);
    imgRGB = repmat(mat2gray(double(AvgImg)),1,1,3);

    % ---- top: peak-position angle ----
    hT   = image(axT, imgRGB); hold(axT,'on'); axis(axT,'image','off');
    hScT = scatter(axT, coord(:,1), coord(:,2), dotSize, theta, 'filled');
    colormap(axT, hsv(256)); caxis(axT,[0 2*pi]);
    cbT = colorbar(axT); cbT.Label.String = 'peak position (angle)';
    cbT.Ticks = [0 pi/2 pi 3*pi/2 2*pi]; cbT.TickLabels = {'0','\pi/2','\pi','3\pi/2','2\pi'};
    hSelT = plot(axT, coord(sel,1), coord(sel,2), 'o','MarkerSize',13,'MarkerEdgeColor','k','LineWidth',2,'HitTest','off');
    title(axT, sprintf('Peak position  (%d synapses) - click one', N));
    set(hT ,'ButtonDownFcn',@(s,e)onClick(axT));
    set(hScT,'ButtonDownFcn',@(s,e)onClick(axT));

    % ---- bottom: max-lag trace correlation to selected ----
    hB   = image(axB, imgRGB); hold(axB,'on'); axis(axB,'image','off');
    hScB = scatter(axB, coord(:,1), coord(:,2), dotSize, C(:,sel), 'filled');
    colormap(axB, 'turbo'); caxis(axB,[0 0.03]);
    cbB = colorbar(axB); cbB.Label.String = sprintf('max-lag(\\pm%d) trace corr to selected', maxLag);
    hSelB = plot(axB, coord(sel,1), coord(sel,2), 'o','MarkerSize',13,'MarkerEdgeColor','k','LineWidth',2,'HitTest','off');
    set(hB ,'ButtonDownFcn',@(s,e)onClick(axB));
    set(hScB,'ButtonDownFcn',@(s,e)onClick(axB));

    updateBottom();

    % ================= nested =================
    function onClick(ax)
        cp = get(ax,'CurrentPoint'); x = cp(1,1); y = cp(1,2);
        [~, sel] = min((coord(:,1)-x).^2 + (coord(:,2)-y).^2);
        set(hSelT,'XData',coord(sel,1),'YData',coord(sel,2));
        set(hSelB,'XData',coord(sel,1),'YData',coord(sel,2));
        updateBottom();
    end
    function updateBottom()
        set(hScB,'CData',C(:,sel));
        title(axB, sprintf('Max-lag trace corr to synapse %d', sel));
        drawnow;
    end
end
