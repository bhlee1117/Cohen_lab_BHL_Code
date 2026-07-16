function interactive_glu_tuningcorr(AvgImg, coord, tuning, dotSize)
% INTERACTIVE_GLU_TUNINGCORR  Peak-position map + tuning-curve correlation map.
%
%   interactive_glu_tuningcorr(AvgImg, coord, tuning, dotSize)
%   AvgImg  : HxW average glutamate image (grayscale background)
%   coord   : N x 2   synapse locations [x y]  (pass ONLY the synapses to show, e.g. SI_p<0.05)
%   tuning  : N x place_bin   place tuning curve per synapse (e.g. mean-over-laps firing rate)
%   dotSize : optional scatter marker size (default 14)
%
% TOP    : synapses colored by PEAK position, expressed as an ANGLE (0..2*pi) since the
%          VR track is circular -> circular (hsv) colormap.
% BOTTOM : synapses colored by the Pearson correlation of their tuning curve with the
%          CLICKED synapse's tuning curve (caxis [-1 1]). Click a synapse in either panel.

    [H,W] = size(AvgImg);
    if nargin<4 || isempty(dotSize), dotSize = 14; end
    N  = size(coord,1);  pb = size(tuning,2);
    if N==0, warning('No synapses to show.'); return; end

    T = ringmovMean(tuning, 3);              % smooth tuning curves along position (circular)

    % peak position -> angle (circular)
    [~, pk] = max(T, [], 2);
    theta   = (pk-1)/pb * 2*pi;              % 0..2*pi

    % pairwise tuning-curve correlation (Pearson) between the shown synapses
    Tf = T; Tf(isnan(Tf)) = 0;
    C  = corr(Tf');                          % N x N
    C(isnan(C)) = 0;                         % flat curves -> undefined -> 0

    sel = 1;
    fig = figure('Name','Glu tuning correlation','NumberTitle','off','Position',[120 60 780 920]);
    tl  = tiledlayout(fig,2,1,'TileSpacing','compact','Padding','compact');
    axT = nexttile(tl,1);  axB = nexttile(tl,2);
    imgRGB = repmat(mat2gray(double(AvgImg)),1,1,3);

    % ---- top: peak-position angle ----
    hT   = image(axT, imgRGB); hold(axT,'on'); axis(axT,'image','off');
    hScT = scatter(axT, coord(:,1), coord(:,2), dotSize, theta, 'filled');
    colormap(axT, hsv(256)); caxis(axT,[0 2*pi]);
    cbT = colorbar(axT); cbT.Label.String = 'peak position (angle)';
    cbT.Ticks = [0 pi/2 pi 3*pi/2 2*pi]; cbT.TickLabels = {'0','\pi/2','\pi','3\pi/2','2\pi'};
    hSelT = plot(axT, coord(sel,1), coord(sel,2), 'o','MarkerSize',13,'MarkerEdgeColor','w','LineWidth',2,'HitTest','off');
    title(axT, sprintf('Peak position  (%d synapses) - click one', N));
    set(hT ,'ButtonDownFcn',@(s,e)onClick(axT));
    set(hScT,'ButtonDownFcn',@(s,e)onClick(axT));

    % ---- bottom: correlation to selected ----
    hB   = image(axB, imgRGB); hold(axB,'on'); axis(axB,'image','off');
    hScB = scatter(axB, coord(:,1), coord(:,2), dotSize, C(:,sel), 'filled');
    colormap(axB, 'turbo'); caxis(axB,[0 0.5]);
    cbB = colorbar(axB); cbB.Label.String = 'tuning corr to selected';
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
        title(axB, sprintf('Tuning corr to synapse %d', sel));
        drawnow;
    end
end
