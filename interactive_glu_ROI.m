function interactive_glu_ROI(AvgImg, GluCoord, cvec, dFF, GluSpike, t_ax, Lap_dFF, Lap_FR, posUnit, dotSize)
% INTERACTIVE_GLU_ROI  Click glutamate synapses to inspect trace & place tuning.
%
%   interactive_glu_ROI(AvgImg, GluCoord, cvec, dFF, GluSpike, t_ax, Lap_dFF, Lap_FR, posUnit, dotSize)
%
%   AvgImg   : HxW average glutamate image (grayscale background)
%   GluCoord : Nsyn x 2   synapse locations [x(col) y(row)]
%   cvec     : Nsyn x 1   value coloring each synapse (e.g. -log10(SI_p))
%   dFF      : Nsyn x nT  raw glutamate traces
%   GluSpike : Nsyn x nT  glutamate events (0/1)
%   t_ax     : 1 x nT     time axis (s)
%   Lap_dFF  : nLap x place_bin x Nsyn   place x lap mean dF/F   (PlaceTrigger_average)
%   Lap_FR   : nLap x place_bin x Nsyn   place x lap event rate
%   posUnit  : optional {label, fullLength} -> x axis in real units; default = bin index
%   dotSize  : optional scatter marker size (default 8)
%
% LEFT  : AvgImg + synapses colored by cvec (click a synapse to select).
% RIGHT (3 stacked): raw trace + events | dF/F map (lap x position) | GluSpike map (lap x position).

    Nsyn = size(dFF,1);
    nT   = size(dFF,2);
    t_ax = t_ax(:)'; t_ax = t_ax(1:nT);
    pb   = size(Lap_dFF,2);  nLap = size(Lap_dFF,1);
    if nargin<9 || isempty(posUnit), xpos = 1:pb; xlab = 'Position bin';
    else, xpos = ((1:pb)-0.5)/pb*posUnit{2}; xlab = posUnit{1}; end
    if nargin<10 || isempty(dotSize), dotSize = 8; end
    cvec = cvec(:); cvec(isnan(cvec)) = 0;

    [~, sel] = max(cvec);                       % start on the most-tuned synapse

    fig = figure('Name','Interactive Glu ROI','NumberTitle','off','Position',[100 100 1350 780]);
    tl  = tiledlayout(fig,2,3,'TileSpacing','compact','Padding','compact');
    axI = nexttile(tl,1,[1 3]);                 % top row, spans all columns (cell image)
    axT = nexttile(tl,4);                        % bottom-left:   trace
    axD = nexttile(tl,5);                        % bottom-middle: dF/F map
    axS = nexttile(tl,6);                        % bottom-right:  GluSpike map

    % ---- left: image + colored synapses ----
    img  = mat2gray(double(AvgImg));
    hImg = image(axI, repmat(img,1,1,3));       % RGB so it ignores the colormap
    hold(axI,'on'); axis(axI,'image','off');
    hSc  = scatter(axI, GluCoord(:,1), GluCoord(:,2), dotSize, cvec, 'filled');
    colormap(axI,'turbo'); caxis(axI,[0 max(prctile(cvec,99.5),eps)]);
    cb = colorbar(axI); cb.Label.String = 'synapse color (e.g. -log_{10} p)';
    hSel = plot(axI, GluCoord(sel,1), GluCoord(sel,2), 'o', ...
                'MarkerSize',11,'MarkerEdgeColor','w','LineWidth',1.5,'HitTest','off');
    title(axI,'Click a synapse');
    % Callback on the OBJECTS (image + dots), not the axes: works even with axis off.
    set(hImg,'ButtonDownFcn',@onClick);
    set(hSc ,'ButtonDownFcn',@onClick);

    updatePanels();

    % ================= nested =================
    function onClick(~,~)
        cp = get(axI,'CurrentPoint'); x = cp(1,1); y = cp(1,2);
        [~, sel] = min((GluCoord(:,1)-x).^2 + (GluCoord(:,2)-y).^2);
        set(hSel,'XData',GluCoord(sel,1),'YData',GluCoord(sel,2));
        updatePanels();
    end

    function updatePanels()
        % --- top-right: raw trace + events ---
        cla(axT); plot(axT, t_ax, dFF(sel,:), 'k'); hold(axT,'on');
        spk = find(GluSpike(sel,:)>0);
        plot(axT, t_ax(spk), dFF(sel,spk), 'r.', 'MarkerSize',9);
        axis(axT,'tight'); xlabel(axT,'Time (s)'); ylabel(axT,'dF/F');
        title(axT, sprintf('Synapse %d   color=%.2f   (%d events)', sel, cvec(sel), numel(spk)));

        % --- middle: dF/F map (lap x position), position-smoothed (ringmovMean 3) ---
        Md = ringmovMean(Lap_dFF(:,:,sel), 3);
        cla(axD); hd = imagesc(axD, xpos, 1:nLap, Md); set(hd,'AlphaData',~isnan(Md));
        axis(axD,'tight'); set(axD,'YDir','reverse'); ylabel(axD,'Lap');
        title(axD,'dF/F  (lap \times position)'); colorbar(axD);

        % --- right: GluSpike rate map (lap x position), position-smoothed ---
        Ms = ringmovMean(Lap_FR(:,:,sel), 3);
        cla(axS); hs = imagesc(axS, xpos, 1:nLap, Ms); set(hs,'AlphaData',~isnan(Ms));
        axis(axS,'tight'); set(axS,'YDir','reverse'); xlabel(axS,xlab); ylabel(axS,'Lap');
        title(axS,'GluSpike rate  (lap \times position)'); colorbar(axS);
        drawnow;
    end
end
