function interactive_voltglu_ROI(AvgImg, G, V, posUnit, dotSize)
% INTERACTIVE_VOLTGLU_ROI  Browse glutamate synapses AND voltage footprints on one image.
%
%   interactive_voltglu_ROI(AvgImg, G, V, posUnit, dotSize)
%
%   AvgImg : HxW average glutamate image (grayscale background)
%   G (glutamate) struct:
%       .coord   Nsyn x 2   synapse locations [x y]
%       .cvec    Nsyn x 1   color per synapse (e.g. -log10 SI_p)
%       .dFF     Nsyn x nTg traces      .spike Nsyn x nTg events
%       .t       1 x nTg    time (s)
%       .Lap_dFF nLap x pb x Nsyn       .Lap_FR nLap x pb x Nsyn
%   V (voltage) struct  (footprints aligned to the glu image; SAME index order as the rest):
%       .ftprnt  H x W x nV  aligned footprints (contours)
%       .cvec    nV x 1      color per footprint (e.g. -log10 V_SI_p)
%       .trace   nV x Tv     traces        .spike nV x Tv events
%       .t       1 x Tv      time (s)
%       .Lap_FR  nLapV x pb x nV          place map (event rate)
%   posUnit : optional {label, fullLength} for the position axis; default = bin index
%   dotSize : optional glu scatter marker size (default 8)
%
% TOP  : glu image, glu synapses (dots, colored by G.cvec) + voltage footprint contours
%        (colored by V.cvec). Click a dot -> glu row; click inside a footprint -> voltage row.
% ROW2 (glutamate): trace | dF/F map (lap x position) | GluSpike map (lap x position)
% ROW3 (voltage)  : trace | place map (lap x position) | mean place tuning

    [H,W] = size(AvgImg);
    if nargin<4, posUnit = []; end
    if nargin<5 || isempty(dotSize), dotSize = 8; end
    pb = size(G.Lap_dFF,2); nLapG = size(G.Lap_dFF,1); nLapV = size(V.Lap_FR,1);
    if isempty(posUnit), xpos = 1:pb; xlab = 'Position bin';
    else, xpos = ((1:pb)-0.5)/pb*posUnit{2}; xlab = posUnit{1}; end
    Gc = G.cvec(:); Gc(isnan(Gc)) = 0;
    Vc = V.cvec(:); Vc(isnan(Vc)) = 0;
    nV = size(V.ftprnt,3);
    Gt = G.t(:)';  Vt = V.t(:)';

    % voltage footprint contour levels, centroids, and a click-hit label image
    Vlev = zeros(1,nV); Vcolor = vec2cmap(Vc,'turbo',[0 max(prctile(Vc,99.5),eps)]);
    Vcent = get_coord(V.ftprnt);                          % nV x 2 [x y]
    VlabelImg = zeros(H,W);
    for k = 1:nV
        fk = V.ftprnt(:,:,k); Vlev(k) = 0.35*max(fk(:));
        VlabelImg(fk > Vlev(k)) = k;                      % last-writer wins on overlap
    end

    selG = argmaxv(Gc);  selV = argmaxv(Vc);

    fig = figure('Name','Interactive Volt+Glu ROI','NumberTitle','off','Position',[80 60 1400 900]);
    tl  = tiledlayout(fig,3,3,'TileSpacing','compact','Padding','compact');
    axI = nexttile(tl,1,[1 3]);                 % top: image
    axGT=nexttile(tl,4); axGD=nexttile(tl,5); axGS=nexttile(tl,6);   % glu row
    axVT=nexttile(tl,7); axVP=nexttile(tl,8); axVM=nexttile(tl,9);   % volt row

    % ---- top image: glu dots + voltage contours ----
    hImg = image(axI, repmat(mat2gray(double(AvgImg)),1,1,3)); hold(axI,'on'); axis(axI,'image','off');
    hSc  = scatter(axI, G.coord(:,1), G.coord(:,2), dotSize, Gc, 'filled');
    colormap(axI,'turbo'); caxis(axI,[0 max(prctile(Gc,99.5),eps)]);
    cb = colorbar(axI); cb.Label.String = 'glu color (-log_{10} p)';
    hC = gobjects(1,nV);
    for k = 1:nV
        [~,hC(k)] = contour(axI, V.ftprnt(:,:,k), [Vlev(k) Vlev(k)], 'LineColor', Vcolor(k,:), 'LineWidth',1);
        set(hC(k),'HitTest','off');
    end
    hSelG = plot(axI, G.coord(selG,1), G.coord(selG,2), 'o','MarkerSize',11,'MarkerEdgeColor','w','LineWidth',1.5,'HitTest','off');
    title(axI,'Left-click: synapse     Right-click: footprint (contour turns white)');
    set(hSc ,'ButtonDownFcn',@(s,e)onClick());   % left -> glu synapse, right -> volt footprint
    set(hImg,'ButtonDownFcn',@(s,e)onClick());

    highlightVolt();
    updateGlu(); updateVolt();

    % ================= nested =================
    function onClick()
        st = get(fig,'SelectionType');
        cp = get(axI,'CurrentPoint'); x = cp(1,1); y = cp(1,2);
        if strcmp(st,'alt')                          % right-click -> voltage footprint (contour)
            r = round(y); c = round(x); vk = 0;
            if r>=1 && r<=H && c>=1 && c<=W, vk = VlabelImg(r,c); end
            if vk==0, [~,vk] = min((Vcent(:,1)-x).^2 + (Vcent(:,2)-y).^2); end  % else nearest footprint
            selV = vk; highlightVolt(); updateVolt();
        else                                         % left-click -> glutamate synapse
            [~,selG] = min((G.coord(:,1)-x).^2 + (G.coord(:,2)-y).^2);
            set(hSelG,'XData',G.coord(selG,1),'YData',G.coord(selG,2)); updateGlu();
        end
    end
    function highlightVolt()
        for k = 1:nV, set(hC(k),'LineColor',Vcolor(k,:),'LineWidth',1); end
        set(hC(selV),'LineColor',[1 1 1],'LineWidth',2.5);   % selected footprint -> white
    end

    function updateGlu()
        cla(axGT); plot(axGT, Gt, G.dFF(selG,:), 'k'); hold(axGT,'on');
        spk = find(G.spike(selG,:)>0); plot(axGT, Gt(spk), G.dFF(selG,spk), 'r.','MarkerSize',9);
        axis(axGT,'tight'); xlabel(axGT,'Time (s)'); ylabel(axGT,'dF/F');
        title(axGT, sprintf('Glu synapse %d  (-log_{10}p=%.2f, %d ev)', selG, Gc(selG), numel(spk)));
        Md = ringmovMean(G.Lap_dFF(:,:,selG),3);
        cla(axGD); hd=imagesc(axGD,xpos,1:nLapG,Md); set(hd,'AlphaData',~isnan(Md));
        axis(axGD,'tight'); set(axGD,'YDir','reverse'); ylabel(axGD,'Lap'); title(axGD,'Glu dF/F (lap \times pos)'); colorbar(axGD);
        Ms = ringmovMean(G.Lap_FR(:,:,selG),3);
        cla(axGS); hs=imagesc(axGS,xpos,1:nLapG,Ms); set(hs,'AlphaData',~isnan(Ms));
        axis(axGS,'tight'); set(axGS,'YDir','reverse'); xlabel(axGS,xlab); ylabel(axGS,'Lap'); title(axGS,'GluSpike rate (lap \times pos)'); colorbar(axGS);
        drawnow;
    end

    function updateVolt()
        cla(axVT); plot(axVT, Vt, V.trace(selV,:), 'k'); hold(axVT,'on');
        if isfield(V,'spike') && ~isempty(V.spike)
            spk = find(V.spike(selV,:)>0); plot(axVT, Vt(spk), V.trace(selV,spk), 'r.','MarkerSize',6);
        end
        axis(axVT,'tight'); xlabel(axVT,'Time (s)'); ylabel(axVT,'V');
        title(axVT, sprintf('Volt footprint %d  (-log_{10}p=%.2f)', selV, Vc(selV)));
        % middle: voltage-trace place map (continuous) if provided, else the spike map
        if isfield(V,'Lap_trace') && ~isempty(V.Lap_trace)
            Mt = ringmovMean(V.Lap_trace(:,:,selV),3); ttl = 'Volt trace (lap \times pos)';
        else
            Mt = ringmovMean(V.Lap_FR(:,:,selV),3);    ttl = 'Volt place map (lap \times pos)';
        end
        cla(axVP); hv=imagesc(axVP,xpos,1:nLapV,Mt); set(hv,'AlphaData',~isnan(Mt));
        axis(axVP,'tight'); set(axVP,'YDir','reverse'); xlabel(axVP,xlab); ylabel(axVP,'Lap');
        title(axVP,ttl); colorbar(axVP);
        % right-bottom: spiking place map (lap x position), like the glu spike map
        Ms = ringmovMean(V.Lap_FR(:,:,selV),3);
        cla(axVM); hs=imagesc(axVM,xpos,1:nLapV,Ms); set(hs,'AlphaData',~isnan(Ms));
        axis(axVM,'tight'); set(axVM,'YDir','reverse'); xlabel(axVM,xlab); ylabel(axVM,'Lap');
        title(axVM,'Volt spike rate (lap \times pos)'); colorbar(axVM);
        drawnow;
    end
end

function i = argmaxv(v)
[~,i] = max(v); if isempty(i), i = 1; end
end
