function scatter_density(x, y, mksz, cmap, rnge, facealpha)
% SCATTER_DENSITY  Plots scatter points color-coded by local density.
%
%   scatter_density(x, y)
%   scatter_density(x, y, mksz)
%   scatter_density(x, y, mksz, cmap)
%   scatter_density(x, y, mksz, cmap, rnge, facealpha)
%
%   Inputs:
%       x, y  : vectors of equal length
%       mksz  : marker size (default = 30)
%       cmap  : colormap (default = parula(256))

    % ---- Input checking & defaults ----
    if nargin < 3 || isempty(mksz)
        mksz = 30;
    end
    if nargin < 4 || isempty(cmap)
        cmap = parula(256);
    end
    
    if nargin <6 || isempty(facealpha)
        facealpha=0.5;
    end
    assert(isvector(x) && isvector(y) && numel(x) == numel(y), ...
        'x and y must be equal-length vectors.');

    % ---- Estimate density efficiently ----
    data = [x(:), y(:)];
    % Use a grid and interpolate for speed if large dataset
    if numel(x) > 2000
        % Grid-based estimation (faster for large data)
        nGrid = 200;
        [xx, yy] = ndgrid(linspace(min(x), max(x), nGrid), ...
                          linspace(min(y), max(y), nGrid));
        [f,~] = ksdensity(data, [xx(:), yy(:)]);
        Fgrid = griddata(xx(:), yy(:), f, data(:,1), data(:,2), 'nearest');
        density = Fgrid;
    else
        % Direct density estimation
        density = ksdensity(data, data);
        % [Prob, ~, ~, binX, binY]= histcounts2(data(:,1), data(:,2), 'BinMethod', 'auto','Normalization','probability');
        % invalidInd=find(binX==0);
        % binX(invalidInd)=1; binY(invalidInd)=1;
        % linInd=sub2ind(size(Prob),binX,binY);
        % density=Prob(linInd);
        % density(invalidInd)=NaN;
    end

    %density = density / sum(density);

    if nargin < 5 || isempty(rnge)
        rnge=[prctile(density,5) prctile(density,95)];
    end

    % ---- Map density to colors ----
    Cmap_density = vec2cmap(density, cmap, rnge);

    % ---- Plot ----
    scatter(data(:,1), data(:,2), mksz, Cmap_density, 'filled','MarkerFaceAlpha',facealpha);
    hold on;
    colorbar; % optional: to show the density scale
end
