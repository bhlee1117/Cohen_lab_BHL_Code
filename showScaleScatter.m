function showScaleScatter(dSpikeRate, str_rot, ftprnt, cmapName, cmapRange)
% showScaleScatter Scatter plot colored by dSpikeRate within footprint
% showScaleScatter(dSpikeRate, str_rot, ftprnt, cmapName, cmapRange)
%
% Inputs:
%   dSpikeRate : [1 x N] vector of spike rates for each ROI
%   str_rot    : [M x 3] matrix of scatter points [X Y size]
%   ftprnt     : [X x Y x N] ROI footprint masks
%   cmapName   : (optional) colormap name, default 'parula'
%   cmapRange  : (optional) [min max] range for colormap mapping

if nargin < 4 || isempty(cmapName)
    cmapName = 'parula';
end

if nargin < 5 || isempty(cmapRange)
    cmapRange = [min(dSpikeRate), max(dSpikeRate)];
end

N = size(ftprnt,3);
sz = [size(ftprnt,1), size(ftprnt,2)];
% if min(str_rot(:,3))<=0
%     str_rot(:,4)=str_rot(:,4)+5;
% end
if size(str_rot,2)<4
    str_rot(:,4)=1;
end
str_rot(str_rot(:,4)==0,4)=1;

% Initialize all points to gray
scatter(str_rot(:,2), str_rot(:,1), str_rot(:,4), [0.5 0.5 0.5], 'filled');
axis equal tight off; hold on;

% Map dSpikeRate to RGB colors using vec2cmap
colors = vec2cmap(dSpikeRate, cmapName, cmapRange);

for r = 1:N
    DMDmask = ftprnt(:,:,r) > 0;

    % Find points inside the current footprint
    validIdx = str_rot(:,1)>0.5 & str_rot(:,1)<(sz(2)-0.5) & str_rot(:,2)>0.5 & str_rot(:,2)<(sz(1)-0.5);
    points = str_rot(validIdx,:);
    ind = sub2ind(sz, round(points(:,2)), round(points(:,1)));

    insideMask = DMDmask(ind);

    scatter(points(insideMask,2), points(insideMask,1), points(insideMask,4), ...
        colors(r,:), 'filled');
end

colormap(cmapName);
caxis(cmapRange);
end
