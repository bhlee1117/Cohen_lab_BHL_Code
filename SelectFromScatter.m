function idx = SelectFromScatter(X, Y, im)
% SELECT_POINTS_ROI  Draw a freehand ROI and return indices of points inside.
%
%   idx = select_points_ROI(X, Y)
%   idx = select_points_ROI(X, Y, im)
%
%   INPUTS
%     X, Y : N x 1 (or 1 x N) vectors of point coordinates
%     im   : (optional) image to display as background
%
%   OUTPUT
%     idx  : indices into X/Y of points that fall inside the drawn ROI
%
%   USAGE
%     Draw a freehand polygon on the figure, double-click to close.
%     The selected points are highlighted and their indices returned.
%
%   EXAMPLE
%     idx = select_points_ROI(synX, synY, mean_image);
%     selected_synapses = synapseData(idx);
 
X = X(:);
Y = Y(:);
 
figure; hold on;
 
% Show background image if provided
if nargin >= 3 && ~isempty(im)
    imagesc(im);
    colormap(gray);
    axis image;
end
 
% Plot all points
scatter(X, Y, 20, 'c', 'filled');
title('Draw ROI (freehand) — double-click to close');
axis ij;   % image coordinate convention (y increases downward)
 
% Draw freehand ROI
h = drawfreehand('Color','r','LineWidth',1.5);
wait(h);   % block until user double-clicks to finish
 
% Get ROI polygon vertices
roi_pos = h.Position;   % [x y] vertices
 
% Find points inside ROI
in = inpolygon(X, Y, roi_pos(:,1), roi_pos(:,2));
idx = find(in);
 
% Highlight selected points
scatter(X(idx), Y(idx), 30, 'r', 'filled');
title(sprintf('%d points selected', numel(idx)));
drawnow;
 
end