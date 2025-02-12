function [counts, x_center, y_center, h]=scatter_heatmap2(x,y,x_edges,y_edges)

% 
% % Define the bin edges
% x_edges = linspace(min(x), max(x), Nbin);  % 50 bins along X-axis
% y_edges = linspace(min(y), max(y), Nbin);  % 50 bins along Y-axis

% Compute the 2D histogram (heatmap data)
Npoints=length(x);
counts = histcounts2(x, y, x_edges, y_edges)/Npoints;
x_center=mean([x_edges(1:end-1); x_edges(2:end)]);
y_center=mean([y_edges(1:end-1); y_edges(2:end)]);


% Plot the heatmap
%imagesc(x_edges, y_edges, counts');
h=pcolor(x_center, y_center, counts');
shading interp
axis xy;  % Ensure correct axis orientation
colorbar;  % Add a colorbar for reference
end