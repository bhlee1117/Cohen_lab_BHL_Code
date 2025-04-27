function [counts, x_edges, y_edges]=hist_heatmap(x,y,Nbin)


% Define the bin edges
x_edges = linspace(min(x), max(x), Nbin);  % 50 bins along X-axis
y_edges = linspace(min(y), max(y), Nbin);  % 50 bins along Y-axis

% Compute the 2D histogram (heatmap data)
Npoints=length(x);
counts = histcounts2(x, y, x_edges, y_edges)/Npoints;

end