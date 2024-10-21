function scatter_heatmap(x,y,Nbin)


% Define the bin edges
x_edges = linspace(min(x), max(x), Nbin);  % 50 bins along X-axis
y_edges = linspace(min(y), max(y), Nbin);  % 50 bins along Y-axis

% Compute the 2D histogram (heatmap data)
Npoints=length(x);
counts = histcounts2(x, y, x_edges, y_edges)/Npoints;

% Plot the heatmap
imagesc(x_edges, y_edges, counts');
axis xy;  % Ensure correct axis orientation
colorbar;  % Add a colorbar for reference
xlabel('X');
ylabel('Y');
end