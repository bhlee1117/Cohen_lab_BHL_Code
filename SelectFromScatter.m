function idx = SelectFromScatter(A, B)

% Use current axes
ax = gca;
hold(ax, 'on')

% Draw polygon interactively
h = drawpolygon(ax, 'LineWidth', 1.5);
polyPos = h.Position;

% Determine which points are inside
in = inpolygon(A, B, polyPos(:,1), polyPos(:,2));

% Highlight selected points
scatter(ax, A(in), B(in), 50, 'r', 'LineWidth',1.5);

idx = find(in);

end