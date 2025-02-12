function l=imshow_label(dataMatrix, labels ,cmap,legend_label)
    % imshow_label displays an N x N matrix as a colored image using pcolor.
    % Labels are shown on the axes as colored dots.
    %
    % Inputs:
    %   - dataMatrix: An NxN matrix to be displayed as an image.
    %   - labels: A 1xN vector indicating the labels for the rows and columns.

    % Validate inputs

    if nargin<3
        cmap=hsv(length(unique(labels)));
    end
    cmap_nan=[0.3 0.3 0.3];
    [rows, cols] = size(dataMatrix);
    if rows ~= cols
        error('dataMatrix must be a square matrix.');
    end
    if length(labels) ~= rows
        error('Length of labels must match the dimensions of dataMatrix.');
    end
    %cmap=[0.2 0.2 0.2; hsv(max(labels))];

    % Create figure and display the matrix using pcolor
    imagesc(dataMatrix);
    colorbar;
    colormap turbo;
    axis square off;

    % Adjust axis labels
    xticks(0.5:(rows + 0.5)); % Position tick marks in the center of cells
    yticks(0.5:(cols + 0.5));
    xticklabels(labels);
    yticklabels(labels);
    set(gca, 'TickLength', [0 0]);

    % Add color-coded dots for labels on axes
    hold on;
    l=[];
    labels_cat=unique(labels(~isnan(labels)));
    g=1;
    if sum(isnan(labels))>0
        l(g) = plot(find(isnan(labels)) , repmat(0,sum(isnan(labels)),1), 'o', 'MarkerFaceColor',cmap_nan,'MarkerEdgeColor',cmap_nan); % Bottom axis dots
        plot(0, find(isnan(labels)), 'o', 'MarkerFaceColor',cmap_nan,'MarkerEdgeColor',cmap_nan); % Left axis dots
g=g+1;
    end
    for i = labels_cat   
        l(g) = plot(find(labels==i) , repmat(0,sum(labels==i),1), 'o', 'MarkerFaceColor',cmap(i,:),'MarkerEdgeColor',cmap(i,:)); % Bottom axis dots
        plot(0, find(labels==i), 'o', 'MarkerFaceColor',cmap(i,:),'MarkerEdgeColor',cmap(i,:)); % Left axis dots
        g=g+1;
    end
    hold off;
    if sum(isnan(labels))>0
    legend(l,[{'Unlabeled'},legend_label(labels_cat)],'location','northwestoutside')
    else
    legend(l,legend_label(labels_cat),'location','northwestoutside')
    end

    % Set axis limits to make room for labels
    xlim([0.5, rows + 0.5]);
    ylim([0.5, cols + 0.5]);

end