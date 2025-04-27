function plotVennBinary(varargin)
    % Ensure at least two or three input binary vectors
    numVectors = length(varargin);
    if numVectors < 2 || numVectors > 3
        error('Function supports only 2 or 3 binary vectors.');
    end
    
    % Ensure all vectors have the same length and are binary
    vecLength = length(varargin{1});
    for i = 1:numVectors
        vec = varargin{i};
        if length(vec) ~= vecLength
            error('All input vectors must have the same length.');
        end
        if ~all(ismember(vec, [0, 1]))
            error('All input vectors must be binary (contain only 0s and 1s).');
        end
    end

    % Compute unique and intersection counts
    A_only = sum(varargin{1} & ~varargin{2} & (numVectors == 2 || ~varargin{3}));
    B_only = sum(varargin{2} & ~varargin{1} & (numVectors == 2 || ~varargin{3}));
    AB_intersect = sum(varargin{1} & varargin{2} & (numVectors == 2 || ~varargin{3}));
    
    if numVectors == 2
        % Define area proportions for the Venn diagram
        A = [A_only + AB_intersect, B_only + AB_intersect]; % Total counts for each set
        I = AB_intersect; % Intersection count
        
        % Plot Venn diagram
        figure;
        [H, S] = venn(A, I, 'FaceColor', {'r', 'b'}, 'FaceAlpha', {0.6, 0.6}, 'EdgeColor', 'black');
        title('Venn Diagram of Two Binary Vectors');

        % Display text labels in the center of each section
        hold on;
        text(S.ZoneCentroid(1,1), S.ZoneCentroid(1,2), sprintf('%d', A_only), 'FontSize', 12, 'Color', 'r', 'HorizontalAlignment', 'center');
        text(S.ZoneCentroid(2,1), S.ZoneCentroid(2,2), sprintf('%d', B_only), 'FontSize', 12, 'Color', 'b', 'HorizontalAlignment', 'center');
        text(S.ZoneCentroid(3,1), S.ZoneCentroid(3,2), sprintf('%d', AB_intersect), 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
        hold off;
        
    elseif numVectors == 3
        C_only = sum(varargin{3} & ~varargin{1} & ~varargin{2});
        AC_intersect = sum(varargin{1} & varargin{3} & ~varargin{2});
        BC_intersect = sum(varargin{2} & varargin{3} & ~varargin{1});
        ABC_intersect = sum(varargin{1} & varargin{2} & varargin{3});
        
        % Define area proportions for the Venn diagram
        A = [A_only + AB_intersect + AC_intersect, ...
             B_only + AB_intersect + BC_intersect, ...
             C_only + AC_intersect + BC_intersect];
        I = [AB_intersect, AC_intersect, BC_intersect, ABC_intersect];
        
        % Plot Venn diagram
        figure;
        [H, S] = venn(A, I, 'FaceColor', {'r', 'b', 'g'}, 'FaceAlpha', {0.6, 0.6, 0.6}, 'EdgeColor', 'black');
        title('Venn Diagram of Three Binary Vectors');
        
        % Display text labels in the center of each section
        hold on;
        text(S.ZoneCentroid(1,1), S.ZoneCentroid(1,2), sprintf('%d', A_only), 'FontSize', 12, 'Color', 'r', 'HorizontalAlignment', 'center');
        text(S.ZoneCentroid(2,1), S.ZoneCentroid(2,2), sprintf('%d', B_only), 'FontSize', 12, 'Color', 'b', 'HorizontalAlignment', 'center');
        text(S.ZoneCentroid(3,1), S.ZoneCentroid(3,2), sprintf('%d', C_only), 'FontSize', 12, 'Color', 'g', 'HorizontalAlignment', 'center');
        text(S.ZoneCentroid(4,1), S.ZoneCentroid(4,2), sprintf('%d', AB_intersect), 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
        text(S.ZoneCentroid(5,1), S.ZoneCentroid(5,2), sprintf('%d', AC_intersect), 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
        text(S.ZoneCentroid(6,1), S.ZoneCentroid(6,2), sprintf('%d', BC_intersect), 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
        text(S.ZoneCentroid(7,1), S.ZoneCentroid(7,2), sprintf('%d', ABC_intersect), 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'k', 'HorizontalAlignment', 'center');
        hold off;
    end
end
