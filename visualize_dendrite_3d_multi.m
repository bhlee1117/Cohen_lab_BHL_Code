function visualize_dendrite_3d_multi(matfile)
    % Visualize dendrite data in 3D with improved node selection for multiple nodes
    
    if nargin < 1
        matfile = 'test_dendrite_data.mat';
    end
    
    % Load the data
    data = load(matfile);
    
    % Convert types to double to avoid issues
    node_positions = double(data.node_positions);
    
    % Extract edges from adjacency matrix
    graph_adj = full(data.graph_adj);
    [row, col] = find(graph_adj);
    edges = [row, col];
    
    voltages = double(data.voltages);
    
    % Check if denoised voltages exist (and are not empty)
    has_denoised = isfield(data, 'denoised_voltages') && ~isempty(data.denoised_voltages);
    if has_denoised
        denoised_voltages = double(data.denoised_voltages);
    else
        denoised_voltages = [];
    end
    
    % Create figure with vertical layout
    fig = figure('Name', 'Dendrite 3D Explorer', 'Position', [100, 100, 1000, 800]);
    
    % Plot the dendrite tree in 3D (top)
    subplot(2, 1, 1);
    ax_dendrite = gca;
    hold(ax_dendrite, 'on');
    
    % Plot edges in 3D
    edge_handles = [];
    for i = 1:size(edges, 1)
        h = plot3(ax_dendrite, ...
             [node_positions(edges(i,1), 1), node_positions(edges(i,2), 1)], ...
             [node_positions(edges(i,1), 2), node_positions(edges(i,2), 2)], ...
             [node_positions(edges(i,1), 3), node_positions(edges(i,2), 3)], ...
             'k-', 'LineWidth', 1.5);
        edge_handles = [edge_handles; h];
    end
    
    % Define default node color (light gray) instead of using Z-depth coloring
    default_node_color = [0.7 0.7 0.7];  % Light gray
    node_colors = repmat(default_node_color, size(node_positions, 1), 1);
    
    % Create individual scatter objects for each node for better selection
    node_handles = gobjects(size(node_positions, 1), 1);
    for i = 1:size(node_positions, 1)
        node_handles(i) = scatter3(ax_dendrite, ...
                                  node_positions(i, 1), ...
                                  node_positions(i, 2), ...
                                  node_positions(i, 3), ...
                                  150, node_colors(i,:), 'filled');
        % Set specific tag for each node to identify it
        set(node_handles(i), 'Tag', ['Node' num2str(i)]);
        % Make nodes respond to clicks
        set(node_handles(i), 'PickableParts', 'all');
        % Assign the same callback to each node
        set(node_handles(i), 'ButtonDownFcn', @(src, event)node_clicked(i));
    end
    
    title(ax_dendrite, 'Dendrite Structure - Click on Nodes (Max 5)');
    xlabel(ax_dendrite, 'X');
    ylabel(ax_dendrite, 'Y');
    zlabel(ax_dendrite, 'Z');
    grid(ax_dendrite, 'on');
    axis(ax_dendrite, 'equal');
    
    % Create a container for voltage plots in the bottom half
    voltage_container = subplot(2, 1, 2);
    % Make this invisible, we'll create individual axes for each trace
    set(voltage_container, 'Visible', 'off');
    
    % Define selection colors for up to 5 nodes - using more vibrant colors
    selection_colors = [
        1 0 0;    % Red
        0 0 1;    % Blue
        0 0.7 0;  % Green
        1 0.6 0;  % Orange
        0.7 0 1   % Purple
    ];
    
    % Store data for callbacks
    setappdata(fig, 'node_positions', node_positions);
    setappdata(fig, 'voltages', voltages);
    setappdata(fig, 'denoised_voltages', denoised_voltages);
    setappdata(fig, 'has_denoised', has_denoised);
    setappdata(fig, 'node_handles', node_handles);
    setappdata(fig, 'voltage_container', voltage_container);
    setappdata(fig, 'selected_nodes', []);
    setappdata(fig, 'original_node_colors', node_colors);
    setappdata(fig, 'selection_colors', selection_colors);
    setappdata(fig, 'voltage_axes', []);
    setappdata(fig, 'default_node_color', default_node_color);
    
    % Enable 3D rotation
    rotate3d(ax_dendrite, 'on');
    
    % Set context menu with options
    c = uicontextmenu;
    uimenu(c, 'Label', 'Reset View', 'Callback', @reset_view);
    uimenu(c, 'Label', 'Top View', 'Callback', @(~,~)view(ax_dendrite, [0 90]));
    uimenu(c, 'Label', 'Side View', 'Callback', @(~,~)view(ax_dendrite, [0 0]));
    uimenu(c, 'Label', 'Front View', 'Callback', @(~,~)view(ax_dendrite, [90 0]));
    ax_dendrite.UIContextMenu = c;
    
    % Add UI controls
    uicontrol('Style', 'pushbutton', 'String', 'Reset View', ...
              'Position', [20, 20, 100, 30], 'Callback', @reset_view);
    uicontrol('Style', 'pushbutton', 'String', 'Clear Selection', ...
              'Position', [130, 20, 100, 30], 'Callback', @clear_selection);
    
    % Add selection counter text
    selection_text = uicontrol('Style', 'text', 'String', '0 of 5 nodes selected', ...
                             'Position', [240, 20, 150, 30], 'HorizontalAlignment', 'left');
    setappdata(fig, 'selection_text', selection_text);
    
    % Show console message to guide user
    disp('3D visualization ready. Click nodes to select/deselect (up to 5).');
    disp('Use click-and-drag to rotate the view, scroll to zoom, right-click for view options.');
    
    % Node click handler function that takes the node index directly
    function node_clicked(node_idx)
        selected_nodes = getappdata(fig, 'selected_nodes');
        
        % Check if node is already selected
        already_selected = ismember(node_idx, selected_nodes);
        
        if already_selected
            % Deselect the node
            selected_nodes = selected_nodes(selected_nodes ~= node_idx);
            disp(['Node ' num2str(node_idx) ' deselected!']);
        else
            % Check if we've reached the selection limit
            if length(selected_nodes) >= 5
                disp('Maximum of 5 nodes already selected. Deselect a node first.');
                return;
            end
            
            % Add the node to selection
            selected_nodes = [selected_nodes, node_idx];
            disp(['Node ' num2str(node_idx) ' selected!']);
        end
        
        % Update the selection
        setappdata(fig, 'selected_nodes', selected_nodes);
        update_selection();
    end
    
    % Function to update visualization based on current selection
    function update_selection()
        selected_nodes = getappdata(fig, 'selected_nodes');
        node_handles = getappdata(fig, 'node_handles');
        default_node_color = getappdata(fig, 'default_node_color');
        selection_colors = getappdata(fig, 'selection_colors');
        
        % Update selection counter
        selection_text = getappdata(fig, 'selection_text');
        set(selection_text, 'String', [num2str(length(selected_nodes)) ' of 5 nodes selected']);
        
        % Reset all nodes to default color
        for i = 1:length(node_handles)
            set(node_handles(i), 'CData', default_node_color);
        end
        
        % Color selected nodes
        for i = 1:length(selected_nodes)
            node_idx = selected_nodes(i);
            color_idx = min(i, size(selection_colors, 1));
            set(node_handles(node_idx), 'CData', selection_colors(color_idx, :));
        end
        
        % Update voltage plots
        update_voltage_plots();
    end
    
    % Function to update voltage plots
    function update_voltage_plots()
        selected_nodes = getappdata(fig, 'selected_nodes');
        voltage_container = getappdata(fig, 'voltage_container');
        voltage_axes = getappdata(fig, 'voltage_axes');
        selection_colors = getappdata(fig, 'selection_colors');
        has_denoised = getappdata(fig, 'has_denoised');
        voltages = getappdata(fig, 'voltages');
        denoised_voltages = getappdata(fig, 'denoised_voltages');
        
        % Clear any existing voltage axes
        if ~isempty(voltage_axes)
            for i = 1:length(voltage_axes)
                if isvalid(voltage_axes(i))
                    delete(voltage_axes(i));
                end
            end
        end
        
        % If no nodes are selected, clear the plots
        if isempty(selected_nodes)
            % Create a message in the voltage container area
            axes('Position', get(voltage_container, 'Position'));
            text(0.5, 0.5, 'Select nodes to view voltage traces', ...
                 'HorizontalAlignment', 'center', 'FontSize', 14);
            axis off;
            setappdata(fig, 'voltage_axes', []);
            return;
        end
        
        % Get the position of the voltage container
        container_pos = get(voltage_container, 'Position');
        
        % Calculate positions for each subplot
        n_plots = length(selected_nodes);
        plot_height = container_pos(4) / n_plots * 0.9; % 90% height, leave some space
        
        % Time data for plots
        time = 1:size(voltages, 2);
        
        % Create a new axis for each selected node
        new_axes = gobjects(n_plots, 1);
        for i = 1:n_plots
            node_idx = selected_nodes(i);
            color_idx = min(i, size(selection_colors, 1));
            
            % Calculate position: [left bottom width height]
            bottom = container_pos(2) + container_pos(4) - i * plot_height;
            ax_pos = [container_pos(1), bottom, container_pos(3), plot_height*0.9];
            
            % Create axis
            ax = axes('Position', ax_pos);
            new_axes(i) = ax;
            
            % Plot voltages
            hold(ax, 'on');
            if has_denoised
                p1 = plot(ax, time, voltages(node_idx, :), 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
                p2 = plot(ax, time, denoised_voltages(node_idx, :), 'Color', selection_colors(color_idx,:), 'LineWidth', 1.5);
                if i == 1 % Add legend only to the first plot
                    legend(ax, [p1, p2], {'Raw', 'Denoised'}, 'Location', 'northeast');
                end
            else
                plot(ax, time, voltages(node_idx, :), 'Color', selection_colors(color_idx,:), 'LineWidth', 1.5);
            end
            
            % Set title and labels
            title(ax, ['Node ' num2str(node_idx)], 'Color', selection_colors(color_idx,:));
            if i == n_plots % Only the bottom plot gets x label
                xlabel(ax, 'Time Steps');
            else
                set(ax, 'XTickLabel', []);
            end
            ylabel(ax, 'Voltage');
            grid(ax, 'on');
            
            % Adjust y limits for better visualization
            if has_denoised
                ydata = [voltages(node_idx, :), denoised_voltages(node_idx, :)];
            else
                ydata = voltages(node_idx, :);
            end
            y_range = [min(ydata), max(ydata)];
            padding = 0.1 * (y_range(2) - y_range(1));
            if padding > 0  % Avoid setting equal limits if there's no variation
                ylim(ax, [y_range(1) - padding, y_range(2) + padding]);
            end
        end
        
        % Link x axes for synchronization
        if n_plots > 1
            linkaxes(new_axes, 'x');
        end
        
        % Store the new axes
        setappdata(fig, 'voltage_axes', new_axes);
    end
    
    % Function to clear all selections
    function clear_selection(~, ~)
        setappdata(fig, 'selected_nodes', []);
        update_selection();
    end
    
    % Function to reset the 3D view
    function reset_view(~, ~)
        view(ax_dendrite, 3);  % Default 3D view
        axis(ax_dendrite, 'equal');
    end
end