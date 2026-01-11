function [scale pos_time]=show_traces_align_Position_wCS(traces_F, traces_CS, Subthreshold, pos, Window, VR_data, loi,scale,offset)
% [scale pos_time]=show_traces_align_Position_wCS(traces_F, traces_CS, Subthreshold, pos, Window, VR_data, loi,scale,offset)
%   Align traces with a specific position and plot with CS overlay.
%
% Inputs:
%   traces_F   : 2D matrix (cells x time)
%   traces_CS  : same size as traces_F, logical (1 = CS event)
%   pos        : target position
%   Window     : vector of time points (relative window)
%   VR_data    : matrix with VR information (channels x time)
%   loi        : laps of interest (optional)

% Figure setup
if nargin<8
    scale = range(Subthreshold);              % vertical offset between laps
end

if nargin<9
   offset=0;
end
lap_ids = unique(VR_data(8,:));   % Lap IDs

if nargin < 7 || isempty(loi)
    loi = 1:numel(lap_ids);
end

if isempty(traces_CS)
    traces_CS=zeros(size(traces_F));
end

if isempty(Subthreshold)
    Subthreshold=NaN(size(traces_F));
end

% Precompute position-aligned times
n_laps = numel(lap_ids);
pos_time = nan(1, n_laps);
lap_mask = VR_data(8,:);

for l = 1:n_laps
    lap_idx = lap_mask == lap_ids(l);
    [dist, frm] = min(abs(VR_data(5, lap_idx) - pos));
    if dist < 1
        lap_start = find(lap_idx, 1, 'first') - 1;
        pos_time(l) = frm + lap_start;
    end
end

% Create axes and plot
ax = gca; %nexttile([1 1]);
hold(ax, 'on');
cmap_base = [0 0 0];
cmap_CS = [1 0 0];
cmap_sub = [0 0 0];
zero_line = [0.7 0.7 0.7];

g = 1;  % counter for vertical offsets

for l = loi
    if ~isnan(pos_time(l))
        idx = pos_time(l) + Window;

        % Bound checking
        if idx(1) < 1 || idx(end) > size(traces_F,2)
            continue;
        end

        trace = traces_F(1, idx);
        cs_mask = traces_CS(1, idx);
        subtrace= Subthreshold(1,idx);

        % Plot baseline trace
        plot(Window, trace - g*scale - offset, 'Color', cmap_base);
        plot(Window, subtrace - g*scale - offset, 'Color', cmap_sub);

        % Overlay CS trace
        plot_CS = trace;
        plot_CS(cs_mask==0) = NaN;
        plot(Window, plot_CS - g*scale - offset, 'Color', cmap_CS, 'LineWidth', 1.2);

        % Draw baseline zero line
        %plot(Window, -g*scale*ones(size(Window)), 'Color', zero_line);

        g = g + 1;
    end
end

% Adjust axes
axis tight;
xlabel('Time (ms)');
ylabel('Laps');
yticks(-scale * ((g-1):-1:1));
yticklabels(flip(loi(1:1:end)));

end
