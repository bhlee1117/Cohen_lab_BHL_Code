function [spikeIdx, threshold1, threshold2] = interactive_spike_finder(trace)
% interactive_spike_finder Visual spike detector with adjustable thresholds.
%   [spikeIdx, threshold1, threshold2] = interactive_spike_finder(trace)
%
%   INPUT:
%       trace : 1D numeric array
%   OUTPUT:
%       spikeIdx   : indices of detected spikes
%       threshold1 : final threshold1 value
%       threshold2 : final threshold2 value

    % Default values
    threshold1 = 5;
    threshold2 = 2;

    % Create figure
    f = figure('Name','Interactive Spike Finder','NumberTitle','off','Position',[200 200 800 400]);

    % Axes
    ax = axes('Parent',f,'Position',[0.1 0.3 0.85 0.65]);
    plotTrace = plot(ax, trace, 'k'); hold on;  grid on;
    plotSpikes = plot(ax, nan, nan, 'ro', 'MarkerSize', 6, 'DisplayName','Spikes');
    ylabel(ax,'Amplitude'); xlabel(ax,'Time');
    legend(ax, 'Trace', 'Spikes');

    % UI: Threshold 1
    uicontrol('Style','text','String','Threshold1:',...
              'Units','normalized','Position',[0.1 0.15 0.1 0.05],'HorizontalAlignment','right');
    th1Edit = uicontrol('Style','edit','String',num2str(threshold1),...
              'Units','normalized','Position',[0.21 0.15 0.1 0.05]);

    % UI: Threshold 2
    uicontrol('Style','text','String','Threshold2:',...
              'Units','normalized','Position',[0.35 0.15 0.1 0.05],'HorizontalAlignment','right');
    th2Edit = uicontrol('Style','edit','String',num2str(threshold2),...
              'Units','normalized','Position',[0.46 0.15 0.1 0.05]);

    % Update button
    uicontrol('Style','pushbutton','String','Update',...
              'Units','normalized','Position',[0.6 0.15 0.1 0.05],...
              'Callback',@updateThresholds);

    % Done button
    uicontrol('Style','pushbutton','String','Done',...
              'Units','normalized','Position',[0.75 0.15 0.1 0.05],...
              'Callback','uiresume(gcbf);');

    % Initial detection
    updateThresholds();

    % Wait for user
    uiwait(f);

    % Return output
    threshold1 = str2double(th1Edit.String);
    threshold2 = str2double(th2Edit.String);
    spikeIdx = get(plotSpikes, 'XData');

    % Close figure
    if isvalid(f)
        close(f);
    end

    % --- Nested update function
    function updateThresholds(~, ~)
        threshold1 = str2double(th1Edit.String);
        threshold2 = str2double(th2Edit.String);
        spikes = find_spike_bh(trace, threshold1, threshold2);
        spikes_frame= find(spikes);
        set(plotSpikes, 'XData', spikes_frame, 'YData', trace(spikes_frame));
        drawnow;
    end
end
