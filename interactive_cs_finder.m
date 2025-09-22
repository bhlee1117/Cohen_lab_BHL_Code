function [CS_trace, CS_thres] = interactive_cs_finder(sp_soma, tr_hi)
% INTERACTIVE_CS_FINDER  Visual interface to tune CS detection thresholds
% and spike parameters.
%   [CS_trace, CS_thres] = interactive_cs_finder(sp_soma, tr_hi)

% Default parameters
CS_thres = [4 1];
N_Spike = 3;
N_Spike2ISI = 3;
Max_AUC = 150;

% Process subthreshold trace
%tr_sub= get_subthreshold(tr_hi,sp_soma,7,17);
tr_sub = movmean(tr_hi, 20, 2);
tr_sub = tr_sub - movmedian(tr_sub, 300, 2);
nTime = length(tr_sub);

% Create figure and main plot
f = figure('Name','Interactive CS Finder','NumberTitle','off','Position',[200 200 900 500]);
ax = axes('Parent', f, 'Position', [0.1 0.35 0.85 0.6]);
plot(ax, tr_hi, 'w'); hold on; grid on;
plot(ax, tr_sub,'color',[1 0.8 0]);
plot(ax, find(sp_soma>0), tr_hi(find(sp_soma>0)),'go')
plotCS = plot(ax, nan, nan, 'r', 'LineWidth', 1.5, 'DisplayName','CS Trace');
ylabel(ax,'Amplitude'); xlabel(ax,'Time');
legend(ax);

%% --- Threshold UI controls ---
uicontrol('Style','text','String','Threshold1:', 'Units','normalized', ...
    'Position',[0.1 0.25 0.1 0.05],'HorizontalAlignment','right');
th1Edit = uicontrol('Style','edit','String',num2str(CS_thres(1)), 'Units','normalized', ...
    'Position',[0.21 0.25 0.1 0.05]);

uicontrol('Style','text','String','Threshold2:', 'Units','normalized', ...
    'Position',[0.35 0.25 0.1 0.05],'HorizontalAlignment','right');
th2Edit = uicontrol('Style','edit','String',num2str(CS_thres(2)), 'Units','normalized', ...
    'Position',[0.46 0.25 0.1 0.05]);

%% --- Additional parameters UI ---
uicontrol('Style','text','String','# Spikes:', 'Units','normalized', ...
    'Position',[0.1 0.18 0.1 0.05],'HorizontalAlignment','right');
nSpikeEdit = uicontrol('Style','edit','String',num2str(N_Spike), 'Units','normalized', ...
    'Position',[0.21 0.18 0.1 0.05]);

uicontrol('Style','text','String','Max Amp Thres:', 'Units','normalized', ...
    'Position',[0.35 0.18 0.1 0.05],'HorizontalAlignment','right');
maxAmpEdit = uicontrol('Style','edit','String',num2str(Max_AUC), 'Units','normalized', ...
    'Position',[0.46 0.18 0.1 0.05]);

uicontrol('Style','text','String','# Spikes for ISI:', 'Units','normalized', ...
    'Position',[0.6 0.18 0.13 0.05],'HorizontalAlignment','right');
nSpike2ISIEdit = uicontrol('Style','edit','String',num2str(N_Spike2ISI), 'Units','normalized', ...
    'Position',[0.74 0.18 0.1 0.05]);

%% --- Buttons ---
uicontrol('Style','pushbutton','String','Update','Units','normalized', ...
    'Position',[0.6 0.08 0.1 0.06],'Callback',@updateThresholds);
uicontrol('Style','pushbutton','String','Done','Units','normalized', ...
    'Position',[0.75 0.08 0.1 0.06],'Callback','uiresume(gcbf);');

% Initial update
updateThresholds();
uiwait(f);

% Final thresholds after closing
CS_thres = [str2double(th1Edit.String), str2double(th2Edit.String)];

if isvalid(f)
    close(f);
end

%% --- Nested update function ---
    function updateThresholds(~,~)
        % Read values from UI
        CS_thres = [str2double(th1Edit.String), str2double(th2Edit.String)];
        N_Spike = str2double(nSpikeEdit.String);
        Max_AUC = str2double(maxAmpEdit.String);
        N_Spike2ISI = str2double(nSpike2ISIEdit.String);

        % Detect transients using user settings
        [trans, tr_trace] = detect_transient2(tr_sub, CS_thres, sp_soma, 20);
        if isempty(trans.amp)
            CS_trace = zeros(1, nTime);
        else
            transcand = cell2mat(cellfun(@(x) length(x) > 1, trans.ISI, 'UniformOutput', false));
            meanISI_frnt = cellfun(@(x) mean(x(1:min(N_Spike2ISI-1,end))), trans.ISI(transcand));
            meanISI_first3 = zeros(1,length(trans.length));
            meanISI_first3(transcand) = meanISI_frnt;

            CS_ind = find(trans.spike_number >= N_Spike & ...
                meanISI_first3 < 20 & trans.int > Max_AUC);
            CS_trace = ismember(tr_trace, CS_ind);
        end

        % Adjust trace
        bwCS = bwlabel(CS_trace);
        CS_spike = sp_soma .* bwCS;
        [~, CS_spike_time] = unique(CS_spike);
        for b = 1:max(bwCS)
            bfrm = find(bwCS == b);
            CS_trace(bfrm(1):CS_spike_time(bwCS(CS_spike_time) == b) - 5) = 0;
        end

        tr_hi_CS = NaN(1,nTime);
        tr_hi_CS(CS_trace>0) = tr_hi(CS_trace>0);

        % Update plot
        set(plotCS, 'XData', 1:nTime, 'YData', tr_hi_CS);
        drawnow;
    end
end
