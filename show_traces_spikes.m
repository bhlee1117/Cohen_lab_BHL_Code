function ax1 = show_traces_spikes(traces, spikes, otherT, t)
% SHOW_TRACES_SPIKES  Plot traces with spike overlays and optional bottom panel.
%
%   ax1 = show_traces_spikes(traces, spikes)
%   ax1 = show_traces_spikes(traces, spikes, otherT)
%   ax1 = show_traces_spikes(traces, spikes, otherT, t)
%
%   INPUTS
%     traces  : [nROI x nTime] trace matrix
%     spikes  : [nROI x nTime] binary spike matrix
%     otherT  : (optional) [nTrace x nTime] traces for bottom panel.
%               Pass [] to skip the bottom panel.
%     t       : (optional) [1 x nTime] time axis vector. Default: frame index.
%
%   OUTPUT
%     ax1 : handle to the main (top) axes

noi     = 1:size(traces, 1);
frmrate = 1;

% --- defaults ---
if nargin < 3,              otherT = [];                        end
if nargin < 4 || isempty(t), t = (1:size(traces,2)) / frmrate; end

showOtherT = ~isempty(otherT);

clf;

% --- main trace panel ---
if showOtherT
    ax1 = subplot(10, 1, 1:8);
else
    ax1 = axes;
end

scale = max(prctile(traces, 99, 2)) * 1.5;
tr    = traces;

lines = plot(t, tr(noi,:)' + (1:numel(noi)) * scale); %#ok<NASGU>
arrayfun(@(l,c) set(l,'Color',c{:}),lines,num2cell(turbo(size(traces,1)),2))
hold all;
drawnow limitrate;

if any(tovec(spikes == 1))
    S = tr;
    S(~(spikes == 1)) = NaN;
    plot(t, S(noi,:)' + (1:numel(noi)) * scale, 'r.');
    set(gca, 'ytick', (1:numel(noi)) * scale, 'yticklabel', noi);
end

% --- optional bottom panel ---
if showOtherT
    ax1 = [ax1 subplot(10, 1, 9:10)];

    if size(otherT, 1) > size(otherT, 2)
        otherT = otherT';
    end

    for j = 1:size(otherT, 1)
        plot(t, otherT(j,:));
        hold all;
    end

    linkaxes([ax1], 'x');
end

end