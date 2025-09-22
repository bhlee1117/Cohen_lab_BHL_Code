function filt_trace = interactive_bandstop_filter(trace, fs)
% INTERACTIVE_BANDSTOP_FILTER - Visual tool for interactive bandstop filtering.
%   filt_trace = interactive_bandstop_filter(trace, fs)
%   Inputs:
%       trace : 1 x T signal
%       fs    : sampling frequency (Hz)
%   Output:
%       filt_trace : filtered version of input trace

    maxBands = 5;
    bandCount = 1;
    bandLimits = zeros(maxBands, 2);
    bandLimits(1,:) = [58 62]; % Default initial band

    T = length(trace);
    t = (0:T-1)/fs;

    % FFT for power spectrum
    NFFT = 2^nextpow2(T);
    fVec = fs/2*linspace(0,1,NFFT/2+1);
    Pxx = abs(fft(trace,NFFT)/T).^2;
    Pxx = Pxx(1:NFFT/2+1);

    % Create figure (non-overlapping)
    screenSize = get(0, 'ScreenSize');
    figWidth = 900; figHeight = 500;
    figLeft = screenSize(3)/2 - figWidth/2 + rand*100 - 50; % small offset to avoid stacking
    figBottom = screenSize(4)/2 - figHeight/2 + rand*100 - 50;
    f = figure('Name','Interactive Bandstop Filter','NumberTitle','off',...
        'Position',[figLeft figBottom figWidth figHeight]);

    ax2 = subplot(2,1,1); % Trace (moved to top)
    hTrace = plot(ax2, t, trace, 'w');
    hold(ax2, 'on'); hFilt = plot(ax2, t, trace, 'r'); hold(ax2, 'off');
    xlabel(ax2, 'Time (s)'); ylabel(ax2, 'Amplitude'); title(ax2, 'Filtered Trace');
    legend(ax2, 'Original', 'Filtered');

    ax1 = subplot(2,1,2); % Power spectrum (moved to bottom)
    hPow = plot(ax1, fVec, 10*log10(Pxx), 'b');
    ylabel(ax1, 'Power (dB)'); title(ax1, 'Power Spectrum');

    % UI for bands
    bandTexts = gobjects(maxBands,2);
    bandEdits = gobjects(maxBands,2);
    for i = 1:maxBands
        for j = 1:2
        lh_str={'Low','High'};
            bandTexts(i,j) = uicontrol('Style','text','String',sprintf('Band %d [%s]:', i, lh_str{j}),...
                'Units','normalized','Position',[0.05+0.2*(j-1) 0.94-0.05*i 0.08 0.03],'HorizontalAlignment','right');
            bandEdits(i,j) = uicontrol('Style','edit','String',num2str(bandLimits(i,j)),...
                'Units','normalized','Position',[0.14+0.2*(j-1) 0.94-0.05*i 0.08 0.03]);
        end
    end

    % Increase/decrease buttons
    uicontrol('Style','pushbutton','String','+ Band','Units','normalized','Position',[0.45 0.94 0.08 0.04],...
              'Callback',@addBand);
    uicontrol('Style','pushbutton','String','- Band','Units','normalized','Position',[0.55 0.94 0.08 0.04],...
              'Callback',@removeBand);

    % Update & Done buttons
    uicontrol('Style','pushbutton','String','Update','Units','normalized','Position',[0.7 0.94 0.1 0.04],...
              'Callback',@updateFilter);
    uicontrol('Style','pushbutton','String','Done','Units','normalized','Position',[0.82 0.94 0.1 0.04],...
              'Callback','uiresume(gcbf);');

    updateFilter();
    uiwait(f);

    if isvalid(f), close(f); end

    function updateFilter(~,~)
        bands = zeros(bandCount,2);
        for i = 1:bandCount
            bands(i,1) = str2double(bandEdits(i,1).String);
            bands(i,2) = str2double(bandEdits(i,2).String);
        end
        d = trace;
        for i = 1:bandCount
            [b,a] = butter(2, bands(i,:)/(fs/2), 'stop');
            d = filtfilt(b,a,d);
        end
        filt_trace = d;
        set(hFilt, 'YData', filt_trace);
        drawnow;
    end

    function addBand(~,~)
        if bandCount < maxBands
            bandCount = bandCount + 1;
            set([bandTexts(bandCount,:) bandEdits(bandCount,:)], 'Visible','on');
        end
    end

    function removeBand(~,~)
        if bandCount > 1
            set([bandTexts(bandCount,:) bandEdits(bandCount,:)], 'Visible','off');
            bandCount = bandCount - 1;
        end
    end
end